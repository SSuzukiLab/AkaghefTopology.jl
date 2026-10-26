"""Global spline fairing with fixed marked crossings and circular cut boundaries.

Optimisation is a proposal, never the acceptance test. Source edge words are
retained as ordered parameter intervals; D uses those same intervals. Reuses
measure_clearance.branches and publish_marked_display.bind's source-edge/D
contract. scipy.interpolate.BSpline supplies the basis and exact derivatives;
check_smooth_marked.py checks the emitted Bezier control hulls independently.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.sparse import block_diag
from scipy.spatial import cKDTree

from measure_clearance import branches, measure
from draw_marked_cut import render


def unit(v):
    return v / np.linalg.norm(v)


def setup(d):
    pos, xs, bs = branches(d)
    ports = {p['vertex']: pos[b['cap']] for b in d['boundaries'] for p in b['ports']}
    rays = {}
    for b in bs:
        q = np.array([pos[v] for v in b['vertices']])
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]
        b['arc'] = arc
        b['length'] = arc[-1]
        b['old'] = q
        for end in (0, -1):
            v = b['vertices'][end]
            if v in xs:
                # The old compact maps make the immediate rays transverse.
                r = q[1] - q[0] if end == 0 else q[-2] - q[-1]
                rays.setdefault((b['curve'], v), []).append((end, unit(r)))
    tangents = {}
    for key, values in rays.items():
        assert len(values) == 2 and {x[0] for x in values} == {0, -1}
        tangents[key] = unit(sum(r if e == 0 else -r for e, r in values))
    controls, mats, d1s, d2s, targets = [], [], [], [], []
    branch_ids, params = [], []
    offset = 0
    fixed = []
    for i, b in enumerate(bs):
        length = b['length']
        n = max(6, int(np.ceil(length / 55)) + 3)
        knots = np.r_[np.zeros(4), np.linspace(0, 1, n - 2)[1:-1], np.ones(4)]
        basis = BSpline(knots, np.eye(n), 3)
        t = np.linspace(0, 1, max(20, int(np.ceil(length / 3)) + 1))
        target = np.column_stack([np.interp(t * length, b['arc'], b['old'][:, k]) for k in (0, 1)])
        B = basis(t)
        c = np.zeros((n, 2))
        c[0], c[-1] = b['old'][0], b['old'][-1]
        for end, inner in ((0, 1), (-1, -2)):
            v = b['vertices'][end]
            if v in xs:
                direction = tangents[b['curve'], v] * (1 if end == 0 else -1)
            elif b['vertices'][0]==b['vertices'][-1]:
                direction=unit(b['old'][1]-b['old'][-2])*(1 if end==0 else -1)
            else:
                direction = unit(pos[v] - ports[v])
            c[inner] = c[end] + direction * min(18., length / 7)
        pinned = [0, 1, n-2, n-1]
        free = list(range(2, n-2))
        c[free] = np.linalg.lstsq(B[:, free], target - B[:, pinned] @ c[pinned], rcond=None)[0]
        fixed.extend(offset + k for k in pinned)
        b.update(knots=knots, count=n, offset=offset, t=t, basis=basis)
        offset += n
        controls.append(c)
        mats.append(B)
        d1s.append(basis(t, nu=1) / length)
        d2s.append(basis(t, nu=2) / length**2)
        targets.append(target)
        branch_ids.extend([i] * len(t))
        params.extend(t * length)
    return pos, xs, bs, np.vstack(controls), np.array(fixed), block_diag(mats).tocsr(), block_diag(d1s).tocsr(), block_diag(d2s).tocsr(), np.vstack(targets), np.array(branch_ids), np.array(params)


def optimise(d, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20., maxiter=1800):
    pos, xs, bs, c, fixed, B, B1, B2, target, bid, arc = setup(d)
    free = np.setdiff1d(np.arange(len(c)), fixed)
    holes = [(pos[b['cap']],b['radius']*600,{p['vertex'] for p in b['ports']}) for b in d['boundaries']]
    target_gap=minimum_separation+2.5
    # At common crossings only the documented radius-20 neighbourhood is exempt.
    def objective(x, collision_weight, tether):
        c[free] = x.reshape(-1, 2)
        q, v, acc = B @ c, B1 @ c, B2 @ c
        delta = q - target
        value = 1600 * np.sum(acc**2) + .08 * np.sum(v**2) + tether * np.sum(delta**2)
        grad = 3200 * (B2.T @ acc) + .16 * (B1.T @ v) + 2*tether*(B.T @ delta)
        qgrad = np.zeros_like(q)
        pairs = cKDTree(q).query_pairs(target_gap+1.5, output_type='ndarray')
        if len(pairs):
            ii, jj = pairs.T
            along=abs(arc[ii]-arc[jj])
            for k,b in enumerate(bs):
                if b['vertices'][0]==b['vertices'][-1]:
                    cyclic=(bid[ii]==k)&(bid[jj]==k)
                    along[cyclic]=np.minimum(along[cyclic],b['length']-along[cyclic])
            mask = (bid[ii] != bid[jj]) | (along > 35)
            # Adjacent samples of one branch are not separate strands.
            ii, jj = ii[mask], jj[mask]
            mask = np.ones(len(ii), dtype=bool)
            for vertex in xs:
                incident = np.array([vertex in b['ends'] for b in bs])
                near = np.linalg.norm(q - pos[vertex], axis=1) < crossing_radius+.1
                mask &= ~(incident[bid[ii]] & incident[bid[jj]] & (near[ii] | near[jj]))
            ii, jj = ii[mask], jj[mask]
            dif = q[ii] - q[jj]
            dist = np.linalg.norm(dif, axis=1)
            deficit = np.maximum(0, target_gap - dist)
            value += collision_weight * np.sum(deficit**2)
            force = (-2*collision_weight*deficit / np.maximum(dist, 1e-8))[:, None] * dif
            np.add.at(qgrad, ii, force)
            np.add.at(qgrad, jj, -force)
        for center,radius,port_ids in holes:
            dif = q - center
            dist = np.linalg.norm(dif, axis=1)
            margin=np.full(len(q),puncture_clearance+1.5)
            margin_gradient=np.zeros_like(q)
            # Only a branch incident to this exact port gets a tapered collar.
            # Every other branch must keep the full puncture clearance.
            for v in port_ids:
                incident=np.array([v in (b['vertices'][0],b['vertices'][-1]) for b in bs])[bid]
                radial=q-pos[v];distance=np.linalg.norm(radial,axis=1)
                taper=.75*distance;active=incident & (taper<margin)
                margin[active]=taper[active]
                margin_gradient[active]=.75*radial[active]/np.maximum(distance[active,None],1e-8)
            deficit = np.maximum(0, radius+margin - dist)
            value += 10*collision_weight*np.sum(deficit**2)
            qgrad += (-20*collision_weight*deficit / np.maximum(dist, 1e-8))[:, None] * dif
            qgrad += (20*collision_weight*deficit)[:,None]*margin_gradient
        grad += B.T @ qgrad
        return value, grad[free].ravel()
    initial = c.copy()
    for weight, tether in ((1., .003), (10., .001), (50., .0001), (200., .00002)):
        result = minimize(objective, c[free].ravel(), args=(weight, tether), jac=True,
                          method='L-BFGS-B', options=dict(maxiter=maxiter, ftol=1e-11, gtol=1e-6, maxls=40))
        c[free] = result.x.reshape(-1, 2)
        print('fairing', weight, result.nit, result.fun, result.message, flush=True)
    return bs, c.copy(), initial


def export(d, bs, controls):
    d = copy.deepcopy(d)
    pos = {int(v): np.array(q)*600 for v, q in d['positions'].items()}
    old_disk = d.get('disk',[])
    edge_paths = {}
    next_id = max(pos) + 1
    lookup = {}
    for b in bs:
        spline = BSpline(b['knots'], controls[b['offset']:b['offset']+b['count']], 3)
        for v, t in zip(b['vertices'], b['arc']/b['length']):
            pos[v] = spline(t)
        # Only source edge endpoints are semantic stops; all other old samples
        # are replaced by a flatness-bounded subdivision of this entire spline.
        for u, v in zip(b['vertices'], b['vertices'][1:]):
            lookup[b['curve'], u, v] = b
        b['spline'] = spline
        b['parameter'] = {v:t for v,t in zip(b['vertices'], b['arc']/b['length'])}
    for curve in d['curves']:
        for s in curve['segments']:
            old = s['vertices']; u, v = old[0], old[-1]
            b = lookup[curve['id'], old[0], old[1]]
            spline = b['spline']; a, z = b['parameter'][u], b['parameter'][v]
            if b['vertices'][0]==b['vertices'][-1]:
                if u==b['vertices'][0]:a=0.
                if v==b['vertices'][0]:z=1.
            # Each original source edge is inside a crossing/port branch.
            assert z > a
            # Source interval boundaries can coincide with a spline knot up to
            # floating roundoff after reloading. Do not emit zero-length pieces.
            knots = sorted(set([a,z] + [float(t) for t in b['knots'] if a+1e-12<t<z-1e-12]))
            bez=[]
            for lo,hi in zip(knots,knots[1:]):
                q0,q3=spline(lo),spline(hi);dt=(hi-lo)/3
                bez.append(np.array([q0,q0+dt*spline(lo,nu=1),q3-dt*spline(hi,nu=1),q3]))
            samples=[(u,pos[u])]
            def flatten(q):
                # Convex-hull bound relative to the chord's linear parameter.
                err=max(np.linalg.norm(q[1]-(2*q[0]+q[3])/3),np.linalg.norm(q[2]-(q[0]+2*q[3])/3))
                if err <= .003:
                    samples.append((None,q[3]));return
                a=(q[:-1]+q[1:])/2;b=(a[:-1]+a[1:])/2;c=(b[0]+b[1])/2
                flatten(np.array([q[0],a[0],b[0],c]));flatten(np.array([c,b[1],a[2],q[3]]))
            for q in bez:flatten(q)
            vertices=[u]
            for _,q in samples[1:-1]:
                pos[next_id]=q;vertices.append(next_id);next_id+=1
            vertices.append(v)
            s['vertices']=vertices
            s['beziers']=[q.tolist() for q in bez]
            # Map each old microscopic edge to its interval for D reconstruction.
            for a0,z0 in zip(old,old[1:]):edge_paths[a0,z0]=(s,old)
    # D travels entire source edges, possibly in reverse. Resolve through old
    # directed incidence and use exactly the same new curve samples/Bezier data.
    steps=list(zip(old_disk,old_disk[1:]+old_disk[:1]));disk=[];disk_bez=[];last=None
    for u,v in steps:
        if u==v:continue
        forward=(u,v) in edge_paths
        s,old=edge_paths[(u,v) if forward else (v,u)]
        key=(id(s),forward)
        if key==last:continue
        vs=s['vertices'] if forward else s['vertices'][::-1]
        be=s['beziers'] if forward else [q[::-1] for q in s['beziers'][::-1]]
        if disk:assert disk[-1]==vs[0]
        disk.extend(vs if not disk else vs[1:]);disk_bez.extend(be);last=key
    if disk:assert disk[0]==disk[-1]
    d['disk']=disk[:-1] if disk else [];d['disk_beziers']=disk_bez
    d['positions']={str(v):(q/600).tolist() for v,q in pos.items()}
    d['spline_fairing']=dict(method='global cubic B-spline bending energy with strand and circular-hole penalties',
        crossing_positions_fixed=True,port_positions_fixed=True,flatness_bound_svg=.003,
        control_count=len(controls),crossing_exemption_radius=20)
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    args=p.parse_args();root=args.directory
    source=json.loads((root/'KodaSeparated.json').read_text(encoding='utf-8'))
    bs,c,initial=optimise(source)
    d=export(source,bs,c)
    (root/'KodaSmoothCandidate.json').write_text(json.dumps(d),encoding='utf-8')
    (root/'KodaSmoothCandidate.svg').write_text(render(d),encoding='utf-8')
    print(json.dumps(measure(d, threshold=13),indent=2),flush=True)
