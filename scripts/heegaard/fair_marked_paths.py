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


def setup(d, *, knot_refinement=None):
    knot_refinement={} if knot_refinement is None else knot_refinement
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
        factor=knot_refinement.get((b['curve'],b['vertices'][0],b['vertices'][-1]),1)
        assert type(factor) is int and 1<=factor<=16
        n = max(6, int(np.ceil(length / (55/factor))) + 3)
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
            c[inner] = c[end] + direction * min(18., length / 7, length/(n-3)/3)
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


def _separator_seed_binding(seed):
    # Source: the SAME original seed passed to the bounded constructor.
    import hashlib
    return hashlib.sha256(json.dumps(seed,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def source_separations_from_failure(seed,candidate,error):
    """Recover actual source intervals, then verify a SOURCE strict separator.

    Source: existing correspondence's optional diagnostic labels and Shapely
    nearest_points. Floating closest points propose a normal; exact rational
    projections of the actual old moving edges must certify a positive gap.
    """
    from fractions import Fraction as Q
    from shapely.geometry import LineString
    from shapely.ops import nearest_points
    from check_fairing_isotopy import correspondence
    payload=error.args[0] if len(error.args)==1 else error.args
    if not isinstance(payload,tuple) or not payload or payload[0]!='Interpolation collision':return []
    if not isinstance(payload[-1],dict) or 'source_parameters' not in payload[-1]:return []
    a,_,_,names,signatures,parameters=correspondence(seed,candidate,with_source_parameters=True)
    i,j=map(int,payload[1:3]);assert [names[i],names[j]]==list(payload[3:5])
    assert [parameters[i],parameters[j]]==payload[-1]['source_parameters']
    assert [signatures[i],signatures[j]]==[tuple(x) for x in payload[-1]['branch_signatures']]
    points=[parameters[k]['source_points_svg'] for k in (i,j)]
    x,y=nearest_points(LineString(points[0]),LineString(points[1]));normal=np.array(y.coords[0])-np.array(x.coords[0])
    length=np.linalg.norm(normal)
    if not np.isfinite(length) or length==0:return []
    normal=tuple(Q(float(v)) for v in normal/length)
    old=[[tuple(Q(float(v)) for v in point) for point in path] for path in points]
    project=lambda point:sum(x*y for x,y in zip(normal,point))
    left=max(map(project,old[0]));right=min(map(project,old[1]))
    if not left<right:return []
    offset=(left+right)/2;margin=(right-left)/4;binding=_separator_seed_binding(seed)
    return [dict(seed_sha256=binding,branch=list(signatures[k]),interval=list(parameters[k]['interval']),
        normal=list(map(str,normal)),offset=str(offset),margin=str(margin),sign=sign,
        source_points_svg=[[str(v) for v in point] for point in points])
        for k,sign,points in ((i,-1,old[0]),(j,1,old[1]))]


def _separator_penalty_rows(bs,records):
    """Existing BSpline identity basis/derivatives give restricted cubic rows."""
    from fractions import Fraction as Q
    result=[];lookup={(b['curve'],b['vertices'][0],b['vertices'][-1]):b for b in bs}
    for record in records:
        b=lookup[tuple(record['branch'])];lo,hi=record['interval'];assert 0<=lo<hi<=1
        knots=sorted(set([lo,hi]+[float(t) for t in b['knots'] if lo<t<hi]))
        basis=b['basis'];normal=np.array([float(Q(x)) for x in record['normal']])
        for x,y in zip(knots,knots[1:]):
            dt=(y-x)/3
            rows=np.array([basis(x),basis(x)+dt*basis(x,nu=1),basis(y)-dt*basis(y,nu=1),basis(y)])
            result.append((slice(b['offset'],b['offset']+b['count']),rows,normal,
                float(Q(record['offset'])),float(Q(record['margin'])),record['sign']))
    return result


def _separator_penalty(controls,rows,weight):
    value=0.;gradient=np.zeros_like(controls)
    for indices,matrix,normal,offset,margin,sign in rows:
        points=matrix@controls[indices];height=sign*(points@normal-offset)
        deficit=np.maximum(0,margin-height)
        value+=weight*np.sum(deficit**2)
        gradient[indices]+=matrix.T@((-2*weight*sign*deficit)[:,None]*normal)
    return value,gradient


def _restrict_exact_cubic(points,lo,hi):
    # Source: verify_bezier_attachment's exact rational de Casteljau blend/split,
    # here at the actual emitted span restriction endpoints rather than 1/2.
    from fractions import Fraction as Q
    def split(q,t):
        rows=[q]
        while len(rows[-1])>1:
            rows.append([tuple(x+t*(y-x) for x,y in zip(a,b)) for a,b in zip(rows[-1],rows[-1][1:])])
        return [row[0] for row in rows],[row[-1] for row in reversed(rows)]
    q=[tuple(Q(x) for x in p) for p in points];lo,hi=Q(lo),Q(hi);assert 0<=lo<hi<=1
    if lo:q=split(q,lo)[1]
    if hi<1:q=split(q,(hi-lo)/(1-lo))[0]
    return q


def verify_source_separations(seed,candidate,records):
    """Construction feasibility only; ordinary smooth/full motion gates follow.

    Actual exported cubic controls/spans, not optimiser flags, certify the
    strict halfplanes. This is sufficient for these local interval pairs only.
    """
    from fractions import Fraction as Q
    if not records:return
    assert all(r['seed_sha256']==_separator_seed_binding(seed) for r in records)
    _,_,bs=branches(seed);edge_branch={(b['curve'],u,v):(b['curve'],b['vertices'][0],b['vertices'][-1])
        for b in bs for u,v in zip(b['vertices'],b['vertices'][1:])}
    old={(c['id'],tuple(s['source'])):s for c in seed['curves'] for s in c['segments']}
    for record in records:
        lo,hi=map(Q,record['interval']);normal=list(map(Q,record['normal']));offset=Q(record['offset']);pieces=[]
        for curve in candidate['curves']:
            for segment in curve['segments']:
                source=old[curve['id'],tuple(segment['source'])]
                if edge_branch[curve['id'],*source['vertices'][:2]]!=tuple(record['branch']):continue
                spans=segment['spline_parameter_intervals'];assert len(spans)==len(segment['beziers'])
                for (a,b),q in zip(spans,segment['beziers']):
                    a,b=Q(a),Q(b);x,y=max(lo,a),min(hi,b)
                    if x<y:pieces.append((x,y,_restrict_exact_cubic(q,(x-a)/(b-a),(y-a)/(b-a))))
        pieces.sort(key=lambda x:x[0]);assert pieces and pieces[0][0]==lo and pieces[-1][1]==hi
        assert all(a[1]==b[0] for a,b in zip(pieces,pieces[1:]))
        assert all(record['sign']*(sum(x*y for x,y in zip(normal,p))-offset)>0 for _,_,q in pieces for p in q), (
            'Source separator construction infeasible',dict(branch_signatures=[tuple(record['branch'])]))


def optimise(d, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20., maxiter=1800, knot_refinement=None, source_separations=None):
    pos, xs, bs, c, fixed, B, B1, B2, target, bid, arc = setup(d,knot_refinement=knot_refinement)
    source_separations=[] if source_separations is None else source_separations
    assert all(r['seed_sha256']==_separator_seed_binding(d) for r in source_separations)
    separator_rows=_separator_penalty_rows(bs,source_separations)
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
            # Source: cyclic source-arclength distances above. Adjacent samples
            # need a tapered gap; omitting all pairs <=35 hides short folds.
            gaps=np.where(bid[ii] != bid[jj],target_gap,np.minimum(target_gap,along/2))
            mask = (bid[ii] != bid[jj]) | (along > 0)
            # Adjacent samples of one branch are not separate strands.
            ii, jj, gaps = ii[mask], jj[mask], gaps[mask]
            mask = np.ones(len(ii), dtype=bool)
            for vertex in xs:
                incident = np.array([vertex in b['ends'] for b in bs])
                near = np.linalg.norm(q - pos[vertex], axis=1) < crossing_radius+.1
                mask &= ~(incident[bid[ii]] & incident[bid[jj]] & (near[ii] & near[jj]))
            ii, jj, gaps = ii[mask], jj[mask], gaps[mask]
            dif = q[ii] - q[jj]
            dist = np.linalg.norm(dif, axis=1)
            deficit = np.maximum(0, gaps - dist)
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
        if separator_rows:
            extra,extra_gradient=_separator_penalty(c,separator_rows,collision_weight)
            value+=extra;grad+=extra_gradient
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
            # Source: the actual refined BSpline knot partition used above.
            # Motion correspondence must use these spans, not recreate /55.
            s['spline_parameter_intervals']=[[float(lo),float(hi)] for lo,hi in zip(knots,knots[1:])]
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




def offending_fairing_branches(seed, candidate, error, *, crossing_radius=20.):
    """Recover original branch signatures from actual existing judge failures.

    Source: export's retained source intervals, measure_clearance.branches,
    verifier edge/control-leaf order and motion check's branch diagnostics.
    Unknown mathematical/source-map failures are never interpreted as fairing.
    """
    _,xs,bs=branches(seed)
    keys={(b['curve'],b['vertices'][0],b['vertices'][-1]):b for b in bs}
    edges={(b['curve'],u,v):k for k,b in keys.items() for u,v in zip(b['vertices'],b['vertices'][1:])}
    oldsegments={(c['id'],tuple(s['source'])):s for c in seed['curves'] for s in c['segments']}
    def segment_key(cid,s):
        old=oldsegments[cid,tuple(s['source'])]
        return edges[cid,*old['vertices'][:2]]
    payload=error.args[0] if len(error.args)==1 else error.args
    selected=set();contexts=[]
    if isinstance(payload,tuple) and payload:
        name=payload[0]
        if name=='Extra intersection':
            for cid,u,v in (payload[3],payload[4]):
                hits=[s for c in candidate['curves'] if c['id']==cid for s in c['segments'] if any((a,b)==(u,v) for a,b in zip(s['vertices'],s['vertices'][1:]))]
                if len(hits)!=1:return set()
                selected.add(segment_key(cid,hits[0]))
            from shapely import wkt
            contexts.append((set(selected),wkt.loads(payload[-1])))
        elif name=='Non-regular spline piece':
            cid,source=payload[1],payload[2]
            hits=[s for c in candidate['curves'] if c['id']==cid for s in c['segments'] if s['source']==source]
            if len(hits)!=1:return set()
            selected.add(segment_key(cid,hits[0]))
        elif name=='Bezier hull collision':
            from check_smooth_marked import leaves
            controls=[]
            for c in candidate['curves']:
                for s in c['segments']:
                    controls.extend([segment_key(c['id'],s)]*sum(sum(1 for _ in leaves(np.array(q))) for q in s['beziers']))
            selected.update(controls[i] for i in (payload[1],payload[2]))
            from shapely import wkt
            contexts.append((set(selected),wkt.loads(payload[-1])))
        elif name in ('Interpolation collision','Incident edges collapse','Swept hull meets hole','Source separator construction infeasible') and isinstance(payload[-1],dict):
            selected.update(tuple(k) for k in payload[-1].get('branch_signatures',[]))
            if payload[-1].get('location_svg') is not None:
                from shapely.geometry import Point
                contexts.append((set(selected),Point(payload[-1]['location_svg'])))
    elif isinstance(payload,dict) and 'violations' in payload and 'pairs' in payload:
        _,_,current=branches(candidate)
        for pair in payload['pairs']:
            group=set()
            for i in pair['branches']:
                b=current[i];key=(b['curve'],b['vertices'][0],b['vertices'][-1]);selected.add(key);group.add(key)
            from shapely.geometry import Point
            contexts.append((group,Point(np.mean(np.array(pair['points']),axis=0))))
    if not selected or not selected<=keys.keys():return set()
    # Source: actual rejected component and marked incidence. An interior
    # self-collision belongs to its one source branch. Only DISTINCT rejected
    # branches sharing a true crossing near this actual location justify
    # expanding that crossing's star; unrelated far endpoints do not.
    from shapely.geometry import Point
    touched=set()
    for group,location in contexts:
        if len(group)<2:continue
        shared=set.intersection(*(set(keys[k]['ends']) for k in group))&xs.keys()
        for v in shared:
            if location.distance(Point(np.array(seed['positions'][str(v)])*600))<=crossing_radius+.1:touched.add(v)
    selected.update(k for k,b in keys.items() if b['ends']&touched)
    return selected




def fairing_option_bindings(seed):
    """Bind options to actual native/prepared data and schema-defined geometry.

    Source: the existing source-map canonical binary64 geometry payload. No
    saved output controls or coordinates can enter this construction option.
    """
    import hashlib
    from marked_display_source_map import _source_map_json
    from marked_display_source_map_pipeline import source_map_geometry_payload
    context=seed.get('_source_map_input')
    if not isinstance(context,dict) or 'prepared' not in context or 'source_document' not in seed:
        raise ValueError('Fairing options require the actual native/preparation context')
    digest=lambda value:hashlib.sha256(_source_map_json(value)).hexdigest()
    auxiliary={k:seed.get(k) for k in ('genus','cuts','boundary_vertex_pairs','port_transport','preparation_representation_sha256','infinity_source_triangle')}
    return dict(native_document_sha256=digest(seed['source_document']),
        preparation_sha256=digest(context['prepared']),seed_geometry_sha256=digest(source_map_geometry_payload(seed)),
        seed_auxiliary_sha256=digest(auxiliary))


def load_fairing_construction_options(seed, path, *, policy):
    """Read ONLY source-bound integer knot factors, then run the normal fitter.

    Source: existing branches and native/prepared/source-geometry bindings.
    Strict fields exclude any supplied output positions, controls or curves.
    This finite option never upgrades the general adaptive default's status.
    """
    import hashlib
    path=Path(path);raw=path.read_bytes();option=json.loads(raw)
    fields={'schema','native_document_sha256','preparation_sha256','seed_geometry_sha256','seed_auxiliary_sha256','producer_sha256','policy','knot_refinement'}
    if not isinstance(option,dict) or set(option)!=fields or option['schema']!='marked-fairing-construction-options/1':
        raise ValueError('Invalid parameters-only fairing option schema')
    actual=fairing_option_bindings(seed)
    if any(option[k]!=v for k,v in actual.items()):raise ValueError('Fairing option native/preparation/seed binding mismatch')
    names={'fair_marked_paths.py','check_fairing_isotopy.py','check_smooth_marked.py','marked_display_pipeline.py'}
    if not isinstance(option['producer_sha256'],dict) or set(option['producer_sha256'])!=names:
        raise ValueError('Fairing option producer pin fields differ')
    for name in names:
        if hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()!=option['producer_sha256'][name]:
            raise ValueError('Fairing option producer changed: '+name)
    if option['policy']!=policy:raise ValueError('Fairing option numerical policy changed')
    _,_,source_branches=branches(seed)
    valid={(b['curve'],b['vertices'][0],b['vertices'][-1]) for b in source_branches}
    factors={}
    if not isinstance(option['knot_refinement'],list):raise ValueError('Invalid branch factor list')
    for record in option['knot_refinement']:
        if not isinstance(record,dict) or set(record)!={'curve','start','end','factor'}:
            raise ValueError('Invalid branch factor record')
        key=(record['curve'],record['start'],record['end']);factor=record['factor']
        if type(key[0]) is not str or type(key[1]) is not int or type(key[2]) is not int or key not in valid or key in factors:
            raise ValueError('Unknown or duplicate actual source branch signature')
        if type(factor) is not int or factor not in (1,2,4,8,16):raise ValueError('Knot factor outside fixed dyadic construction bounds')
        factors[key]=factor
    return factors,dict(schema=option['schema'],option_sha256=hashlib.sha256(raw).hexdigest(),**actual,
        producer_sha256=option['producer_sha256'],policy=policy,knot_refinement=option['knot_refinement'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    args=p.parse_args();root=args.directory
    source=json.loads((root/'KodaSeparated.json').read_text(encoding='utf-8'))
    bs,c,initial=optimise(source)
    d=export(source,bs,c)
    (root/'KodaSmoothCandidate.json').write_text(json.dumps(d),encoding='utf-8')
    (root/'KodaSmoothCandidate.svg').write_text(render(d),encoding='utf-8')
    print(json.dumps(measure(d, threshold=13),indent=2),flush=True)
