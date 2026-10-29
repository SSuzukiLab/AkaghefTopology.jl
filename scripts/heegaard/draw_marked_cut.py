"""Draw a source-marked surface after explicit auxiliary cuts, without alpha-cut projection.

The corner quotient follows marked_reference.regions; every displayed edge retains
its source triangle occurrence. This first entry point uses the supplied Koda data.
"""
import argparse
import hashlib
import json
import math
import sys
from html import escape
from collections import defaultdict
from pathlib import Path
from fractions import Fraction

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import networkx as nx


def spread(pos, triangles, pins, caps, *, preserve_linear_orientation=False):
    ids=sorted(pos);ix={v:i for i,v in enumerate(ids)}
    P=np.array([pos[v] for v in ids]);T=np.array([[ix[v] for v in t] for t in triangles])
    E=np.array(sorted({tuple(sorted((ix[u],ix[v]))) for t in triangles for u,v in zip(t,t[1:]+t[:1])}))
    free=np.array([ix[v] for v in ids if v not in pins])
    W=np.array([8. if set(t)&caps else 1. for t in triangles]);W*=30/W.sum()
    def fun(x):
        Q=P.copy();Q[free]=x.reshape((-1,2));a,b,c=(Q[T[:,j]] for j in range(3))
        ab=b-a;ac=c-a;areas=ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0]
        if min(areas)<=0:return 1e100,np.zeros_like(x)
        delta=Q[E[:,0]]-Q[E[:,1]];value=float((delta*delta).sum()-np.dot(W,np.log(areas)))
        grad=np.zeros_like(Q);np.add.at(grad,E[:,0],2*delta);np.add.at(grad,E[:,1],-2*delta)
        for j,u,v in ((0,b,c),(1,c,a),(2,a,b)):
            g=np.column_stack((u[:,1]-v[:,1],v[:,0]-u[:,0]))*(-W/areas)[:,None]
            np.add.at(grad,T[:,j],g)
        return value,grad[free].ravel()
    x=P[free].ravel().copy();value,grad=fun(x);history=[]
    for iteration in range(2500):
        q=grad.copy();coeff=[]
        for s,y,r in reversed(history):
            alpha=r*np.dot(s,q);coeff.append(alpha);q-=alpha*y
        if history:
            s,y,_=history[-1];q*=np.dot(s,y)/np.dot(y,y)
        else:q/=max(1.,np.linalg.norm(q))
        for (s,y,r),alpha in zip(history,reversed(coeff)):q+=s*(alpha-r*np.dot(y,q))
        direction=-q;slope=float(np.dot(grad,direction));step=1.
        if slope>=0:direction=-grad/max(1.,np.linalg.norm(grad));slope=float(np.dot(grad,direction));history=[]
        for trial in range(70):
            candidate=x+step*direction;nv,ng=fun(candidate)
            acceptable=nv<value+1e-4*step*slope
            if acceptable and preserve_linear_orientation:
                # Source: check_fairing_isotopy's polynomial interval method.
                # A line-search endpoint alone cannot exclude an intermediate
                # triangle flip. This numerical guard is a proposal filter;
                # the retained full chart is independently checked exactly.
                start=P.copy();end=P.copy()
                start[free]=x.reshape((-1,2));end[free]=candidate.reshape((-1,2))
                a=start[T];delta=end[T]-a
                u=a[:,1]-a[:,0];v=a[:,2]-a[:,0]
                du=delta[:,1]-delta[:,0];dv=delta[:,2]-delta[:,0]
                determinant=lambda s,t:s[:,0]*t[:,1]-s[:,1]*t[:,0]
                c=determinant(u,v);b=determinant(du,v)+determinant(u,dv);q=determinant(du,dv)
                t=np.where(q>0,np.clip(-b/np.where(q!=0,2*q,1),0,1),0)
                acceptable=float(np.min(np.minimum(np.minimum(c,c+b+q),c+t*b+t*t*q)))>0
            if acceptable:break
            step*=.5
        else:break
        s=candidate-x;y=ng-grad;sy=float(np.dot(s,y))
        if sy>1e-16:history.append((s,y,1/sy));history=history[-20:]
        x,value,grad=candidate,nv,ng
        if np.linalg.norm(grad)<1e-6:break
    P[free]=x.reshape((-1,2))
    print('Coarse layout:',len(ids),'vertices;',iteration,'iterations',flush=True)
    return {v:P[i] for v,i in ix.items()}


def build(source, kernel, output, prepared_path):
    sys.path.insert(0, str(kernel))
    import marked_reference as m
    doc = m.read_document(source)
    st = m.validate(doc)
    prepared=json.loads(prepared_path.read_text(encoding='utf-8'))
    assert prepared['schema']=='marked-display-complex/1'
    assert prepared['source_file_sha256']==hashlib.sha256(source.read_bytes()).hexdigest(), 'Stale prepared source'
    cuts = set(prepared['cut_curves'])
    blocked = {e for e, owners in st.used.items() if set(owners) & cuts}
    audit = m.regions(st, blocked)
    assert len(audit)==1 and audit[0]['genus']==0 and audit[0]['boundary_components']==2*prepared['genus']
    rs=[dict(boundary_cycles=[b['source'] for b in prepared['boundaries']])]
    print('Prepared cut verified: connected genus 0,',2*prepared['genus'],'boundaries', flush=True)
    ep = m.fine_edge_occurrences(st)
    corners = [(f, i, j) for f, i in st.triangles for j in range(3)]
    vd = m.DSU(corners)
    ed = m.DSU(list(ep))
    def cut(a):
        return (a[0], m.edge(a[2], a[3])) in blocked
    def corner(a, v):
        return (a[0], a[1], st.triangles[a[:2]].index(v))
    for a, b in ep.items():
        if not cut(a) and not cut(b):
            ed.union(a,b)
            vd.union(corner(a, a[2]), corner(b, b[3]))
            vd.union(corner(a, a[3]), corner(b, b[2]))
    roots = sorted({vd.find(c) for c in corners})
    index = {r: i for i, r in enumerate(roots)}
    cv = {c: index[vd.find(c)] for c in corners}
    native_cv={tuple([parts[0],int(parts[1]),int(parts[2])]):v for key,v in prepared['source_cut_corner_classes'].items() for parts in [key.split(',')]}
    assert cv==native_cv, 'Julia/Python source-corner quotient disagreement'
    cv=native_cv
    tris = {tuple(t['source']):tuple(t['vertices']) for t in prepared['triangles']}
    def occurrence(a):
        return [cv[corner(a, a[2])], cv[corner(a, a[3])]]
    def newnode(label):
        v=len(roots);roots.append(label);return v
    edge_nodes={}
    def midpoint(a):
        if a not in ep:a=(*a[:2],a[3],a[2])
        k=ed.find(a)
        if k not in edge_nodes:edge_nodes[k]=newnode(('edge',k))
        return edge_nodes[k]
    def full_occurrence(a):
        u,v=occurrence(a)
        return [u,midpoint(a),v]
    boundaries = []
    for cyc in rs[0]['boundary_cycles']:
        ids = [v for a in cyc for v in full_occurrence(tuple(a))[:2]]
        labels = {o for a in cyc for o in st.used[(a[0], m.edge(a[2], a[3]))]}
        assert len(labels) == 1
        metadata=next(b for b in prepared['boundaries'] if b['source']==cyc)
        boundaries.append(dict(cut=next(iter(labels)), vertices=ids, source=cyc,
                               side=metadata['side'],row=metadata['row'],column=metadata['column']))
    # An unoccupied original face becomes the chart's outer face, not a puncture.
    marked_faces={t[0] for mark in prepared['region_marks'] for t in mark['source_region']['triangles']}
    occupied_faces={d.face for _,ds in st.curves.values() for d in ds}
    candidates = [f for f in st.patches if f not in marked_faces and f not in occupied_faces]
    assert candidates, 'No unused source patch for the exterior chart; add an explicit auxiliary chart disk'
    outer = max(candidates, key=lambda f: len(st.boundaries[f]))
    local_nodes = defaultdict(set)
    for (f, i), t in st.triangles.items():
        for j, v in enumerate(t):
            local_nodes[(f, v)].add(cv[(f, i, j)])
    outer_ids=[]
    ob=st.boundaries[outer]
    for u,v in zip(ob,ob[1:]+ob[:1]):
        oc=st.incidence[outer][m.edge(u,v)][0]
        outer_ids.extend(full_occurrence((outer,oc[0],u,v))[:2])
    # Cap the six cut boundaries temporarily. These triangles are never drawn.
    # Flag subdivision avoids degeneracy from two-edge boundary circles and
    # parallel source edges. Midpoints are keyed by edge occurrence, not endpoints.
    mesh=[];weights=[]
    for (f,i),t in tris.items():
        if f==outer:continue
        center=newnode(('triangle',f,i))
        original=st.triangles[(f,i)]
        for u,v in m.directed_triangle(original):
            a,mid,b=full_occurrence((f,i,u,v))
            mesh.extend(((a,mid,center),(mid,b,center)))
            weights.extend((1.,1.))
    for b in boundaries:
        center = newnode(('cap',b['cut']))
        ids = b['vertices']
        b['cap'] = center
        mesh.extend((ids[(i+1) % len(ids)], ids[i], center) for i in range(len(ids)))
        weights.extend([len(st.triangles)/len(ids)]*len(ids))
    adj = defaultdict(set)
    for t in mesh:
        for a, b in zip(t, t[1:]+t[:1]):
            adj[a].add(b); adj[b].add(a)
    # Boundary order is reversed: exterior chart disk is on the other side.
    fixed = {v: np.array([math.cos(-2*math.pi*i/len(outer_ids)), math.sin(-2*math.pi*i/len(outer_ids))]) for i,v in enumerate(outer_ids)}
    inside = sorted(set(adj)-set(fixed)); ii = {v:i for i,v in enumerate(inside)}
    A = lil_matrix((len(inside), len(inside))); B = np.zeros((len(inside),2))
    for v,i in ii.items():
        A[i,i] = len(adj[v])
        for w in adj[v]:
            if w in ii: A[i,ii[w]] = -1
            else: B[i] += fixed[w]
    sol = spsolve(A.tocsr(), B)
    pos = {**fixed, **{v:sol[i] for v,i in ii.items()}}
    # Integer grid realization of the SAME oriented rotation system avoids
    # exponential crowding of the harmonic embedding. No planar re-embedding.
    rotations=defaultdict(dict)
    for t in mesh:
        for j,u in enumerate(t):rotations[u][t[(j+1)%3]]=t[(j+2)%3]
    # Complete the exterior face around the removed patch.
    for j,u in enumerate(outer_ids):
        rotations[u][outer_ids[(j+1)%len(outer_ids)]]=outer_ids[(j-1)%len(outer_ids)]
    cyclic={}
    for u,successor in rotations.items():
        v=min(successor);ns=[]
        while v not in ns:ns.append(v);v=successor[v]
        assert len(ns)==len(successor)
        cyclic[u]=ns
    embedding=nx.PlanarEmbedding();embedding.set_data(cyclic);embedding.check_structure()
    grid=nx.combinatorial_embedding_to_pos(embedding,fully_triangulate=False)
    vals=np.array(list(grid.values()),dtype=float);lo=vals.min(axis=0);hi=vals.max(axis=0)
    pos={v:(np.array(p)-(lo+hi)/2)/(max(hi-lo)/2) for v,p in grid.items()}
    def area(t):
        a,b,c = (pos[v] for v in t)
        return float(np.linalg.det(np.array([b-a,c-a])))
    areas = [area(t) for t in mesh]
    if max(areas)<0:
        pos={v:np.array([p[0],-p[1]]) for v,p in pos.items()}
        areas=[area(t) for t in mesh]
    assert min(areas) > 1e-15, (min(areas), max(areas), outer)
    # Every source curve dart maps to a unique edge occurrence after the cut.
    curves = []
    for c in prepared['curves']:
        cid=c['id']
        segments = []
        for segment in c['segments']:
            f,u,v=segment['source']
            occ = st.incidence[f][m.edge(u,v)][0]
            a = (f,occ[0],u,v)
            # The occurrence has to follow the curve direction.
            path=full_occurrence(a)
            assert [path[0],path[-1]]==segment['vertices']
            segments.append(dict(vertices=path, source=segment['source']))
        curves.append(dict(id=cid,role=c['role'],segments=segments))
    crossings = []
    for x in doc['crossings']:
        if all(st.curves[c][0]['role'] in ('alpha','beta') for c in x['curves']):
            # crossing_at is the already validated local address index.
            addr = next(k for k,v in st.crossing_at.items() if v['id']==x['id'])
            ns=local_nodes[addr]; assert len(ns)==1
            crossings.append(dict(x, vertex=next(iter(ns))))
    assert len(crossings)==len(prepared['crossings'])
    assert len(prepared['region_marks'])==1, 'This SVG view currently displays one joining disk'
    disk = prepared['region_marks'][0]['source_region']
    assert len(prepared['region_marks'][0]['cut_pieces'])==1, 'This geometry adapter does not yet draw split D; prepared data retains all pieces'
    assert disk['embedded_closed_disk']
    disk_boundary = [v for a in disk['boundary_cycles'][0] for v in full_occurrence(tuple(a))[:2]]
    # Remove only unmarked triangulation scaffolding. Keep alpha, beta and O
    # together with the exact cyclic orders induced by the source cut mesh.
    angles=defaultdict(dict)
    def graph_edge(u,v,at_u,at_v):
        assert v not in angles[u]
        angles[u][v]=math.atan2(*(at_u-pos[u])[::-1])
        angles[v][u]=math.atan2(*(at_v-pos[v])[::-1])
    for cid,(c,ds) in st.curves.items():
        if cid in cuts:continue
        for dart in ds:
            oc=st.incidence[dart.face][m.edge(dart.u,dart.v)][0]
            u,mid,v=full_occurrence((dart.face,oc[0],dart.u,dart.v))
            graph_edge(u,v,pos[mid],pos[mid])
    ports=[]
    for b in boundaries:
        for v in b['vertices']:
            if v not in angles:continue
            assert len(angles[v])==1
            graph_edge(v,b['cap'],pos[b['cap']],pos[v])
            ports.append((v,b['cap']))
    reduced=nx.PlanarEmbedding()
    reduced.set_data({u:sorted(ns,key=ns.get,reverse=True) for u,ns in angles.items()})
    reduced.check_structure()
    from networkx.algorithms.planar_drawing import triangulate_embedding
    triangulated,outer_cycle=triangulate_embedding(reduced,fully_triangulate=True)
    gp=nx.combinatorial_embedding_to_pos(triangulated,fully_triangulate=True)
    allfaces=[];seen=set()
    for u,v in triangulated.edges():
        if (u,v) not in seen:allfaces.append(tuple(triangulated.traverse_face(u,v,seen)))
    forbidden=set(disk_boundary)|{b['cap'] for b in boundaries}|{c['vertex'] for c in crossings}
    eligible=[t for t in allfaces if not set(t)&forbidden]
    assert eligible,'No unmarked exterior triangle'
    capset={b['cap'] for b in boundaries}
    base_edges={tuple(sorted((u,v))) for u in angles for v in angles[u]}
    def dist_seg(p,a,b):
        ab=b-a;t=float(np.dot(p-a,ab)/np.dot(ab,ab));return float(np.linalg.norm(p-(a+np.clip(t,0,1)*ab)))
    def quality(pp):
        radii=[min(dist_seg(pp[cap],pp[u],pp[v]) for u,v in base_edges if cap not in (u,v)) for cap in capset]
        cr=[pp[c['vertex']] for c in crossings]
        gaps=[np.linalg.norm(a-b) for j,a in enumerate(cr) for b in cr[j+1:]]
        return min(radii)+.2*min(gaps)+.003*sum(math.log(max(r,1e-20)) for r in radii)
    def harmonic(outer_cycle):
        ids=sorted(gp);free_ids=[v for v in ids if v not in outer_cycle];ix={v:i for i,v in enumerate(free_ids)}
        pin={v:np.array([math.cos(-2*math.pi*j/3),math.sin(-2*math.pi*j/3)]) for j,v in enumerate(outer_cycle)}
        mat=lil_matrix((len(ix),len(ix)));rhs=np.zeros((len(ix),2))
        for v,i in ix.items():
            ns=list(triangulated.neighbors(v));mat[i,i]=len(ns)
            for w in ns:
                if w in ix:mat[i,ix[w]]=-1
                else:rhs[i]+=pin[w]
        coords=spsolve(mat.tocsr(),rhs)
        pp={**pin,**{v:coords[i] for v,i in ix.items()}}
        ff=[t for t in allfaces if set(t)!=set(outer_cycle)]
        a,b,c=(pp[v] for v in ff[0])
        if np.linalg.det(np.array([b-a,c-a]))<0:pp={v:np.array([p[0],-p[1]]) for v,p in pp.items()}
        return pp,ff
    candidates=[]
    for oc in eligible:
        pp,ff=harmonic(oc);candidates.append((quality(pp),oc,pp,ff))
    best=None
    for _,oc,pp,ff in sorted(candidates,key=lambda x:x[0],reverse=True)[:5]:
        pp=spread(pp,ff,set(oc),capset);score=quality(pp)
        if best is None or score>best[0]:best=(score,oc,pp,ff)
    _,outer_cycle,pos,faces=best
    print('Selected readability score',best[0],flush=True)
    assert min(area(t) for t in faces)>0
    for c in curves:
        for s in c['segments']:s['vertices']=s['vertices'][::2]
    disk_boundary=disk_boundary[::2]
    # Extend the cut collars radially into the temporary caps, then excise a
    # round disk. This does not move or remove any alpha-beta intersection.
    reduced_edges={tuple(sorted((u,v))) for u in angles for v in angles[u]}
    def distance(p,a,b):
        ab=b-a;t=float(np.dot(p-a,ab)/np.dot(ab,ab));return float(np.linalg.norm(p-(a+np.clip(t,0,1)*ab)))
    for b in boundaries:
        center=pos[b['cap']]
        clearance=min(distance(center,pos[u],pos[v]) for u,v in reduced_edges if b['cap'] not in (u,v))
        b['radius']=.38*clearance
        b['ports']=[]
        for v,cap in ports:
            if cap!=b['cap']:continue
            p=center+(pos[v]-center)*b['radius']/np.linalg.norm(pos[v]-center)
            new=max(pos)+1;pos[new]=p;b['ports'].append(dict(vertex=new,source_vertex=v))
            for c in curves:
                for s in c['segments']:
                    if s['vertices'][0]==v:s['vertices'].insert(0,new)
                    if s['vertices'][-1]==v:s['vertices'].append(new)
    def crossing_sign(x):
        rays={};v=x['vertex']
        for c in curves:
            if c['id'] not in x['curves']:continue
            for s in c['segments']:
                vs=s['vertices']
                if vs[0]==v:rays[c['id']+'+']=pos[vs[1]]-pos[v]
                if vs[-1]==v:rays[c['id']+'-']=pos[vs[-2]]-pos[v]
        order=sorted(rays,key=lambda k:math.atan2(rays[k][1],rays[k][0]))
        a,b=x['curves'];k=order.index(a+'+');order=order[k:]+order[:k]
        assert order[2]==a+'-' and order[1].startswith(b) and order[3].startswith(b)
        return 1 if order[1]==b+'+' else -1
    observed=[crossing_sign(x) for x in crossings]
    if observed==[-x['sign'] for x in crossings]:
        pos={v:np.array([p[0],-p[1]]) for v,p in pos.items()}
        faces=[tuple(reversed(t)) for t in faces]
    assert [crossing_sign(x) for x in crossings]==[x['sign'] for x in crossings]
    # Check the actual rounded SVG coordinates, not just the unrounded solver.
    qp={v:(Fraction(f'{650+540*p[0]:.8f}'),Fraction(f'{620-540*p[1]:.8f}')) for v,p in pos.items()}
    segs=[(c['id'],u,v) for c in curves for s in c['segments'] for u,v in zip(s['vertices'],s['vertices'][1:])]
    for i,(cid,u,v) in enumerate(segs):
        for did,a,b in segs[i+1:]:
            if {u,v}&{a,b}:continue
            assert not m.segments_meet(qp[u],qp[v],qp[a],qp[b]),('Extra curve crossing',cid,did,u,v,a,b)
    positive=[area(t) for t in faces];assert min(positive)>0
    print('Geometry verified:',len(crossings),'signed crossings; no additional curve intersections after SVG rounding',flush=True)
    data = dict(schema='marked-cut-view/1',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                source=str(source),cuts=prepared['cut_curves'],genus=prepared['genus'],cut_genus=0,boundary_components=2*prepared['genus'],
                preparation_sha256=hashlib.sha256(prepared_path.read_bytes()).hexdigest(),grid=prepared['grid'],
                cut_reference_outer_face=outer,chart_outer_triangle=list(outer_cycle),positions={str(v):p.tolist() for v,p in pos.items()},
                boundaries=boundaries,curves=curves,crossings=crossings,disk=disk_boundary,
                source_cut_corner_classes={','.join(map(str,c)):v for c,v in cv.items()},
                triangles=[list(t) for t in faces],
                minimum_double_area=min(positive),verified_signed_crossings=len(crossings),extra_curve_intersections=0,
                joining_disk=m.check_joining_disk(st,prepared['region_marks'][0]['seed']))
    output.mkdir(parents=True,exist_ok=True)
    (output/'KodaMarkedCut.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    svg = render(data)
    (output/'KodaMarkedCut.svg').write_text(svg,encoding='utf-8')
    print(json.dumps({k:data[k] for k in ('chart_outer_triangle','minimum_double_area','joining_disk')}),flush=True)


def attachment_boundary_path(d, part):
    """Read the complete source-bound cut interval in the literal chart domain.

    Source: source-map full native-lambda boundary trace and exact CDT constraint
    subdivisions. This path is consumed by that SAME polygon-to-circle chart;
    it is not a new endpoint route or an independently centered SVG arc.
    """
    layer=d['source_map_boundary_paths']
    assert layer['schema']=='marked-display-attachment-boundaries/1'
    b=next(b for b in layer['boundaries'] if (b['cut'],b['side'])==(part['cut'],part['side']))
    u,v=part['from_source_vertex'],part['to_source_vertex'];ports=b['source_ports']
    start,end=ports[str(u)],ports[str(v)];points=b['points_svg'];indices=[start]
    while indices[-1]!=end or len(indices)==1:
        indices.append((indices[-1]+1)%len(points));assert len(indices)<=len(points)+1
    # Verify the actual native edge word, including a possible full boundary
    # traversal. Endpoint equality alone does not specify which cut interval.
    native=b['native_boundary_vertices'];i=native.index(u);word=[u];darts=[]
    while word[-1]!=v or len(word)==1:
        darts.append(b['native_boundary_darts'][i])
        i=(i+1)%len(native);word.append(native[i]);assert len(word)<=len(native)+1
    assert part['source_darts']==darts, 'D cut sweep differs from retained native interval'
    intervals={(x['from_source_vertex'],x['to_source_vertex']):x for x in b['consecutive_port_intervals']}
    native_ports=[x for x in word if str(x) in ports]
    certified=[]
    for a,z in zip(native_ports,native_ports[1:]):
        chain=intervals[a,z]['polygon_vertices']
        assert not certified or certified[-1]==chain[0]
        certified.extend(chain if not certified else chain[1:])
    assert certified==indices, 'D polygon interval differs from certified full boundary trace'
    return np.asarray([points[i] for i in indices],dtype=float)


def disk_boundary_cycles(d, piece):
    """Choose a rendering start without changing a directed cyclic boundary."""
    cycles=piece['boundary_cycles']
    if d.get('svg_coordinate_encoding')!='float64-roundtrip':return cycles
    out=[]
    for cycle in cycles:
        start=next((i for i,part in enumerate(cycle) if part['kind']=='curve'),None)
        if start is None:raise ValueError('A joining-disk boundary containing only auxiliary cuts is unsupported by the current marked preparation')
        out.append(cycle[start:]+cycle[:start])
    return out


def disk_piece_cycles(d):
    """Resolve source-bound D pieces from shared curve paths and cut-side arcs.

    Source: marked_display.jl::prepare_marked_display display_boundary_cycles.
    Surface-region darts have the interior on the left, hence cut-side arcs run
    clockwise in the oriented plane. Their end ports and wrap are fixed by the
    saved side order; no straight chord or alternate endpoint route is chosen.
    Arc samples have sagitta <= .003 SVG units, matching spline flattening.
    """
    pos={int(v):np.asarray(q,dtype=float)*600 for v,q in d['positions'].items()}
    if d.get('svg_coordinate_encoding')=='float64-roundtrip':
        for v,q in d['source_endpoint_normalization']['shared_svg_endpoints'].items():pos[int(v)]=np.asarray(q,dtype=float)
    segments={(c['id'],tuple(s['source'])):s for c in d['curves'] for s in c['segments']}
    boundaries={(b['cut'],b['side']):b for b in d['boundaries']}
    out=[]
    for piece in d.get('disk_pieces',[]):
        cycles=[]
        for cycle in disk_boundary_cycles(d,piece):
            points=[]
            for part in cycle:
                if part['kind']=='curve':
                    s=segments[part['curve'],tuple(part['source'])]
                    vs=s['vertices'] if part['direction']==1 else s['vertices'][::-1]
                    q=[pos[v] for v in vs]
                else:
                    assert part['kind']=='cut','Unknown fill boundary interval'
                    b=boundaries[part['cut'],part['side']]
                    ports={p['source_vertex']:p['vertex'] for p in b['ports']}
                    u=ports[part['from_source_vertex']];v=ports[part['to_source_vertex']]
                    if d.get('svg_coordinate_encoding')=='float64-roundtrip':
                        q=list(attachment_boundary_path(d,part))
                        assert np.array_equal(q[0],pos[u]) and np.array_equal(q[-1],pos[v]), 'D attachment endpoints differ from actual curves'
                    else:
                        center=pos[b['cap']];r=b['radius']*600
                        a=math.atan2(*(pos[u]-center)[::-1]);z=math.atan2(*(pos[v]-center)[::-1])
                        sweep=-((a-z)%(2*math.pi))
                        assert -2*math.pi+1e-10<sweep<-1e-10,'Degenerate cut-circle interval'
                        n=max(1,math.ceil(abs(sweep)/(2*math.acos(max(-1.,1-.003/r)))))
                        q=[center+r*np.array([math.cos(t),math.sin(t)]) for t in np.linspace(a,a+sweep,n+1)]
                        q[0]=pos[u];q[-1]=pos[v]
                assert not points or np.linalg.norm(points[-1]-q[0])<1e-6,'Disconnected D piece boundary'
                points.extend(q if not points else q[1:])
            assert len(points)>3 and np.linalg.norm(points[-1]-points[0])<1e-6,'Open D piece boundary'
            cycles.append(np.asarray(points[:-1]))
        out.append((piece,cycles))
    return out


def orientation_glyph_layer(d, svg):
    """Display-only arrows from the saved directed route, never edge sorting.

    Source: layout_marked_document.coarse_view expands edge_sign into the
    segment's ordered vertices; fair_marked_paths.export keeps that direction
    in each cubic's control points. Shapely hull/label exclusion selects an
    uncluttered glyph position without modifying any mathematical geometry.
    """
    import xml.etree.ElementTree as ET
    from shapely.geometry import Point, Polygon, MultiPoint, box
    from shapely.strtree import STRtree
    ns={'s':'http://www.w3.org/2000/svg'};root=ET.fromstring(svg)
    obstacles=[]
    for text in root.findall('s:text',ns):
        x,y=float(text.get('x')),float(text.get('y'));size=float(text.get('font-size',13))
        width=max(size*.65*len(text.text or ''),size)
        x0=x-width/2 if text.get('text-anchor')=='middle' else x
        obstacles.append(box(x0-4,y-size-4,x0+width+4,y+5))
    for b in d['boundaries']:
        x,y=d['positions'][str(b['cap'])]
        obstacles.append(Point(600*x,-600*y).buffer(600*b['radius']+8))
    for c in d['crossings']:
        x,y=d['positions'][str(c['vertex'])];obstacles.append(Point(600*x,-600*y).buffer(18))
    hulls=[];owners=[];candidates={}
    for ci,c in enumerate(d['curves']):
        choices=[]
        for s in c['segments']:
            if s.get('beziers'):
                pieces=[np.asarray(q)*[1,-1] for q in s['beziers']]
            else:
                points=[np.asarray(d['positions'][str(v)])*[600,-600] for v in s['vertices']]
                pieces=[np.array([a,a+(b-a)/3,a+2*(b-a)/3,b]) for a,b in zip(points,points[1:])]
            for q in pieces:
                hulls.append(MultiPoint(q).convex_hull);owners.append(ci)
                p=(q[0]+3*q[1]+3*q[2]+q[3])/8
                tangent=3*(q[3]+q[2]-q[1]-q[0])/4;length=np.linalg.norm(tangent)
                if length>1e-10:choices.append((float(np.linalg.norm(q[-1]-q[0])),p,tangent/length))
        candidates[ci]=sorted(choices,key=lambda x:-x[0])[:128]
    tree=STRtree(hulls) if hulls else None;out=[]
    for ci,c in enumerate(d['curves']):
        best=None
        for _,p,tangent in candidates[ci]:
            normal=np.array([-tangent[1],tangent[0]])
            q=np.array([p+5*tangent,p-4*tangent+3*normal,p-4*tangent-3*normal]);glyph=Polygon(q)
            if any(glyph.intersects(o) for o in obstacles):continue
            near=tree.query(glyph.buffer(60)) if tree is not None else []
            distances=[glyph.distance(hulls[j]) for j in near if owners[j]!=ci]
            clearance=min(distances,default=60.)
            if clearance<3:continue
            score=min(clearance,min((glyph.distance(o) for o in obstacles),default=60.))
            if best is None or score>best[0]:best=(score,q,glyph)
        if best is None:raise ValueError('No unobstructed orientation glyph position for '+c['id'])
        _,q,glyph=best;obstacles.append(glyph.buffer(5))
        points=' '.join(f'{x:.8f},{y:.8f}' for x,y in q)
        color=escape(c.get('color',{'alpha':'#ba433f','beta':'#255dba','aux':'#357e64'}[c['role']]),quote=True)
        out.append(f'<g data-curve="{escape(c["id"],quote=True)}" data-orientation="source-directed"><polygon points="{points}" fill="{color}" stroke="{color}" stroke-width=".5"/></g>')
    return ''.join(out)


def render(d, *, orientation_glyphs=False):
    scale=600
    pos={int(v):(scale*p[0],-scale*p[1]) for v,p in d['positions'].items()}
    roundtrip=d.get('svg_coordinate_encoding')=='float64-roundtrip'
    if roundtrip:
        # Source: marked_display_pipeline.normalize_shared_endpoints. These
        # are the same actual SVG-unit endpoints used by the certified cubics;
        # avoid a second divide/multiply rounding at D/crossing/port readback.
        for v,q in d['source_endpoint_normalization']['shared_svg_endpoints'].items():
            pos[int(v)]=(float(q[0]),-float(q[1]))
    def coordinate(x):return repr(float(x)) if roundtrip else f'{x:.8f}'
    visible={v for c in d['curves'] for s in c['segments'] for v in s['vertices']}
    bounds=[pos[v] for v in visible]
    for b in d['boundaries']:
        x,y=pos[b['cap']];r=scale*b['radius']
        bounds.extend([(x-r,y-r),(x+r,y+r)])
    if not bounds:
        # The empty sphere still has its explicitly selected chart at infinity.
        bounds=[pos[v] for v in d['chart_outer_boundary']]
    left=min(p[0] for p in bounds)-45;top=min(p[1] for p in bounds)-75
    width=max(p[0] for p in bounds)-left+45;height=max(p[1] for p in bounds)-top+45
    def points(vs): return ' '.join(f'{coordinate(pos[v][0])},{coordinate(pos[v][1])}' for v in vs)
    def bezier_path(pieces):
        def xy(p):return f'{coordinate(p[0])},{coordinate(-p[1])}'
        if roundtrip:
            assert all(a[-1]==b[0] for a,b in zip(pieces,pieces[1:])), 'Unnormalized SVG cubic join'
        return 'M'+xy(pieces[0][0])+' '+' '.join('C'+' '.join(xy(p) for p in q[1:]) for q in pieces)
    disk_shape=(f'<path data-disk="D" d="{bezier_path(d["disk_beziers"])} Z" fill="#f1c85a" fill-opacity=".5"/>'
                if d.get('disk_beziers') else
                f'<polygon data-disk="D" points="{points(d["disk"])}" fill="#f1c85a" fill-opacity=".5"/>' if d.get('disk') else '')
    # Source: disk_piece_cycles above. SVG's even-odd rule also represents the
    # infinity-containing component: the viewport is a chart window, not a new
    # physical boundary. The old single-D branch above remains byte-identical.
    piece_cycles=disk_piece_cycles(d)
    segments={(c['id'],tuple(s['source'])):s for c in d['curves'] for s in c['segments']}
    cut_sides={(b['cut'],b['side']):b for b in d['boundaries']}
    for piece,cycles in piece_cycles:
        path=(f'M{left},{top} h{width} v{height} h{-width} Z ' if piece['contains_infinity'] else '')
        for parts,cycle in zip(disk_boundary_cycles(d,piece),cycles):
            if roundtrip:
                first=parts[0]
                assert first['kind']=='curve','A D piece must start on its retained curve interval'
                first_segment=segments[first['curve'],tuple(first['source'])]
                first_vertex=first_segment['vertices'][0 if first['direction']==1 else -1]
                path+=f'M{coordinate(pos[first_vertex][0])},{coordinate(pos[first_vertex][1])} '
            else:path+=f'M{cycle[0][0]:.8f},{-cycle[0][1]:.8f} '
            for part in parts:
                if part['kind']=='curve':
                    s=segments[part['curve'],tuple(part['source'])]
                    vs=s['vertices'] if part['direction']==1 else s['vertices'][::-1]
                    if s.get('beziers'):
                        be=s['beziers'] if part['direction']==1 else [q[::-1] for q in s['beziers'][::-1]]
                        path+=bezier_path(be).split(' ',1)[1]+' '
                    else:path+=' L'+points(vs[1:])+' '
                else:
                    if roundtrip:
                        chain=attachment_boundary_path(d,part)
                        path+=' '.join(f'L{coordinate(x)},{coordinate(-y)}' for x,y in chain[1:])+' '
                        continue
                    b=cut_sides[part['cut'],part['side']];ports={p['source_vertex']:p['vertex'] for p in b['ports']}
                    u=pos[ports[part['from_source_vertex']]];v=pos[ports[part['to_source_vertex']]];c=pos[b['cap']]
                    a=math.atan2(u[1]-c[1],u[0]-c[0]);z=math.atan2(v[1]-c[1],v[0]-c[0]);sweep=(z-a)%(2*math.pi)
                    r=600*b['radius'];path+=f'A{coordinate(r)},{coordinate(r)} 0 {int(sweep>math.pi)} 1 {coordinate(v[0])},{coordinate(v[1])} '
            path+='Z '
        disk_shape+=f'<path data-disk="{escape(str(piece["mark"]),quote=True)}" d="{path}" fill="#f1c85a" fill-opacity=".5" fill-rule="evenodd"/>'
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left} {top} {width} {height}" role="img" aria-label="Marked surface: genus {d["genus"]}, {len(d["boundaries"])} punctures, {len(d["crossings"])} interior crossings">',
         f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="#fcfbf7"/>',
         '<style>text{font-family:Arial,sans-serif;fill:#253247}.label{paint-order:stroke;stroke:#fcfbf7;stroke-width:4;stroke-linejoin:round}</style>',
         f'<text x="{left+20}" y="{top+35}" font-size="23">{escape(str(d.get("title", "Heegaard diagram · g = "+str(d["genus"]))))}</text>',
         disk_shape]
    # Source: native periodic 2-chain/compiler, actual closed cubic boundaries.
    from periodic_annotations import layer
    out.append(layer(d))
    colors={'alpha':'#ba433f','beta':'#255dba','aux':'#357e64'}
    for c in d['curves']:
        out.append(f'<g data-curve="{escape(c["id"],quote=True)}" fill="none" stroke="{escape(c.get("color",colors[c["role"]]),quote=True)}" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">')
        for s in c['segments']:
            if s.get('beziers'):out.append(f'<path d="{bezier_path(s["beziers"])}"/>')
            else:out.append(f'<polyline points="{points(s["vertices"])}"/>')
        out.append('</g>')
        if c.get('label_at'):
            x,y=c['label_at'];out.append(f'<text class="label" x="{600*x}" y="{-600*y}" text-anchor="middle" font-size="18">{escape(c["id"])}</text>')
    counts=defaultdict(int)
    for b in d['boundaries']:
        counts[b['cut']]+=1; name=b.get('label',b['cut']+('⁺' if b.get('row',counts[b['cut']]-1)==0 else '⁻'))
        x,y=pos[b['cap']]
        out.append(f'<circle data-puncture="{escape(name,quote=True)}" cx="{x}" cy="{y}" r="{scale*b["radius"]}" fill="white" stroke="#39444f" stroke-width="1.5"/>')
        out.append(f'<text class="label" x="{x}" y="{y}" text-anchor="middle" font-size="13">{escape(name)}</text>')
    for c in d['crossings']:
        x,y=pos[c['vertex']]
        out.append(f'<circle data-crossing="{escape(c["id"],quote=True)}" cx="{x}" cy="{y}" r="2.4" fill="#182538"/><text class="label" x="{x+6}" y="{y-6}" font-size="13">{escape(c["id"])}</text>')
    from shapely.geometry import Point, Polygon
    from shapely.ops import polylabel
    for p in d.get('port_labels',[]):
        x,y=pos[p['vertex']];out.append(f'<circle cx="{x}" cy="{y}" r="3" fill="#253247"/><text class="label" x="{x+p.get("dx",7)}" y="{y+p.get("dy",-8)}" font-size="12">{escape(p["label"])}</text>')
    if d.get('disk'):
        label=polylabel(Polygon([pos[v] for v in d['disk']]),tolerance=2)
        x,y=label.x,label.y
        out.append(f'<text x="{x}" y="{y}" font-size="20">D</text>')
    if piece_cycles:
        from shapely.geometry import box
        for piece,cycles in piece_cycles:
            fill=box(left,-top-height,left+width,-top) if piece['contains_infinity'] else Polygon()
            for cycle in cycles:fill=fill.symmetric_difference(Polygon(cycle))
            for b in d['boundaries']:
                x,y=pos[b['cap']];fill=fill.difference(Point(x,-y).buffer(600*b['radius'],quad_segs=128))
            if fill.geom_type=='MultiPolygon':fill=max(fill.geoms,key=lambda p:p.area)
            label=polylabel(fill,tolerance=2)
            out.append(f'<text x="{label.x}" y="{-label.y}" font-size="20">D</text>')
    out.append('</svg>')
    svg=''.join(out)
    return svg[:-6]+orientation_glyph_layer(d,svg)+'</svg>' if orientation_glyphs else svg


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--kernel',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--prepared',type=Path,required=True,help='Julia-produced marked-display-complex/1; no inferred marking')
    a=p.parse_args();build(a.source,a.kernel,a.output,a.prepared)
