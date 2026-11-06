"""Realize the native cut cellulation; never route curves independently.

Input comes from Heegaard.prepare_marked_display. The caller names the source
triangle containing infinity. A collar inserted INSIDE that triangle gives a
convex exterior frame, without replacing the exterior region. SciPy realizes
the whole mesh; marked_cut_layout then moves its boundaries to the fixed grid.
The same source-edge paths carry disconnected curves, crossings and D.

Source: src/heegaard/marked_display.jl::prepare_marked_display supplies the
oriented quotient; draw_marked_cut.py::build supplies the capped harmonic mesh
pattern. Flag subdivision and the collar inside the named infinity triangle
are this adapter's construction, not a new route inferred from curve endpoints.
"""
import argparse
import copy
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from marked_cut_layout import fixed_mesh
from verify_fixed_grid import verify
from draw_marked_cut import render


def rotation_grid_embedding(triangles,outer):
    """Reuse the existing whole-mesh integer drawing with its exact outer face.

    Source: draw_marked_cut.py::build, its PlanarEmbedding/set_data and
    networkx.combinatorial_embedding_to_pos calls; the existing backend runtime
    is selected by marked_display_plot.jl::layout_marked.
    Here the unique nontriangular face is the caller's retained infinity frame;
    it is verified before the library draws it, never selected for readability.
    The source triangle orientation and all edge paths stay unchanged.
    """
    import networkx as nx
    rotation=defaultdict(dict)
    for t in triangles:
        for j,u in enumerate(t):
            v,w=t[(j+1)%3],t[(j+2)%3]
            assert v not in rotation[u] or rotation[u][v]==w,'Inconsistent source rotation'
            rotation[u][v]=w
    for j,u in enumerate(outer):rotation[u][outer[(j+1)%len(outer)]]=outer[j-1]
    cyclic={}
    for u,successor in rotation.items():
        v=min(successor);neighbors=[]
        while v not in neighbors:
            neighbors.append(v);v=successor[v]
        assert v==neighbors[0] and len(neighbors)==len(successor),'Disconnected source vertex link'
        cyclic[u]=neighbors
    embedding=nx.PlanarEmbedding();embedding.set_data(cyclic);embedding.check_structure()
    seen=set();faces=[]
    for u,v in embedding.edges():
        if (u,v) not in seen:faces.append(embedding.traverse_face(u,v,seen))
    exterior=[face for face in faces if len(face)>3]
    def same_cycle(a,b):
        return len(a)==len(b) and any(a==b[k:]+b[:k] for k in range(len(b)))
    assert len(exterior)==1 and (same_cycle(exterior[0],outer) or same_cycle(exterior[0],outer[::-1])), 'Specified infinity frame is not the unique outer face'
    grid=nx.combinatorial_embedding_to_pos(embedding,fully_triangulate=False)
    assert set(grid)==set(rotation),'Integer embedding changed source vertices'
    values=np.asarray(list(grid.values()),dtype=float);lo=values.min(axis=0);hi=values.max(axis=0)
    scale=max(hi-lo)/2;assert scale>0
    pos={v:(np.asarray(q,dtype=float)-(lo+hi)/2)/scale for v,q in grid.items()}
    # The frame's given source order bounds the deleted disk, hence is clockwise
    # in the positive source chart. Use the entire cycle, not three possibly
    # collinear outer grid points, to select the single global orientation.
    signed=sum(pos[u][0]*pos[v][1]-pos[u][1]*pos[v][0] for u,v in zip(outer,outer[1:]+outer[:1]))
    assert abs(signed)>1e-13,'Degenerate explicit infinity frame'
    if signed>0:pos={v:q*np.array([1.,-1.]) for v,q in pos.items()}
    return pos


def subdivide(faces,edges,paths,new):
    """Flag subdivision with edge OCCURRENCES, including loops/parallel edges.

Each saved edge path is expanded through its own midpoint. Two rounds suffice
for the source delta-complexes used here; acceptance still checks simplicity.
"""
    outedges={};halves={};mids={}
    def edge(u,v):
        e=len(outedges);outedges[e]=(u,v);return e
    for e,(u,v) in edges.items():
        m=new();mids[e]=m;halves[e]=(edge(u,m),edge(m,v))
    out=[]
    for f in faces:
        center=new();vs=f['vertices'];spokes=[edge(v,center) for v in vs]
        bary=f.get('source_barycentric')
        center_bary=[sum(q[k] for q in bary)/3 for k in range(3)] if bary else None
        for j,u in enumerate(vs):
            v=vs[(j+1)%3];e=f['edges'][j];s=f['edge_signs'][j];m=mids[e]
            a,b=halves[e] if s>0 else halves[e][::-1];radial=edge(m,center)
            mid_bary=[(bary[j][k]+bary[(j+1)%3][k])/2 for k in range(3)] if bary else None
            out.append(dict(vertices=[u,m,center],edges=[a,radial,spokes[j]],edge_signs=[s,1,-1],source=f['source'],
                source_barycentric=[bary[j],mid_bary,center_bary] if bary else None))
            out.append(dict(vertices=[m,v,center],edges=[b,spokes[(j+1)%3],radial],edge_signs=[s,1,-1],source=f['source'],
                source_barycentric=[mid_bary,bary[(j+1)%3],center_bary] if bary else None))
    expanded={key:[pair for e,s in word for pair in
        ([(halves[e][0],1),(halves[e][1],1)] if s>0 else [(halves[e][1],-1),(halves[e][0],-1)])]
        for key,word in paths.items()}
    return out,outedges,expanded


def coarse_view(p):
    assert p['schema']=='marked-display-complex/1'
    doc=p['source_document'];outer=tuple(p['outer_triangle'])
    faces=[dict(t,source_barycentric=[[Fraction(int(i==j)) for i in range(3)] for j in range(3)]) for t in p['triangles']]
    assert any(tuple(t['source'])==outer for t in faces),'Unknown infinity triangle'
    assert len(p['region_marks'])<=1,'Multiple joining-disk fill is not yet supported'
    ids={v for t in faces for v in t['vertices']};nextid=max(ids)+1
    def new():
        nonlocal nextid
        n=nextid;nextid+=1;return n
    # Keep edge occurrences distinct through subdivision; endpoints alone do
    # not identify an edge on a glued surface.
    edge_pairs={}
    for t in faces:
        for i,e in enumerate(t['edges']):
            uv=(t['vertices'][i],t['vertices'][(i+1)%3])
            uv=uv if t['edge_signs'][i]>0 else uv[::-1]
            assert e not in edge_pairs or edge_pairs[e]==uv,'Native edge involution mismatch'
            edge_pairs[e]=uv
    boundaries=copy.deepcopy(p['boundaries']);spoke_ids={}
    def edge(u,v):
        e=max(edge_pairs,default=-1)+1;edge_pairs[e]=(u,v);return e
    for b in boundaries:
        b['cap']=new();b['ports']=[];b['source_vertices']=b['vertices'][:]
        vs=b['vertices'];spokes=[edge(v,b['cap']) for v in vs]
        for j,u in enumerate(vs):
            v=vs[(j+1)%len(vs)];spoke_ids[b['cap'],u]=spokes[j]
            faces.append(dict(vertices=[v,u,b['cap']],edges=[b['edges'][j],spokes[j],spokes[(j+1)%len(vs)]],
                edge_signs=[-b['edge_signs'][j],1,-1],source=None))
    paths={e:[(e,1)] for e in edge_pairs}
    def simple():
        return (all(u!=v for u,v in edge_pairs.values()) and
            len({tuple(sorted(uv)) for uv in edge_pairs.values()})==len(edge_pairs) and
            all(len(set(t['vertices']))==3 for t in faces))
    subdivisions=0
    while not simple():
        assert subdivisions<2,'Flag subdivision did not yield a simplicial cut mesh'
        faces,edge_pairs,paths=subdivide(faces,edge_pairs,paths,new);subdivisions+=1
    def edge_path(e,sign=1):
        word=paths[e] if sign>0 else [(k,-s) for k,s in paths[e][::-1]]
        out=[]
        for k,s in word:
            u,v=edge_pairs[k] if s>0 else edge_pairs[k][::-1]
            assert not out or out[-1]==u
            if not out:out.append(u)
            out.append(v)
        return out
    chosen=next(t for t in faces if t['source'] is not None and tuple(t['source'])==outer)
    triangles=[tuple(t['vertices']) for t in faces if t is not chosen]
    # Source: subdivide preserves each native triangle's source label. Keep
    # this whole-cell correspondence for later common-refinement transport,
    # including unmarked scaffold; None denotes an artificial cut cap.
    triangle_sources=[copy.deepcopy(t['source']) for t in faces if t is not chosen]
    triangle_barycentric=[t.get('source_barycentric') for t in faces if t is not chosen]
    # Twelve sides make this the unique non-triangular outer face. All inserted
    # vertices lie in the SAME specified source triangle, away from its edges.
    frame=[new() for _ in range(12)]
    # A rational inner triangle in the SAME native face supplies the second
    # chart at infinity. Keep exact local coordinates of every flag child and
    # frame vertex, not only its parent-face label.
    cb=chosen['source_barycentric'];infinity_barycentric=[sum(q[k] for q in cb)/3 for k in range(3)]
    inner=[[(q[k]+infinity_barycentric[k])/2 for k in range(3)] for q in cb]
    native_bary={v:q for v,q in zip(chosen['vertices'],cb)}
    for j in range(3):
        for k in range(4):native_bary[frame[4*j+k]]=[Fraction(4-k,4)*inner[j][l]+Fraction(k,4)*inner[(j+1)%3][l] for l in range(3)]
    for j,u in enumerate(chosen['vertices']):
        v=chosen['vertices'][(j+1)%3];arc=[frame[(4*j+k)%12] for k in range(5)]
        triangles.append((u,v,arc[-1]))
        triangles.extend((u,b,a) for a,b in zip(arc,arc[1:]))
        triangle_sources.extend([copy.deepcopy(chosen['source']) for _ in range(5)])
        triangle_barycentric.extend([[native_bary[v] for v in t] for t in triangles[-5:]])
    adj=defaultdict(set);incidence=defaultdict(list)
    for t in triangles:
        assert len(set(t))==3,'Non-simplicial source triangle requires subdivision'
        for u,v in zip(t,t[1:]+t[:1]):
            adj[u].add(v);adj[v].add(u);incidence[tuple(sorted((u,v)))].append((u,v))
    assert all(len(x)==1 or (len(x)==2 and x[0]==x[1][::-1]) for x in incidence.values()),'Invalid oriented mesh'
    assert len(adj)-len(incidence)+len(triangles)==1,'Capped cut mesh is not a disk'
    boundary_edges={e for e,x in incidence.items() if len(x)==1}
    assert boundary_edges=={tuple(sorted((u,v))) for u,v in zip(frame,frame[1:]+frame[:1])}
    # Barycentric coordinates with a fixed convex boundary. This is only a
    # candidate until every face orientation is checked below.
    fixed={v:np.array([np.cos(-2*np.pi*j/len(frame)),np.sin(-2*np.pi*j/len(frame))]) for j,v in enumerate(frame)}
    free=sorted(set(adj)-set(fixed));ix={v:i for i,v in enumerate(free)}
    mat=lil_matrix((len(free),len(free)));rhs=np.zeros((len(free),2))
    for v,i in ix.items():
        mat[i,i]=len(adj[v])
        for w in adj[v]:
            if w in ix:mat[i,ix[w]]=-1
            else:rhs[i]+=fixed[w]
    solution=spsolve(mat.tocsr(),rhs)
    pos={**fixed,**{v:solution[i] for v,i in ix.items()}}
    def area(t):
        a,b,c=(pos[v] for v in t);x=b-a;y=c-a;return x[0]*y[1]-x[1]*y[0]
    harmonic_minimum=float(min(area(t) for t in triangles));candidate_method='harmonic convex-frame coordinates'
    if harmonic_minimum<=1e-13:
        # Same acceptance criterion after the proposal changes: harmonic
        # crowding is a failed candidate, not permission to publish tiny faces.
        print('Harmonic candidate rejected; minimum double area',harmonic_minimum,flush=True)
        pos=rotation_grid_embedding(triangles,frame)
        candidate_method='NetworkX integer grid from retained rotation and explicit infinity frame'
    assert min(area(t) for t in triangles)>1e-13,'Whole-mesh candidate has collapsed or reversed faces'
    curves=copy.deepcopy(p['curves'])
    for c in curves:
        c['color']={'alpha':'#ba433f','beta':'#255dba','aux':'#73538f'}[c['role']]
        for s in c['segments']:
            expanded=edge_path(s['edge'],s['edge_sign'])
            assert [expanded[0],expanded[-1]]==s['vertices']
            s['vertices']=expanded
    patch={x['face']:x for x in doc['patches']}
    crossings=[]
    for x in p['crossings']:
        f,v=x['at'];nodes={p['source_cut_corner_classes'][f'{f},{i},{j}']
            for i,t in enumerate(patch[f]['triangles']) for j,w in enumerate(t) if w==v}
        assert len(nodes)==1,'An interior marked crossing was cut'
        crossings.append(dict(x,vertex=nodes.pop()))
    ports_at={s['vertices'][j] for c in curves for s in c['segments'] for j in (0,-1)}
    def distance(q,a,b):
        delta=b-a;t=np.dot(q-a,delta)/np.dot(delta,delta)
        return np.linalg.norm(q-(a+np.clip(t,0,1)*delta))
    for b in boundaries:
        cap=b['cap'];center=pos[cap]
        clearance=min(distance(center,pos[u],pos[v]) for u,v in incidence if cap not in (u,v))
        b['radius']=.2*clearance
        for v in b['vertices']:
            if v not in ports_at:continue
            spoke=edge_path(spoke_ids[cap,v],-1);neighbor=spoke[1]
            n=new();pos[n]=center+(pos[neighbor]-center)*b['radius']/np.linalg.norm(pos[neighbor]-center)
            b['ports'].append(dict(vertex=n,source_vertex=v,cap_neighbor=neighbor))
            for c in curves:
                for s in c['segments']:
                    if s['vertices'][0]==v:s['vertices'][:0]=[n]+spoke[1:-1]
                    if s['vertices'][-1]==v:s['vertices'].extend(spoke[1:-1][::-1]+[n])
    for b in boundaries:
        # Source: the explicit cap fan added before flag subdivision. Retain
        # ALL radial bridges and boundary edges, not only visible curve ports;
        # they fix the collar correspondence and its marked gluing class.
        b['source_cap_spokes']={str(v):edge_path(spoke_ids[b['cap'],v],-1) for v in b['vertices']}
        b['source_boundary_paths']=[edge_path(e,s) for e,s in zip(b['edges'],b['edge_signs'])]
    native_edge_paths={str(e):edge_path(e) for e in sorted({e for t in p['triangles'] for e in t['edges']})}
    paths={tuple(s['source']):s['vertices'] for c in curves for s in c['segments']}
    disk=[]
    # Source: prepare_marked_display::display_boundary_cycles. Resolve fill
    # intervals through the very same curve paths, including after fairing.
    pieces=[dict(mark=mark['id'],source_triangles=piece['triangles'],
        boundary_cycles=piece['display_boundary_cycles'],contains_infinity=list(outer) in piece['triangles'])
        for mark in p['region_marks'] for piece in mark['cut_pieces']]
    legacy_disk=(len(pieces)==1 and not pieces[0]['contains_infinity'] and
        all(part['kind']=='curve' for cycle in pieces[0]['boundary_cycles'] for part in cycle))
    if legacy_disk:
        for f,_,u,v in p['region_marks'][0]['source_region']['boundary_cycles'][0]:
            vs=paths[(f,u,v)] if (f,u,v) in paths else paths[(f,v,u)][::-1]
            disk.extend(vs[:-1])
    result=dict(schema='marked-cut-view/1',genus=p['genus'],source_document=doc,
        preparation_representation_sha256=p['source_representation_sha256'],
        marking_id=p['marking_id'],cuts=p['cut_curves'],visible_curves=p['visible_curves'],
        grid=p['grid'],positions={str(v):q.tolist() for v,q in pos.items()},
        boundaries=boundaries,boundary_vertex_pairs=p['boundary_vertex_pairs'],
        curves=curves,crossings=crossings,disk=disk,triangles=[list(t) for t in triangles],
        source_triangle_labels=triangle_sources,
        source_triangle_barycentric=[[[str(x) for x in q] for q in t] if t else None for t in triangle_barycentric],
        infinity_chart=dict(source_triangle=list(outer),point_barycentric=[str(x) for x in infinity_barycentric],
            frame_barycentric=[[str(x) for x in native_bary[v]] for v in frame]),
        source_edge_paths=native_edge_paths,
        source_cut_corner_classes=p['source_cut_corner_classes'],chart_outer_boundary=frame,
        infinity_source_triangle=list(outer),
        display_request=p.get('display_request',dict(cut_curves=p['cut_curves'],visible_curves=p['visible_curves'],infinity_triangle=list(outer))),
        disk_region_marks=p['region_marks'],
        initial_embedding_check=dict(method='whole cut cellulation with explicit infinity triangle',
            flag_subdivisions=subdivisions,
            candidate_method=candidate_method,harmonic_candidate_minimum_double_area=harmonic_minimum,
            specified_outer_frame_verified=True,
            oriented_faces=len(triangles),minimum_double_area=min(area(t) for t in triangles),
            selected_exterior_retained=True,source_edge_words_retained=True))
    if 'exterior_face_marker' in p:result['exterior_face_marker']=p['exterior_face_marker']
    if not legacy_disk:result['disk_pieces']=pieces
    for key in ('source_figure','source_status'):
        if key in doc.get('metadata',{}):result[key]=doc['metadata'][key]
    return result


def source_map_layout_candidate(p):
    """Prepare a retained whole-chart candidate using the existing optimizers.

    Sources: draw_marked_cut.spread and Triangle's p/C/e constrained planar
    triangulation switches. No source point or marked route is deleted from
    the payload. Only the numerical scaffold is replaced; the two complete
    charts must subsequently pass the exact same-domain verifier before this
    candidate can be accepted as a source map.
    """
    import triangle
    from draw_marked_cut import spread
    coarse=coarse_view(p)
    positions={int(v):np.array(q) for v,q in coarse['positions'].items()}
    original=[list(t) for t in coarse['triangles']]
    caps={b['cap'] for b in coarse['boundaries']}
    moved=spread(positions,original,set(coarse['chart_outer_boundary']),caps,preserve_linear_orientation=True)
    assert all(np.array_equal(positions[v],moved[v]) for v in coarse['chart_outer_boundary']), 'Specified infinity frame moved'
    coarse['positions']={str(v):q.tolist() for v,q in moved.items()}
    edges=lambda t:{tuple(sorted((int(u),int(v)))) for u,v in zip(t,t[1:]+t[:1])}
    oldedges={edge for t in original for edge in edges(t)}
    keep=edges(coarse['chart_outer_boundary'])
    for t in original:
        if caps.intersection(t):keep.update(edges(t))
    portcaps={port['vertex']:b['cap'] for b in coarse['boundaries'] for port in b['ports']}
    for curve in coarse['curves']:
        for segment in curve['segments']:
            word=[portcaps.get(v,v) for v in segment['vertices']]
            keep.update(tuple(sorted((u,v))) for u,v in zip(word,word[1:]))
    if coarse['disk']:keep.update(edges(coarse['disk']))
    assert keep<=oldedges, 'A retained path is not an original source edge'
    ids=sorted({v for edge in keep for v in edge});index={v:i for i,v in enumerate(ids)}
    vertices=np.array([moved[v] for v in ids]);constraints=np.array([[index[u],index[v]] for u,v in sorted(keep)])
    candidate=triangle.triangulate(dict(vertices=vertices,segments=constraints),'pCe')
    assert np.array_equal(candidate['vertices'],vertices), 'Source reduction inserted or moved vertices'
    reduced=[[ids[int(i)] for i in t] for t in candidate['triangles']]
    reducededges={edge for t in reduced for edge in edges(t)}
    assert keep<=reducededges, 'Source reduction lost a constrained edge'
    for cap in caps:
        assert {tuple(sorted(t)) for t in original if cap in t}=={tuple(sorted(t)) for t in reduced if cap in t}, 'Source cap star changed'
    def area(t):
        a,b,c=(moved[v] for v in t);u=b-a;v=c-a
        return float(u[0]*v[1]-u[1]*v[0])
    assert all(area(t)>0 for t in original+reduced), 'Source reduction is not positively oriented'
    for source,target in (('source_triangle_labels','reduction_source_triangle_labels'),('source_triangle_barycentric','reduction_source_triangle_barycentric')):
        if source in coarse:coarse[target]=coarse.pop(source)
    coarse['reduction_source_triangles']=original;coarse['triangles']=reduced
    coarse['source_reduction_candidate']=dict(method='Triangle pCe on retained source constraints',
        triangle_version=triangle.__version__,constraints=[list(e) for e in sorted(keep)],
        exact_full_domain_validation_required=True)
    from marked_display_source_map import verify_native_chart
    coarse['native_chart_check']=verify_native_chart(coarse,p)
    return coarse


def layout(p, *, _source_map_geometry=False):
    coarse=source_map_layout_candidate(p) if _source_map_geometry else coarse_view(p)
    result=fixed_mesh(coarse,retain_source_map=_source_map_geometry)
    result['geometry_check']=verify(result,coarse)
    if _source_map_geometry:
        # Private handoff to the SAME common plotting pipeline. This context
        # is removed from the saved view after compact reconstruction data is
        # assembled; it is not a second application or a replacement source.
        result['_source_map_input']=dict(prepared=p,coarse=coarse)
    result['svg']=render(result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--source-map-geometry',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();p=json.loads(args.input.read_text(encoding='utf-8'))
    d=layout(p,_source_map_geometry=args.source_map_geometry)
    args.output.write_text(json.dumps(d,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    args.output.with_suffix('.svg').write_text(d['svg'],encoding='utf-8')
    print(json.dumps(d['geometry_check']),flush=True)
