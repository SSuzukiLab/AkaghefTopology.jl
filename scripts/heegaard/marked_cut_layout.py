"""Realize an oriented marked cut mesh on the standard two-by-genus circle grid.

Source curve-edge references, crossing rotations and joining disk are retained.
Only verified geometry is exported. This numerical solver is not a completeness
claim for arbitrary triangulations or general genus.
"""
import copy
import math
from collections import defaultdict

import numpy as np


def cross(a,b): return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def source_mesh(data):
    """Reconstruct the exact retained source mesh without any target solve.

    Source: fixed_mesh's affine parent construction. Replay uses the SAME
    producer, rather than trusting saved source coordinates or invoking a
    second harmonic/untangling optimization merely to recover its domain.
    """
    return fixed_mesh(data,retain_source_map=True,_source_mesh_only=True)


def fixed_mesh(data, *, retain_source_map=False, _source_mesh_only=False):
    """Realize the same cut triangulation with prescribed boundary collars.

    A positive oriented realization gives the explicit facewise PL map; no arc
    is independently selected from endpoint data. Invalid iterates are never
    published. The regularized area objective is only a numerical search tool.
    """
    from scipy.optimize import minimize
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import spsolve
    if _source_mesh_only and not retain_source_map:
        raise ValueError('Exact source reconstruction requires retained parent coordinates')
    d=copy.deepcopy(data);old={int(v):np.array(p) for v,p in data['positions'].items()}
    capids={b['cap'] for b in d['boundaries']};tris=[tuple(t) for t in d['triangles']]
    # Source: coarse_view's exact native-face/edge correspondence. Bind these
    # arrays to their OWN source triangulation and coordinates before replacing
    # the display mesh; they must not appear to index the refined final faces.
    if 'source_triangle_labels' in data or 'reduction_source_triangle_labels' in data:
        d['source_chart']=dict(positions=copy.deepcopy(data['positions']),
            triangles=copy.deepcopy(data.get('reduction_source_triangles',data['triangles'])),
            triangle_labels=copy.deepcopy(data.get('reduction_source_triangle_labels',data.get('source_triangle_labels'))),
            triangle_native_barycentric=copy.deepcopy(data.get('reduction_source_triangle_barycentric',data.get('source_triangle_barycentric'))),
            native_edge_paths=copy.deepcopy(data.get('source_edge_paths',{})),
            infinity_chart=copy.deepcopy(data.get('infinity_chart')))
        for key in ('source_triangle_labels','source_triangle_barycentric','source_edge_paths','infinity_chart'):
            d.pop(key,None)
    # Source: layout_marked_document.coarse_view's oriented disk chart. A closed
    # sphere has no auxiliary circles to prescribe; retain that same chart.
    if d['genus']==0:
        assert not d['boundaries'] and not d['cuts']
        d['layout_method']='whole sphere cellulation with explicit infinity triangle'
        d['minimum_double_area']=float(min(cross(old[b]-old[a],old[c]-old[a]) for a,b,c in tris))
        assert d['minimum_double_area']>1e-13
        if retain_source_map:
            from fractions import Fraction
            exact={v:tuple(Fraction(float(x)) for x in q) for v,q in old.items()}
            assert all((exact[b][0]-exact[a][0])*(exact[c][1]-exact[a][1])-(exact[b][1]-exact[a][1])*(exact[c][0]-exact[a][0])>0 for a,b,c in tris)
            outer=list(d['chart_outer_boundary'])
            if sum(exact[u][0]*exact[v][1]-exact[u][1]*exact[v][0] for u,v in zip(outer,outer[1:]+outer[:1]))<0:outer.reverse()
            ids=sorted({v for t in tris for v in t})
            d['source_map_mesh']=dict(vertex_ids=ids,source_positions={str(v):list(map(str,exact[v])) for v in ids},
                triangles=[list(t) for t in tris],cap_positions={},boundary_cycles=[],outer=outer)
        return d['source_map_mesh'] if _source_mesh_only else d
    positions=dict(old);pins={};boundary={};cutnodes={};nextid=max(old)+1
    # Source: the affine parent constructions below. Keep their exact weights
    # for full-cell reconstruction instead of interpreting separately rounded
    # midpoint coordinates as a conforming source triangulation.
    from fractions import Fraction
    exact={v:tuple(Fraction(float(x)) for x in q) for v,q in old.items()} if retain_source_map else None
    def new(p, weights=None, exact_override=None):
        nonlocal nextid
        i=nextid;nextid+=1;positions[i]=np.array(p)
        if retain_source_map:
            assert weights is not None or exact_override is not None, 'Source-map vertex has no retained construction'
            exact[i]=tuple(exact_override) if exact_override is not None else tuple(sum(w*exact[v][j] for v,w in weights) for j in (0,1))
        return i
    adj=defaultdict(set)
    for t in tris:
        for u,v in zip(t,t[1:]+t[:1]):adj[u].add(v);adj[v].add(u)
    genus=d['genus']
    assert d.get('grid')==dict(rows=2,columns=genus,equal_radii=True), 'Use the prepared display complex, not inferred hole-side order'
    for b in d['boundaries']:
        cap=b['cap'];center=old[cap];neighbors=sorted(adj[cap],key=lambda v:math.atan2(*(old[v]-center)[::-1]))
        # Position is a display convention; side/column identities are supplied
        # by the oriented cut quotient, never guessed from list order or labels.
        target=np.array([(b['column']-(genus-1)/2)*.65,(1 if b['row']==0 else -1)*5/12])
        ns=[]
        theta=np.array([math.atan2(*(old[v]-center)[::-1]) for v in neighbors])
        uniform=2*math.pi*np.arange(len(neighbors))/len(neighbors);uniform+=np.mean(theta-uniform)
        ring_radius=.12
        # Dense ring, with original cap spokes retained as labelled ports.
        for j,v in enumerate(neighbors):
            theta0=uniform[j];theta1=uniform[(j+1)%len(neighbors)]
            if j==len(neighbors)-1:theta1+=2*math.pi
            node=new(center+.05*(old[v]-center), [(cap,Fraction(19,20)),(v,Fraction(1,20))]);cutnodes[(cap,v)]=node;ns.append(node)
            pins[node]=target+ring_radius*np.array([math.cos(theta0),math.sin(theta0)])
            for k in range(1,4):
                a=theta0+(theta1-theta0)*k/4
                # Source: the cap-corner clipping below. Interpolate the actual
                # adjacent spoke cuts, so every inserted node stays in its
                # source wedge, regardless of unequal angles and edge lengths.
                left=center+.05*(old[v]-center)
                right=center+.05*(old[neighbors[(j+1)%len(neighbors)]]-center)
                n=new((1-k/4)*left+(k/4)*right, [(cap,Fraction(19,20)),(v,Fraction(4-k,80)),(neighbors[(j+1)%len(neighbors)],Fraction(k,80))]);ns.append(n)
                pins[n]=target+ring_radius*np.array([math.cos(a),math.sin(a)])
        boundary[cap]=ns;b['grid_center']=target.tolist();b['mesh_boundary']=ns
    # Clip the cap corners in each source triangle, retaining the boundary arcs.
    mesh=[]
    for t in tris:
        polygon=[]
        for j,v in enumerate(t):
            if v not in capids:polygon.append(v);continue
            prev,nxt=t[j-1],t[(j+1)%3];a=cutnodes[(v,prev)];b=cutnodes[(v,nxt)];ring=boundary[v]
            k=ring.index(a);sequence=[a]
            # Triangle interior is on the left; holes have clockwise boundary.
            while sequence[-1]!=b:k=(k-1)%len(ring);sequence.append(ring[k])
            polygon.extend(sequence)
        center=new(np.mean([positions[v] for v in polygon],axis=0), [(v,Fraction(1,len(polygon))) for v in polygon])
        mesh.extend((center,u,v) for u,v in zip(polygon,polygon[1:]+polygon[:1]))
    # Add an exterior collar so no mathematical curve vertex is pinned at infinity.
    outer=list(d.get('chart_outer_boundary',d.get('chart_outer_triangle',[])))
    assert len(outer)>=3,'An explicit retained exterior chart boundary is required'
    # Source: layout_marked_document.rotation_grid_embedding. Its outer grid
    # can contain collinear consecutive vertices; the full cycle determines
    # orientation even when the first three points do not.
    outer_area=sum(cross(old[u],old[v]) for u,v in zip(outer,outer[1:]+outer[:1]))
    assert abs(outer_area)>1e-13,'Degenerate exterior chart boundary'
    if outer_area<0:outer.reverse()
    # Source: scipy.optimize.linprog (HiGHS), a Chebyshev center of the
    # polygon kernel. Homothety about this verified kernel point makes a
    # positive source collar even for collinear outer vertices. Source
    # coordinates must not be replaced with independently phased final pins.
    from scipy.optimize import linprog
    oq=np.array([old[v] for v in outer]);edges=np.roll(oq,-1,axis=0)-oq
    lengths=np.linalg.norm(edges,axis=1)
    normals=np.column_stack((-edges[:,1],edges[:,0]))/lengths[:,None]
    kernel=linprog([0.,0.,-1.],A_ub=np.column_stack((-normals,np.ones(len(outer)))),
        b_ub=-np.sum(normals*oq,axis=1),bounds=[(None,None),(None,None),(0,None)],method='highs')
    assert kernel.success and kernel.x[2]>1e-13,'Source exterior has no verified positive polygon kernel'
    origin=kernel.x[:2];ext=[]
    for j,v in enumerate(outer):
        a=2*math.pi*j/len(outer)-math.pi/2;r=max(2.8,.65*genus+1)
        exact_point=[2*exact[v][j]-Fraction(float(origin[j])) for j in (0,1)] if retain_source_map else None
        n=new(origin+2*(old[v]-origin), exact_override=exact_point);pins[n]=np.array([r*math.cos(a),r*math.sin(a)]);ext.append(n)
    for j,u in enumerate(outer):
        v=outer[(j+1)%len(outer)];a,b=ext[j],ext[(j+1)%len(outer)]
        mesh.extend(((a,b,v),(a,v,u)))
    # Interior diagonals are routing scaffolding, not straight-line constraints
    # between the prescribed circles. Subdivide them so they can bend.
    edge_count=defaultdict(int)
    for t in mesh:
        for u,v in zip(t,t[1:]+t[:1]):edge_count[tuple(sorted((u,v)))]+=1
    mids={e:new((positions[e[0]]+positions[e[1]])/2, [(e[0],Fraction(1,2)),(e[1],Fraction(1,2))]) for e,count in edge_count.items() if count==2}
    refined=[]
    for t in mesh:
        c=new(np.mean([positions[v] for v in t],axis=0), [(v,Fraction(1,len(t))) for v in t])
        for u,v in zip(t,t[1:]+t[:1]):
            mid=mids.get(tuple(sorted((u,v))))
            refined.extend([(c,u,v)] if mid is None else [(c,u,mid),(c,mid,v)])
    mesh=refined
    ids=sorted({v for t in mesh for v in t});ix={v:i for i,v in enumerate(ids)}
    incidence=defaultdict(list)
    for t in mesh:
        for a,b in zip(t,t[1:]+t[:1]):incidence[tuple(sorted((a,b)))].append((a,b))
    assert all(len(occ)==1 or (len(occ)==2 and occ[0]==occ[1][::-1]) for occ in incidence.values()),'Inconsistent mesh orientation'
    assert len(ids)-len(incidence)+len(mesh)==1-2*genus,('Bad Euler',len(ids)-len(incidence)+len(mesh))
    T=np.array([[ix[v] for v in t] for t in mesh]);P=np.array([positions[v] for v in ids]);free=[v for v in ids if v not in pins];fi={v:i for i,v in enumerate(free)}
    # A later harmonic solve is only a proposal. Check the retained source
    # realization first; invalid collars must never masquerade as a warm start.
    source_areas=cross(P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]])
    assert np.all(source_areas>0),('Invalid source collar realization',float(min(source_areas)),int(sum(source_areas<=0)))
    d['source_mesh_check']=dict(vertices=len(ids),faces=len(T),minimum_double_area=float(min(source_areas)),positive=True)
    if retain_source_map:
        def exact_determinant(t):
            a,b,c=(exact[v] for v in t)
            return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
        assert all(exact_determinant(t)>0 for t in mesh), 'Exact parent source mesh is not positive'
        d['source_map_mesh']=dict(vertex_ids=ids,source_positions={str(v):list(map(str,exact[v])) for v in ids},
            triangles=[list(t) for t in mesh],cap_positions={str(b['cap']):list(map(str,exact[b['cap']])) for b in data['boundaries']},
            boundary_cycles=[b['mesh_boundary'] for b in d['boundaries']],outer=ext)
    print('Source collar mesh',d['source_mesh_check'],flush=True)
    if _source_mesh_only:return d['source_map_mesh']
    adjacency=defaultdict(set)
    for t in mesh:
        for u,v in zip(t,t[1:]+t[:1]):adjacency[u].add(v);adjacency[v].add(u)
    A=lil_matrix((len(free),len(free)));rhs=np.zeros((len(free),2))
    for v,i in fi.items():
        A[i,i]=len(adjacency[v])
        for w in adjacency[v]:
            if w in fi:A[i,fi[w]]=-1
            else:rhs[i]+=pins[w]
    solution=spsolve(A.tocsr(),rhs)
    for v in ids:P[ix[v]]=pins[v] if v in pins else solution[fi[v]]
    freeidx=np.array([ix[v] for v in free]);target_area=17./len(T)
    def objective(x,eps):
        Q=P.copy();Q[freeidx]=x.reshape((-1,2));a,b,c=(Q[T[:,j]] for j in range(3));D=cross(b-a,c-a)
        S=np.sum((a-b)**2+(b-c)**2+(c-a)**2,axis=1)
        root=np.sqrt(D*D+eps*eps);chi=np.where(D>=0,(D+root)/2,eps*eps/(2*(root-D)))
        derivative=.5*(1+D/root);beta=1/target_area
        num=S+beta*(D*D+target_area*target_area);value=np.sum(num/chi)
        dd=2*beta*D/chi-num*derivative/(chi*chi);grad=np.zeros_like(Q)
        for j,u,v,w in ((0,a,b,c),(1,b,c,a),(2,c,a,b)):
            ds=2*(2*u-v-w)/chi[:,None]
            da=np.column_stack((v[:,1]-w[:,1],w[:,0]-v[:,0]))*dd[:,None]
            np.add.at(grad,T[:,j],ds+da)
        return float(value),grad[freeidx].ravel()
    x=P[freeidx].ravel()
    for eps in (.03,.01,.003,.001,.0003,.0001,.00003):
        opt=minimize(lambda x:objective(x,eps),x,jac=True,method='L-BFGS-B',options={'maxiter':4000,'ftol':1e-11,'maxls':50,'maxcor':25})
        x=opt.x;P[freeidx]=x.reshape((-1,2));areas=cross(P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]])
        print('Fixed-boundary mesh',eps,'negative faces',sum(areas<=0),'min',min(areas),'iterations',opt.nit,opt.message,flush=True)
        if min(areas)>1e-7:break
    assert min(areas)>1e-8,('Could not realize fixed boundary mesh',min(areas))
    pp={v:P[ix[v]] for v in ids}
    for b in d['boundaries']:
        cap=b['cap'];target=np.array(b.pop('grid_center'));pp[cap]=target;b['radius']=55/600
        for port in b['ports']:
            spoke=cutnodes[(cap,port.get('cap_neighbor',port['source_vertex']))];q=pp[spoke];newport=max(pp)+1
            pp[newport]=target+(q-target)/np.linalg.norm(q-target)*b['radius']
            previous=port['vertex'];port['vertex']=newport
            for curve in d['curves']:
                for s in curve['segments']:
                    if s['vertices'][0]==previous:s['vertices'][:1]=[newport,spoke]
                    if s['vertices'][-1]==previous:s['vertices'][-1:]=[spoke,newport]
    def refine_path(vs,closed=False):
        out=[]
        for u,v in zip(vs,vs[1:]+(vs[:1] if closed else [])):
            out.append(u);mid=mids.get(tuple(sorted((u,v))))
            if mid is not None:out.append(mid)
        if not closed:out.append(vs[-1])
        return out
    for curve in d['curves']:
        for s in curve['segments']:s['vertices']=refine_path(s['vertices'])
    d['disk']=refine_path(d['disk'],closed=True)
    d['positions']={str(v):p.tolist() for v,p in pp.items()};d['triangles']=[list(t) for t in mesh]
    d['layout_method']='fixed 2-by-g boundary cellulation';d['minimum_double_area']=float(min(areas))
    d['source_triangle_refinement']=[list(t) for t in tris]
    return d


if __name__=='__main__':
    import argparse,hashlib,json
    from pathlib import Path
    from draw_marked_cut import render
    from verify_fixed_grid import verify
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=json.loads(args.input.read_text(encoding='utf-8'))
    assert hashlib.sha256(Path(source['source']).read_bytes()).hexdigest()==source['source_sha256']
    result=fixed_mesh(source)
    result['geometry_check']=verify(result,source)
    output=args.output or args.input.with_name('KodaFixedMesh.json')
    output.write_text(json.dumps(result,ensure_ascii=False),encoding='utf-8')
    output.with_suffix('.svg').write_text(render(result),encoding='utf-8')
    print(json.dumps(result['geometry_check']),flush=True)
