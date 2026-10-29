"""Fresh exact graph-isotopy gate using existing attachment predicates.

Source: marked_display_source_map.verify_bezier_attachment: exact dyadic
control hulls, positive directed projection, shared-joint cones and strict
nonshared separation. Every candidate pair, including nonneighbors in the same
cubic, is checked at synchronized time. No sampled/floating pass is accepted.
"""
import math
from fractions import Fraction as Q
from periodic_point_transport import need,hull,separation,sub,dot,split_cubic,serial

MAX_DEPTH=24

def blend_rows(a,b):
    return [tuple((x+y)/2 for x,y in zip(p,q)) for p,q in zip(a,b)]

def time_split(pair):
    middle=blend_rows(pair['start'],pair['end'])
    a=dict(pair,end=middle);b=dict(pair,start=middle)
    return a,b

def param_split(pair):
    a,b=split_cubic(pair['start'],Q(1,2));c,d=split_cubic(pair['end'],Q(1,2))
    return dict(pair,start=a,end=c),dict(pair,start=b,end=d)

def projection_witness(pair):
    a=pair['start'];b=pair['end'];x=sub(a[-1],a[0]);y=sub(b[-1],b[0])
    differences=[sub(v,u) for row in (a,b) for u,v in zip(row,row[1:])]
    for n in (x,y,tuple(u+v for u,v in zip(x,y))):
        values=[dot(n,v) for v in differences]
        if min(values)>0:return dict(normal=serial(n),margin=str(min(values)))
    return None

def regular_certificate(pair,depth=0):
    w=projection_witness(pair)
    if w is not None:return dict(separator=w)
    need(depth<MAX_DEPTH,'graph leaf derivative projection cannot be certified')
    a,b=time_split(pair)
    return dict(time_children=[regular_certificate(a,depth+1),regular_certificate(b,depth+1)])

def relative_certificate(a,b,depth=0):
    # ALL independent spatial parameter control differences at SAME time.
    for k in range(2):
        lower=min(min(p[k] for p in a[row])-max(p[k] for p in b[row]) for row in ('start','end'))
        upper=max(max(p[k] for p in a[row])-min(p[k] for p in b[row]) for row in ('start','end'))
        if lower>0 or upper<0:
            n=[0,0];n[k]=1 if lower>0 else -1
            return dict(separator=dict(normal=list(map(str,n)),gap=str(lower if lower>0 else -upper)))
    points=[sub(x,y) for row in ('start','end') for x in a[row] for y in b[row]]
    w=separation(hull(points),[(Q(0),Q(0))])
    if w is not None:return dict(separator=w)
    need(depth<MAX_DEPTH,'nonincident graph family contact or refinement bound')
    direction=depth%3
    if direction==0:
        a0,a1=time_split(a);b0,b1=time_split(b)
        return dict(time_children=[relative_certificate(a0,b0,depth+1),relative_certificate(a1,b1,depth+1)])
    if direction==1:
        a0,a1=param_split(a)
        return dict(first_parameter_children=[relative_certificate(a0,b,depth+1),relative_certificate(a1,b,depth+1)])
    b0,b1=param_split(b)
    return dict(second_parameter_children=[relative_certificate(a,b0,depth+1),relative_certificate(a,b1,depth+1)])

def joint_witness(a,b,key):
    i=0 if a['left']==key else -1;j=0 if b['left']==key else -1
    need(all(a[row][i]==b[row][j] for row in ('start','end')),'source joint endpoints disagree')
    X=[sub(p,a[row][i]) for row in ('start','end') for k,p in enumerate(a[row]) if k!=i%4]
    Y=[sub(p,b[row][j]) for row in ('start','end') for k,p in enumerate(b[row]) if k!=j%4]
    # Source: existing verify_bezier_attachment normalized-ray proposal;
    # every resulting halfspace sign remains exact rational below.
    import numpy as np
    x=np.array([float(v) for v in X[-1]]);y=np.array([float(v) for v in Y[-1]])
    guesses=[tuple(Q(float(v)) for v in x/np.linalg.norm(x)-y/np.linalg.norm(y))]
    guesses.append(tuple(x-y for x,y in zip(X[-1],Y[-1])))
    guesses.extend(((-v[1],v[0]) for v in X+Y))
    guesses.extend(X+Y)
    for n in guesses:
        x=[dot(n,v) for v in X];y=[dot(n,v) for v in Y]
        if min(x)>0 and max(y)<0 or max(x)<0 and min(y)>0:
            return dict(normal=serial(n),first_min=str(min(x)),first_max=str(max(x)),second_min=str(min(y)),second_max=str(max(y)))
    return None

def joint_certificate(a,b,key,depth=0):
    w=joint_witness(a,b,key)
    if w is not None:return dict(separator=w)
    need(depth<MAX_DEPTH,'incident graph leaves lack exact moving-joint cone separation')
    a0,a1=time_split(a);b0,b1=time_split(b)
    return dict(time_children=[joint_certificate(a0,b0,key,depth+1),joint_certificate(a1,b1,key,depth+1)])

def graph_hole_certificate(pair,polygon,depth=0):
    w=separation(hull(pair['start']+pair['end']),polygon)
    if w is not None:return dict(separator=w)
    need(depth<MAX_DEPTH,'graph family meets hole or hole-hull refinement bound')
    if depth%2:
        a,b=time_split(pair)
    else:a,b=param_split(pair)
    return dict(children=[graph_hole_certificate(a,polygon,depth+1),graph_hole_certificate(b,polygon,depth+1)])

def spatial_leaves(pair,depth=0):
    # Source: existing .003 SVG flatness proposal. Exact regularity/separation
    # below decide validity; this bound only controls finite hull granularity.
    def small(q):
        return max(math.hypot(*(float(p[k]-((3-i)*q[0][k]+i*q[-1][k])/3) for k in range(2))) for i,p in enumerate(q[1:3],1))<=.003
    if small(pair['start']) and small(pair['end']):return [pair]
    need(depth<25,'graph spatial subdivision bound')
    a,b=param_split(pair);mid=(pair['left'],'dyadic',depth,tuple(a['start'][-1]),tuple(a['end'][-1]))
    a['right']=mid;b['left']=mid
    return spatial_leaves(a,depth+1)+spatial_leaves(b,depth+1)

def normalized_outer_planes(outer):
    """Exact positive-scalar deduplication of the SAME oriented inequalities.

    Source: orient(u,v,q)>0; no hull, tolerance or reversed-normal merging.
    Degenerate zero edges retain the original strict predicate rejection.
    """
    unique={};mapping=[]
    for i,(u,v) in enumerate(zip(outer,outer[1:]+outer[:1])):
        dx,dy=sub(v,u);raw=(-dy,dx,dy*u[0]-dx*u[1])
        need(raw[0]!=0 or raw[1]!=0,'degenerate outer zero edge')
        scale=abs(raw[0] if raw[0]!=0 else raw[1]);plane=tuple(x/scale for x in raw)
        need(scale>0 and tuple(scale*x for x in plane)==raw,'outer positive-multiple identity failed')
        if plane not in unique:unique[plane]=len(unique)
        mapping.append(dict(original_edge=int(i),unique_plane=int(unique[plane]),positive_scale=str(scale)))
    return list(unique),dict(logical_edges=len(outer),unique_planes=len(unique),mapping=mapping)


def strict_outer_margin(planes,q):
    return min(nx*q[0]+ny*q[1]+c for nx,ny,c in planes)


def certify(family,emit=None):
    import math
    from shapely.geometry import box
    from shapely.strtree import STRtree
    leaves=[leaf for pair in family['pairs'] for leaf in spatial_leaves(pair)]
    regular=[];counts=dict(regular=0,pairs=0,incident=0,holes=0)
    def retain(kind,key,value):
        if emit is not None:emit(kind,key,value)
        return value
    for i,p in enumerate(leaves):
        proof=retain("regular",i,regular_certificate(p));counts["regular"]+=1
        if emit is None:regular.append(proof)
    # Complete directed source-labelled join incidence must be four cycles.
    incidences={}
    for i,p in enumerate(leaves):
        for side in ('left','right'):incidences.setdefault(p[side],[]).append((i,side))
    need(all(len(x)==2 and {side for _,side in x}=={'left','right'} for x in incidences.values()),'graph join incidence is not directed cycles')
    boxes=[]
    for p in leaves:
        points=p['start']+p['end'];lo=[min(q[k] for q in points) for k in range(2)];hi=[max(q[k] for q in points) for k in range(2)]
        boxes.append(box(*(math.nextafter(float(x),-math.inf) for x in lo),*(math.nextafter(float(x),math.inf) for x in hi)))
    tree=STRtree(boxes);pairs=[]
    for i,p in enumerate(leaves):
        for j in tree.query(boxes[i]):
            j=int(j)
            if j<=i:continue
            q=leaves[j];common={p['left'],p['right']}&{q['left'],q['right']}
            if common:
                need(len(common)==1,'unexpected doubled graph adjacency')
                proof=joint_certificate(p,q,next(iter(common)))
            else:proof=relative_certificate(p,q)
            record=dict(first=i,second=j,shared=bool(common),proof=proof)
            retain("pair",(i,j),record);counts["pairs"]+=1;counts["incident"]+=bool(common)
            if emit is None:pairs.append(record)
    holes=[]
    for i,p in enumerate(leaves):
        for cut,side,shape in family['holes']:
            record=dict(leaf=i,cut=cut,side=side,proof=graph_hole_certificate(p,shape))
            retain('hole',(i,cut,side),record);counts['holes']+=1
            if emit is None:holes.append(record)
    outer=family['outer'];planes,outer_coverage=normalized_outer_planes(outer)
    for p in leaves:
        for row in ('start','end'):
            for q in p[row]:
                need(strict_outer_margin(planes,q)>0, 'graph controls leave fixed convex exterior')
    return dict(schema='periodic-exact-graph-family/1',source_document_sha256=family['source_document_sha256'],
                leaves=len(leaves),counts=counts,outer_halfspaces=outer_coverage,regular=regular,pairs=pairs,hole_certificates=holes,
                scope='fresh exact old-R-PL to actual cubic embedded-family certificate, no serialized pass reused')
