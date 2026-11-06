"""Five coupled point-track proposal and independent exact relative checker.

Private ordinary-call contract: fresh common source-map/attachment gates and
the Julia inverse compiler precede this module. Numerical optimization proposes
only; no candidate becomes an image until the complete exact checker passes.
Sources: existing relative convex-hull predicates, de Casteljau subdivision,
exact quadratic minimum and SciPy/Shapely proposal libraries already in use.
"""
from fractions import Fraction as Q
import copy,hashlib,json
from periodic_point_transport import need,rat,serial,hull,separation,sub,dot,split_cubic,minimum
from periodic_graph_certificate import param_split,time_split,normalized_outer_planes,strict_outer_margin

MAX_DEPTH=24

def inverse_points(compiled,family,recipe):
    """Check actual barycentric coordinates, not a serialized pass flag."""
    T=family['triangles'];A=family['A'];R=family['R'];vertices={int(v) for t in T for v in t}
    def meshhash(P):
        payload=dict(triangles=T,positions={str(v):serial(P[v]) for v in vertices})
        return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
    m=compiled['inverse_anchor_mesh']
    need(m['A_sha256']==meshhash(A) and m['R_sha256']==meshhash(R) and m['triangles']==len(T),'inverse stage mesh binding differs')
    def exact_strings(values,length,what):
        import re
        need(type(values) is list and len(values)==length and all(type(v) is str and re.fullmatch(r'-?[0-9]+/[1-9][0-9]*',v) for v in values),'invalid '+what+' rational grammar')
        parsed=list(map(rat,values))
        need(serial(parsed)==values,'noncanonical '+what+' rational spelling')
        return parsed
    starts=[]
    for p in compiled['points']:
        inv=p['inverse_anchor_point'];index=inv['triangle_index']
        need(type(index) is int and 1<=index<=len(T),'invalid inverse triangle index type/range')
        vertices_record=inv['triangle_vertices']
        need(type(vertices_record) is list and len(vertices_record)==3 and all(type(v) is int for v in vertices_record),'invalid inverse triangle vertex types/shape')
        triangle=T[index-1]
        need(list(triangle)==vertices_record,'inverse triangle index mismatch')
        weights=exact_strings(inv['lambda'],3,'inverse lambda');need(len(weights)==3 and sum(weights)==1 and min(weights)>=0,'invalid inverse barycentric weights')
        qA=tuple(exact_strings(p['point_svg_exact'],2,'source qA'));qR=tuple(exact_strings(inv['qR'],2,'inverse qR'))
        need(qA==tuple(exact_strings(inv['qA'],2,'inverse qA')),'inverse initial point changed')
        need(qA==tuple(sum(weights[j]*A[int(triangle[j])][k] for j in range(3)) for k in range(2)),'inverse A barycentric identity failed')
        need(qR==tuple(sum(weights[j]*R[int(triangle[j])][k] for j in range(3)) for k in range(2)),'inverse R barycentric identity failed')
        starts.append(qR)
    need(len(starts)==5,'requires all five authored points')
    return starts

def inverse_anchor_tracks(compiled,family,recipe):
    """Reverse the already freshly certified R-to-A mesh in the SAME API call.

    Fixed lambda gives q(t)=sum lambda_j*((1-t)A_j+tR_j) on its
    actual moving face. Global mesh/graph compatibility is a private caller
    prerequisite, never an input Boolean. Visible graph is fixed here, checked
    by compile_family; source qA is retained separately from qR.
    """
    ends=inverse_points(compiled,family,recipe)
    starts=[tuple(map(rat,p['point_svg_exact'])) for p in compiled['points']]
    pairs=[]
    for i in range(5):
        for j in range(i):
            value=point_pair_minimum(starts[i],ends[i],starts[j],ends[j])
            need(value>0,'inverse anchor point tracks collide')
            pairs.append(dict(first=i,second=j,minimum_squared_distance=str(value)))
    return ends,dict(stage='reverse actual certified A-to-R mesh',
        points=[dict(qA=serial(a),qR=serial(b),barycentric=p['inverse_anchor_point'])
            for a,b,p in zip(starts,ends,compiled['points'])],point_pairs=pairs)


def point_pair_minimum(a0,a1,b0,b1):
    d=sub(a0,b0);v=sub(sub(a1,a0),sub(b1,b0))
    return minimum(dot(d,d),2*dot(d,v),dot(v,v))

def relative_point_certificate(pair,q0,q1,depth=0):
    points=[sub(p,q) for row,q in ((pair['start'],q0),(pair['end'],q1)) for p in row]
    w=separation(hull(points),[(Q(0),Q(0))])
    if w is not None:return dict(separator=w)
    need(depth<MAX_DEPTH,'moving point meets graph or exact refinement bound')
    if depth%2==0:
        a,b=time_split(pair);mid=tuple((x+y)/2 for x,y in zip(q0,q1))
        return dict(time_children=[relative_point_certificate(a,q0,mid,depth+1),relative_point_certificate(b,mid,q1,depth+1)])
    a,b=param_split(pair)
    return dict(parameter_children=[relative_point_certificate(a,q0,q1,depth+1),relative_point_certificate(b,q0,q1,depth+1)])

def fixed_polygon_certificate(shape,q0,q1):
    w=separation(shape,list(dict.fromkeys((q0,q1))))
    need(w is not None,'point track meets fixed hole/support hull')
    return w

def exact_track_check(family,tracks,emit=None):
    """Full synchronized time/source rectangle and all ten point pairs."""
    need(len(tracks)==5,'all five tracks required')
    knots=tracks[0]['times'];need(knots[0]==0 and knots[-1]==1 and all(a<b for a,b in zip(knots,knots[1:])),'track time coverage incomplete')
    need(all(t['times']==knots and len(t['points'])==len(knots) for t in tracks),'tracks are not synchronized')
    planes,outer_coverage=normalized_outer_planes(family['outer'])
    proof=[];counts=dict(slabs=0,graph_rectangles=0,hole_segments=0,point_pairs=0)
    def retain(kind,key,value):
        if emit is not None:emit(kind,key,value)
        return value if emit is None else None
    for slab,(lo,hi) in enumerate(zip(knots,knots[1:])):
        # Restrict the ACTUAL affine-time graph family to the common slab.
        pairs=[]
        for p in family['pairs']:
            start=[tuple((1-lo)*x+lo*y for x,y in zip(a,b)) for a,b in zip(p['start'],p['end'])]
            end=[tuple((1-hi)*x+hi*y for x,y in zip(a,b)) for a,b in zip(p['start'],p['end'])]
            pairs.append(dict(p,start=start,end=end))
        pointproof=[]
        for i,track in enumerate(tracks):
            q0,q1=track['points'][slab:slab+2]
            graph=[]
            for cell,p in enumerate(pairs):
                value=retain('graph',(slab,i,cell),relative_point_certificate(p,q0,q1));counts['graph_rectangles']+=1
                if emit is None:graph.append(value)
            holes=[]
            for hole,(_,_,shape) in enumerate(family['holes']):
                value=retain('hole',(slab,i,hole),fixed_polygon_certificate(shape,q0,q1));counts['hole_segments']+=1
                if emit is None:holes.append(value)
            outer=family['outer']
            margins=[strict_outer_margin(planes,q) for q in (q0,q1)]
            need(min(margins)>0,'moving point leaves fixed convex exterior')
            pointproof.append(dict(point=i,graph=graph,holes=holes,outer_minimum=str(min(margins))))
        distinct=[]
        for i in range(5):
            for j in range(i):
                value=point_pair_minimum(tracks[i]['points'][slab],tracks[i]['points'][slab+1],tracks[j]['points'][slab],tracks[j]['points'][slab+1])
                need(value>0,'point tracks collide at synchronized time')
                record=dict(first=i,second=j,minimum_squared_distance=str(value));retain('point_pair',(slab,i,j),record);counts['point_pairs']+=1
                if emit is None:distinct.append(record)
        counts['slabs']+=1
        if emit is None:proof.append(dict(interval=[str(lo),str(hi)],points=pointproof,pairs=distinct))
        else:retain('slab',slab,dict(interval=[str(lo),str(hi)],outer_minima=[p['outer_minimum'] for p in pointproof]))
    return proof if emit is None else counts

def propose_tracks(family,starts,steps):
    """SciPy constructs endpoints jointly; no manually chosen endpoint target.

    Source-labelled moving complementary regions are numerical barriers ONLY.
    The independent relative checker, not region samples, accepts the result.
    """
    import numpy as np
    from scipy.optimize import minimize
    from shapely.geometry import Polygon,Point
    from shapely.ops import nearest_points
    from periodic_graph_certificate import spatial_leaves
    groups={}
    for pair in family['pairs']:groups.setdefault(pair['left'][0],[]).extend(spatial_leaves(pair))
    rows=list(groups.values());outer=Polygon([(float(x),float(y)) for x,y in family['outer']])
    base=outer
    for _,_,shape in family['holes']:base=base.difference(Polygon([(float(x),float(y)) for x,y in shape]))
    def regions(t):
        polygons=[Polygon([tuple(float((1-t)*a+t*b) for a,b in zip(p['start'][0],p['end'][0])) for p in group]) for group in rows]
        result=[]
        for signature in signatures:
            region=base
            for inside,polygon in zip(signature,polygons):region=region.intersection(polygon) if inside else region.difference(polygon)
            need(not region.is_empty,'numerical moving region became empty')
            result.append(region)
        return result
    initial=[Polygon([(float(p['start'][0][0]),float(p['start'][0][1])) for p in group]) for group in rows]
    signatures=[[polygon.contains(Point(tuple(map(float,p)))) for polygon in initial] for p in starts]
    paths=[[p] for p in starts];previous=np.array([[float(x),float(y)] for x,y in starts])
    for k in range(1,steps+1):
        t=k/steps;allowed=regions(t);guess=previous.copy()
        margins=[]
        for i,region in enumerate(allowed):
            margin=min(.02,max(region.area,1e-16)**.5/1000)
            inner=region.buffer(-margin)
            if inner.is_empty:inner=region;margin=0.
            if not inner.contains(Point(guess[i])):
                guess[i]=np.asarray(nearest_points(inner,Point(guess[i]))[0].coords[0])
            margins.append(margin)
        def constraints(z):
            points=z.reshape((5,2));out=[]
            for point,region,margin in zip(points,allowed,margins):
                p=Point(point);distance=p.distance(region.boundary)
                out.append((distance if region.contains(p) else -distance)-margin/2)
            for i in range(5):
                for j in range(i):out.append(np.sum((points[i]-points[j])**2)-1e-12)
            return np.asarray(out)
        result=minimize(lambda z:float(np.sum((z.reshape((5,2))-previous)**2)),guess.ravel(),method='SLSQP',
            constraints=[dict(type='ineq',fun=constraints)],options=dict(maxiter=160,ftol=1e-11))
        need(np.min(constraints(result.x))>=-1e-7,'numerical coupled point proposal failed')
        previous=result.x.reshape((5,2))
        for i,p in enumerate(previous):paths[i].append(tuple(Q(float(x)) for x in p))
    times=[Q(k,steps) for k in range(steps+1)]
    return [dict(times=times,points=path) for path in paths]

def rounding_check(family,endpoints):
    qhat=[tuple(Q(float(x)) for x in q) for q in endpoints];proof=[]
    planes,outer_coverage=normalized_outer_planes(family['outer'])
    for i,(q,d) in enumerate(zip(endpoints,qhat)):
        segment=list(dict.fromkeys((q,d)));graph=[]
        for p in family['pairs']:
            fixed=dict(p,start=p['end'],end=p['end'])
            graph.append(relative_point_certificate(fixed,q,d))
        supports=[fixed_polygon_certificate(shape,q,d) for _,_,shape in family['final_supports']]
        holes=[fixed_polygon_certificate(shape,q,d) for _,_,shape in family['holes']]
        outer=family['outer'];minimum_outer=min(strict_outer_margin(planes,p) for p in segment)
        need(minimum_outer>0,'decoded marker leaves exterior')
        proof.append(dict(point=i,graph=graph,supports=supports,holes=holes,outer_minimum=str(minimum_outer)))
    pairs=[]
    for i in range(5):
        for j in range(i):
            value=point_pair_minimum(endpoints[i],qhat[i],endpoints[j],qhat[j]);need(value>0,'decoded markers collide')
            pairs.append(dict(first=i,second=j,minimum_squared_distance=str(value)))
    return qhat,dict(points=proof,pairs=pairs)


def endpoint_wire_payload(original_qA,moved,qhat):
    """Preserve rational definitions; actual renderer uses builtin Float64 JSON.

    Source: existing native JSON writer rejects foreign NumPy scalar types.
    Decoded coordinates are exact binary rationals of the emitted builtin float.
    """
    need(len(original_qA)==len(moved)==len(qhat)==5,'all five endpoint records required')
    output=[]
    for index,(qA,q,d) in enumerate(zip(original_qA,moved,qhat)):
        numeric=[float(x) for x in d]
        need(all(type(x) is float for x in numeric),'foreign renderer numeric scalar')
        need(tuple(Q(x) for x in numeric)==tuple(d),'decoded glyph is not emitted binary rational')
        output.append(dict(point_index=int(index),point_svg_exact=serial(qA),
            point_svg_moved_exact=serial(q),point_svg_decoded=serial(d),point_svg=numeric))
    return output


def transport_all_points(view,recipe,compiled,diagnostic_dir=None):
    """Private ordinary-call endpoint; never accepts a caller pass flag.

    SAME ordinary call must precede this with fresh source-map cells and
    attachment validation. Fresh graph certificate and every track inequality
    are regenerated here, using actual source-bound recipe and curve controls.
    Source: recorded exact graph/moving-point component plus ambient extension
    construction N_A o J o H o M_AR o phi_A; J is not the interior anchor PL map.
    """
    from periodic_graph_family import compile_family
    from periodic_graph_certificate import certify
    import marked_display_source_map as maps
    from marked_display_source_map_pipeline import source_map_geometry_payload
    from pathlib import Path
    import time
    def encode(obj):
        if isinstance(obj,Q):return str(obj.numerator)+'/'+str(obj.denominator)
        if isinstance(obj,dict):return {k:encode(v) for k,v in obj.items()}
        if isinstance(obj,(list,tuple)):return [encode(v) for v in obj]
        return obj
    def keep(name,data):
        if diagnostic_dir is None:return
        target=Path(diagnostic_dir)/name
        need(not target.exists(),'refusing transport diagnostic overwrite')
        target.write_bytes(maps._source_map_json(encode(data)))
    source=copy.deepcopy(view['source_document'])
    original=[tuple(map(rat,p['point_svg_exact'])) for p in compiled['points']]
    family=compile_family(view,recipe)
    geometry_hash=hashlib.sha256(maps._source_map_json(source_map_geometry_payload(view))).hexdigest()
    recipe_hash=hashlib.sha256(maps._source_map_json(recipe)).hexdigest()
    def stream():
        digest=hashlib.sha256();counter=[0]
        def emit(kind,key,value):
            digest.update(json.dumps([kind,key,value],sort_keys=True,separators=(',',':')).encode());counter[0]+=1
        return digest,counter,emit
    graph_digest,graph_count,graph_emit=stream()
    graph=certify(family,emit=graph_emit)
    graph_summary=dict(schema=graph['schema'],source_document_sha256=graph['source_document_sha256'],
        leaves=graph['leaves'],counts=graph['counts'],outer_halfspaces=graph['outer_halfspaces'],inequality_records=graph_count[0],certificate_stream_sha256=graph_digest.hexdigest())
    keep('periodic-fresh-graph-summary.json',graph_summary)
    starts,inverse=inverse_anchor_tracks(compiled,family,recipe)
    keep('periodic-fresh-inverse-anchor.json',inverse)
    attempts=[]
    for steps in (16,32,64):
        before=time.monotonic();attempt=dict(slabs=int(steps));attempts.append(attempt)
        try:
            tracks=propose_tracks(family,starts,steps)
            need(all(t['points'][0]==q for t,q in zip(tracks,starts)),'inverse start changed')
            keep('periodic-candidate-tracks-'+str(steps)+'.json',tracks)
            digest,count,emit=stream();coverage=exact_track_check(family,tracks,emit=emit)
            moved=[t['points'][-1] for t in tracks]
            decoded,rounding=rounding_check(family,moved)
            # Every rounding inequality ran; its full diagnostic is optional,
            # compact logical report binds regenerated coverage and digest only.
            rounding_raw=maps._source_map_json(encode(rounding))
            keep('periodic-rounding-proof-'+str(steps)+'.json',rounding)
            points=endpoint_wire_payload(original,moved,decoded)
            need(view['source_document']==source,'point transport mutated source')
            attempt.update(status='exact all-time transport passed',seconds=time.monotonic()-before)
            result=dict(schema='periodic-moving-point-transport/1',source_document_sha256=compiled['source_document_sha256'],
                source_representation_sha256=compiled['source_representation_sha256'],base_geometry_sha256=geometry_hash,
                base_recipe_sha256=recipe_hash,graph=graph_summary,inverse_anchor=inverse,
                stages=family['stages']+['relative graph/boundary extension J fixes moved points; not interior anchor PL map','final analytic normalization identity on certified endpoint/rounding segment'],
                tracks=encode(tracks),points=points,coverage=coverage,inequality_records=count[0],
                certificate_stream_sha256=digest.hexdigest(),rounding=dict(points=5,point_pairs=10,
                    graph_rectangles=5*len(family['pairs']),proof_sha256=hashlib.sha256(rounding_raw).hexdigest()),
                selected_map='N_A compose J compose H compose M_AR compose phi_A',
                literal_F1_point_track_claim=False,source_points_relocated=False,attempts=attempts)
            keep('periodic-fresh-point-transport.json',result)
            return result
        except Exception as error:
            attempt.update(status='bounded candidate rejected',error=str(error),seconds=time.monotonic()-before)
            keep('periodic-candidate-failure-'+str(steps)+'.json',attempt)
    raise ValueError('periodic moving-point transport: all bounded candidates rejected: '+json.dumps(attempts))
