"""Exact same-region certificate for disjoint closed zero-port Jordan loops.

Sources: verify_bezier_attachment exact hull/separating-axis predicates and
build_rendered_mesh; port_domain.directed_seed source interval correspondence.
Parameter-matched local interval hulls conservatively contain the FULL final-A-PL to
actual emitted cubic homotopy, including fairing and endpoint correction.
No claim that this arbitrary graph homotopy is embedded is needed: avoiding
each labelled Jordan loop preserves its winding number about the point.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json


def need(ok, message):
    if not ok:
        raise ValueError('periodic point transport: '+message)


def rat(x):
    # Native Julia Rational strings may use //; canonical output always n/d.
    return Q(x.replace('//','/')) if isinstance(x,str) else Q(x)


def serial(q):
    return [str(x.numerator)+'/'+str(x.denominator) for x in q]


def canonical(x):
    return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()


def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def orient(a,b,c):return cross(sub(b,a),sub(c,a))


# Source: verify_bezier_attachment::hull/separate, exact rational arithmetic.
def hull(points):
    points=sorted(set(points))
    if len(points)<3:return points
    lower=[];upper=[]
    for chain,seq in ((lower,points),(upper,reversed(points))):
        for q in seq:
            while len(chain)>1 and orient(chain[-2],chain[-1],q)<=0:chain.pop()
            chain.append(q)
    return lower[:-1]+upper[:-1]


def separation(A,B):
    """Return an exact strict separating-axis witness, fail closed on contact."""
    for polygon in (A,B):
        for u,v in zip(polygon,polygon[1:]+polygon[:1]):
            e=sub(v,u);n=(-e[1],e[0])
            aa=[dot(n,q) for q in A];bb=[dot(n,q) for q in B]
            if max(aa)<min(bb):
                return dict(normal=serial(n),gap=str(min(bb)-max(aa)))
            if max(bb)<min(aa):
                return dict(normal=serial(n),gap=str(min(aa)-max(bb)))
    return None


def split_cubic(q,t):
    # Source: verify_bezier_attachment::split, exact de Casteljau at retained t.
    blend=lambda a,b:tuple((1-t)*x+t*y for x,y in zip(a,b))
    ab=[blend(a,b) for a,b in zip(q,q[1:])]
    bc=[blend(a,b) for a,b in zip(ab,ab[1:])];z=blend(*bc)
    return [q[0],ab[0],bc[0],z],[z,bc[1],ab[2],q[-1]]


def restrict_cubic(q,lo,hi):
    need(0<=lo<hi<=1,'cubic restriction outside actual span')
    right=q if lo==0 else split_cubic(q,lo)[1]
    return right if hi==1 else split_cubic(right,(hi-lo)/(1-lo))[0]


def separated_transport_leaves(line,cubic,segment,depth=0,path=()):
    # Source: verify_bezier_attachment dyadic de Casteljau hull subdivision.
    # Strict separation is unchanged; only conservative hull granularity changes.
    shape=hull(line+cubic);witness=separation(shape,segment)
    if witness is not None:return [dict(dyadic_subinterval=list(path),**witness)]
    need(depth<25,'point/rounding segment meets refined full fairing+anchor sweep or depth bound')
    a,b=split_cubic(line,Q(1,2));c,d=split_cubic(cubic,Q(1,2))
    return (separated_transport_leaves(a,c,segment,depth+1,path+(0,))+
            separated_transport_leaves(b,d,segment,depth+1,path+(1,)))


def minimum(c,b,a):
    values=[c,a+b+c]
    if a>0 and 0 < -b/(2*a) < 1:values.append(c-b*b/(4*a))
    return min(values)


def outer_margin(u,v,U,V,p):
    # Source: verify_bezier_attachment's moving_det quadratic minimum.
    e=sub(v,u);de=sub(sub(V,U),e);w=sub(p,u);dw=sub(u,U)
    return minimum(cross(e,w),cross(de,w)+cross(e,dw),cross(de,dw))


def require_distinct_decoded_marks(annotation):
    marks=[tuple(map(rat,p['point_svg_decoded'])) for p in annotation['points'] if p['kind']=='basepoint']
    need(len(marks)==2 and marks[0]!=marks[1],
         'distinct native basepoints collapse to same decoded SVG point')
    points=annotation['points']
    for i,a in enumerate(points):
        for b in points[:i]:
            exact_a=tuple(map(rat,a['point_svg_exact']));exact_b=tuple(map(rat,b['point_svg_exact']))
            decoded_a=tuple(map(rat,a['point_svg_decoded']));decoded_b=tuple(map(rat,b['point_svg_decoded']))
            need(exact_a==exact_b or decoded_a!=decoded_b,
                 'distinct exact labelled points collapse to same decoded SVG point')


def literal_cut_chart_cycle(curve):
    """Reassemble a native-patch word as a literal planar ID cycle.

    Source: vendor seam_word records surface.gluings, not display-cut ports.
    No seam flag or quotient identification is used to close this cycle.
    """
    need(curve['closed'] and curve['segments'],'open/empty curve unsupported')
    path=[];seen_sources=set()
    for segment in curve['segments']:
        source=tuple(segment['source']);need(source not in seen_sources,'repeated semantic source interval');seen_sources.add(source)
        vertices=segment['vertices']
        need(len(vertices)>=2,'empty directed semantic interval')
        if path:need(path[-1]==vertices[0],'curve requires a quotient-bank jump to join')
        path.extend(vertices if not path else vertices[1:])
    need(path[0]==path[-1],'curve requires a quotient-bank jump to close')
    need(len(path)>=4 and len(set(path[:-1]))==len(path)-1,
         'literal cut-chart cycle repeats a vertex')
    return path


def require_injective_positions(points,seen,message):
    need(len(set(points))==len(points) and not (set(points)&seen),message)
    seen.update(points)


def require_actual_interval_controls(actual,shared):
    cubics=[[tuple(map(rat,p)) for p in q] for q in actual['beziers']]
    need(cubics and all(len(q)==4 for q in cubics),'incomplete actual cubic interval')
    need(all(a[-1]==b[0] for a,b in zip(cubics,cubics[1:])), 'actual cubic interval gap')
    need(cubics[0][0]==tuple(map(rat,shared[str(actual['vertices'][0])])) and
         cubics[-1][-1]==tuple(map(rat,shared[str(actual['vertices'][-1])])),
         'actual cubic endpoints disagree with bound source interval vertices')
    return cubics


def prove(view,recipe,annotation):
    import marked_display_source_map as maps
    need(not view['crossings'] and not recipe['source_chart']['crossings'],
         'only zero-crossing closed-loop case is supported')
    need(not any(b['ports'] for b in recipe['source_chart']['boundaries']),
         'visible ports are unsupported')
    require_distinct_decoded_marks(annotation)
    # H.plot freshly ran this exact endpoint graph checker on these bound
    # cubic controls. This fresh recipe component remains byte-for-byte bound
    # by the calling Julia finalizer; saved check flags never bypass rebuilding.
    need(recipe['checks']['attachment']['literal_normalized_beziers_to_exact_anchor_isotopy'] is True,
         'missing fresh exact endpoint graph gate')
    pm=recipe['port_map'];ids=pm['original_vertex_ids'];index={int(v):i for i,v in enumerate(ids)}
    _,R,A=maps.build_rendered_mesh(pm,recipe['exact_circle_anchors'])
    seed=recipe['port_domain']['directed_seed']
    need(not seed['crossings'],'seed has crossings')
    need(len(seed['curves'])==len(view['curves']),'visible curve count changed')
    # Source: measure_clearance.branches and check_fairing_isotopy's OLD-R
    # Float64 arc labels. They are parameter labels, not exact length claims.
    import copy
    import numpy as np
    from measure_clearance import branches
    oldseed=copy.deepcopy(seed)
    oldseed['positions']={str(v):[float(rat(x)) for x in pm['target_positions'][str(i)]]
                          for i,v in enumerate(ids)}
    positions,_,bs=branches(oldseed);lookup={}
    for branch in bs:
        q=np.array([positions[v] for v in branch['vertices']])
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))]
        need(arc[-1]>0,'zero source branch length')
        arc/=arc[-1]
        labels={v:rat(float(t)) for v,t in zip(branch['vertices'],arc)}
        closed=branch['vertices'][0]==branch['vertices'][-1]
        if closed:labels[branch['vertices'][0]]=Q(0)
        for u,v in zip(branch['vertices'],branch['vertices'][1:]):
            lookup[branch['curve'],u,v]=(branch,labels,closed)
    carriers=[];carrier_payload=[];source_loop_vertices=set();actual_loop_vertices=set();source_loop_positions=set();actual_loop_positions=set();loop_bindings=[]
    hole_vertices={int(v) for h in recipe['exact_circle_anchors']['holes'] for v in h['vertices']}
    for old,new in zip(seed['curves'],view['curves'],strict=True):
        for key in ('id','role','closed','word','seam_word'):
            need(old[key]==new[key],'source curve identity/order changed: '+key)
        source_cycle=literal_cut_chart_cycle(old);actual_cycle=literal_cut_chart_cycle(new)
        mapped_source={index[int(v)] for v in source_cycle[:-1]}
        need(not (mapped_source & hole_vertices),'visible curve meets display-cut hole bank')
        need(not (set(source_cycle[:-1]) & source_loop_vertices),'source loops share a cut-chart vertex')
        need(not (set(actual_cycle[:-1]) & actual_loop_vertices),'actual loops share a cut-chart vertex')
        source_positions=[A[index[int(v)]] for v in source_cycle[:-1]]
        require_injective_positions(source_positions,source_loop_positions,'distinct source cut-chart vertices alias exact positions')
        source_loop_vertices.update(source_cycle[:-1]);actual_loop_vertices.update(actual_cycle[:-1])
        loop_bindings.append(dict(curve=old['id'],semantic_segments=len(old['segments']),
            source_cycle_vertices=len(source_cycle)-1,actual_cycle_vertices=len(actual_cycle)-1,
            native_patch_seam_word=old['seam_word'],literal_cycle_closed_without_quotient_identification=True))
        need(len(old['segments'])==len(new['segments']),'semantic interval count changed')
        old_end=None;new_end=None;firstA=None;firstC=None;actual_break_positions=[]
        for i,(source,actual) in enumerate(zip(old['segments'],new['segments'],strict=True)):
            need(source['vertices'][0]==actual['vertices'][0] and source['vertices'][-1]==actual['vertices'][-1], 'semantic interval endpoint IDs changed')
            for key in ('source','edge','edge_sign'):
                need(source[key]==actual[key],'directed semantic interval changed: '+key)
            vertices=source['vertices'];poly=[A[index[int(v)]] for v in vertices]
            need(len(poly)>=2,'incomplete source PL interval')
            branch,labels,closed=lookup[old['id'],vertices[0],vertices[1]]
            params=[labels[v] for v in vertices]
            if closed and vertices[-1]==branch['vertices'][0]:params[-1]=Q(1)
            need(all(a<b for a,b in zip(params,params[1:])),'nonincreasing old-R parameter labels')
            cubics=require_actual_interval_controls(actual,view['source_endpoint_normalization']['shared_svg_endpoints'])
            spans=[[rat(x) for x in span] for span in actual.get('spline_parameter_intervals',[])]
            need(cubics and len(spans)==len(cubics) and all(len(q)==4 for q in cubics),
                 'actual per-cubic spline_parameter_intervals are required')
            need(all(len(span)==2 and span[0]<span[1] for span in spans),'invalid actual spline span')
            need(spans[0][0]==params[0] and spans[-1][1]==params[-1] and
                 all(a[1]==b[0] for a,b in zip(spans,spans[1:])),
                 'actual spline spans do not exactly cover original source interval')
            actual_break_positions.extend(q[0] for q in cubics)
            need(old_end is None or old_end==poly[0],'A-polyline directed interval join changed')
            need(new_end is None or new_end==cubics[0][0],'actual cubic directed interval join changed')
            if firstA is None:firstA=poly[0];firstC=cubics[0][0]
            old_end=poly[-1];new_end=cubics[-1][-1]
            carrier_payload.append(dict(curve=old['id'],source=source['source'],interval=i,
                A_polyline=[serial(p) for p in poly],old_R_parameter_labels=[str(t) for t in params],
                actual_spans=[[str(t) for t in span] for span in spans],
                actual_cubics=[[serial(p) for p in q] for q in cubics]))
            # Split at BOTH actual cubic spans and every original A-PL corner,
            # retaining the OLD-R labels. No A arclength recomputation occurs.
            for bi,(q,(lo,hi)) in enumerate(zip(cubics,spans,strict=True)):
                knots=sorted({lo,hi}|{t for t in params if lo<t<hi})
                for left,right in zip(knots,knots[1:]):
                    j=next(j for j in range(len(params)-1) if params[j]<=left<right<=params[j+1])
                    def pl(t):
                        z=(t-params[j])/(params[j+1]-params[j])
                        return tuple(x+z*(y-x) for x,y in zip(poly[j],poly[j+1]))
                    a,b=pl(left),pl(right)
                    line=[a,tuple((2*x+y)/3 for x,y in zip(a,b)),
                          tuple((x+2*y)/3 for x,y in zip(a,b)),b]
                    restricted=restrict_cubic(q,(left-lo)/(hi-lo),(right-lo)/(hi-lo))
                    shape=hull(line+restricted)
                    need(len(shape)>=2,'degenerate local interval hull unsupported')
                    carriers.append(dict(curve=old['id'],source=source['source'],interval=i,
                        spline_piece=bi,parameter_interval=[str(left),str(right)],hull=shape,line_controls=line,actual_controls=restricted))
        need(old_end==firstA and new_end==firstC,'loop endpoint closure changed')
        require_injective_positions(actual_break_positions,actual_loop_positions,'distinct actual cubic joins alias exact positions')
    holes=[]
    for hole in recipe['exact_circle_anchors']['holes']:
        ring=list(map(int,hole['vertices']))
        holes.append((hole['cut'],hole['side'],hull([chart[v] for chart in (R,A) for v in ring])))
    final_supports=[]
    for c in recipe['graph_free_collars']['collars']:
        final_supports.append((c['cut'],c['side'],hull([tuple(600*rat(x) for x in p) for p in c['outer_positions']])))
    proofs=[];decoded_marks=[]
    for point in annotation['points']:
        exact=tuple(map(rat,point['point_svg_exact']))
        decoded=tuple(map(rat,point['point_svg_decoded']))
        need(decoded==tuple(rat(float(x)) for x in point['point_svg']),
             'decoded marker differs from actual binary64 coordinate')
        segment=list(dict.fromkeys((exact,decoded)))
        witnesses=[]
        for c in carriers:
            try:
                leaves=separated_transport_leaves(c['line_controls'],c['actual_controls'],segment)
            except ValueError as error:
                raise ValueError(str(error)+' at '+str((point['name'],c['curve'],c['source'],c['parameter_interval']))) from error
            witnesses.append(dict(curve=c['curve'],source=c['source'],interval=c['interval'],spline_piece=c['spline_piece'],parameter_interval=c['parameter_interval'],separated_dyadic_leaves=leaves))
        hole_witness=[]
        for cut,side,shape in holes:
            witness=separation(shape,segment)
            need(witness is not None,'point/rounding segment meets moving hole hull: '+str((point['name'],cut,side)))
            hole_witness.append(dict(cut=cut,side=side,**witness))
        support_witness=[]
        for cut,side,shape in final_supports:
            witness=separation(shape,segment)
            need(witness is not None,'point/rounding segment meets final analytic support: '+str((point['name'],cut,side)))
            support_witness.append(dict(cut=cut,side=side,**witness))
        outer=list(map(int,recipe['exact_circle_anchors']['outer']))
        margin=min(outer_margin(R[u],R[v],A[u],A[v],p)
                   for u,v in zip(outer,outer[1:]+outer[:1]) for p in segment)
        need(margin>0,'point/rounding segment leaves moving exterior frame')
        if point['kind']=='basepoint':decoded_marks.append(decoded)
        proofs.append(dict(name=point['name'],kind=point['kind'],point_svg_exact=point['point_svg_exact'],
                           point_svg_decoded=point['point_svg_decoded'],curve_hull_separators=witnesses,
                           moving_hole_separators=hole_witness,final_support_separators=support_witness,
                           moving_outer_minimum_orientation=str(margin)))
    return dict(schema='marked-periodic-point-transport/1',
                source_document_sha256=annotation['source_document_sha256'],
                exact_carriers_sha256=sha256(canonical(carrier_payload)).hexdigest(),
                source_map_sha256_before_annotations=sha256(maps._source_map_json(recipe)).hexdigest(),
                proof='labelled Jordan-loop winding signatures preserved by full semantic-interval homotopies avoiding closed exact-to-decoded point segments; finite labelled points corrected within corresponding final complementary regions',
                scope='disjoint literal cut-chart closed Jordan loops, zero crossings/display ports; native patch seams reassembled without quotient closure; conservative local parameter-matched hulls',
                literal_loop_reassembly=loop_bindings,
                literal_F1_point_track_claim=False,points=proofs)
