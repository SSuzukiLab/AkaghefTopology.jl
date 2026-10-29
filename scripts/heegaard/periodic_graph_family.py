"""Source-bound old-R-PL/actual-cubic family, no point acceptance yet.
Source: existing v8 actual parameter correspondence and literal-loop reassembly.
Private ordinary-call preconditions: fresh common source-map cells+attachment.
"""
from periodic_point_transport import *

def compile_family(view,recipe):
    import marked_display_source_map as maps
    need(not view['crossings'] and not recipe['source_chart']['crossings'], 'requires zero crossings')
    need(not any(b['ports'] for b in recipe['source_chart']['boundaries']), 'display ports unsupported')
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
    visible={index[int(v)] for c in seed['curves'] for s in c['segments'] for v in s['vertices']}
    need(all(A[v]==R[v] for v in visible),'reverse anchor map moves visible graph; stage generalization required')
    need(not recipe['exact_bezier_replacements'],'nonidentity Bezier anchor correction requires explicit stage')
    pairs=[];groups={}
    for c in carriers:
        groups.setdefault(c['curve'],[]).append(c)
    for cid,parts in groups.items():
        for i,c in enumerate(parts):
            nextc=parts[(i+1)%len(parts)]
            need(c['line_controls'][-1]==nextc['line_controls'][0] and c['actual_controls'][-1]==nextc['actual_controls'][0], 'family source incidence mismatch')
            pairs.append(dict(start=c['line_controls'],end=c['actual_controls'],
                left=(cid,i),right=(cid,(i+1)%len(parts)),source=c['source'],interval=c['parameter_interval']))
    return dict(pairs=pairs,R=R,A=A,triangles=pm['triangles'],holes=holes,
        outer=[R[int(v)] for v in recipe['exact_circle_anchors']['outer']],
        final_supports=final_supports,loop_bindings=loop_bindings,
        source_document_sha256=recipe['source_document_sha256'],
        stages=['reverse certified A-to-R mesh; visible graph fixed',
                'fresh exact old-R-PL to actual cubic family',
                'Bezier endpoint-anchor correction identity verified from actual data'])
