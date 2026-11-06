"""Common marked-cut layout fairing; no named-example or genus-specific routing.

Reuses fair_marked_paths, check_smooth_marked, and the existing SVG renderer.
The input is an explicit feasible embedded layout, not endpoints to be routed.
"""
import argparse,copy,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from draw_marked_cut import render
from fair_marked_paths import optimise,export


def mapped_port_position(p,maps):
    """Replay pair_ports' shared ambient map, also for source-map witnesses."""
    for center,r,support,src,f in maps:
        q=p-center;rho=np.linalg.norm(q)
        if rho<1e-10 or rho>=support:continue
        a=src[0]+(math.atan2(q[1],q[0])-src[0])%(2*math.pi)
        t=np.clip((rho-r)/(support-r),0,1);weight=1-t*t*(3-2*t)
        angle=a+weight*(float(f(a))-a)
        return center+rho*np.array([math.cos(angle),math.sin(angle)])
    return p.copy()


def normalize_shared_endpoints(result,seed):
    """Normalize the new source-map branch by retained vertex identity.

    Source: fair_marked_paths.export and its shared source interval endpoints.
    Floating spline evaluation may produce distinct last bits for one crossing.
    Use one stored endpoint per vertex, shift its adjacent control by the same
    amount, then rerun the ordinary geometry/isotopy gates. Legacy callers do
    not invoke this helper. No equality tolerance is promoted to exact equality.
    """
    from check_smooth_marked import leaves
    fixed={x['vertex'] for x in seed['crossings']}
    fixed.update(p['vertex'] for b in seed['boundaries'] for p in b['ports'])
    endpoints={v for c in result['curves'] for s in c['segments'] for v in (s['vertices'][0],s['vertices'][-1])}
    canonical={v:np.array((seed if v in fixed else result)['positions'][str(v)])*600 for v in endpoints}
    changed=0;internal_changed=0;maximum=0.
    for c in result['curves']:
        for s in c['segments']:
            for first,v in ((True,s['vertices'][0]),(False,s['vertices'][-1])):
                q=s['beziers'][0 if first else -1];i,j=(0,1) if first else (3,2)
                delta=canonical[v]-q[i]
                if np.any(delta):changed+=1;maximum=max(maximum,float(np.max(np.abs(delta))))
                q[i]=canonical[v].tolist();q[j]=(np.array(q[j])+delta).tolist()
            for left,right in zip(s['beziers'],s['beziers'][1:]):
                delta=np.array(left[-1])-right[0]
                if np.any(delta):internal_changed+=1;maximum=max(maximum,float(np.max(np.abs(delta))))
                right[0]=left[-1][:];right[1]=(np.array(right[1])+delta).tolist()
            pieces=[np.array(q) for q in s['beziers']]
            points=[pieces[0][0]]+[leaf[-1] for q in pieces for leaf in leaves(q)]
            assert len(points)==len(s['vertices']), 'Endpoint normalization changed Bezier sampling topology'
            for v,q in zip(s['vertices'],points):result['positions'][str(v)]=(q/600).tolist()
    # Source: fair_marked_paths.export reconstructs D from these SAME directed
    # source intervals. Rebind after endpoint replacement: reversed copies may
    # otherwise retain the old point-list objects at a shared endpoint.
    disk=result.get('disk',[])
    if disk:
        edges={}
        for curve in result['curves']:
            for segment in curve['segments']:
                for u,v in zip(segment['vertices'],segment['vertices'][1:]):
                    edges[u,v]=(segment,True);edges[v,u]=(segment,False)
        beziers=[];previous=None
        for u,v in zip(disk,disk[1:]+disk[:1]):
            segment,forward=edges[u,v];key=(id(segment),forward)
            if key==previous:continue
            pieces=segment['beziers']
            beziers.extend(pieces if forward else [q[::-1] for q in pieces[::-1]])
            previous=key
        result['disk_beziers']=beziers
    result['source_endpoint_normalization']=dict(changed_occurrences=changed,internal_piece_joins_changed=internal_changed,maximum_coordinate_change_svg=maximum,
        shared_svg_endpoints={str(v):q.tolist() for v,q in canonical.items()})
    # The certified coordinates are the numeric Float64 controls, not their
    # eight-decimal approximation. The renderer's matching branch round-trips
    # them through SVG text and verifies that decoder contract separately.
    result['svg_coordinate_encoding']='float64-roundtrip'
    return result


def pair_ports(data):
    """One orientation-reversing circle parameter for each paired cut boundary.

Transport the entire drawing through disjoint annular angular homeomorphisms,
not just its endpoints. Outside the collars, crossings and geometry are fixed.
"""
    d=copy.deepcopy(data);pos={int(v):np.array(q)*600 for v,q in d['positions'].items()}
    pairing={a:b for u,v in d['boundary_vertex_pairs'] for a,b in ((u,v),(v,u))}
    for b in d['boundaries']:
        for p in b['ports']:
            assert abs(np.linalg.norm(pos[p['vertex']]-pos[b['cap']])-600*b['radius'])<1e-6,'Port is not on its hole boundary'
    maps=[];targets={}
    for cut in d['cuts']:
        pair=sorted([b for b in d['boundaries'] if b['cut']==cut],key=lambda b:b['row'])
        assert len(pair)==2 and [b['row'] for b in pair]==[0,1]
        top,bottom=pair
        bysource={p['source_vertex']:p for p in top['ports']}
        order=[v for v in top['vertices'] if v in bysource]
        # Old Koda geometry's boundary vertices include cap-subdivision nodes;
        # retain the native source boundary order explicitly when available.
        order=[v for v in top.get('source_vertices',top['vertices']) if v in bysource]
        assert len(order)==len(top['ports']) and len(order)==len(bottom['ports'])
        if not order:continue
        bp={p['source_vertex']:p for p in bottom['ports']}
        n=len(order);phases=[]
        for j,v in enumerate(order):
            for b,p,sgn in ((top,bysource[v],1),(bottom,bp[pairing[v]],-1)):
                q=pos[p['vertex']]-pos[b['cap']]
                phases.append(sgn*math.atan2(q[1],q[0])+2*math.pi*j/n)
        phase=math.atan2(sum(math.sin(a) for a in phases),sum(math.cos(a) for a in phases))
        for j,v in enumerate(order):
            angle=phase-2*math.pi*j/n
            targets[bysource[v]['vertex']]=angle;targets[bp[pairing[v]]['vertex']]=-angle
        already_paired=all(np.linalg.norm((pos[bysource[v]['vertex']]-pos[top['cap']])-(pos[bp[pairing[v]]['vertex']]-pos[bottom['cap']])*[1,-1])<1e-7 for v in order)
        if already_paired:
            for b in pair:
                for p in b['ports']:
                    q=pos[p['vertex']]-pos[b['cap']];targets[p['vertex']]=math.atan2(q[1],q[0])
        for b in pair:
            center=pos[b['cap']];r=b['radius']*600
            angle=sorted((math.atan2(*(pos[p['vertex']]-center)[::-1]),p['vertex']) for p in b['ports'])
            src=np.array([a for a,v in angle]);dst=np.array([targets[v] for a,v in angle])
            for j in range(1,len(dst)):
                dst[j]+=2*math.pi*math.ceil((dst[j-1]-dst[j]+1e-10)/(2*math.pi))
            assert dst[-1]-dst[0]<2*math.pi-1e-8,'Boundary port cyclic order mismatch'
            dst+=round(np.mean(src-dst)/(2*math.pi))*2*math.pi
            f=PchipInterpolator(np.r_[src-2*math.pi,src,src+2*math.pi],np.r_[dst-2*math.pi,dst,dst+2*math.pi])
            cross_dist=min((np.linalg.norm(pos[x['vertex']]-center) for x in d['crossings']),default=r+150)
            hole_dist=min((np.linalg.norm(pos[c['cap']]-center)-c['radius']*600 for c in d['boundaries'] if c is not b),default=r+150)
            support=min(r+100,cross_dist-5,(r+hole_dist)/2-2)
            assert support>r+3,'No collar room for changing ports without moving a crossing'
            maps.append((center,r,support,src,f))
    def mapped(p):
        return mapped_port_position(p,maps)
    out={v:mapped(p) for v,p in pos.items()};nextid=max(pos)+1;cache={};vertex_sources={}
    def edge(u,v):
        nonlocal nextid
        if (u,v) in cache:return cache[u,v]
        a,b=pos[u],pos[v];delta=b-a
        def sample(t0,t1,depth=0):
            p,q=mapped(a+t0*delta),mapped(a+t1*delta)
            error=max(np.linalg.norm(mapped(a+(t0+(t1-t0)*f)*delta)-((1-f)*p+f*q)) for f in (.25,.5,.75))
            if error<.002 and np.linalg.norm(delta)*(t1-t0)<8:return [(t1,q)]
            assert depth<25,'Port transport refinement failed'
            mid=(t0+t1)/2;return sample(t0,mid,depth+1)+sample(mid,t1,depth+1)
        points=sample(0,1);vs=[u]
        for t,q in points[:-1]:
            out[nextid]=q;vs.append(nextid)
            # Source: the dyadic subdivision in sample(). Retain the actual
            # source-edge parameter so a whole-mesh PL witness can include
            # every emitted curve corner, independently of rounded PCHIP.
            vertex_sources[str(nextid)]=dict(edge=[u,v],parameter=str(Fraction(t)))
            nextid+=1
        vs.append(v);cache[u,v]=vs;cache[v,u]=vs[::-1];return vs
    for c in d['curves']:
        for s in c['segments']:
            vs=s['vertices'];s['vertices']=[w for u,v in zip(vs,vs[1:]) for w in edge(u,v)[:-1]]+[vs[-1]]
            s.pop('beziers',None)
    vs=d.get('disk',[]);d['disk']=[w for u,v in zip(vs,vs[1:]+vs[:1]) for w in edge(u,v)[:-1]]
    d.pop('disk_beziers',None);d['positions']={str(v):(p/600).tolist() for v,p in out.items()}
    d['port_convention']='paired reflection (x,y) -> (x,-y), cyclic order retained, normal exits'
    # Source: the explicit mapped() construction above and SciPy PCHIP. These
    # coefficients retain the whole ambient collar map, so omitted source
    # cells can be transported too, rather than retaining only curve samples.
    d['port_transport']=dict(method='monotone angular maps in disjoint circular collars',crossings_fixed=True,
        vertex_sources=vertex_sources,
        maps=[dict(center_svg=center.tolist(),inner_radius_svg=float(r),outer_radius_svg=float(support),
            angle_origin=float(src[0]),spline_breaks=f.x.tolist(),spline_coefficients=f.c.tolist(),
            radial_blend='1 - t*t*(3 - 2*t), t=clamp((rho-inner)/(outer-inner),0,1)')
            for center,r,support,src,f in maps])
    return d


def _plot_once(data, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20., correct_ports=True, maxiter=1800, _source_map_geometry=False, _source_map_proof_output=None, _fairing_knot_refinement=None, _fairing_capture=None, _fairing_options_path=None, _source_map_context=None, _fairing_source_separations=None):
    from check_smooth_marked import check
    if data.get('schema')!='marked-cut-view/1':raise ValueError('Explicit marked-cut-view/1 embedded layout required')
    if 'source_map' in data:
        raise ValueError('Restore the source-bound reconstruction base before plotting a saved source-map view')
    if not all(np.isfinite(x) and x>0 for x in (minimum_separation,puncture_clearance,crossing_radius)):raise ValueError('Positive finite drawing clearances required')
    if data.get('disks'):raise ValueError('Use source-bound disk_pieces for joining-disk fills')
    if data.get('genus',0)>0 and not data.get('boundary_vertex_pairs'):raise ValueError('Source boundary vertex pairing is required; do not infer it from endpoints')
    seed=pair_ports(data) if correct_ports else copy.deepcopy(data)
    returned_seed=seed
    option_binding=None
    if _fairing_options_path:
        from fair_marked_paths import load_fairing_construction_options
        if _fairing_knot_refinement:raise ValueError('Cannot combine finite options with adaptive branch factors')
        _fairing_knot_refinement,option_binding=load_fairing_construction_options(seed,_fairing_options_path,
            policy=dict(minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,maxiter=maxiter))

    # Preserve the existing legacy empty-view result. Native/saved full-map
    # inputs continue through the mandatory reconstruction chain below.
    if not seed['curves'] and not _source_map_geometry:
        from verify_fixed_grid import verify
        result=copy.deepcopy(seed);result['geometry_check']=verify(result,seed)
        result['smooth_check']=dict(crossing_coordinates_unchanged=0,source_disk_boundary_unchanged=True,
            holes_and_ports_unchanged=True,regular_spline_pieces=True,no_extra_bezier_intersections=True,
            normal_port_tangents=True,empty_visible_family=True)
        result['layout_isotopy_check']=dict(method='identity on empty visible family',collisions=0,hole_entries=0)
        result['svg']=render(result)
        return result,seed
    # Source: fair_marked_paths.setup requires a nonempty branch matrix.
    # Empty graph is an identity proposal, not an exemption from source/mesh,
    # D, infinity, gluing, whole-step or output-readback checks below.
    if not seed['curves']:
        assert not seed['crossings'] and not any(b['ports'] for b in seed['boundaries'])
        assert not seed.get('disk') and not seed.get('disk_pieces') and not seed.get('disk_region_marks'), 'Empty visibility cannot hide a required D boundary'
        result=copy.deepcopy(seed)
    else:
        bs,c,_=optimise(seed,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,maxiter=maxiter,knot_refinement=_fairing_knot_refinement,source_separations=_fairing_source_separations)
        result=export(seed,bs,c)
        if option_binding is not None:result['spline_fairing']['construction_options']=option_binding
    # Internal full-source reconstruction branch only. Repair the ACTUAL data
    # before ordinary geometry/isotopy checks and SVG emission; normalizing a
    # separate proof copy would not bind the output shown or saved by the app.
    if _source_map_geometry:normalize_shared_endpoints(result,seed)
    # Failure capture deliberately retains the pre-finish source seed; the
    # current candidate carries any literal paths established before failure.
    if _fairing_capture is not None:_fairing_capture.update(seed=seed,candidate=result)
    if _fairing_source_separations:
        from fair_marked_paths import verify_source_separations
        verify_source_separations(seed,result,_fairing_source_separations)
    from check_fairing_isotopy import check as check_motion
    def check_output(proposed,original):
        checked=check(proposed,original,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius)
        checked['layout_isotopy_check']=check_motion(original,checked)
        if not original['curves']:
            checked['smooth_check']['empty_visible_family']=True
            checked['geometry_provenance']='identity proposal on empty visible family; full source-map checks retained'
        return checked
    if not _source_map_geometry:result=check_output(result,seed)
    result['drawing_policy']=dict(minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,port_collar_radius=20)
    if _source_map_geometry:
        from marked_display_source_map_pipeline import construct_source_map
        from marked_display_source_map import encode_source_map,_source_map_json
        context=data.get('_source_map_input')
        if not context:raise ValueError('A freshly source-bound reconstruction context is required')
        fixed=copy.deepcopy(data);fixed.pop('_source_map_input',None)
        policy=dict(minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius)
        if _source_map_proof_output:
            stage=Path(_source_map_proof_output)
            stage.with_suffix('.proposed.json').write_bytes(_source_map_json(result))
            stage.with_suffix('.seed.json').write_bytes(_source_map_json(seed))
        # Build the source-bound cut chains BEFORE rendering split-D fills.
        # The new decoder branch consumes those chains, never temporary arcs.
        recipe=construct_source_map(context['prepared'],context['coarse'],fixed,result,seed,policy=policy,check_output=check_output,_call_context=_source_map_context,_proof_output=_source_map_proof_output)
        # Only a fully checked finish enriches the returned/CLI sidecar seed.
        # The original seed remains the immutable context key/failure capture.
        returned_seed=copy.deepcopy(seed)
        returned_seed['source_map_boundary_paths']=copy.deepcopy(recipe['attachment_boundaries'])
        # Keep the completed mathematical stages before a later envelope/body
        # assertion can fail. This is the caller's existing temporary directory.
        if _source_map_proof_output:
            Path(_source_map_proof_output).write_bytes(_source_map_json(recipe))
        result['source_map']=encode_source_map(recipe)
        result['source_map_verification_pending']=True
        for key in ('_source_map_input','source_map_mesh','source_chart'):result.pop(key,None)
    else:
        result['svg']=render(result)
    # Source: draw_marked_cut.orientation_glyph_layer; optional display layer
    # uses the same directed source routes. The accepted default SVG is intact.
    try:
        result['orientation_svg']=render(result,orientation_glyphs=True)
        result.pop('orientation_display_error',None)
    except ValueError as exc:
        result.pop('orientation_svg',None)
        result['orientation_display_error']=str(exc)
    if _source_map_geometry:
        from marked_display_source_map import _source_map_json
        if _source_map_proof_output:
            Path(_source_map_proof_output).with_suffix('.candidate.json').write_bytes(_source_map_json(result))
        # HTTP body limits belong to serve.py; this result is an in-process plot.
    return result,returned_seed



def plot(data, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20., correct_ports=True, maxiter=1800, _source_map_geometry=False, _source_map_proof_output=None):
    """Bounded knot construction on the SAME original feasible source.

    Source: existing fairer/checkers; only rejected branch knots are refined.
    Native preparation/fixed mesh are already in data and are never rerun here.
    Unknown failures and exhausted bounds remain failures with actual artifacts.
    """
    import os
    finite_options=os.environ.get('HEEGAARD_FAIRING_CONSTRUCTION_OPTIONS')
    import hashlib,tempfile
    from fair_marked_paths import offending_fairing_branches,source_separations_from_failure
    source_map_context=None
    if _source_map_geometry:
        from marked_display_source_map_pipeline import _SourceMapCallContext
        source_map_context=_SourceMapCallContext()
    refinement={};source_separations=[];history=[];folder=None
    # Fixed construction bounds, never acceptance thresholds.
    maximum_attempts=8;maximum_factor=16
    def retain(capture,error,attempt):
        nonlocal folder
        if folder is None:
            folder=(Path(_source_map_proof_output).parent/'fairing-refinement') if _source_map_proof_output else Path(tempfile.mkdtemp(prefix='marked-fairing-refinement-'))
            folder.mkdir(exist_ok=True)
            (folder/'source-seed.json').write_text(json.dumps(capture['seed'],allow_nan=False),encoding='utf-8')
        path=folder/f'attempt-{attempt}.json'
        path.write_text(json.dumps(capture['candidate'],allow_nan=False),encoding='utf-8')
        history.append(dict(attempt=attempt,candidate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),error=str(error),source_separations=copy.deepcopy(source_separations),knot_refinement=[dict(curve=k[0],start=k[1],end=k[2],factor=v) for k,v in sorted(refinement.items())]))
        (folder/'history.json').write_text(json.dumps(dict(maximum_attempts=maximum_attempts,maximum_factor=maximum_factor,attempts=history),indent=2),encoding='utf-8')
    def retain_original(capture,error,attempt):
        # Diagnostic I/O must never replace the rejected mathematical result.
        try:
            retain(capture,error,attempt)
        except Exception as retention_error:
            error.add_note('Fairing retention failed: '+type(retention_error).__name__+': '+str(retention_error))
            return False
        error.add_note('Actual fairing-stage graph retained at '+str(folder))
        return True
    for attempt in range(1 if finite_options else maximum_attempts):
        capture={}
        try:
            if finite_options:
                # Same finite numerical call and gates; capture is diagnostic.
                result=_plot_once(data,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,correct_ports=correct_ports,maxiter=maxiter,_source_map_geometry=_source_map_geometry,_source_map_proof_output=_source_map_proof_output,_fairing_options_path=finite_options,_fairing_capture=capture,_source_map_context=source_map_context)
            else:
                result=_plot_once(data,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,correct_ports=correct_ports,maxiter=maxiter,_source_map_geometry=_source_map_geometry,_source_map_proof_output=_source_map_proof_output,_fairing_knot_refinement=refinement,_fairing_capture=capture,_source_map_context=source_map_context,_fairing_source_separations=source_separations)
            return result
        except Exception as error:
            if not capture:raise
            # Retain before any decoder can fail. Only a recognized ordinary
            # AssertionError may change knots; every other error stops here.
            if not retain_original(capture,error,attempt):raise
            if finite_options or not isinstance(error,AssertionError):raise
            try:
                selected=offending_fairing_branches(capture['seed'],capture['candidate'],error,crossing_radius=crossing_radius)
                proposed=source_separations_from_failure(capture['seed'],capture['candidate'],error)
                for record in proposed:
                    if record not in source_separations:source_separations.append(record)
            except Exception as diagnostic_error:
                error.add_note('Fairing diagnostic decoder failed: '+type(diagnostic_error).__name__+': '+str(diagnostic_error))
                selected=None
            if not selected:raise
            if attempt+1==maximum_attempts or any(refinement.get(k,1)>=maximum_factor for k in selected):
                raise AssertionError(('Bounded fairing refinement failed',str(error),'actual stages retained at '+str(folder))) from error
            for k in selected:refinement[k]=2*refinement.get(k,1)
            print('Fairing knot refinement',attempt+1,'branches',len(selected),'maximum factor',max(refinement.values()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--maxiter',type=int,default=1800)
    p.add_argument('--prepared',type=Path);p.add_argument('--minimum-separation',type=float,default=12.)
    p.add_argument('--source-map-geometry',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--source-map-proof-output',type=Path,help=argparse.SUPPRESS)
    p.add_argument('--recheck',action='store_true',help='Recheck existing spline output against its saved seed; do not optimise again')
    p.add_argument('--puncture-clearance',type=float,default=12.);p.add_argument('--crossing-radius',type=float,default=20.)
    a=p.parse_args();d=json.loads(a.input.read_text(encoding='utf-8'))
    if a.recheck:
        if 'source_map' in d:
            raise ValueError('Spline-only recheck cannot revalidate a saved full source map; use bound reconstruction replay')
        from check_smooth_marked import check
        from check_fairing_isotopy import check as check_motion
        seed=json.loads(a.input.with_suffix('.seed.json').read_text(encoding='utf-8'))
        d=check(d,seed,minimum_separation=a.minimum_separation,puncture_clearance=a.puncture_clearance,crossing_radius=a.crossing_radius)
        d['layout_isotopy_check']=check_motion(seed,d);d['svg']=render(d)
        a.output.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8');a.output.with_suffix('.svg').write_text(d['svg'],encoding='utf-8')
        print(json.dumps(d['smooth_check']),flush=True);raise SystemExit(0)
    if a.prepared:
        prepared=json.loads(a.prepared.read_text(encoding='utf-8'))
        if 'source_map' in d:
            if not a.source_map_geometry:raise ValueError('Saved source-map input requires full reconstruction replay')
            from marked_display_source_map_pipeline import restore_source_map_view
            d=restore_source_map_view(d,prepared)
        assert prepared['source_cut_corner_classes']==d['source_cut_corner_classes']
        d['boundary_vertex_pairs']=prepared['boundary_vertex_pairs']
        for b in d['boundaries']:
            b['source_vertices']=next(c['vertices'] for c in prepared['boundaries'] if c['cut']==b['cut'] and c['side']==b['side'])
    result,seed=plot(d,maxiter=a.maxiter,minimum_separation=a.minimum_separation,puncture_clearance=a.puncture_clearance,crossing_radius=a.crossing_radius,_source_map_geometry=a.source_map_geometry,_source_map_proof_output=a.source_map_proof_output)
    if a.source_map_proof_output:
        if not a.source_map_geometry:raise ValueError('Source-map proof output requires the complete reconstruction branch')
        # Source: bounded decode_source_map; the SAME calling Julia process
        # reads these exact canonical bytes, validates every cell/collar, and
        # only then removes source_map_verification_pending. No nested Julia.
        from marked_display_source_map import decode_source_map,_source_map_json
        a.source_map_proof_output.write_bytes(_source_map_json(decode_source_map(result['source_map'])))
    a.output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    a.output.with_suffix('.svg').write_text(result.pop('svg'),encoding='utf-8')
    a.output.with_suffix('.seed.json').write_text(json.dumps(seed,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result['smooth_check']),flush=True)
