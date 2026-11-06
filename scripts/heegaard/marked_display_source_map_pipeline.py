"""Whole-source reconstruction used by the existing marked-display pipeline.

Source: marked_display_source_map builders/validators, fixed_mesh's retained
exact parent geometry and prepare_marked_display's native source intervals.
This module creates a verification request for the SAME calling Julia process;
Python transport success alone never authorizes publishing the marked view.
"""
import copy
from fractions import Fraction as F


GEOMETRY_FIELDS=('positions','boundaries','curves','crossings','disk','disk_beziers','disk_pieces',
                 'disk_region_marks','display_request','svg_coordinate_encoding',
                 'source_endpoint_normalization','source_map_boundary_paths')


def source_map_geometry_payload(view):
    """Canonicalize only schema-declared binary64 geometry, retaining all data.

    Source: fixed_mesh/export/normalize_shared_endpoints and the renderer's
    float64-roundtrip numeric decoder. JSON3 may serialize a coordinate 0.0 as
    0. Both decode to the same binary64 point; IDs, orientation signs, source
    words and every other scalar retain their original JSON type and value.
    This is a hash representation only: it never mutates the saved view.
    """
    import math
    geometry=copy.deepcopy({k:view.get(k) for k in GEOMETRY_FIELDS})
    def scalar(x):
        if type(x) not in (int,float):raise ValueError('Geometry coordinate must be a JSON number, not a Boolean or string')
        y=float(x)
        if not math.isfinite(y):raise ValueError('Geometry coordinate must be finite binary64')
        return y
    def point(q):
        if not isinstance(q,list) or len(q)!=2:raise ValueError('Geometry point must have two coordinates')
        return [scalar(x) for x in q]
    def cubics(rows):
        for q in rows:
            if not isinstance(q,list) or len(q)!=4:raise ValueError('Cubic must have four controls')
        return [[point(p) for p in q] for q in rows]
    geometry['positions']={v:point(q) for v,q in geometry['positions'].items()}
    for boundary in geometry['boundaries']:boundary['radius']=scalar(boundary['radius'])
    for curve in geometry['curves']:
        for segment in curve['segments']:
            if 'beziers' in segment:segment['beziers']=cubics(segment['beziers'])
            # Source: adaptive fairer's actual normalized branch spans.
            if 'spline_parameter_intervals' in segment:
                segment['spline_parameter_intervals']=[point(q) for q in segment['spline_parameter_intervals']]
    if geometry['disk_beziers'] is not None:geometry['disk_beziers']=cubics(geometry['disk_beziers'])
    normalization=geometry['source_endpoint_normalization']
    if normalization is not None:
        normalization['shared_svg_endpoints']={v:point(q) for v,q in normalization['shared_svg_endpoints'].items()}
        # This bound is also explicitly produced as a binary64 SVG distance;
        # counts alongside it remain integers and are not normalized.
        normalization['maximum_coordinate_change_svg']=scalar(normalization['maximum_coordinate_change_svg'])
    paths=geometry['source_map_boundary_paths']
    if paths is not None:
        for boundary in paths['boundaries']:boundary['points_svg']=[point(q) for q in boundary['points_svg']]
    # Source: same native points/chain consumed by ordinary SVG renderer.
    if 'periodic_annotations' in view:
        from periodic_annotations import geometry_payload
        geometry['periodic_annotations']=geometry_payload(view['periodic_annotations'])
    for key in ('exterior_face_marker', 'exterior_face_glyph'):
        if key in view:
            geometry[key] = copy.deepcopy(view[key])
    return geometry


def verify_source_map_geometry(view,recipe):
    """Check the existing geometry digest after the declared numeric decoder.

    Legacy recipes retain their exact original digest. They may match either
    their original type-sensitive JSON or the producer's binary64 coordinate
    representation restored at the named schema fields. No digest is rewritten
    and no field is omitted. New recipes explicitly select the latter encoding.
    """
    import hashlib
    import marked_display_source_map as maps
    mode=recipe.get('output_geometry_numeric_encoding')
    if mode not in (None,'schema-binary64/1'):raise ValueError('Unknown source-map geometry numeric encoding')
    canonical=hashlib.sha256(maps._source_map_json(source_map_geometry_payload(view))).hexdigest()
    expected=recipe['output_geometry_sha256']
    if canonical==expected:return 'schema-binary64/1'
    if mode is None:
        legacy=hashlib.sha256(maps._source_map_json({k:view.get(k) for k in GEOMETRY_FIELDS})).hexdigest()
        if legacy==expected:return 'legacy-type-sensitive-json'
    raise ValueError('Source-map displayed geometry changed')


def restore_source_map_view(view, prepared):
    """Bind a saved recipe to its actual source/output, then rebuild its base.

    Source: bounded source-map decoder and verify_native_chart. A faired saved
    position table is not reused as the old triangulation's feasible base.
    Every mathematical gate is rerun by construct_source_map and Julia; saved
    check booleans and the validation string are never authorization tokens.
    """
    import hashlib
    import marked_display_source_map as maps
    recipe=maps.decode_source_map(view['source_map'])
    if 'periodic_annotations' in view or 'periodic_annotations' in recipe:
        assert view.get('periodic_annotations')==recipe.get('periodic_annotations'), 'Saved periodic annotations changed'
        from periodic_annotations import verify_numeric
        verify_numeric(view)
    authored_marker = 'exterior_face_marker' in view['source_document'].get('metadata', {})
    if authored_marker or 'exterior_face_marker' in view or 'exterior_face_glyph' in recipe:
        if not authored_marker or 'exterior_face_marker' not in view or 'exterior_face_glyph' not in recipe:
            raise ValueError('Saved exterior marker source/payload presence disagrees')
        if view.get('exterior_face_glyph') != recipe.get('exterior_face_glyph'):
            raise ValueError('Saved exterior glyph changed')
        from exterior_face_annotations import verify
        verify(view, recipe)
        maps.verify_svg_coordinates(view)
    digest=lambda value:hashlib.sha256(maps._source_map_json(value)).hexdigest()
    assert digest(view['source_document'])==recipe['source_document_sha256'], 'Source-map native document changed'
    verify_source_map_geometry(view,recipe)
    assert hashlib.sha256(view['svg'].encode()).hexdigest()==recipe['output_svg_sha256'], 'Source-map SVG changed'
    coarse=copy.deepcopy(recipe['source_chart']);coarse['source_document']=view['source_document']
    maps.verify_native_chart(coarse,prepared)
    fixed=copy.deepcopy(recipe['fixed_replay_view']);fixed['source_document']=view['source_document']
    fixed['source_map_mesh']=copy.deepcopy(recipe['fixed_map']['source'])
    assert all(fixed['positions'][v]==q for v,q in recipe['fixed_map']['target_positions'].items())
    assert fixed['display_request']==coarse['display_request']==recipe['display_request']
    fixed['_source_map_input']=dict(prepared=prepared,coarse=coarse)
    return fixed


def verify_fixed_source_paths(coarse, fixed, mesh, native_collars, circle_collars):
    """Bind each directed fixed path to its native path and both endpoint collars.

    Source: the explicit cap-spoke construction in fixed_mesh and the retained
    sector radial bridges. Equal ordered exact source polylines on actual mesh
    edges have the saved target polyline as their facewise PL image. Only
    forward collinear subdivisions are removed; no endpoint route is inferred.
    """
    old={int(v):tuple(map(F,q)) for v,q in coarse['positions'].items()}
    source={int(v):tuple(map(F,q)) for v,q in mesh['source_positions'].items()}
    target={int(v):tuple(map(F,q)) for v,q in fixed['positions'].items()}
    edges={tuple(sorted((u,v))) for t in mesh['triangles'] for u,v in zip(t,t[1:]+t[:1])}
    ports={p['vertex']:(i,p) for i,b in enumerate(coarse['boundaries']) for p in b['ports']}
    def simplified(points):
        out=[]
        for point in points:
            if out and out[-1]==point:continue
            while len(out)>1:
                a,b=out[-2:];u=(b[0]-a[0],b[1]-a[1]);v=(point[0]-b[0],point[1]-b[1])
                if u[0]*v[1]!=u[1]*v[0] or u[0]*v[0]+u[1]*v[1]<=0:break
                out.pop()
            out.append(point)
        return out
    checked=endpoints=0
    for oldcurve,newcurve in zip(coarse['curves'],fixed['curves'],strict=True):
        for key in ('id','role','closed','word','seam_word'):assert oldcurve[key]==newcurve[key]
        for oldsegment,newsegment in zip(oldcurve['segments'],newcurve['segments'],strict=True):
            for key in ('source','edge','edge_sign'):assert oldsegment[key]==newsegment[key]
            vs=oldsegment['vertices'];ns=newsegment['vertices']
            expected=[old[v] for v in vs];trim_left=int(vs[0] in ports);trim_right=int(vs[-1] in ports)
            for first in (True,False):
                port_id=vs[0 if first else -1]
                if port_id not in ports:continue
                bi,port=ports[port_id];vword=vs[::-1] if first else vs;nword=ns[::-1] if first else ns
                spoke,display_port=nword[-2:];native=port['source_vertex']
                assert spoke in mesh['boundary_cycles'][bi]
                boundary=coarse['boundaries'][bi]
                nc=native_collars['boundaries'][bi];cc=circle_collars['boundaries'][bi]
                assert (nc['cut'],nc['side'])==(cc['cut'],cc['side'])==(boundary['cut'],boundary['side'])
                sector=next(s for s in nc['sectors'] if s['boundary_edge'][0]==native)
                sb=[tuple(map(F,q)) for q in sector['source_chains'][0]]
                tb=[tuple(map(F,q)) for q in sector['target_chains'][0]]
                assert sb[-1]==old[native]
                assert tb==[sb[0]]+[old[v] for v in boundary['source_cap_spokes'][str(native)][:0:-1]]+[source[spoke]]
                j=vword.index(native);assert j>0
                assert simplified([old[vword[j-1]],sb[0],old[native]])==[old[vword[j-1]],old[native]]
                sector=next(s for s in cc['sectors'] if s['boundary_edge'][0]==spoke)
                sb=[tuple(map(F,q)) for q in sector['source_chains'][0]]
                tb=[tuple(map(F,q)) for q in sector['target_chains'][0]]
                assert sb[-1]==target[spoke] and tb==[sb[0],target[spoke],target[display_port]]
                assert simplified([target[nword[-3]],sb[0],target[spoke]])==[target[nword[-3]],target[spoke]]
                expected[0 if first else -1]=source[spoke];endpoints+=1
            clipped=ns[trim_left:len(ns)-trim_right if trim_right else len(ns)]
            assert simplified(expected)==simplified([source[v] for v in clipped])
            assert all(tuple(sorted((u,v))) in edges for u,v in zip(clipped,clipped[1:]))
            checked+=1
    for key in ('crossings','disk_pieces','disk_region_marks','source_document','display_request'):
        assert coarse.get(key)==fixed.get(key),('Fixed mesh changed a marked feature',key)
    for crossing in coarse['crossings']:
        v=crossing['vertex'];assert source[v]==old[v]
    disk=fixed.get('disk',[])
    if disk:
        assert simplified([old[v] for v in coarse['disk']]+[old[coarse['disk'][0]]])==simplified([source[v] for v in disk]+[source[disk[0]]])
        assert all(tuple(sorted((u,v))) in edges for u,v in zip(disk,disk[1:]+disk[:1]))
    return dict(directed_source_intervals=checked,source_bound_port_endpoints=endpoints,
                fixed_paths_bound_without_expanded_overlay=True)


def _prepare_source_map_prefix(prepared, coarse, fixed, seed):
    """Build the full source-only data from actual inputs, never a pass flag."""
    import marked_display_source_map as maps
    from marked_cut_layout import source_mesh
    native=maps.verify_native_chart(coarse,prepared)
    mesh=source_mesh(coarse)
    assert mesh==fixed['source_map_mesh']
    collars=maps.build_collar_sectors(coarse,mesh)
    circle_input,circle_mesh=maps.circle_collar_input(fixed)
    circle_collars=maps.build_collar_sectors(circle_input,circle_mesh)
    path_binding=verify_fixed_source_paths(coarse,fixed,mesh,collars,circle_collars)
    domain,args=maps.build_port_domain(fixed,circle_input,seed)
    port_map=maps.build_port_map(*args)
    anchors=maps.build_exact_circle_anchors(port_map,fixed,boundary_cycles=domain['holes'])
    support=maps.build_graph_free_collars(port_map,anchors)
    trace=maps.build_boundary_trace(coarse,mesh,fixed,collars,circle_collars)
    gluing=maps.build_boundary_gluing(coarse,fixed,trace)
    family=maps.build_fixed_quotient(port_map,anchors,gluing,support,fixed)
    literal=maps.build_literal_boundary_paths(port_map,anchors,fixed,family)
    return dict(native=native,mesh=mesh,collars=collars,circle_collars=circle_collars,path_binding=path_binding,domain=domain,port_map=port_map,anchors=anchors,support=support,family=family,literal=literal)


class _SourceMapCallContext:
    """Private, empty-created, in-memory prefix owner for ONE plot invocation.

    No dictionary/file/receipt can supply a prepared context. Full source data
    is freshly built on first use. Copies isolate every candidate from storage;
    exact input/producer and stored-data digests reject mutation between uses.
    """
    __slots__=('_bindings','_components','_integrity')
    def __init__(self):
        self._bindings=None;self._components=None;self._integrity=None
    @staticmethod
    def _digest(value):
        import hashlib
        from marked_display_source_map import _source_map_json
        return hashlib.sha256(_source_map_json(value)).hexdigest()
    @classmethod
    def _bind(cls,prepared,coarse,fixed,seed,policy):
        import hashlib
        from pathlib import Path
        import marked_display_source_map as maps
        base=Path(maps.__file__).parent
        producers=[Path(__file__)]+[base/name for name in (
            'marked_display_source_map.py','marked_cut_layout.py',
            'marked_display_pipeline.py','draw_marked_cut.py','periodic_annotations.py',
            'fair_marked_paths.py','check_smooth_marked.py','check_fairing_isotopy.py')]
        return dict(inputs={name:cls._digest(value) for name,value in (
            ('prepared',prepared),('coarse',coarse),('fixed',fixed),('seed',seed),('policy',policy))},
            producers={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in producers})
    def _get(self,prepared,coarse,fixed,seed,policy):
        import copy
        binding=self._bind(prepared,coarse,fixed,seed,policy)
        if self._components is None:
            components=_prepare_source_map_prefix(prepared,coarse,fixed,seed)
            if binding!=self._bind(prepared,coarse,fixed,seed,policy):
                raise ValueError('Source-map inputs/producers changed during preparation')
            self._components=copy.deepcopy(components)
            self._integrity=self._digest(self._components);self._bindings=binding
        elif binding!=self._bindings:
            raise ValueError('Source-map per-call context input/producer mismatch')
        if self._digest(self._components)!=self._integrity:
            raise ValueError('Source-map per-call context data mutated')
        return copy.deepcopy(self._components)
    def _assert_current(self,prepared,coarse,fixed,seed,policy):
        if self._components is None or self._bind(prepared,coarse,fixed,seed,policy)!=self._bindings:
            raise ValueError('Source-map per-call context input/producer mismatch at finish')
        if self._digest(self._components)!=self._integrity:
            raise ValueError('Source-map per-call context data mutated at finish')


def construct_source_map(prepared, coarse, fixed, result, seed, *, policy, check_output=None, _call_context=None, _proof_output=None):
    """Reconstruct the complete current-source recipe from actual common output.

    The returned request still requires exact whole-mesh and collar validation
    by the calling Julia process. A saved recipe is never used as a pass flag.
    """
    if _call_context is not None and type(_call_context) is not _SourceMapCallContext:
        raise TypeError('Only the private same-plot-call context is supported')
    import hashlib
    import marked_display_source_map as maps
    from draw_marked_cut import render
    if _call_context is None:
        components=_prepare_source_map_prefix(prepared,coarse,fixed,seed)
    else:
        components=_call_context._get(prepared,coarse,fixed,seed,policy)
    if _proof_output is not None:
        # Source: actual prefix built above. Diagnostic evidence only: never loaded.
        from pathlib import Path
        binding=_SourceMapCallContext._bind(prepared,coarse,fixed,seed,policy)
        Path(_proof_output).with_suffix('.prefix.json').write_bytes(maps._source_map_json(dict(
            diagnostic_only=True,input_producer_binding=binding,components=components)))
    native=components['native']
    mesh=components['mesh']
    collars=components['collars']
    circle_collars=components['circle_collars']
    path_binding=components['path_binding']
    domain=components['domain']
    port_map=components['port_map']
    anchors=components['anchors']
    support=components['support']
    family=components['family']
    literal=components['literal']
    result['source_map_boundary_paths']=copy.deepcopy(literal)
    seed_for_checks=copy.deepcopy(seed)
    seed_for_checks['source_map_boundary_paths']=copy.deepcopy(literal)
    # Source: disk_piece_cycles consumes the certified polygonal cut chains in
    # the new decoder branch. They must exist before ordinary D geometry checks.
    if check_output is not None:
        checked=check_output(result,seed_for_checks)
        result.clear();result.update(checked)
    replacement=maps.build_exact_anchor_beziers(result,seed_for_checks,port_map,anchors,**policy,_diagnostic_output=(str(Path(_proof_output).with_suffix('.anchor-diagnostic.json')) if _proof_output is not None else None))
    attachment=maps.verify_bezier_attachment(result,seed_for_checks,port_map,anchors,replacement,support)
    result['svg']=render(result)
    numeric=maps.verify_svg_coordinates(result)
    if _call_context is not None:_call_context._assert_current(prepared,coarse,fixed,seed,policy)
    return pack_source_map(coarse,fixed,result,mesh,collars,circle_collars,domain,
        port_map,anchors,support,family,replacement,literal,
        dict(native=native,source_paths=path_binding,attachment=attachment,numeric=numeric))


def pack_source_map(coarse,fixed,result,mesh,collars,circle_collars,domain,
                    port_map,anchors,support,family,replacement,literal,checks):
    """Serialize retained components without recomputing their geometry.

    Source: construct_source_map above. This assembly does not validate or
    authorize components; the common caller freshly executes every check.
    Keeping transport separate also permits a failed carrier check to be
    repaired against retained, explicitly scoped computational evidence.
    """
    import hashlib
    import marked_display_source_map as maps
    # Keep full exact cells reconstructible. The native source is already
    # losslessly retained on result; do not duplicate its JSON in every layer.
    source_chart={k:v for k,v in coarse.items() if k not in ('source_document','svg','orientation_svg')}
    replay={k:v for k,v in fixed.items() if k not in ('source_document','svg','orientation_svg','_source_map_input','source_map_mesh','source_chart')}
    retained_domain={k:domain[k] for k in ('source_positions','triangles','holes')}
    retained_domain['directed_seed']={k:domain['seed'].get(k) for k in ('curves','crossings','disk','disk_pieces','disk_region_marks','port_transport')}
    recipe=dict(schema='marked-display-source-map/1',
        source_representation_sha256=coarse['preparation_representation_sha256'],
        display_request=coarse['display_request'],source_chart=source_chart,
        fixed_map=dict(source=mesh,target_positions={k:fixed['positions'][k] for k in mesh['source_positions']}),
        collars=[collars,circle_collars],port_domain=retained_domain,port_map=port_map,
        exact_circle_anchors=anchors,graph_free_collars=support,fixed_quotient=family,
        exact_bezier_replacements=replacement['replacements'],attachment_boundaries=literal,
        fixed_replay_view=replay,normalization=result['source_endpoint_normalization'],
        checks=checks,
        geometry_contract='encoded polygonal cut-surface graph mapped by the retained attachment chart to the decoded-circle quotient; exact anchors are witness-only, not literal circle-glyph incidence',
        validation='awaiting same-process exact whole-mesh/collar verification; saved checks never bypass replay')
    recipe['infinity_extension']=maps.build_infinity_extension(coarse,mesh,fixed)
    geometry=source_map_geometry_payload(result)
    recipe['output_geometry_numeric_encoding']='schema-binary64/1'
    recipe['source_document_sha256']=hashlib.sha256(maps._source_map_json(result['source_document'])).hexdigest()
    recipe['output_geometry_sha256']=hashlib.sha256(maps._source_map_json(geometry)).hexdigest()
    recipe['output_svg_sha256']=hashlib.sha256(result['svg'].encode()).hexdigest()
    recipe,transport=maps.canonical_source_map_recipe(recipe)
    recipe['transport_normalization']=transport
    return recipe
