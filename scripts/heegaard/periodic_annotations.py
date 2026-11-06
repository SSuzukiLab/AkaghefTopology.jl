"""Narrow executable ordinary-render hook for the source-bound annulus case.

Source: draw_marked_cut.render's closed cubic path encoder and compound D
cycles; native validation/finite points come from periodic_annotations.proposed.jl.
No point placement, substitute geometry, source mutation, or new map engine.
"""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from html import escape
import json
import xml.etree.ElementTree as ET


def need(ok, message):
    if not ok:
        raise ValueError('periodic annotations: ' + message)


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(',', ':')).encode()


def bind(view, compiled, recipe, diagnostic_dir=None):
    """Consume same-process compiler output before render and final pack hashes."""
    need(compiled['schema'] == 'marked-periodic-annotations/1', 'unknown schema')
    need(compiled['source_document_sha256'] == recipe['source_document_sha256'],
         'wrong source recipe')
    need(compiled['source_representation_sha256'] == recipe['source_representation_sha256'],
         'wrong representation')
    # Same canonical JSON encoder as existing source-map hashing.
    import marked_display_source_map as maps
    need(sha256(maps._source_map_json(view['source_document'])).hexdigest() ==
         compiled['source_document_sha256'], 'source changed')
    need(compiled['boundary_checked'] is True, 'missing same-process chain check')
    result = deepcopy(compiled)
    from periodic_moving_points import transport_all_points
    transport=transport_all_points(view,recipe,compiled,diagnostic_dir=diagnostic_dir)
    need(len(transport['points'])==len(result['points'])==5,'missing moved endpoints')
    for p,endpoint in zip(result['points'],transport['points']):
        need(p['point_svg_exact']==endpoint['point_svg_exact'],'source qA changed')
        p.update({k:v for k,v in endpoint.items() if k!='point_index'})
    result['point_transport']=transport
    view['periodic_annotations'] = result
    return result


def closed_curve_pieces(view, cid):
    curves = [c for c in view['curves'] if c['id'] == cid]
    need(len(curves) == 1 and curves[0]['closed'], 'missing actual closed curve')
    pieces = [q for s in curves[0]['segments'] for q in s.get('beziers', [])]
    need(pieces and all(len(q) == 4 for q in pieces), 'missing actual cubic geometry')
    need(all(a[-1] == b[0] for a, b in zip(pieces, pieces[1:] + pieces[:1])),
         'actual closed cubic curve has a gap')
    return pieces


def layer(view):
    """Inserted by draw_marked_cut.render before curves, using its SVG units."""
    annotation = view.get('periodic_annotations')
    if annotation is None:
        return ''
    def xy(p):
        return repr(float(p[0])) + ',' + repr(-float(p[1]))
    def path(pieces):
        return 'M' + xy(pieces[0][0]) + ' ' + ' '.join(
            'C' + ' '.join(xy(p) for p in q[1:]) for q in pieces) + ' Z'
    out = ['<g data-periodic-annotations="P">']
    for i, region in enumerate(annotation['regions'], 1):
        n = region['coefficient']
        need(type(n) is int, 'coefficient must be an integer')
        if n == 0:
            continue
        cycles = region['boundary_curve_cycles']
        need(len(cycles) == 2, 'first case supports annulus fills only')
        compound = ' '.join(path(closed_curve_pieces(view, c)) for c in cycles)
        color = '#edbb58' if n > 0 else '#82b5da'
        out.append(f'<path data-periodic-region="{i}" d="{compound}" '
                   f'fill="{color}" fill-opacity=".28" fill-rule="evenodd" stroke="none"/>')
    for i, p in enumerate(annotation['points']):
        x, y = p['point_svg']; y = -y
        name = str(p['name'])
        if p['kind'] == 'coefficient':
            name = f"{p['coefficient']:+d}" if p['coefficient'] else '0'
        else:
            need(name in ('w', 'z'), 'unknown basepoint')
            out.append(f'<circle data-periodic-point="{i}" cx="{x!r}" cy="{y!r}" '
                       'r="3" fill="#182538"/>')
        # Text anchor is the actual point: no silent geometric displacement.
        out.append(f'<text data-periodic-label="{i}" class="label" x="{x!r}" '
                   f'y="{y!r}" text-anchor="middle" font-size="18">{escape(name)}</text>')
    out.append('</g>')
    return ''.join(out)


def verify_numeric(view):
    """Recover every new native-point image from the actual emitted SVG."""
    annotation = view.get('periodic_annotations')
    if annotation is None:
        return None
    tree = ET.fromstring(view['svg'])
    labels = {int(e.attrib['data-periodic-label']): e for e in tree.iter()
              if 'data-periodic-label' in e.attrib}
    glyphs = {int(e.attrib['data-periodic-point']): e for e in tree.iter()
              if 'data-periodic-point' in e.attrib}
    need(set(labels) == set(range(len(annotation['points']))), 'missing/extra labels')
    need(set(glyphs) == {i for i, p in enumerate(annotation['points'])
                         if p['kind'] == 'basepoint'}, 'missing/extra basepoint glyphs')
    for i, p in enumerate(annotation['points']):
        x, y = p['point_svg']
        from periodic_point_transport import rat
        exact = tuple(map(rat, p['point_svg_moved_exact']))
        decoded = tuple(map(rat, p['point_svg_decoded']))
        need(tuple(rat(float(q)) for q in exact) == decoded == (rat(x),rat(y)),
             'transported endpoint/decoded binary-rational coordinates disagree')
        # Exact binary-rational check against actual decoded SVG circle glyphs.
        for circle in tree.iter():
            if 'data-puncture' not in circle.attrib:
                continue
            center=(rat(float(circle.attrib['cx'])),-rat(float(circle.attrib['cy'])))
            radius=rat(float(circle.attrib['r']))
            need(radius>0 and sum((a-b)**2 for a,b in zip(decoded,center))>radius**2,
                 'decoded marked point lies in/on actual SVG puncture circle')
        need((float(labels[i].attrib['x']), -float(labels[i].attrib['y'])) == (x, y),
             'label coordinate differs from mapped native point')
        name = str(p['name']) if p['kind'] == 'basepoint' else (
            f"{p['coefficient']:+d}" if p['coefficient'] else '0')
        need(labels[i].text == name, 'label content changed')
        if i in glyphs:
            need((float(glyphs[i].attrib['cx']), -float(glyphs[i].attrib['cy'])) == (x, y),
                 'basepoint glyph differs from mapped native point')
            need((rat(float(glyphs[i].attrib['cx'])), -rat(float(glyphs[i].attrib['cy']))) == decoded,
                 'basepoint glyph binary-rational point changed')
    paths = [e for e in tree.iter() if 'data-periodic-region' in e.attrib]
    need(len(paths) == sum(r['coefficient'] != 0 for r in annotation['regions']),
         'missing/extra annulus fill')
    # Re-evaluate whole layer from the bound actual cubics; all controls and
    # cycle closure must match, not merely a label or region count.
    expected = ET.fromstring('<svg>' + layer(view) + '</svg>')
    expected_paths = [e.attrib for e in expected.iter() if 'data-periodic-region' in e.attrib]
    need([e.attrib for e in paths] == expected_paths, 'fill boundary differs from actual cubics')
    return dict(points=len(labels), basepoints=len(glyphs), filled_annuli=len(paths),
                svg_sha256=sha256(view['svg'].encode()).hexdigest())


def geometry_payload(annotation):
    """Append to existing source_map_geometry_payload; exact records preserved."""
    if annotation is None:
        return None
    result = deepcopy(annotation)
    for p in result['points']:
        p['point_svg'] = [float(x) for x in p['point_svg']]
    return result
