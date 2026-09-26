"""Bind the supplied front/back PL witness to the exact marked document.

This is an audited source projection, not the certified rational disk atlas.
No expected o-data participates in construction. Original vendor files are read-only.
"""
import base64
import hashlib
import html
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'vendor/koda_genus3'
sys.path.insert(0, str(SOURCE / 'src'))
import build_case as B
import marked_reference as M
import trace_input as T
import verify_case as V


def build():
    doc = json.loads((SOURCE / 'data/figure6.heegaard.marked.json').read_text())
    full = json.loads((SOURCE / 'data/figure6.fullmarked.json').read_text())
    cycles, pairs, signs, sheets, witnesses, segments = B.input_data()
    rotations, twin, faces = B.ribbon(cycles, pairs, signs)
    reconstructed = B.make_document(cycles, pairs, signs, rotations, twin, faces, sheets)
    assert M.representation_digest(reconstructed) == M.representation_digest(full)
    reconstructed['curves'] = [c for c in reconstructed['curves'] if c['role'] != 'aux']
    reconstructed['crossings'] = [x for x in reconstructed['crossings'] if x['id'].startswith(('c', 'x'))]
    digest = M.representation_digest(doc)
    assert M.representation_digest(reconstructed) == digest
    audit = V.pixel_intersection_audit()
    source_bytes = (SOURCE / 'data/source.png').read_bytes()
    assert hashlib.sha256(source_bytes).hexdigest() == doc['metadata']['source_sha256']
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1700 860" role="img" aria-label="種数3の前面背面Heegaard投影とjoining disk">',
           '<style>.trace{fill:none;stroke-width:5;stroke-linejoin:round;stroke-linecap:round}.alpha{stroke:#b63837}.beta{stroke:#275db2}.trace[data-sheet="B"]{stroke-dasharray:12 10}.dim{opacity:.12}.point-label{font:22px sans-serif;fill:#152435;stroke:white;stroke-width:6;paint-order:stroke}.source-image{opacity:.24}</style>',
           '<defs><marker id="trace-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="3" markerHeight="3" orient="auto"><path d="M0 0 L10 5 L0 10z" fill="context-stroke"/></marker></defs>',
           '<rect width="1700" height="860" fill="white"/>',
           '<image class="source-image" width="1700" height="840" href="data:image/png;base64,' + base64.b64encode(source_bytes).decode() + '"/>']
    # Trace D along the declared source curves, retaining all intermediate PL vertices.
    boundary = []
    for cid, start, end in T.D_RIGHT:
        path = [v for _, vs in T.PIECES[cid] for v in vs[:-1]]
        i = path.index(start)
        while True:
            boundary.append(B.point(path[i]))
            i = (i + 1) % len(path)
            if path[i] == end:
                break
    pts = lambda vs: ' '.join(f'{x},{y}' for x, y in vs)
    svg.append(f'<polygon class="region" points="{pts(boundary)}" fill="#efd27a" fill-opacity=".65"/>')
    svg.append('<text class="region" x="555" y="574" font-size="38" fill="#775b13">D</text>')
    for cid, pieces in T.PIECES.items():
        for sheet, vs in pieces:
            role = 'alpha' if cid.startswith('alpha') else 'beta'
            svg.append(f'<polyline class="trace {role}" data-curve="{cid}" data-sheet="{sheet}" points="{pts(map(B.point,vs))}" marker-mid="url(#trace-arrow)"/>')
    offsets = {'c1':(-48,24),'c2':(4,-16),'c3':(18,8),'c4':(14,30),'c5':(13,28),'c6':(-47,28),
               'x1':(-52,5),'x2':(-145,5),'x3':(18,-8),'x4':(-65,-18),'x5':(10,-18)}
    for x in doc['crossings']:
        name = x['id']; px, py = T.POINTS[name]; dx, dy = offsets[name]
        label = name + (' (back)' if name == 'x2' else '')
        svg.append(f'<g data-point="{name}"><circle cx="{px}" cy="{py}" r="5" fill="white" stroke="#333" stroke-width="2"/><text x="{px+dx}" y="{py+dy}" class="point-label">{html.escape(label)}</text></g>')
    svg.append('</svg>')
    result = {'representation_sha256': digest, 'svg': ''.join(svg), 'audit': audit,
              'source_sha256': hashlib.sha256(source_bytes).hexdigest(),
              'scope': 'Recorded PL trace checked on separate front/back sheets and reconstructed to the identical marked cellulation. Raster underlay is source material, not a certified layout.'}
    out = ROOT / 'artifacts/heegaard/cellulated/projections'
    out.mkdir(parents=True, exist_ok=True)
    (out / (digest + '.json')).write_text(json.dumps(result, ensure_ascii=False), encoding='utf-8')
    print(json.dumps({'projection': digest, 'trace_audit': audit}))


if __name__ == '__main__':
    build()
