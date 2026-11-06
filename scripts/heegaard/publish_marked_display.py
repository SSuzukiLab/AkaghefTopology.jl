"""Bind a saved geometric realization to the independently prepared Julia complex.

This does not change the source marking or call an endpoint router. Publication
fails on a stale source, a changed curve edge, hole side, crossing, or D boundary.
Reusing this realization is explicitly recorded, not called fresh generation.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from shapely.geometry import LineString, Polygon
from draw_marked_cut import render
from verify_fixed_grid import verify


def bind(d, p, base):
    assert p['schema']=='marked-display-complex/1'
    assert p['source_file_sha256']==d['source_sha256']==hashlib.sha256(Path(d['source']).read_bytes()).hexdigest()
    assert d['source_cut_corner_classes']==p['source_cut_corner_classes']
    assert d['genus']==p['genus'] and d['cuts']==p['cut_curves']
    assert len(d['curves'])==len(p['curves'])
    for c,expected in zip(d['curves'],p['curves']):
        assert (c['id'],c['role'])==(expected['id'],expected['role'])
        assert [s['source'] for s in c['segments']]==[s['source'] for s in expected['segments']]
    assert [{k:x[k] for k in ('id','at','curves','sign')} for x in d['crossings']]==p['crossings']
    d=copy.deepcopy(d);d['grid']=p['grid']
    for b in d['boundaries']:
        occurrences={tuple(a) for a in b['source']}
        expected=next(x for x in p['boundaries'] if {tuple(a) for a in x['source']}==occurrences)
        assert b['cut']==expected['cut']
        b.update({k:expected[k] for k in ('side','row','column')})
        center=np.array(d['positions'][str(b['cap'])])*600
        target=[390*(b['column']-(d['genus']-1)/2),250*(1-2*b['row'])]
        assert np.max(abs(center-target))<1e-6, 'Saved geometry has different marked hole sides'
    assert len(p['region_marks'])==1 and len(p['region_marks'][0]['cut_pieces'])==1
    # Reconstruct D's entire boundary from the SOURCE triangle darts, not the
    # saved yellow polygon. Compare their realised curves after SVG rounding.
    positions={int(v):np.round(np.array(x)*600,7) for v,x in d['positions'].items()}
    paths={tuple(s['source']):s['vertices'] for c in d['curves'] for s in c['segments']}
    boundary=[]
    for f,_,u,v in p['region_marks'][0]['source_region']['boundary_cycles'][0]:
        vs=paths.get((f,u,v))
        if vs is None: vs=paths[(f,v,u)][::-1]
        points=[positions[w] for w in vs]
        if boundary: assert np.linalg.norm(boundary[-1]-points[0])<1e-5
        boundary.extend(points if not boundary else points[1:])
    expected=LineString(boundary)
    actual=Polygon([positions[v] for v in d['disk']]).boundary
    assert expected.hausdorff_distance(actual)<1e-5, 'D is not the original marked region'
    d['geometry_check']=verify(d,base)
    d['prepared_check']=dict(source_corner_quotient=True,curve_edge_words=True,hole_sides=True,
        crossing_labels_and_signs=True,joining_disk_boundary=True)
    d['preparation_representation_sha256']=p['source_representation_sha256']
    d['marking_id']=p['marking_id']
    d['geometry_provenance']='saved fixed-grid PL realization, independently revalidated against the Julia-prepared complex'
    d['source_document']=json.loads(Path(d['source']).read_text(encoding='utf-8'))
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path);args=parser.parse_args();directory=args.directory
    read=lambda name:json.loads((directory/name).read_text(encoding='utf-8'))
    p=read('KodaPrepared.json');d=bind(read('KodaGrid.json'),p,read('KodaMarkedCut.json'))
    d['preparation_sha256']=hashlib.sha256((directory/'KodaPrepared.json').read_bytes()).hexdigest()
    (directory/'KodaDisplay.json').write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8')
    (directory/'KodaDisplay.svg').write_text(render(d),encoding='utf-8')
    print(json.dumps(dict(prepared=d['prepared_check'],geometry=d['geometry_check']),ensure_ascii=False))
