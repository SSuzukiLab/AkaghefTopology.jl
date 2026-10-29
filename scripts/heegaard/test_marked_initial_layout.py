"""Focused whole-cellulation initial-layout and shared-fairing checks.

Source: test/heegaard/marked_initial_layout.jl emits native prepared inputs.
These are finite numerical checks, not a completeness theorem for the solver.
"""
import argparse
import copy
import json
from pathlib import Path

from layout_marked_document import layout
from marked_display_pipeline import plot


def check_case(prepared):
    before=copy.deepcopy(prepared)
    seed=layout(prepared)
    assert prepared==before,'Initial layout mutated mathematical source'
    assert seed['source_document']==prepared['source_document']
    assert seed['infinity_source_triangle']==prepared['outer_triangle']
    assert seed['initial_embedding_check']['selected_exterior_retained']
    assert seed['geometry_check']['extra_intersections']==0
    pieces=seed.get('disk_pieces',[])
    if pieces:
        expected={tuple(t) for mark in prepared['region_marks'] for t in mark['source_region']['triangles']}
        assert {tuple(t) for piece in pieces for t in piece['source_triangles']}==expected
        assert sum(piece['contains_infinity'] for piece in pieces)==int(tuple(prepared['outer_triangle']) in expected)
    result,_=plot(seed)
    assert result.get('disk_pieces')==seed.get('disk_pieces')
    assert result['source_document']==prepared['source_document']
    assert result['infinity_source_triangle']==prepared['outer_triangle']
    assert result['smooth_check']['source_disk_boundary_unchanged']
    assert result['layout_isotopy_check']['collisions']==0
    assert '<svg' in result['svg']
    # Strict JSON is the actual Julia/server interchange, including empty cases.
    json.dumps(result,allow_nan=False)
    return result,seed


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result,seed=check_case(json.loads(args.prepared.read_text(encoding='utf-8')))
    args.output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    args.output.with_suffix('.seed.json').write_text(json.dumps(seed,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    args.output.with_suffix('.svg').write_text(result['svg'],encoding='utf-8')
    print(json.dumps(dict(initial=result['initial_embedding_check'],geometry=result['geometry_check'],smooth=result['smooth_check']),allow_nan=False),flush=True)
