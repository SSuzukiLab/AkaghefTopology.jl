"""Executable post-cell-check proposal called FROM ordinary Julia H.plot.

Source: existing draw_marked_cut.render, maps.canonical_source_map_recipe,
pipeline.source_map_geometry_payload. Renderer requires the narrow hook in
HOOKS.md; refusing its absence prevents success using an opaque sidecar.
"""
import argparse
import hashlib
import json
from pathlib import Path
from periodic_annotations import bind, verify_numeric, need


def finish(view, recipe, compiled, diagnostic_dir=None):
    import marked_display_source_map as maps
    from marked_display_source_map_pipeline import source_map_geometry_payload
    from draw_marked_cut import render
    bind(view, compiled, recipe, diagnostic_dir=diagnostic_dir)
    view['svg'] = render(view)
    check = verify_numeric(view)  # fails unless ordinary render consumed data
    numeric = maps.verify_svg_coordinates(view)  # retain all original numeric gates too
    recipe['periodic_annotations'] = view['periodic_annotations']
    recipe['checks']['periodic_annotations'] = check
    recipe['checks']['numeric'] = numeric
    geometry = source_map_geometry_payload(view)
    need(geometry.get('periodic_annotations') == view['periodic_annotations'],
         'ordinary geometry digest does not consume annotations')
    recipe['output_geometry_sha256'] = hashlib.sha256(maps._source_map_json(geometry)).hexdigest()
    recipe['output_svg_sha256'] = hashlib.sha256(view['svg'].encode()).hexdigest()
    recipe, transport = maps.canonical_source_map_recipe(recipe)
    recipe['transport_normalization'] = transport
    view['source_map'] = maps.encode_source_map(recipe)
    return view, recipe


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('view', 'recipe', 'compiled', 'output', 'proof'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    need(not a.output.exists() and not a.proof.exists(), 'refusing overwrite')
    view, recipe = finish(*(json.loads(getattr(a, name).read_text(encoding='utf-8'))
                             for name in ('view', 'recipe', 'compiled')), diagnostic_dir=a.output.parent)
    a.output.write_bytes(__import__('marked_display_source_map')._source_map_json(view))
    a.proof.write_bytes(__import__('marked_display_source_map')._source_map_json(recipe))
