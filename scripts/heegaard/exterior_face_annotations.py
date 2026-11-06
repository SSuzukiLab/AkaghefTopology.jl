"""Standalone exterior-face glyph, using the existing ordinary SVG renderer."""
import hashlib
from html import escape
import xml.etree.ElementTree as ET

from exterior_face_boundary import source_bound_cycles
from exterior_face_membership import certify_exterior_marker
from periodic_point_transport import rat


def layer(view):
    record = view.get("exterior_face_glyph")
    if record is None:
        return ""
    x, y = map(lambda q: float(rat(q)), record["point"])
    name = escape(record["name"])
    return (f'<circle data-exterior-face-glyph="point" cx="{x!r}" cy="{-y!r}" r="3" fill="#182538"/>'
            f'<text data-exterior-face-glyph="name" class="label" x="{x!r}" y="{-y!r}">{name}</text>')


def viewport_and_holes(svg):
    root = ET.fromstring(svg)
    x, top, width, height = map(rat, root.attrib["viewBox"].split())
    fence = [(x, -top), (x + width, -top), (x + width, -top - height), (x, -top - height)]
    holes = [
        ((float(e.attrib["cx"]), -float(e.attrib["cy"])), float(e.attrib["r"]))
        for e in root.iter() if "data-puncture" in e.attrib
    ]
    return root, fence, holes


def verify(view, recipe):
    record = view["exterior_face_glyph"]
    name = view["exterior_face_marker"]["selector"]["name"]
    if record["name"] != name:
        raise ValueError("exterior glyph source name changed")
    root, fence, holes = viewport_and_holes(view["svg"])
    items = [e for e in root.iter() if "data-exterior-face-glyph" in e.attrib]
    if [e.attrib["data-exterior-face-glyph"] for e in items] != ["point", "name"]:
        raise ValueError("ordinary renderer must emit one exterior glyph and name")
    dot, text = items
    svg_namespace = "{http://www.w3.org/2000/svg}"
    if (dot.tag, text.tag) != (svg_namespace + "circle", svg_namespace + "text"):
        raise ValueError("wrong exterior SVG glyph element types")
    if any("transform" in element.attrib for element in root.iter()):
        raise ValueError("transformed exterior SVG is outside the coordinate contract")
    point = [rat(float(dot.attrib["cx"])), -rat(float(dot.attrib["cy"]))]
    label = [rat(float(text.attrib["x"])), -rat(float(text.attrib["y"]))]
    if point != label or point != list(map(rat, record["point"])) or text.text != name:
        raise ValueError("exterior glyph SVG readback changed")
    check = certify_exterior_marker(point, source_bound_cycles(view, recipe), fence, holes)
    if check["boundary_windings"] != record["boundary_windings"]:
        raise ValueError("exterior winding readback changed")
    return check


def finish(view, recipe):
    import marked_display_source_map as maps
    from marked_display_source_map_pipeline import source_map_geometry_payload
    from draw_marked_cut import render

    _, fence, holes = viewport_and_holes(view["svg"])
    cycles = source_bound_cycles(view, recipe)
    name = view["exterior_face_marker"]["selector"]["name"]
    left, top = fence[0]
    right, bottom = fence[2]
    # Bounded deterministic viewport search. These are proposals, independent
    # of the selected paper. Exact membership is the sole geometric judge.
    for n in (8, 16, 32):
        for i in range(1, n):
            for j in range(1, n):
                candidate = (float(left + (right - left) * i / n),
                             float(bottom + (top - bottom) * j / n))
                try:
                    result = certify_exterior_marker(candidate, cycles, fence, holes)
                except ValueError:
                    continue
                view["exterior_face_glyph"] = {"name": name, **result}
                break
            else:
                continue
            break
        else:
            continue
        break
    else:
        raise ValueError("no certified exterior glyph position in bounded viewport search")
    view["svg"] = render(view)
    recipe["checks"]["exterior_face_glyph"] = verify(view, recipe)
    recipe["checks"]["numeric"] = maps.verify_svg_coordinates(view)
    recipe["exterior_face_glyph"] = view["exterior_face_glyph"]
    recipe["output_geometry_sha256"] = hashlib.sha256(
        maps._source_map_json(source_map_geometry_payload(view))).hexdigest()
    recipe["output_svg_sha256"] = hashlib.sha256(view["svg"].encode()).hexdigest()
    recipe, transport = maps.canonical_source_map_recipe(recipe)
    recipe["transport_normalization"] = transport
    view["source_map"] = maps.encode_source_map(recipe)
    return view, recipe


if __name__ == "__main__":
    import json
    from pathlib import Path
    import sys
    from marked_display_source_map import _source_map_json

    view_path, recipe_path, output, proof = map(Path, sys.argv[1:])
    view, recipe = finish(json.loads(view_path.read_bytes()), json.loads(recipe_path.read_bytes()))
    for path, result in ((output, view), (proof, recipe)):
        with path.open("xb") as stream:
            stream.write(_source_map_json(result))
