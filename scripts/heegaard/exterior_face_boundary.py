"""Detached source-to-boundary adapter for the OS04 region marker.

Source: draw_marked_cut.disk_piece_cycles and attachment_boundary_path.
Run only after ordinary source-map readback. No geometry is inferred here.
"""
from draw_marked_cut import attachment_boundary_path
from exterior_face_membership import exact_point
from periodic_point_transport import require_actual_interval_controls


def source_bound_cycles(view, recipe):
    marker = view["exterior_face_marker"]
    if marker != recipe["source_chart"]["exterior_face_marker"]:
        raise ValueError("exterior face disagrees with source-bound partition")
    if marker["selector"]["infinity"] != view["display_request"]["infinity_triangle"]:
        raise ValueError("exterior face infinity changed")
    indexed = {
        (curve["id"], tuple(segment["source"])): segment
        for curve in view["curves"]
        for segment in curve["segments"]
    }
    shared = view["source_endpoint_normalization"]["shared_svg_endpoints"]
    cycles = []
    for cycle in marker["display_boundary_cycles"]:
        controls = []
        for part in cycle:
            if part["kind"] == "curve":
                segment = indexed[(part["curve"], tuple(part["source"]))]
                cubics = require_actual_interval_controls(segment, shared)
                direction = part["direction"]
                if direction == -1:
                    cubics = [list(reversed(c)) for c in reversed(cubics)]
                elif direction != 1:
                    raise ValueError("bad directed source dart")
                controls.extend(cubics)
            elif part["kind"] == "cut":

                path = [exact_point(q) for q in attachment_boundary_path(view, part)]
                controls.extend(zip(path, path[1:]))
            else:
                raise ValueError("unknown exterior face boundary dart")
        if not controls or any(
            first[-1] != second[0]
            for first, second in zip(controls, controls[1:] + controls[:1])
        ):
            raise ValueError("decoded exterior boundary does not close exactly")
        cycles.append(controls)
    return cycles
