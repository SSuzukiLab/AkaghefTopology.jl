"""Exact decoded point membership in the OS04 exterior face.

Sources: draw_marked_cut.attachment_boundary_path supplies directed bank
intervals; periodic_point_transport supplies exact rational cubic splitting
and hull separation. Input cycles must already be bound to native face darts.
"""
from fractions import Fraction

from periodic_point_transport import hull, rat, restrict_cubic, separation


def exact_point(point):
    if len(point) != 2:
        raise ValueError("expected a planar point")
    return tuple(map(rat, point))


def cubic_chords_avoiding_point(controls, point, depth=0):
    outside_box = any(
        point[k] < min(p[k] for p in controls) or point[k] > max(p[k] for p in controls)
        for k in range(2)
    )
    if outside_box or separation(hull(controls), [point]) is not None:
        return [(controls[0], controls[-1])]
    if depth == 20:
        raise ValueError("strict separation from actual face cubic unproved")
    half = Fraction(1, 2)
    left = restrict_cubic(controls, Fraction(0), half)
    right = restrict_cubic(controls, half, Fraction(1))
    return (cubic_chords_avoiding_point(left, point, depth + 1)
            + cubic_chords_avoiding_point(right, point, depth + 1))


def winding(point, edges):
    number = 0
    if not edges or any(v != a for (_, v), (a, _) in zip(edges, edges[1:] + edges[:1])):
        raise ValueError("directed face boundary is not closed")
    px, py = point
    for (ax, ay), (bx, by) in edges:
        det = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
        if det == 0 and min(ax, bx) <= px <= max(ax, bx) and min(ay, by) <= py <= max(ay, by):
            raise ValueError("point lies on directed face boundary")
        if ay <= py < by and det > 0:
            number += 1
        if by <= py < ay and det < 0:
            number -= 1
    return number


def certify_exterior_marker(candidate, face_cycles, enclosure, svg_holes):
    """Verify a decoded interior representative of the selected exterior face.

    Each cycle item is a directed line segment or cubic control quadruple in
    positive-y plot units. All bank lines are the certified literal attachment
    polygon, while svg_holes independently supplies actual decoded circles.
    Enclosure is a finite oriented polygon in those same units.
    """
    if not face_cycles or len(enclosure) < 3:
        raise ValueError("complete face cycles and enclosure required")
    point = exact_point(candidate)
    numbers = []
    for cycle in face_cycles:
        edges = []
        for controls in cycle:
            q = [exact_point(p) for p in controls]
            if len(q) == 2:
                edges.append(tuple(q))
            elif len(q) == 4:
                edges.extend(cubic_chords_avoiding_point(q, point))
            else:
                raise ValueError("expected directed line or cubic")
        numbers.append(winding(point, edges))
    if sum(numbers) != 0:
        raise ValueError("point is outside the source-selected exterior face")
    fence = list(map(exact_point, enclosure))
    if abs(winding(point, list(zip(fence, fence[1:] + fence[:1])))) != 1:
        raise ValueError("point outside the finite exterior enclosure")
    for center, radius in svg_holes:
        center = exact_point(center)
        radius = rat(radius)
        if radius <= 0 or sum((a - b) ** 2 for a, b in zip(point, center)) <= radius ** 2:
            raise ValueError("point is in or on an actual SVG hole")
    return {"point": [str(q) for q in point], "boundary_windings": numbers}
