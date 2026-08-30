"""Deterministic corridor realization for the Poincare paired-puncture fixture.

This is a fixture generator, not a generic layout algorithm.  It realizes the
certified terminal DCEL by routing a planar tetrahedral skeleton through fixed
corridors and matching the prescribed port order inside disjoint annuli.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


CENTERS = {1: (1.0, 1.0), 2: (1.0, -1.0), 3: (2.0, 1.0), 4: (2.0, -1.0)}
RADIUS = 0.28
COLLAR = 0.42

# (start port, end port), in beta-successor order.
EDGES = [
    (14, 2), (13, 3), (12, 4), (11, 24), (15, 9), (6, 22),
    (17, 1), (7, 18), (21, 19), (20, 10), (5, 23), (16, 8),
]

PORT_TERMINAL = {
    **{port: 1 for port in range(1, 8)},
    **{port: 2 for port in range(8, 15)},
    **{port: 3 for port in range(15, 20)},
    **{port: 4 for port in range(20, 25)},
}
PORT_ANGLE = {
    **{port: (2 * (port - 1) + 1) * math.pi / 7 for port in range(1, 8)},
    **{port: (2 * (port - 8) + 1) * math.pi / 7 for port in range(8, 15)},
    **{port: (2 * (port - 15) + 1) * math.pi / 5 for port in range(15, 20)},
    **{port: (2 * (port - 20) + 1) * math.pi / 5 for port in range(20, 25)},
}

# Collar slot angles preserve each terminal's cyclic port order.  Values are
# deliberately unwrapped; only sin/cos sees them modulo 2*pi.
SLOT_ANGLE = {
    1: math.radians(20), 2: math.radians(120), 3: math.radians(160),
    4: math.radians(200), 5: math.radians(250), 6: math.radians(300),
    7: math.radians(340),
    8: math.radians(-120), 9: math.radians(-80), 10: math.radians(-40),
    11: math.radians(0), 12: math.radians(80), 13: math.radians(120),
    14: math.radians(160),
    15: math.radians(-60), 16: math.radians(-20), 17: math.radians(180),
    18: math.radians(200), 19: math.radians(270),
    20: math.radians(-80), 21: math.radians(-20), 22: math.radians(60),
    23: math.radians(100), 24: math.radians(180),
}

# Corridor points between the two collar slots.  The face structure is:
# outer face (1,12,7), with terminal 4 inside; duplicate edge pairs occupy
# disjoint narrow bands between the four triangular faces.
CORRIDORS = {
    1: [(0.45, -1.25), (0.20, -0.55), (0.20, 0.55), (0.45, 1.35)],
    2: [(0.54, -0.55), (0.48, 0.00), (0.54, 0.55)],
    3: [(0.66, -0.45), (0.62, 0.00), (0.66, 0.45)],
    4: [(1.38, -0.78), (1.62, -0.78)],
    5: [(2.58, 0.55), (2.72, -0.35), (2.65, -1.52), (1.42, -1.55)],
    6: [(1.36, 0.58), (1.58, 0.10), (1.78, -0.42)],
    7: [(1.65, 1.52), (1.35, 1.52)],
    8: [(1.40, 0.82), (1.60, 0.82)],
    9: [(2.44, -1.16), (2.44, 0.55)],
    10: [(1.62, -1.42), (1.38, -1.42)],
    11: [(1.25, 0.55), (1.43, 0.02), (1.67, -0.58)],
    12: [(2.82, 0.72), (2.95, -0.20), (2.95, -1.92), (0.72, -1.92)],
}


def point(center: tuple[float, float], radius: float, angle: float) -> tuple[float, float]:
    return center[0] + radius * math.cos(angle), center[1] + radius * math.sin(angle)


def port_point(port: int) -> tuple[float, float]:
    return point(CENTERS[PORT_TERMINAL[port]], RADIUS, PORT_ANGLE[port])


def collar_path(port: int, reverse: bool = False, samples: int = 12) -> list[tuple[float, float]]:
    center = CENTERS[PORT_TERMINAL[port]]
    start, finish = PORT_ANGLE[port], SLOT_ANGLE[port]
    # Choose the unwrapping compatible with the terminal-wide order above.
    while finish - start > math.pi:
        finish -= 2 * math.pi
    while finish - start < -math.pi:
        finish += 2 * math.pi
    path = []
    for index in range(samples + 1):
        fraction = index / samples
        radius = RADIUS + (COLLAR - RADIUS) * fraction
        angle = start + (finish - start) * fraction
        path.append(point(center, radius, angle))
    return list(reversed(path)) if reverse else path


def routes() -> list[list[tuple[float, float]]]:
    result = []
    for edge_id, (start, finish) in enumerate(EDGES, 1):
        first = collar_path(start)
        last = collar_path(finish, reverse=True)
        result.append(first + CORRIDORS[edge_id] + last)
    return result


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a, b, p):
    return orient(a, b, p) == 0.0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def intersects(a, b, c, d):
    if max(min(a[0], b[0]), min(c[0], d[0])) > min(max(a[0], b[0]), max(c[0], d[0])):
        return False
    if max(min(a[1], b[1]), min(c[1], d[1])) > min(max(a[1], b[1]), max(c[1], d[1])):
        return False
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 * o2 < 0.0 and o3 * o4 < 0.0:
        return True
    return on_segment(a, b, c) or on_segment(a, b, d) or on_segment(c, d, a) or on_segment(c, d, b)


def segment_distance(center, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    denominator = dx * dx + dy * dy
    if denominator == 0.0:
        return math.hypot(center[0] - a[0], center[1] - a[1])
    parameter = max(0.0, min(1.0, ((center[0] - a[0]) * dx + (center[1] - a[1]) * dy) / denominator))
    return math.hypot(center[0] - (a[0] + parameter * dx), center[1] - (a[1] + parameter * dy))


def audit(data):
    crossings = []
    self_crossings = []
    intrusions = []
    for first in range(len(data)):
        for second in range(first + 1, len(data)):
            for i, (a, b) in enumerate(zip(data[first], data[first][1:])):
                if any(intersects(a, b, c, d) for c, d in zip(data[second], data[second][1:])):
                    crossings.append((first + 1, second + 1, i + 1))
                    break
        for i, (a, b) in enumerate(zip(data[first], data[first][1:])):
            for j in range(i + 2, len(data[first]) - 1):
                if intersects(a, b, data[first][j], data[first][j + 1]):
                    self_crossings.append((first + 1, i + 1, j + 1))
            for terminal, center in CENTERS.items():
                if segment_distance(center, a, b) < RADIUS - 1e-12:
                    intrusions.append((first + 1, i + 1, terminal))
    return crossings, self_crossings, intrusions


def main():
    data = routes()
    crossings, self_crossings, intrusions = audit(data)
    print("crossings", crossings)
    print("self_crossings", self_crossings)
    print("intrusions", intrusions)
    output = Path(__file__).with_name("manual_routes.json")
    output.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    raise SystemExit(1 if crossings or self_crossings or intrusions else 0)


if __name__ == "__main__":
    main()
