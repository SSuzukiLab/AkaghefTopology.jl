"""Regression for production fixed-frame outer support probe selection."""

from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marked_display_source_map as maps


def _orient(x, y, z):
    return maps._planar_orient(x, y, z)


def _triangle():
    outer = (0, 1, 2)
    chart = {0: (0.0, 0.0), 1: (400.0, 0.0), 2: (200.0, 346.0)}
    return chart, outer


def test_production_probes_accept_aabb_protrudes():
    chart, outer = _triangle()
    controls = [(120.0, 80.0), (280.0, 80.0), (200.0, 250.0)]
    probes = maps.fixed_frame_outer_support_probe_points(controls)
    aabb = maps.aabb_corner_points(controls)
    assert probes == maps.convex_hull_points(controls)
    assert aabb != probes
    assert any(_orient(chart[u], chart[v], q) <= 0 for u, v in zip(outer, outer[1:] + outer[:1]) for q in aabb)
    maps._fixed_attachment_outer_support_probes(chart, outer, probes, _orient)
    try:
        maps._fixed_attachment_outer_support_probes(chart, outer, aabb, _orient)
    except AssertionError:
        return
    raise AssertionError('AABB corners must fail when protruding from convex outer')


def test_production_probes_reject_outside_control():
    chart, outer = _triangle()
    controls = [(100.0, 50.0), (400.0, 50.0)]
    probes = maps.fixed_frame_outer_support_probe_points(controls)
    try:
        maps._fixed_attachment_outer_support_probes(chart, outer, probes, _orient)
    except AssertionError:
        return
    raise AssertionError('production probes must reject outside control')


if __name__ == '__main__':
    test_production_probes_accept_aabb_protrudes()
    test_production_probes_reject_outside_control()
    print('OK attachment chart support regressions')
