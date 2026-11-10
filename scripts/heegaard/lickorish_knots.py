"""Carry LC crossing visits through placement without changing coordinates.

cf:docs/heegaard/261109_cycle_knots.md 配置から渡す交点 — declared crossing vertices for the smoothing inlet
"""
from collections import defaultdict


def attach_crossing_knots(document, arcs, owners, map_point):
    if not any('event_indices' in arc for arc in arcs): return
    visits = defaultdict(list)
    locations = defaultdict(list)
    for name, curve in document['curves'].items():
        for index, point in enumerate(curve['points']):
            locations[(owners[name], tuple(point))].append(dict(curve=name,index=index))
    used = defaultdict(int)
    for arc in arcs:
        if 'event_indices' not in arc: continue
        for label, point in zip(arc['mids'],arc['events_model'][1:-1]):
            key = (arc['curve'],tuple(map_point(arc['face'],point)))
            visits[label].append(locations[key][used[key]])
            used[key] += 1
    if visits:
        document['crossing_knots'] = [dict(label=label,visits=visits[label]) for label in sorted(visits)]
