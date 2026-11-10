"""Combinatorial subdivision and one harmonic solve for an LC cut disk.

cf:docs/heegaard/261109_tutte_construction.md 組合せ的な細分 — incidence subdivision and fixed-boundary solve; reconstructed
"""
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

DEBUG_MODE = False


def triangulate_face(arcs, words, corners, *, exterior=False):
    """Keep edge occurrences, including loops and parallel edges, distinct."""
    fixed = []
    def vertex(point=None):
        v = len(fixed)
        fixed.append(tuple(point) if point is not None else None)
        return v
    endpoints = [[vertex(a['mP']), vertex(a['mQ'])] for a in arcs]
    boundary = [v for pair in endpoints for v in pair]
    for p in corners:
        if not any(math.dist(p, fixed[v]) < 1e-12 for v in boundary):
            boundary.append(vertex(p))
    center = np.mean(corners, axis=0)
    boundary.sort(key=lambda v: math.atan2(fixed[v][1]-center[1], fixed[v][0]-center[0]))
    crossings = {m: vertex() for m in sorted({m for a in arcs for m in a['mids']})}
    edges, rotation = [], {}
    def edge(u, v):
        d = 2*len(edges)
        edges.append((u, v))
        rotation.setdefault(u, []).append(d)
        rotation.setdefault(v, []).append(d+1)
        return d
    boundary_darts = [edge(u, v) for u, v in zip(boundary, boundary[1:]+boundary[:1])]
    visits, chains, events = {}, [], []
    for a, (s, t) in zip(arcs, endpoints):
        original = [s]+[crossings[m] for m in a['mids']]+[t]
        word = words[a['curve']]
        start_index = next(i for i,x in enumerate(word) if abs(x)==a['start'])
        chain, event_indices = [s], [0]
        for u, v in zip(original, original[1:]):
            x, y = vertex(), vertex()
            edge(u, x); edge(x, y); edge(y, v)
            chain.extend((x, y, v)); event_indices.append(len(chain)-1)
        chains.append(chain); events.append(event_indices)
        for k, m in enumerate(a['mids'], 1):
            # The incoming/outgoing darts at this occurrence are recorded
            # before assigning the cyclic rotation at the shared vertex.
            v = crossings[m]
            incoming_neighbor, outgoing_neighbor = chain[3*k-1], chain[3*k+1]
            darts = rotation[v]
            target = lambda d: edges[d//2][1-(d%2)]
            incoming = next(d for d in darts if target(d)==incoming_neighbor)
            outgoing = next(d for d in darts if target(d)==outgoing_neighbor)
            signed = word[(start_index+k)%len(word)]
            visits[(m, signed>0)] = (outgoing, incoming)
    for m, v in crossings.items():
        po, pi = visits[(m, True)]; no, ni = visits[(m, False)]
        rotation[v] = [po, ni, pi, no] if exterior else [po, no, pi, ni]
    for i, v in enumerate(boundary):
        forward, backward = boundary_darts[i], boundary_darts[i-1]^1
        interior = [d for d in rotation[v] if d not in (forward, backward)]
        rotation[v] = [forward]+interior+[backward]
    successor = {}
    for darts in rotation.values():
        for i, d in enumerate(darts): successor[d^1] = darts[i-1]
    cells, seen = [], set()
    outer = boundary_darts[0]^1
    for start in range(2*len(edges)):
        if start in seen: continue
        walk, d = [], start
        while d not in seen:
            seen.add(d); walk.append(d); d = successor[d]
        if outer not in walk: cells.append(walk)
    original_vertex_count = len(fixed)
    triangles = []
    for walk in cells:
        ring = [edges[d//2][d%2] for d in walk]
        inner = [vertex() for _ in ring]; hub = vertex()
        for i, a in enumerate(ring):
            j = (i+1)%len(ring); b = ring[j]; p, q = inner[i], inner[j]
            triangles.extend(((a,b,q), (a,q,p), (p,q,hub)))
    return dict(fixed=fixed, boundary=boundary, triangles=triangles,
                chains=chains, events=events,
                disk_euler=original_vertex_count-len(edges)+len(cells),
                cell_walks=cells)


def embed_face(arcs, words, corners, *, exterior=False):
    mesh = triangulate_face(arcs, words, corners, exterior=exterior)
    fixed, triangles = mesh['fixed'], mesh['triangles']
    neighbors = [set() for _ in fixed]
    for triangle in triangles:
        for a,b in zip(triangle, triangle[1:]+triangle[:1]):
            neighbors[a].add(b); neighbors[b].add(a)
    free = [v for v,p in enumerate(fixed) if p is None]
    indices = {v:i for i,v in enumerate(free)}
    rows, cols, values = [], [], []
    rhs = np.zeros((len(free),2))
    for v in free:
        i = indices[v]
        rows.append(i); cols.append(i); values.append(len(neighbors[v]))
        for w in sorted(neighbors[v]):
            if fixed[w] is None:
                rows.append(i); cols.append(indices[w]); values.append(-1)
            else: rhs[i] += fixed[w]
    matrix = coo_matrix((values,(rows,cols)), shape=(len(free),len(free))).tocsc()
    solution = spsolve(matrix, rhs)
    points = [p if p is not None else tuple(solution[indices[v]]) for v,p in enumerate(fixed)]
    mesh['points'] = points
    for a, chain, events in zip(arcs, mesh['chains'], mesh['events']):
        a['poly_model'] = [points[v] for v in chain]
        a['event_indices'] = events
        a['events_model'] = [a['poly_model'][k] for k in events]
        a['event_paths_model'] = [a['poly_model'][i:j+1] for i,j in zip(events,events[1:])]
    if DEBUG_MODE:
        problems = audit_mesh(mesh)
        if problems: raise AssertionError(f"{arcs[0]['face']}: {problems}; words={words}")
    return mesh


def audit_mesh(mesh):
    """Construction diagnostics, never called by the default drawing path."""
    problems = []
    if mesh['disk_euler'] != 1:
        problems.append(f"fixed-boundary rotation has Euler characteristic {mesh['disk_euler']}, not a disk")
    triangles, points = mesh['triangles'], mesh['points']
    boundary = set(mesh['boundary'])
    boundary_edges = {frozenset((u,v)) for u,v in zip(mesh['boundary'],mesh['boundary'][1:]+mesh['boundary'][:1])}
    seen, counts = set(), {}
    for triangle in triangles:
        key = frozenset(triangle)
        if len(key)!=3 or key in seen: problems.append('not a simplicial triangle')
        seen.add(key)
        p,q,r = (points[v] for v in triangle)
        area = (q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])
        if not math.isfinite(area) or area<=0: problems.append('nonpositive or nonfinite triangle')
        for u,v in zip(triangle,triangle[1:]+triangle[:1]):
            e = frozenset((u,v)); counts[e] = counts.get(e,0)+1
            if u in boundary and v in boundary and e not in boundary_edges:
                problems.append('interior edge between boundary vertices')
    if any(n != (1 if e in boundary_edges else 2) for e,n in counts.items()):
        problems.append('triangle incidence is not a disk triangulation')
    return list(dict.fromkeys(problems))
