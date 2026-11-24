#!/usr/bin/env python3
"""Explicit finite face pairings to SnapPy JSON.

cf: snappy.snap.t3mlite.Mcomplex.snappy_triangulation and Mcomplex.isosig
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — schema and evidence limits
"""

from __future__ import annotations
import argparse
import contextlib
import io
import itertools
import json
import math
import subprocess
import sys


class SnapPyUnavailable(RuntimeError):
    pass


class InvalidTriangulation(ValueError):
    pass


class HomSearchLimit(RuntimeError):
    pass


def attempt(function):
    try:
        return {"status": "ok", "value": function()}
    except Exception as exc:
        return {"status": "error", "error": {"type": type(exc).__name__, "message": str(exc)}}


def backend():
    try:
        import snappy
        from snappy.snap.t3mlite import Mcomplex, Tetrahedron
        from snappy.snap.t3mlite.simplex import TwoSubsimplices
    except (ImportError, OSError) as exc:
        raise SnapPyUnavailable(str(exc)) from exc
    return snappy, Mcomplex, Tetrahedron, TwoSubsimplices


def validate_pairings(data):
    """Check total face use, permutations, connectedness and supplied orientation."""
    n = data["tetrahedra"]
    if type(n) is not int or n <= 0:
        raise InvalidTriangulation("tetrahedra must be a positive integer")
    seen = set()
    neighbors = [set() for _ in range(n)]
    signs = data.get("signs")
    if signs is not None and (len(signs) != n or any(s not in (-1, 1) for s in signs)):
        raise InvalidTriangulation("signs must contain one +/-1 per tetrahedron")
    for g in data["gluings"]:
        s, t, a, b = (g[k] for k in ("tet", "to_tet", "omitted", "to_omitted"))
        if any(type(x) is not int for x in (s, t, a, b)) or not (
            0 <= s < n and 0 <= t < n and 0 <= a < 4 and 0 <= b < 4
        ):
            raise InvalidTriangulation("face index out of range")
        p = g["permutation"]
        if len(p) != 4 or any(type(x) is not int for x in p) or sorted(p) != list(range(4)) or p[a] != b:
            raise InvalidTriangulation("incompatible vertex permutation")
        for face in ((s, a), (t, b)):
            if face in seen:
                raise InvalidTriangulation(f"face used more than once: {face}")
            seen.add(face)
        if signs is not None:
            parity = (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
            if signs[s] != -parity * signs[t]:
                raise InvalidTriangulation("incompatible tetrahedron orientations")
        neighbors[s].add(t)
        neighbors[t].add(s)
    if len(seen) != 4 * n:
        raise InvalidTriangulation("unpaired face")
    reached, pending = {0}, [0]
    while pending:
        for t in neighbors[pending.pop()]:
            if t not in reached:
                reached.add(t)
                pending.append(t)
    if len(reached) != n:
        raise InvalidTriangulation("disconnected triangulation")


def complex_topology(complex_):
    v, e, f, t = (len(getattr(complex_, name)) for name in ("Vertices", "Edges", "Faces", "Tetrahedra"))
    return {"vertices": v, "edges": e, "faces": f, "tetrahedra": t,
            "euler_characteristic": v-e+f-t,
            "vertex_link_genera": [int(x.link_genus()) for x in complex_.Vertices]}


def _orient_gluings(data):
    """Legacy reader for a sign array that is not already a positive frame.

    ``hc_triangulation`` writes every sign as +1 and conjugates each negative
    tetrahedron into the vertex numbering before the file is stored. SnapPy
    does not correct that output. A record whose first sign is still -1 is a
    file written before that change; conjugating its negative tetrahedra by
    the odd swap (0 1) recovers the orientation the sign array named.
    """
    signs = data.get("signs")
    if not signs or signs[0] != -1:
        return data
    swap = (1, 0, 2, 3)

    def image(vertex, flip):
        return swap[vertex] if flip else vertex

    flips = [sign == -1 for sign in signs]
    gluings = []
    for gluing in data["gluings"]:
        flip_source = flips[gluing["tet"]]
        flip_target = flips[gluing["to_tet"]]
        permutation = gluing["permutation"]
        row = dict(gluing)
        row["omitted"] = image(gluing["omitted"], flip_source)
        row["to_omitted"] = image(gluing["to_omitted"], flip_target)
        row["permutation"] = [
            image(permutation[image(vertex, flip_source)], flip_target) for vertex in range(4)
        ]
        gluings.append(row)
    oriented = dict(data)
    oriented["gluings"] = gluings
    oriented["signs"] = [1] * len(signs)
    return oriented


def make_complex(data, api):
    validate_pairings(data)
    data = _orient_gluings(data)
    validate_pairings(data)
    _, Mcomplex, Tetrahedron, faces = api
    tets = [Tetrahedron() for _ in range(data["tetrahedra"])]
    for g in data["gluings"]:
        tets[g["tet"]].attach(faces[g["omitted"]], tets[g["to_tet"]], tuple(g["permutation"]))
    complex_ = Mcomplex(tets)
    topology = complex_topology(complex_)
    if topology["euler_characteristic"] != 0 or any(x != 0 for x in topology["vertex_link_genera"]):
        raise InvalidTriangulation("closed finite manifold requires sphere vertex links and Euler characteristic zero")
    return complex_, topology


def presentation(group):
    names, words = list(group.generators()), list(group.relators())
    if any(len(x) != 1 or not x.islower() or not x.isascii() for x in names):
        raise ValueError("SnapPy generator syntax outside single ASCII lower-case letters")
    indices = {name: i for i, name in enumerate(names)}
    relators = []
    for word in words:
        letters = []
        for char in word:
            if char.lower() not in indices:
                raise ValueError(f"unsupported SnapPy relator token: {char!r}")
            letters.append((indices[char.lower()], 1 if char.islower() else -1))
        relators.append(letters)
    return {"generators": names, "relators": words, "letters": relators, "display": str(group)}


def hom_count(pres, target, limit=2_000_000):
    """Count assignments, not conjugacy classes; early checks prune satisfied prefixes.

    cf: lib/topology.jl/src/heegaard/hc_compare.jl::hc_hom_count — independent
    Python evaluation of the SnapPy presentation, never of HC data.
    """
    n, relations = len(pres["generators"]), pres["letters"]
    if target in ("Z2", "Z3", "Z5"):
        q, identity = int(target[1:]), 0
        multiply = lambda a, b: (a + b) % q
        inverse = lambda a: (-a) % q
    elif target == "S3":
        permutations = list(itertools.permutations(range(3)))
        q, identity = 6, 0
        index = {p: i for i, p in enumerate(permutations)}
        table = [[index[tuple(p[r[k]] for k in range(3))] for r in permutations] for p in permutations]
        inverses = [next(j for j in range(q) if table[i][j] == identity) for i in range(q)]
        multiply = lambda a, b: table[a][b]
        inverse = lambda a: inverses[a]
    else:
        raise ValueError(target)
    active = sorted({i for word in relations for i, _ in word})
    remap = {old: new for new, old in enumerate(active)}
    ready = [[] for _ in active]
    for word in relations:
        if word:
            mapped = [(remap[i], exp) for i, exp in word]
            ready[max(i for i, _ in mapped)].append(mapped)
    assignment = [identity] * len(active)
    visited = 0

    def visit(depth):
        nonlocal visited
        if depth == len(active):
            return 1
        total = 0
        for value in range(q):
            visited += 1
            if visited > limit:
                raise HomSearchLimit(f"{target}: exceeded {limit} assignment nodes; no count returned")
            assignment[depth] = value
            valid = True
            for word in ready[depth]:
                result = identity
                for i, exp in word:
                    result = multiply(result, assignment[i] if exp == 1 else inverse(assignment[i]))
                if result != identity:
                    valid = False
                    break
            if valid:
                total += visit(depth + 1)
        return total

    return str(visit(0) * q ** (n - len(active)))


def finite_float(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"nonfinite numeric result: {value}")
    return number


def analyze(data, references=(), isosig_only=False, hom_limit=2_000_000):
    api = backend()
    snappy = api[0]
    complex_, topology = make_complex(data, api)
    result = {"ok": True, "schema": "heegaard.snappy-result/1", "snappy_version": snappy.__version__,
              "topology": topology, "isosig": complex_.isosig(),
              "orientation_scope": "face signs retained in input; isoSig and isometry comparisons are unoriented"}
    if isosig_only:
        return result
    tri = complex_.snappy_triangulation(remove_finite_vertices=True)
    manifold = tri.with_hyperbolic_structure()
    fillings = [list(map(float, cusp["filling"])) for cusp in manifold.cusp_info()]
    result["fillings"] = fillings
    result["all_cusps_filled"] = all(any(x != 0 for x in pair) for pair in fillings)
    if not result["all_cusps_filled"]:
        raise InvalidTriangulation("finite-vertex removal left an unfilled cusp; refusing to analyze a different manifold")
    result["filled_triangulation"] = attempt(lambda: manifold.filled_triangulation().triangulation_isosig(decorated=False))
    result["solution_type"] = attempt(manifold.solution_type)
    result["volume"] = attempt(lambda: finite_float(manifold.volume()))

    def homology():
        h1 = manifold.homology()
        divisors = [int(x) for x in h1.elementary_divisors()]
        return {"display": str(h1), "elementary_divisors": divisors,
                "free_rank": divisors.count(0), "torsion": sorted(x for x in divisors if x > 1)}

    result["homology"] = attempt(homology)
    result["fundamental_group"] = attempt(lambda: presentation(manifold.fundamental_group()))
    if result["fundamental_group"]["status"] == "ok":
        pres = result["fundamental_group"]["value"]
        result["hom_counts"] = {g: attempt(lambda g=g: hom_count(pres, g, hom_limit)) for g in ("Z2", "Z3", "Z5", "S3")}
    else:
        result["hom_counts"] = {g: {"status": "not_computed", "reason": "fundamental_group_unavailable"}
                                for g in ("Z2", "Z3", "Z5", "S3")}
    result["identify"] = attempt(lambda: [str(m) for m in manifold.identify()])
    result["references"] = {name: attempt(lambda name=name: bool(manifold.is_isometric_to(snappy.Manifold(name))))
                            for name in references}
    result["hyperbolicity"] = "not_certified; solution_type and volume are raw numerical outputs"
    return result


def serialize_complex(complex_, faces):
    """Serialize t3mlite's actual maps; no separate local gluing convention."""
    gluings, used = [], set()
    for tet in complex_.Tetrahedra:
        for a, face in enumerate(faces):
            if (tet.Index, a) in used:
                continue
            neighbor, p = tet.Neighbor[face], tet.Gluing[face].tuple()
            if neighbor is None:
                raise InvalidTriangulation("boundary face in filled triangulation")
            b = p[a]
            used.update(((tet.Index, a), (neighbor.Index, b)))
            gluings.append({"tet": tet.Index, "omitted": a, "to_tet": neighbor.Index,
                            "to_omitted": b, "permutation": list(p)})
    data = {"schema": "heegaard.closed-triangulation/1", "tetrahedra": len(complex_), "gluings": gluings}
    validate_pairings(data)
    return data


def export_catalog(name):
    """Reuse SnapPy's filled triangulation and spine exporter; stop before branching."""
    snappy, Mcomplex, _, faces = backend()
    manifold = snappy.Manifold(name)
    if any(cusp["is_complete"] for cusp in manifold.cusp_info()):
        raise InvalidTriangulation("catalog entry must be closed (all cusps filled)")
    complex_ = Mcomplex(manifold.filled_triangulation())
    topology = complex_topology(complex_)
    if any(g != 0 for g in topology["vertex_link_genera"]):
        raise InvalidTriangulation("filled triangulation still has non-sphere links")
    data = serialize_complex(complex_, faces)
    spine = io.StringIO()
    complex_.save(spine, format="spine")
    return {"ok": True, "snappy_version": snappy.__version__, "catalog_name": name,
            "triangulation": data, "topology": topology, "isosig": complex_.isosig(),
            "dual_spine": spine.getvalue(), "removed_vertex_balls": topology["vertices"],
            "spine_scope": "dual spine after removing one open ball per finite vertex",
            "branching": {"status": "not_computed",
            "reason": "unbranched dual spine does not supply a normal flow o-graph"},
            "heegaard_diagram": {"status": "not_computed", "reason": "branching_and_ograph_encoding_required"}}


def pachner_neighbors(data):
    """Deterministic safe 2->3 neighbors, checked by a legal inverse 3->2.

    cf: snappy.snap.t3mlite.Mcomplex.two_to_three/three_to_two/rebuild (3.3.2).
    isoSig evidence and state deduplication here are unoriented.
    """
    api = backend()
    _, _, _, faces = api
    original, original_topology = make_complex(data, api)
    key = original.isosig()
    states, seen, rejected = [], {key}, []
    for g in data['gluings']:
        tet_index, omitted = g['tet'], g['omitted']
        child = original.copy()
        child.rebuild()
        objects = {id(t): t.Index for t in child.Tetrahedra}
        if not child.two_to_three(faces[omitted], child.Tetrahedra[tet_index], unsafe_mode=False):
            continue
        child.rebuild()  # native move does not rebuild finite vertex classes
        child_data = serialize_complex(child, faces)
        checked, topology = make_complex(child_data, api)
        if topology['vertices'] != original_topology['vertices'] or topology['euler_characteristic'] != 0:
            raise InvalidTriangulation('legal 2->3 changed the finite vertex count or Euler characteristic')
        child_key = checked.isosig()
        if child_key in seen:
            continue
        inverse = None
        for edge_index, edge in enumerate(checked.Edges):
            if edge.valence() != 3 or not edge.distinct():
                continue
            trial = checked.copy()
            trial.rebuild()
            if trial.three_to_two(trial.Edges[edge_index], unsafe_mode=False):
                trial.rebuild()
                if trial.isosig() == key:
                    inverse = edge_index
                    break
        if inverse is None:
            rejected.append({'tet': tet_index, 'omitted': omitted, 'reason': 'inverse_move_not_verified'})
            continue
        seen.add(child_key)
        surviving = {str(objects[id(t)]): t.Index for t in child.Tetrahedra if id(t) in objects}
        new_tetrahedra = [t.Index for t in child.Tetrahedra if id(t) not in objects]
        states.append({'triangulation': child_data, 'isosig': child_key, 'topology': topology,
            'move': {'kind': 'two_to_three', 'tet': tet_index, 'omitted': omitted,
                     'before_isosig': key, 'after_isosig': child_key,
                     'surviving_tetrahedra': surviving, 'new_tetrahedra': new_tetrahedra,
                     'inverse_kind': 'three_to_two', 'inverse_edge': inverse,
                     'inverse_restores_isosig': True, 'unsafe_mode': False}})
    return {'ok': True, 'isosig': key, 'states': states, 'rejected': rejected,
            'orientation_scope': 'unoriented isomorphism; native safe local moves',
            'snappy_version': api[0].__version__}


def compare_triangulations(left, right):
    api = backend()
    a, _ = make_complex(left, api)
    b, _ = make_complex(right, api)
    if a.isosig() == b.isosig():
        return {"ok": True, "same_manifold": True, "evidence": "combinatorial_isomorphism",
                "orientation": "unoriented"}
    ma, mb = (c.snappy_triangulation().with_hyperbolic_structure() for c in (a, b))
    comparison = attempt(lambda: bool(ma.is_isometric_to(mb)))
    return {"ok": True,
            "same_manifold": True if comparison.get("value") is True else None,
            "evidence": "snappy_numerical_isometry" if comparison.get("value") is True else "unresolved",
            "orientation": "unoriented", "is_isometric_to": comparison}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="JSON file, or - for stdin")
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument("--isosig-only", action="store_true")
    parser.add_argument("--export-catalog", action="store_true")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--pachner-neighbors", action="store_true")
    parser.add_argument("--hom-limit", type=int, default=2_000_000)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.worker:
        try:
            if args.timeout <= 0:
                raise ValueError("timeout must be positive")
            # A native kernel call cannot be interrupted by a Python thread.
            # subprocess.run kills and reaps this exact child on timeout.
            result = subprocess.run([sys.executable, __file__, *sys.argv[1:], "--worker"],
                input=sys.stdin.read() if args.input == "-" else None,
                capture_output=True, text=True, encoding="utf-8", timeout=args.timeout)
            sys.stdout.write(result.stdout)
            sys.stderr.write(result.stderr)
            return result.returncode
        except subprocess.TimeoutExpired:
            print(json.dumps({"ok": False, "error": {"type": "SnapPyTimeout", "message": f"deadline {args.timeout}s exceeded"}}))
            return 2
        except Exception as exc:
            print(json.dumps({"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}}))
            return 2
    try:
        if args.hom_limit <= 0:
            raise ValueError("hom-limit must be positive")
        if args.input == "-":
            data = json.load(sys.stdin)
        else:
            with open(args.input, encoding="utf-8-sig") as stream:
                data = json.load(stream)
        with contextlib.redirect_stdout(sys.stderr):
            if args.pachner_neighbors:
                result = pachner_neighbors(data)
            elif args.compare:
                result = compare_triangulations(data["left"], data["right"])
            elif args.export_catalog:
                result = export_catalog(data["catalog_name"])
            else:
                result = analyze(data, args.reference, args.isosig_only, args.hom_limit)
    except Exception as exc:
        result = {"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}}
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    return 0 if result["ok"] else 2


if __name__ == "__main__":
    sys.exit(main())
