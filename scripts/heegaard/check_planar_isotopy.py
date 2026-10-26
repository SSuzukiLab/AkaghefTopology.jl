"""Independent check that an alpha-cut drawing realizes its combinatorial data.

Equivalence: the drawing is a plane graph G (beta routes + puncture circles) in
S^2.  If every beta route is simple, disjoint from the other routes and from the
open disks, and the ports sit on each circle in the order the data prescribes,
then G is determined up to orientation-preserving homeomorphism of S^2 by its
rotation system.  This script checks those hypotheses geometrically and then
compares the faces traced from the drawing with the complementary regions of
alpha+beta traced from the signed intersection data alone.  It shares no code
with Topology.Heegaard; only the documented sign convention is reused:
positive crossing = CCW (alpha_out, beta_out, alpha_in, beta_in).

usage: python check_planar_isotopy.py DIAGRAM.json LAYOUT.json [...pairs] [--out report.json]
"""
import json, math, sys
from collections import Counter, defaultdict
import numpy as np


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def diagram_of(doc):
    while "curves" not in doc and "diagram" in doc:
        doc = doc["diagram"]
    return doc


def layout_of(doc):
    return doc["layout"] if "layout" in doc and "routes" not in doc else doc


# ---------- combinatorial side ----------

def combinatorial_faces(d):
    curves = [("alpha", c) for c in d["curves"]["alpha"]] + [("beta", c) for c in d["curves"]["beta"]]
    sign = {x["id"]: x["sign"] for x in d["intersections"]}
    darts = []  # (origin, origin_port, target, target_port, edge_key)
    for fam, c in curves:
        v = c["visits"]
        for k in range(len(v)):
            a, b = v[k], v[(k + 1) % len(v)]
            key = (fam, c["id"], k)
            darts.append((a, fam + "_out", b, fam + "_in", key))
            darts.append((b, fam + "_in", a, fam + "_out", key))
    out = {(o, p): i for i, (o, p, *_rest) in enumerate(darts)}
    rot = {x: (["alpha_out", "beta_out", "alpha_in", "beta_in"] if s == 1 else
               ["alpha_out", "beta_in", "alpha_in", "beta_out"]) for x, s in sign.items()}
    nxt = []
    for (_o, _op, t, tp, _k) in darts:
        r = rot[t]
        nxt.append(out[(t, r[(r.index(tp) - 1) % 4])])
    seen, faces = set(), []
    for i in range(len(darts)):
        if i in seen:
            continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j); cyc.append(darts[j][4]); j = nxt[j]
        faces.append(cyc)
    V, E = len(sign), len(darts) // 2
    genus = (2 - (V - E + len(faces))) // 2
    return faces, genus


# ---------- geometric side ----------

def bezier_samples(route, per_cubic=12):
    c = route.get("curve")
    if isinstance(c, dict) and c.get("type") == "bezier_chain":
        pts = []
        for seg in c["segments"]:
            P = np.array(seg["points"], float)
            if seg.get("type") == "cubic" and len(P) == 4:
                t = np.linspace(0, 1, per_cubic + 1)[:, None]
                S = ((1 - t) ** 3) * P[0] + 3 * ((1 - t) ** 2) * t * P[1] + 3 * (1 - t) * t * t * P[2] + t ** 3 * P[3]
            else:
                S = P
            pts.extend(S if not pts else S[1:])
        return np.array(pts)
    return np.array(route["points"], float)


def segments_cross(A, B, eps=1e-12):
    """Proper or touching intersection between segment arrays A (n,2,2) and B (m,2,2)."""
    p, r = A[:, None, 0, :], A[:, None, 1, :] - A[:, None, 0, :]
    q, s = B[None, :, 0, :], B[None, :, 1, :] - B[None, :, 0, :]
    cross = lambda u, v: u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]
    rxs = cross(r, s)
    qp = q - p
    with np.errstate(divide="ignore", invalid="ignore"):
        t = cross(qp, s) / rxs
        u = cross(qp, r) / rxs
    hit = (np.abs(rxs) > eps) & (t >= -eps) & (t <= 1 + eps) & (u >= -eps) & (u <= 1 + eps)
    return hit


def crossing_report(polys, punct):
    segs, owner, idx = [], [], []
    for r, P in enumerate(polys):
        for i in range(len(P) - 1):
            segs.append((P[i], P[i + 1])); owner.append(r); idx.append(i)
    S = np.array(segs); owner = np.array(owner); idx = np.array(idx)
    n = len(S)
    lo, hi = S.min(axis=1), S.max(axis=1)
    bad = []
    step = 400
    for a in range(0, n, step):
        A = S[a:a + step]
        H = segments_cross(A, S)
        ov = (lo[a:a + step, None, 0] <= hi[None, :, 0]) & (lo[None, :, 0] <= hi[a:a + step, None, 0]) & \
             (lo[a:a + step, None, 1] <= hi[None, :, 1]) & (lo[None, :, 1] <= hi[a:a + step, None, 1])
        H &= ov
        I, J = np.nonzero(H)
        I = I + a
        keep = I < J
        for i, j in zip(I[keep], J[keep]):
            if owner[i] == owner[j] and abs(idx[i] - idx[j]) <= 1:
                continue  # consecutive segments of one route share a vertex
            bad.append((int(owner[i]), int(idx[i]), int(owner[j]), int(idx[j])))
    inside = []
    for r, P in enumerate(polys):
        for c in punct:
            dist = np.hypot(P[1:-1, 0] - c["center"][0], P[1:-1, 1] - c["center"][1])
            k = np.nonzero(dist < c["radius"] * (1 - 1e-9))[0]
            if len(k):
                inside.append((r, c["id"], int(len(k)), float(dist.min() / c["radius"])))
    return bad, inside


def geometric_faces(d, L, polys):
    alpha = d["curves"]["alpha"]; beta = {c["id"]: c for c in d["curves"]["beta"]}
    punct = {p["id"]: p for p in L["punctures"]}
    ports = {p["id"]: p for p in L["ports"]}
    issues = []
    # beta edge of each route
    route_key = []
    for r in L["routes"]:
        v = beta[r["curve_id"]]["visits"]
        a, b = ports[r["start_port"]]["intersection_id"], ports[r["end_port"]]["intersection_id"]
        ks = [k for k in range(len(v)) if v[k] == a and v[(k + 1) % len(v)] == b]
        k = r["segment_id"] - 1 if (r["segment_id"] - 1) in ks else (ks[0] if len(ks) == 1 else None)
        if k is None:
            issues.append(f"route {r['segment_id']} of beta {r['curve_id']} matches no unique beta edge")
        route_key.append(("beta", r["curve_id"], k))
    if Counter(route_key) != Counter(("beta", c["id"], k) for c in beta.values() for k in range(len(c["visits"]))):
        issues.append("routes are not in bijection with beta edges")
    # ports per circle, CCW
    halfedges = defaultdict(list)  # port -> [(angle, edge_id, dir_tag)]
    edges = {}  # edge_id -> (port_u, port_v, key, sample polyline u->v)
    eid = 0
    for r, key, P in zip(L["routes"], route_key, polys):
        u, v = r["start_port"], r["end_port"]
        edges[eid] = (u, v, key, P)
        du, dv = P[1] - P[0], P[-2] - P[-1]
        halfedges[u].append((math.atan2(du[1], du[0]), eid, +1))
        halfedges[v].append((math.atan2(dv[1], dv[0]), eid, -1))
        eid += 1
    for bid, c in punct.items():
        cx, cy = c["center"]
        on = [p for p in ports.values() if p["boundary_id"] == bid]
        ang = {p["id"]: math.atan2(p["point"][1] - cy, p["point"][0] - cx) % (2 * math.pi) for p in on}
        for p in on:
            rad = math.hypot(p["point"][0] - cx, p["point"][1] - cy)
            if abs(rad - c["radius"]) > 1e-6 * max(1, c["radius"]):
                issues.append(f"port {p['id']} is off its circle")
        order = sorted(on, key=lambda p: ang[p["id"]])
        a_curve = alpha[c["pair"] - 1]; v = a_curve["visits"]
        want = v if c["side"] == "plus" else list(reversed(v))
        got = [p["intersection_id"] for p in order]
        if not any(got == want[i:] + want[:i] for i in range(len(want))):
            issues.append(f"circle {bid}: CCW port order {got} differs from data {want}")
        m = len(order)
        for i, p in enumerate(order):
            q = order[(i + 1) % m]
            a0, a1 = ang[p["id"]], ang[q["id"]]
            if m == 1 or a1 <= a0:
                a1 += 2 * math.pi
            ip, iq = p["intersection_id"], q["intersection_id"]
            if c["side"] == "plus":
                ks = [k for k in range(len(v)) if v[k] == ip and v[(k + 1) % len(v)] == iq]
            else:
                ks = [k for k in range(len(v)) if v[k] == iq and v[(k + 1) % len(v)] == ip]
            key = ("alpha", a_curve["id"], ks[0] if len(ks) == 1 else None)
            if len(ks) != 1:
                issues.append(f"circle {bid}: arc {ip}->{iq} matches {len(ks)} alpha edges")
            t = np.linspace(a0, a1, 24)
            P = np.stack([cx + c["radius"] * np.cos(t), cy + c["radius"] * np.sin(t)], 1)
            edges[eid] = (p["id"], q["id"], key, P)
            halfedges[p["id"]].append((a0 + math.pi / 2, eid, +1))       # CCW tangent
            halfedges[q["id"]].append((a1 - math.pi / 2 + 0.0, eid, -1))  # CW tangent at q
            eid += 1
    rot = {p: sorted(((math.atan2(math.sin(a), math.cos(a)), e, s) for a, e, s in hs)) for p, hs in halfedges.items()}
    pos = {(p, e, s): i for p, hs in rot.items() for i, (_a, e, s) in enumerate(hs)}
    def head(e, s):
        u, v = edges[e][0], edges[e][1]
        return (v, u) if s == +1 else (u, v)
    seen, faces = set(), []
    for p, hs in rot.items():
        for _a, e, s in hs:
            if (e, s) in seen:
                continue
            cyc, cur = [], (e, s)
            while cur not in seen:
                seen.add(cur); cyc.append(cur)
                v, _u = head(*cur)
                twin = (cur[0], -cur[1])
                i = pos[(v, twin[0], twin[1])]
                _a2, e2, s2 = rot[v][(i - 1) % len(rot[v])]
                cur = (e2, s2)
            faces.append(cyc)
    def area(cyc):
        pts = np.concatenate([edges[e][3] if s == 1 else edges[e][3][::-1] for e, s in cyc])
        x, y = pts[:, 0], pts[:, 1]
        return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    disks, regions = [], []
    for cyc in faces:
        keys = [edges[e][2] for e, _s in cyc]
        is_disk = all(k[0] == "alpha" for k in keys) and all(s == 1 for _e, s in cyc) and \
            len({edges[e][0] for e, _s in cyc} | {edges[e][1] for e, _s in cyc}) and \
            len({ports[edges[e][0]]["boundary_id"] for e, _s in cyc}) == 1
        (disks if is_disk else regions).append((keys, area(cyc)))
    V = len(ports); E = len(edges)
    return regions, disks, issues, V - E + len(faces)


def cyc_eq(a, b):
    return len(a) == len(b) and any(a == b[i:] + b[:i] for i in range(len(b) or 1))


def check(diagram_path, layout_path, mutate=None):
    d = diagram_of(load(diagram_path)); L = layout_of(load(layout_path))
    if mutate:
        mutate(d, L)
    comb, genus = combinatorial_faces(d)
    polys = [bezier_samples(r) for r in L["routes"]]
    bad, inside = crossing_report(polys, L["punctures"])
    regions, disks, issues, euler = geometric_faces(d, L, polys)
    fwd = [any(cyc_eq(k, c) for c in comb) for k, _a in regions]
    rev = [any(cyc_eq(list(reversed(k)), c) for c in comb) for k, _a in regions]
    matched_fwd = all(fwd) and len(regions) == len(comb)
    matched_rev = all(rev) and len(regions) == len(comb)
    outer = [i for i, (_k, a) in enumerate(regions) if a < 0]
    outer_region = None
    if len(outer) == 1:
        k = regions[outer[0]][0]
        k2 = k if matched_fwd else list(reversed(k))
        outer_region = next((i + 1 for i, c in enumerate(comb) if cyc_eq(k2, c)), None)
    ok = (not bad and not inside and not issues and euler == 2 and len(disks) == len(L["punctures"])
          and (matched_fwd or matched_rev) and len(outer) == 1)
    return {"diagram": diagram_path, "layout": layout_path, "genus": genus,
            "combinatorial_regions": len(comb), "drawn_regions": len(regions), "disk_faces": len(disks),
            "euler_characteristic_of_plane_graph": euler,
            "route_crossings": len(bad), "route_crossing_examples": bad[:5],
            "routes_entering_disks": inside[:5], "data_issues": issues,
            "faces_match": "same orientation" if matched_fwd else ("mirror orientation" if matched_rev else "NO"),
            "outer_face_combinatorial_region": outer_region, "unbounded_faces": len(outer),
            "samples": int(sum(len(p) for p in polys)), "pass": ok}


def _mirror(d, L):
    """Reflect the drawing in x: same combinatorics on the page, opposite orientation."""
    for p in L["punctures"]:
        p["center"][0] *= -1
    for p in L["ports"]:
        p["point"][0] *= -1
    for r in L["routes"]:
        r.pop("curve", None)
        for q in r["points"]:
            q[0] *= -1


def _swap_ports(d, L):
    """Exchange which intersections two ports on one circle carry (data/drawing mismatch)."""
    by = defaultdict(list)
    for p in L["ports"]:
        by[p["boundary_id"]].append(p)
    ps = next(v for v in by.values() if len(v) >= 3)
    a, b = ps[0], ps[1]
    a["point"], b["point"] = b["point"], a["point"]


def _cross_routes(d, L):
    """Push the midpoint of one route through another route."""
    r0, r1 = L["routes"][0], L["routes"][1]
    for r in (r0, r1):
        r.pop("curve", None)
    m = r1["points"][len(r1["points"]) // 2]
    r0["points"].insert(len(r0["points"]) // 2, [m[0] + 1e-3, m[1] + 1e-3])
    r0["points"].insert(len(r0["points"]) // 2 + 1, [m[0] - 1e-3, m[1] - 1e-3])


def negative_controls(diagram_path, layout_path):
    res = {}
    for name, f in (("mirror", _mirror), ("swap_ports", _swap_ports), ("cross_routes", _cross_routes)):
        r = check(diagram_path, layout_path, f)
        res[name] = {"pass": r["pass"], "faces_match": r["faces_match"], "route_crossings": r["route_crossings"],
                     "data_issues": r["data_issues"][:2]}
    return res


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--negative-controls":
        print(json.dumps(negative_controls(args[1], args[2]), indent=1, ensure_ascii=False)); sys.exit(0)
    out = None
    if "--out" in args:
        i = args.index("--out"); out = args[i + 1]; del args[i:i + 2]
    results = [check(args[i], args[i + 1]) for i in range(0, len(args), 2)]
    for r in results:
        print(("PASS " if r["pass"] else "FAIL ") + r["layout"],
              {k: r[k] for k in ("genus", "combinatorial_regions", "drawn_regions", "route_crossings",
                                 "faces_match", "outer_face_combinatorial_region", "data_issues")})
    if out:
        with open(out, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=1, ensure_ascii=False)
    sys.exit(0 if all(r["pass"] for r in results) else 1)
