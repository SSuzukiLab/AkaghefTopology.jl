"""Regina exact TV consumer; stdin/stdout are JSON, never numerical evaluation.

cf: Regina 7.4.1, Triangulation3.turaevViro and Cyclotomic.__getitem__ —
exact field/coefficients; REPORT-3.md records the checked API conventions.
"""
import argparse
import importlib.metadata
import json
import subprocess
import sys


def compute(data):
    import regina
    if "jobs" in data:
        return dict(ok=True, schema="heegaard.regina-tv-jobs/1",
                    jobs=[compute(job) for job in data["jobs"]])
    if "requests" in data:
        request_data = dict(data)
        requests = request_data.pop("requests")
        if not requests:
            raise ValueError("requests must be nonempty")
        return dict(ok=True, schema="heegaard.regina-tv-batch/1",
                    results=[compute(dict(request_data, r=r, parity=p)) for r, p in requests])
    n = data["tetrahedra"]
    if type(n) is not int or n <= 0:
        raise ValueError("tetrahedra must be a positive integer")
    r, parity = data["r"], data["parity"]
    if type(r) is not int or r < 3 or type(parity) is not bool:
        raise ValueError("r >= 3 and boolean parity are required")
    t = regina.Triangulation3()
    t.newTetrahedra(n)
    seen = set()
    for g in data["gluings"]:
        a, b, f, h = (g[k] for k in ("tet", "to_tet", "omitted", "to_omitted"))
        p = g["permutation"]
        if (any(type(x) is not int for x in (a, b, f, h)) or
                not (0 <= a < n and 0 <= b < n and 0 <= f < 4 and 0 <= h < 4) or
                len(p) != 4 or any(type(x) is not int for x in p) or
                sorted(p) != [0, 1, 2, 3] or p[f] != h):
            raise ValueError("invalid face pairing")
        if (a, f) == (b, h) or (a, f) in seen or (b, h) in seen:
            raise ValueError("a face must be paired exactly once with another face")
        seen.update(((a, f), (b, h)))
        t.tetrahedron(a).join(f, t.tetrahedron(b), regina.Perm4(*p))
    if len(seen) != 4 * n:
        raise ValueError("unpaired face")
    if not (t.isValid() and t.isClosed() and t.isConnected() and t.isOrientable()):
        raise ValueError("expected a valid connected orientable closed finite triangulation")
    if t.isIdeal():
        raise ValueError("ideal vertices are outside the finite closed input contract")
    if data.get("simplify", True):
        t.simplify()
    value = t.turaevViro(r, parity)
    expected = r if r % 2 and not parity else 2 * r
    if value.field() != expected:
        raise RuntimeError("Regina returned an unexpected or uninitialised cyclotomic field")
    return dict(ok=True, schema="heegaard.regina-tv/1", r=r, parity=parity,
                field=value.field(), coefficients=[str(value[i]) for i in range(value.degree())],
                regina_distribution=importlib.metadata.version("regina"),
                regina_engine=regina.versionString(), input_tetrahedra=n,
                computed_tetrahedra=t.size(), isosig=t.isoSig(),
                homology=str(t.homology()), normalization="Turaev-Viro-1992")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--input")
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args()
    try:
        if args.worker:
            result = compute(json.load(sys.stdin))
        else:
            # The native parent kills and reaps its worker on timeout, including
            # when Windows calls this script through wsl.exe.
            if args.input:
                with open(args.input, encoding="utf-8-sig") as source:
                    payload = source.read()
            else:
                payload = sys.stdin.read()
            proc = subprocess.run([sys.executable, __file__, "--worker"],
                                  input=payload, text=True, capture_output=True,
                                  timeout=args.timeout, check=False)
            if proc.returncode:
                raise RuntimeError(proc.stderr or "worker exited without a reply")
            result = json.loads(proc.stdout)
    except Exception as error:
        result = dict(ok=False, error=dict(type=type(error).__name__, message=str(error)))
    print(json.dumps(result, ensure_ascii=True))


if __name__ == "__main__":
    main()
