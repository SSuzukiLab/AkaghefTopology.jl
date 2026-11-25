"""Select the test files a change can affect (policy: test/TESTING.md).

Usage (from lib/topology.jl):
    python test/impacted.py                     # changed files from git (worktree vs HEAD, plus untracked)
    python test/impacted.py src/a.jl src/b.jl   # explicit changed files
    python test/impacted.py --tier integration  # allow tiers up to integration (default unit)
    python test/impacted.py --depth 2           # follow callers two levels (default 1)
    python test/impacted.py --why               # print the reason for each selected file
    python test/impacted.py --files-only        # only the file list, for the runner

Run the selection on the warm server:
    sh test/client.sh $(python test/impacted.py --files-only)

The dependency map is mechanical:
  * a source file depends on another when it uses a name the other defines
    (function, struct, const, macro, abstract type). Names defined in more than
    two files (overload families) and names shorter than four characters are not
    dependencies; export/import/using lines are ignored;
  * a test file depends on every source file whose names it uses, and on the
    manifest's `requires`;
  * a change selects every test that depends on a changed file or on a caller of
    it up to --depth levels. The full closure covers most of the library, so the
    daily acceptance run covers deeper effects; a miss found there is fixed in
    this map, not by running everything on every change.
"""
import os
import re
import subprocess
import sys
import tomllib
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # lib/topology.jl
TEST = os.path.join(ROOT, "test")
ORDER = {"unit": 1, "integration": 2, "acceptance": 3, "manual": 4}

DEF_RE = re.compile(
    r"^\s*(?:function\s+(?:[\w.]+\.)?([A-Za-z_]\w*!?)"
    r"|(?:mutable\s+)?struct\s+([A-Za-z_]\w*)"
    r"|abstract\s+type\s+([A-Za-z_]\w*)"
    r"|const\s+([A-Za-z_]\w*)"
    r"|macro\s+([A-Za-z_]\w*)"
    r"|([A-Za-z_]\w*!?)\s*\([^=\n]*\)\s*(?:where\b[^=\n]*)?=(?!=))",
    re.M)
IDENT_RE = re.compile(r"[A-Za-z_]\w*!?")
COMMENT_RE = re.compile(r"#=.*?=#|#[^\n]*", re.S)
STRING_RE = re.compile(r'"(?:\\.|[^"\\])*"')
EXPORT_RE = re.compile(r"^\s*(?:export|public|import|using)\b.*(?:,[ \t]*\n.*)*", re.M)
SKIP = {"show", "length", "getindex", "setindex!", "iterate", "hash", "isequal",
        "convert", "copy", "size", "eltype", "push!", "print", "println", "string",
        "main", "run", "check", "apply", "validate", "test", "eval", "include"}


def jl_files(*dirs):
    out = []
    for d in dirs:
        for base, _, files in os.walk(os.path.join(ROOT, d)):
            if os.sep + "." in base:
                continue
            out += [os.path.relpath(os.path.join(base, f), ROOT).replace(os.sep, "/")
                    for f in files if f.endswith(".jl")]
    return sorted(out)


def read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8", errors="replace") as io:
        return io.read().replace("\r\n", "\n")


def code_only(text):
    t = STRING_RE.sub('""', COMMENT_RE.sub("", text))
    return EXPORT_RE.sub("", t)


def build_map():
    sources = jl_files("src", "examples")
    defs = defaultdict(set)
    text = {}
    for f in sources:
        t = code_only(read(f))
        text[f] = t
        for m in DEF_RE.finditer(t):
            name = next(g for g in m.groups() if g)
            if name not in SKIP:
                defs[name].add(f)
    defs = {n: fs for n, fs in defs.items() if len(fs) <= 2 and len(n) >= 4}
    return sources, defs, text


def uses(defs, t, own):
    found = set()
    for name in set(IDENT_RE.findall(t)):
        for g in defs.get(name, ()):
            if g != own:
                found.add(g)
    return found


def main(argv):
    tier, depth = "unit", 1
    why = files_only = False
    changed = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("--tier", "--depth"):
            if a == "--tier":
                tier = argv[i + 1]
            else:
                depth = int(argv[i + 1])
            i += 2
            continue
        if a == "--why":
            why = True
        elif a == "--files-only":
            files_only = True
        else:
            changed.append(a.replace("\\", "/"))
        i += 1
    if not changed:
        diff = subprocess.run(["git", "-C", ROOT, "diff", "--name-only", "HEAD"],
                              capture_output=True, text=True).stdout
        new = subprocess.run(["git", "-C", ROOT, "ls-files", "--others", "--exclude-standard"],
                             capture_output=True, text=True).stdout
        changed = [p for p in (diff + new).split() if p.endswith(".jl") or p.endswith(".toml")]

    sources, defs, text = build_map()
    rdeps = defaultdict(set)
    for f in sources:
        for d in uses(defs, text[f], f):
            rdeps[d].add(f)

    impacted = set(c for c in changed if c in text)
    frontier = list(impacted)
    for _ in range(depth):
        nxt = []
        for f in frontier:
            for g in rdeps[f]:
                if g not in impacted:
                    impacted.add(g)
                    nxt.append(g)
        frontier = nxt

    with open(os.path.join(TEST, "manifest.toml"), "rb") as io:
        manifest = tomllib.load(io)
    timings = {}
    tpath = os.path.join(TEST, "timings.tsv")
    if os.path.exists(tpath):
        for line in open(tpath, encoding="utf-8", errors="replace"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                try:
                    timings[parts[2]] = float(parts[3])
                except ValueError:
                    pass

    changed_tests = {c[len("test/"):] for c in changed if c.startswith("test/")}
    selected = []
    deferred = []   # higher-tier tests that use a changed file directly
    direct = set(c for c in changed if c in text)
    for e in manifest["test"]:
        tf = e["file"]
        if ORDER[e["tier"]] > ORDER[tier]:
            path = os.path.join(TEST, tf)
            if e["tier"] != "manual" and os.path.exists(path) and uses(defs, code_only(read("test/" + tf)), None) & direct:
                deferred.append((tf, e["tier"]))
            continue
        reason = None
        if tf in changed_tests:
            reason = "test file changed"
        else:
            path = os.path.join(TEST, tf)
            tt = code_only(read("test/" + tf)) if os.path.exists(path) else ""
            hit = sorted(uses(defs, tt, None) & impacted)
            req = [r for r in e.get("requires", [])
                   if os.path.normpath(os.path.join("test", r)).replace(os.sep, "/") in impacted]
            if hit:
                reason = "uses " + ", ".join(hit[:3]) + (" ..." if len(hit) > 3 else "")
            elif req:
                reason = "requires " + ", ".join(req)
        if reason:
            selected.append((tf, e["tier"], reason))

    total = sum(timings.get(tf, 0.0) for tf, _, _ in selected)
    unknown = sum(1 for tf, _, _ in selected if tf not in timings)
    for tf, t, reason in selected:
        if why and not files_only:
            print(f"{tf}\t{t}\t{timings.get(tf, float('nan')):.1f}s\t{reason}")
        else:
            print(tf)
    if not files_only:
        for tf, t in deferred:
            print(f"DEFERRED {tf} ({t}, about {timings.get(tf, float('nan')):.0f} s): uses a changed file; runs with --tier {t} before commit")
        note =", manifest changed: run --check too" if "test/manifest.toml" in changed else ""
        print(f"IMPACTED {len(selected)} files, about {total:.0f} s ({unknown} without a timing); "
              f"changed {len(changed)}, impacted sources {len(impacted)}, tier <= {tier}, depth {depth}{note}")


if __name__ == "__main__":
    main(sys.argv[1:])
