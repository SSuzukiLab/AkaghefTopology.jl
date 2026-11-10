# Source: research threads/PJ260922/code/lc_check.py (2026-11-09).
"""Check a Lickorish Code file against docs/specs/lickorish-code.md section 4.

The specification lives in the canonical lib/topology.jl checkout. This is the
deterministic input check only. Passing it does not show that a word transcribes
its source figure, and it says nothing about the drawing.

File format: a first line `genus: <g>`, then one `name : word` line per component,
K first, as in section 3 of the specification. Lines `place`, `disk`, `point`,
and `weight` follow section 3.1. `#` starts a comment.

Exit codes: 0 pass; 1 the data fails a check; 2 the file cannot be read as Lickorish
Data, including a form the specification does not cover.
"""
import json
import re
import sys
from pathlib import Path

META = frozenset({'place', 'disk', 'point', 'weight'})
ROT = [('-', 'out'), ('+', 'out'), ('-', 'in'), ('+', 'in')]
SIDE_RE = re.compile(r'^(.+)([+-])$')
ARC_RE = re.compile(r'^(.+)/([+-])(\d+)$')


class Unsupported(Exception):
    """The input needs something the specification does not define."""


def sign_char(x):
    return '+' if x > 0 else '-'


def rot_next(sign, way):
    return ROT[(ROT.index((sign, way)) + 1) % 4]


def dedupe(items):
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def parse(text):
    genus = None
    words = {}
    places = {}
    disks = {}
    points = {}
    weights = {}
    for raw in text.splitlines():
        line = raw.split('#', 1)[0].strip()
        if not line:
            continue
        if line.startswith('genus:'):
            genus = int(line.split(':', 1)[1])
            continue
        left, sep, rest = line.partition(':')
        if not sep:
            raise Unsupported(f'not a component line: {raw!r}')
        head = left.split()[0]
        if head in META:
            kind = head
            parts = left.split()
            if len(parts) != 2:
                raise Unsupported(f'bad {kind} line: {raw!r}')
            name = parts[1]
            if kind == 'place':
                ref = rest.strip()
                if not ref:
                    raise Unsupported(f'place without region: {raw!r}')
                if name in places:
                    raise Unsupported(f'duplicate place for {name}')
                places[name] = ref
            elif kind == 'disk':
                refs = rest.split()
                if not refs:
                    raise Unsupported(f'disk without regions: {raw!r}')
                if name in disks:
                    raise Unsupported(f'duplicate disk {name}')
                disks[name] = refs
            elif kind == 'point':
                target = rest.strip()
                if not target:
                    raise Unsupported(f'point without target: {raw!r}')
                if name in points:
                    raise Unsupported(f'duplicate point {name}')
                points[name] = target
            elif kind == 'weight':
                try:
                    value = int(rest.strip())
                except ValueError as error:
                    raise Unsupported(f'weight value: {error}') from error
                if name in weights:
                    raise Unsupported(f'duplicate weight on {name}')
                weights[name] = value
            continue
        name = left.strip()
        if name in META:
            raise Unsupported(f'reserved component name {name}')
        if name in words:
            raise Unsupported(f'duplicate component {name}')
        words[name] = [int(x) for x in rest.split()]
    if genus is None:
        raise Unsupported('missing genus line')
    return genus, words, places, disks, points, weights


def k_names(g):
    return [f'l{i}' for i in range(1, g + 1)] + [f'm{i}' for i in range(1, g + 1)] + [f'n{i}' for i in range(1, g)]


def k_constants(g):
    """Constant K words (spec section 4.2, lib/topology.jl reference_surface.jl)."""
    lm = {i: i for i in range(1, g + 1)}
    ln = {}
    label = g
    for i in range(1, g):
        ln[i, i] = label + 1
        ln[i + 1, i] = label + 2
        label += 2
    words = {}
    for i in range(1, g + 1):
        word = [-lm[i]]
        if i < g:
            word.append(-ln[i, i])
        if i > 1:
            word.append(ln[i, i - 1])
        words[f'l{i}'] = word
        words[f'm{i}'] = [lm[i]]
    for i in range(1, g):
        words[f'n{i}'] = [ln[i, i], -ln[i + 1, i]]
    return words


def canon(word):
    return min((tuple(word[i:] + word[:i]) for i in range(len(word))), default=())


def nonempty_words(words):
    return {name: word for name, word in words.items() if word}


def scan_labels(word_items):
    problems = []
    counts = {}
    for name, word in word_items:
        for x in word:
            if x == 0:
                problems.append(f'{name}: label 0 is not allowed')
            else:
                counts.setdefault(abs(x), {'+': 0, '-': 0})[sign_char(x)] += 1
    for label, c in sorted(counts.items()):
        if c != {'+': 1, '-': 1}:
            problems.append(f'label {label} appears {c}')
    return problems, counts


def rotation_occ(words):
    occ = {}
    for name, word in words.items():
        for pos, x in enumerate(word):
            occ.setdefault(abs(x), {})[sign_char(x)] = (name, pos)
    return occ


def along(words, occ, label, sign, way):
    name, pos = occ[label][sign]
    word = words[name]
    x = word[(pos + (1 if way == 'out' else -1)) % len(word)]
    return abs(x), sign_char(x), ('in' if way == 'out' else 'out')


def surface(words):
    """Vertices, edges, faces and connected components of the rotation system."""
    if not words:
        return 0, 0, 0, 0
    occ = rotation_occ(words)
    parent = {k: k for k in occ}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    seen = set()
    faces = 0
    for label in occ:
        for sign, way in ROT:
            z = (label, sign, way)
            if z in seen:
                continue
            faces += 1
            while z not in seen:
                seen.add(z)
                other, s2, w2 = along(words, occ, *z)
                parent[root(z[0])] = root(other)
                s3, w3 = rot_next(s2, w2)
                z = (other, s3, w3)
    return len(occ), 2 * len(occ), faces, len({root(k) for k in occ})


def face_names(words):
    """Face ids F1..Fn in spec section 3.1 order (rotation corners, not pd_code slot order)."""
    if not words:
        return {}, 0
    occ = rotation_occ(words)
    seen = set()
    names = {}
    counter = 1
    for label in sorted(occ):
        for sign, way in ROT:
            start = (label, sign, way)
            if start in seen:
                continue
            fid = f'F{counter}'
            counter += 1
            z = start
            while z not in seen:
                seen.add(z)
                names.setdefault(fid, []).append(z)
                other, s2, w2 = along(words, occ, *z)
                s3, w3 = rot_next(s2, w2)
                z = (other, s3, w3)
    return names, counter - 1


def has_crossings(words):
    return any(words.values())


def curve_sides(words):
    return {f'{name}{ch}' for name in words for ch in '+-'}


def region_set(genus, words, places, disks, face_ids):
    regions = set(face_ids)
    regions |= curve_sides(words)
    regions |= set(disks)
    if genus == 0 and not has_crossings(words):
        regions.add('outer')
    return regions


def resolve_place_chain(comp, places, seen):
    """Follow place nesting through side regions; detect cycles."""
    if comp in seen:
        return True
    seen = seen | {comp}
    parent = places.get(comp)
    if parent is None:
        return False
    m = SIDE_RE.match(parent)
    if not m:
        return False
    return resolve_place_chain(m.group(1), places, seen)


def place_tree_ok(genus, words, places):
    problems = []
    empty = [n for n, w in words.items() if not w]
    for name in words:
        if words[name] and name in places:
            problems.append(f'place on nonempty component {name}')
    for name in empty:
        if name not in places:
            problems.append(f'empty component {name} without place')
    for name in places:
        if name not in words:
            problems.append(f'place names unknown component {name}')
        elif words[name]:
            problems.append(f'place on nonempty component {name}')
    if genus == 0 and not has_crossings(words):
        if not empty:
            problems.append('g=0 without crossings requires empty components')
        for name in empty:
            seen = set()
            if resolve_place_chain(name, places, seen):
                problems.append(f'place containment cycle at {name}')
        for ref in places.values():
            if ref != 'outer' and not SIDE_RE.match(ref):
                problems.append(f'g=0 no crossings: place {ref!r} must reach outer tree')
    else:
        for ref in places.values():
            if ref == 'outer':
                problems.append('outer only when g=0 and no crossings')
    if genus >= 1:
        for ref in places.values():
            if ref == 'outer':
                problems.append('outer not allowed when g>=1')
    for name in places:
        seen = set()
        if resolve_place_chain(name, places, seen):
            problems.append(f'place containment cycle at {name}')
    return problems


def label_in_word(words, component, signed):
    label = abs(signed)
    sign = sign_char(signed)
    if component not in words:
        return False
    for x in words[component]:
        if abs(x) == label and sign_char(x) == sign:
            return True
    return False


def arc_ok(words, target):
    m = ARC_RE.match(target)
    if not m:
        return False, f'not an arc target: {target!r}'
    comp, sign, digits = m.group(1), m.group(2), int(m.group(3))
    if comp not in words or not words[comp]:
        return False, f'arc component {comp} missing or empty'
    signed = digits if sign == '+' else -digits
    if not label_in_word(words, comp, signed):
        return False, f'arc {target}: label {sign}{digits} not in {comp}'
    return True, None


def point_ok(words, target, regions):
    if target in regions:
        return True, None
    if target in words and not words[target]:
        return True, None
    return arc_ok(words, target)


def region_ok(regions, ref):
    return ref in regions


def check_extensions(genus, words, places, disks, points, weights):
    problems = []
    active = nonempty_words(words)
    face_map, nf = face_names(active)
    face_ids = {f'F{i}' for i in range(1, nf + 1)}
    crossings = has_crossings(active)
    regions = region_set(genus, words, places, disks, face_ids)

    problems.extend(place_tree_ok(genus, words, places))

    for comp, ref in places.items():
        if genus >= 1 or crossings:
            if ref == 'outer':
                problems.append(f'place {comp}: outer not allowed')
            elif not region_ok(regions, ref):
                problems.append(f'place {comp}: unknown region {ref!r}')
        elif ref != 'outer' and not region_ok(regions, ref):
            problems.append(f'place {comp}: unknown region {ref!r}')

    if genus == 0 and crossings:
        for comp, ref in places.items():
            if ref == 'outer':
                problems.append('g=0 with crossings: outer not allowed')

    for dname, refs in disks.items():
        for ref in refs:
            if not region_ok(regions, ref):
                problems.append(f'disk {dname}: unknown region {ref!r}')

    for pname, target in points.items():
        ok, reason = point_ok(words, target, regions)
        if not ok:
            problems.append(f'point {pname}: {reason}')

    for wname in weights:
        if not region_ok(regions, wname):
            problems.append(f'weight: unknown region {wname!r}')

    empty_count = sum(1 for w in words.values() if not w)
    if empty_count != len(places):
        problems.append('each empty component needs exactly one place line')

    return problems


def check(genus, words, places, disks, points, weights):
    problems = []
    active = nonempty_words(words)
    crossings = bool(active)

    if genus == 0:
        if crossings:
            problems.extend(check_g0_crossings(words))
        else:
            problems.extend(check_g0_nocross(words, places))
        problems.extend(check_extensions(genus, words, places, disks, points, weights))
        return dedupe(problems)

    empty = [name for name, word in words.items() if not word]
    if not crossings and not empty:
        raise Unsupported('no curve data')

    names = k_names(genus)
    if list(words)[: len(names)] != names:
        problems.append(f'K components must come first in order {names}')

    label_problems, counts = scan_labels(words.items())
    problems.extend(label_problems)

    if problems:
        problems.extend(check_extensions(genus, words, places, disks, points, weights))
        return dedupe(problems)

    limit = 3 * genus - 2
    constants = k_constants(genus)
    for name, expected in constants.items():
        if name not in words:
            continue
        restricted = [x for x in words[name] if abs(x) <= limit]
        if canon(restricted) != canon(expected):
            problems.append(f'{name}: K labels {restricted}, expected cyclic {expected}')
    for label in range(1, limit + 1):
        if crossings and label not in counts:
            problems.append(f'K label {label} is missing')

    if problems:
        problems.extend(check_extensions(genus, words, places, disks, points, weights))
        return dedupe(problems)

    v, e, f, components = surface(active if active else words)
    if active and components != 1:
        problems.append(f'not connected: {components} components')
    if active:
        chi = v - e + f
        if chi != 2 - 2 * genus:
            problems.append(f'Euler characteristic {chi}, expected {2 - 2 * genus}')
        only_k = {name: [x for x in words[name] if abs(x) <= limit] for name in names if name in words}
        if all(only_k.values()):
            fk = surface({k: v for k, v in only_k.items() if v})[2]
            if fk != max(genus, 1):
                problems.append(f'K cuts into {fk} faces, expected {max(genus, 1)}')

    problems.extend(check_extensions(genus, words, places, disks, points, weights))
    return dedupe(problems)


def check_g0_crossings(words):
    problems = []
    active = nonempty_words(words)
    label_problems, _ = scan_labels(active.items())
    problems.extend(label_problems)
    v, e, f, components = surface(active)
    if components != 1:
        problems.append(f'not connected: {components} components')
    chi = v - e + f
    if chi != 2:
        problems.append(f'Euler characteristic {chi}, expected 2')
    return problems


def check_g0_nocross(words, places):
    problems = []
    if has_crossings(words):
        problems.append('g=0 without crossings: all components must be empty')
    return problems


def main(argv):
    if len(argv) != 2:
        print('usage: lc_check.py <file.lc>', file=sys.stderr)
        return 2
    path = Path(argv[1])
    try:
        genus, words, places, disks, points, weights = parse(path.read_text(encoding='utf-8'))
        problems = check(genus, words, places, disks, points, weights)
    except (Unsupported, ValueError) as error:
        print(json.dumps(dict(file=path.name, status='unsupported', reason=str(error)), ensure_ascii=False))
        return 2
    active = nonempty_words(words)
    v, e, f, _ = surface(active) if active else (0, 0, 0, 0)
    print(json.dumps(dict(file=path.name, status='fail' if problems else 'pass', genus=genus,
                          vertices=v, edges=e, faces=f, problems=problems), ensure_ascii=False))
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
