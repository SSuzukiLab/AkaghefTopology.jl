"""Draw one Lickorish configuration by calling the existing placement.

cf: lib/topology.jl/scripts/heegaard/lc_placement.py place — 穴・Tutte・放射配置。ここは座標を計算しない
"""
import argparse
import json
import os
import re
from pathlib import Path

from lc_draw import read_lc
from lc_placement import place
from lc_render import write_svg

_DIGITS = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
_PATH_TOKEN = re.compile(r'[MLHVZmlhvz]|[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')


def _bank_label(index, sign):
    return 'l' + str(index).translate(_DIGITS) + sign


def present(doc, roles):
    """Attach roles and ask the renderer for the cut already stored on doc."""
    roles = {str(key): str(value) for key, value in (roles or {}).items()}
    doc['show_reference'] = True
    doc['roles'] = roles
    labels = dict(doc.get('labels') or {})
    for name in doc.get('holes', {}):
        if name.startswith('top'):
            labels[name] = _bank_label(name[3:], '⁺')
        elif name.startswith('bot'):
            labels[name] = _bank_label(name[3:], '⁻')
    doc['labels'] = labels
    for key, curve in doc['curves'].items():
        support = key.split('#', 1)[0]
        curve['support_id'] = support
        if support in roles:
            curve['role'] = roles[support]
    return doc


def draw_diagram(lc_path, svg_path, roles=None, png=False, title='Heegaard diagram'):
    doc, _, _ = place(read_lc(lc_path))
    present(doc, roles)
    svg_path = Path(svg_path)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    write_svg(doc, svg_path, title=title)
    if png:
        rasterize_svg(svg_path, svg_path.with_suffix('.png'))
    return svg_path


def _path_commands(d):
    tokens = _PATH_TOKEN.findall(d)
    index = 0
    command = None
    cx = cy = 0.0
    start = (0.0, 0.0)
    while index < len(tokens):
        token = tokens[index]
        if token in 'MLHVZmlhvz':
            command = token
            index += 1
        if command in 'Zz':
            yield ('L', start[0], start[1])
            cx, cy = start
            continue
        if command in 'Mm':
            x, y = float(tokens[index]), float(tokens[index + 1])
            index += 2
            if command == 'm':
                x += cx
                y += cy
            cx, cy = x, y
            start = (cx, cy)
            yield ('M', cx, cy)
        elif command in 'Ll':
            x, y = float(tokens[index]), float(tokens[index + 1])
            index += 2
            if command == 'l':
                x += cx
                y += cy
            cx, cy = x, y
            yield ('L', cx, cy)
        elif command in 'Hh':
            x = float(tokens[index])
            index += 1
            if command == 'h':
                x += cx
            cx = x
            yield ('L', cx, cy)
        elif command in 'Vv':
            y = float(tokens[index])
            index += 1
            if command == 'v':
                y += cy
            cy = y
            yield ('L', cx, cy)
        else:
            raise ValueError(f'unsupported SVG path command {command!r}')


def _attr(text, name):
    match = re.search(rf'\b{name}="([^"]*)"', text)
    return match.group(1) if match else None


def _font(size):
    from PIL import ImageFont
    for path in (r'C:\Windows\Fonts\segoeui.ttf', r'C:\Windows\Fonts\calibri.ttf', r'C:\Windows\Fonts\arial.ttf'):
        if os.path.exists(path):
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def _draw_chain(draw, points, fill, width, pattern):
    if len(points) < 2:
        return
    if not pattern:
        draw.line(points, fill=fill, width=width, joint='curve')
        return
    drawing = True
    remain = pattern[0]
    slot = 0
    cursor = points[0]
    for target in points[1:]:
        while True:
            dx = target[0] - cursor[0]
            dy = target[1] - cursor[1]
            dist = (dx * dx + dy * dy) ** 0.5
            if dist < 1e-6:
                cursor = target
                break
            step = min(remain, dist)
            nxt = (cursor[0] + dx / dist * step, cursor[1] + dy / dist * step)
            if drawing:
                draw.line([cursor, nxt], fill=fill, width=width)
            cursor = nxt
            remain -= step
            if remain <= 1e-6:
                slot = (slot + 1) % len(pattern)
                remain = pattern[slot]
                drawing = slot % 2 == 0
            if step >= dist - 1e-6:
                cursor = target
                break


def rasterize_svg(svg_path, png_path, width=1400):
    """Rasterize this renderer's polylines with Pillow. No new coordinates."""
    from PIL import Image, ImageDraw
    text = Path(svg_path).read_text(encoding='utf-8')
    vx, vy, vw, vh = (float(part) for part in _attr(text, 'viewBox').split())
    oversample = 2
    pixel_w = width * oversample
    pixel_h = max(1, round(pixel_w * vh / vw))
    image = Image.new('RGB', (pixel_w, pixel_h), '#fcfbf7')
    draw = ImageDraw.Draw(image)
    scale = pixel_w / vw

    def xy(x, y):
        return ((x - vx) / vw * (pixel_w - 1), (y - vy) / vh * (pixel_h - 1))

    def paint_path(attrs):
        stroke = _attr(attrs, 'stroke')
        path = _attr(attrs, 'd')
        if not stroke or stroke == 'none' or not path:
            return
        stroke_width = float(_attr(attrs, 'stroke-width') or 2)
        dash = _attr(attrs, 'stroke-dasharray')
        pattern = [float(part) * scale for part in dash.split()] if dash else None
        chains = []
        current = []
        for op, x, y in _path_commands(path):
            point = xy(x, y)
            if op == 'M' and current:
                chains.append(current)
                current = [point]
            else:
                current.append(point)
        if current:
            chains.append(current)
        pen = max(1, round(stroke_width * scale))
        for chain in chains:
            _draw_chain(draw, chain, stroke, pen, pattern)

    def paint_circle(attrs):
        if _attr(attrs, 'r') is None:
            return
        cx = float(_attr(attrs, 'cx'))
        cy = float(_attr(attrs, 'cy'))
        radius = float(_attr(attrs, 'r'))
        stroke = _attr(attrs, 'stroke') or '#333333'
        fill = _attr(attrs, 'fill')
        x0, y0 = xy(cx - radius, cy - radius)
        x1, y1 = xy(cx + radius, cy + radius)
        box = [x0, y0, x1, y1]
        if fill and fill != 'none':
            draw.ellipse(box, fill=fill)
        if stroke and stroke != 'none':
            pen = max(1, round(float(_attr(attrs, 'stroke-width') or 1.5) * scale))
            draw.ellipse(box, outline=stroke, width=pen)

    font = _font(max(12, round(18 * scale)))

    def paint_text(attrs, body):
        x, y = xy(float(_attr(attrs, 'x')), float(_attr(attrs, 'y')))
        fill = _attr(attrs, 'fill') or '#253247'
        if _attr(attrs, 'text-anchor') == 'middle':
            box = font.getbbox(body)
            x -= (box[2] - box[0]) / 2
        draw.text((x, y), body, fill=fill, font=font, anchor='ls')

    for match in re.finditer(r'<(path|circle)\b([^>]*)/>|<text\b([^>]*)>([^<]*)</text>', text):
        if match.group(1) == 'path':
            paint_path(match.group(2))
        elif match.group(1) == 'circle':
            paint_circle(match.group(2))
        else:
            paint_text(match.group(3), match.group(4))
    image = image.resize((width, max(1, round(pixel_h / oversample))), Image.Resampling.LANCZOS)
    image.save(png_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('lc')
    parser.add_argument('svg')
    parser.add_argument('--roles', type=Path)
    parser.add_argument('--png', action='store_true')
    parser.add_argument('--title', default='Heegaard diagram')
    args = parser.parse_args()
    roles = json.loads(args.roles.read_text(encoding='utf-8')) if args.roles else {}
    draw_diagram(args.lc, args.svg, roles=roles, png=args.png, title=args.title)


if __name__ == '__main__':
    main()
