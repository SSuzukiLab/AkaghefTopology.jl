"""Polyline-only SVG output for mock-compatible y-up placement documents.

Source: owner-approved app/index.html::draw; mock3.py::SVG export.
The rendering bounds and colors follow the reference; no curve interpolation.
"""
from html import escape
from pathlib import Path
import subprocess


def svg(data, title='Lickorish Data', dots=False):
    w,h=data['W'],data['H']
    # Export the entire computed polyline. The interactive mock's fixed canvas
    # bounds can crop an outer arc; enlarging only the viewport preserves data.
    pts=[(x,h-y) for c in data['curves'].values() for x,y in c['points']]
    pts += [(x-r,h-y-r) for (x,y),r in data['holes'].values()]
    pts += [(x+r,h-y+r) for (x,y),r in data['holes'].values()]
    xmin=min([-(w*.25)]+[p[0] for p in pts])-40
    ymin=min([-(h*.25)]+[p[1] for p in pts])-40
    xmax=max([w*1.25]+[p[0] for p in pts])+40
    ymax=max([h*1.25]+[p[1] for p in pts])+40
    vb=(xmin,ymin,xmax-xmin,ymax-ymin)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb[0]} {vb[1]} {vb[2]} {vb[3]}">',
         f'<rect x="{vb[0]}" y="{vb[1]}" width="{vb[2]}" height="{vb[3]}" fill="#fcfbf7"/>',
         f'<title>{escape(title)}</title>']
    if data.get('square'):
        out.append(f'<path d="M0,0 H{w} V{h} H0 Z" fill="none" stroke="#555" stroke-width="1.5"/>')
        for x,y,label in [(w/2,-12,'l₁ →'),(w/2,h+24,'l₁ →'),(-28,h/2,'m₁ ↓'),(w+28,h/2,'m₁ ↓')]:
            out.append(f'<text x="{x}" y="{y}" text-anchor="middle" font-family="Arial" font-size="14">{label}</text>')
    for name,(c,r) in data['holes'].items():
        out.append(f'<circle cx="{c[0]}" cy="{h-c[1]}" r="{r}" fill="#fcfbf7" stroke="#333" stroke-width="1.5"/>')
    for name,curve in data['curves'].items():
        color='#ba433f' if name.startswith(('a','α')) else '#255dba'
        path='M'+' L'.join(f'{x},{h-y}' for x,y in curve['points'])+(' Z' if curve['closed'] else '')
        out.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        if dots:
            out.extend(f'<circle cx="{x}" cy="{h-y}" r="1.8" fill="{color}"/>' for x,y in curve['points'])
    for name,(c,r) in data['holes'].items():
        label=escape(data.get('labels',{}).get(name,name))
        out.append(f'<text x="{c[0]}" y="{h-c[1]+4}" text-anchor="middle" font-family="Arial" font-size="12" fill="#253247">{label}</text>')
    out.append('</svg>')
    return '\n'.join(out)


def write_svg(data,path,**kwargs):
    Path(path).write_text(svg(data,**kwargs),encoding='utf-8')


def export_png(svg_path,png_path,*,inkscape='inkscape',width=1500):
    """Use installed Inkscape for SVG rasterization; no image synthesis."""
    subprocess.run([inkscape,str(svg_path),'--export-type=png',f'--export-filename={png_path}',
                    f'--export-width={width}'],check=True,capture_output=True)
