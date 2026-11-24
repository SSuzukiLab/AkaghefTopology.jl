"""Polyline-only SVG output for Lickorish Code placement documents.

Source: owner-approved app/index.html::draw; mock3.py::SVG export.
The rendering bounds and colors follow the reference; no curve interpolation.
"""
# cf: lib/topology.jl/docs/conventions/LickorishCodes.md — placement is displayed in screen coordinates
from html import escape
from pathlib import Path
import subprocess


def svg(data, title='Lickorish Code', dots=False):
    # cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/docs/DESIGN.md 現行ソースを見て必要になった修正 — show_reference では配置済みの K を描き、色は role
    w,h=data['W'],data['H']
    # Export the entire computed polyline. The interactive mock's fixed canvas
    # bounds can crop an outer arc; enlarging only the viewport preserves data.
    pts=[(x,y) for c in data['curves'].values() for x,y in c['points']]
    pts += [(x-r,y-r) for (x,y),r in data['holes'].values()]
    pts += [(x+r,y+r) for (x,y),r in data['holes'].values()]
    if data.get('show_reference'):
        for pieces in data.get('reference_curves', {}).values():
            for piece in pieces:
                pts.extend(piece['points'])
    xmin=min([-(w*.25)]+[p[0] for p in pts])-40
    ymin=min([-(h*.25)]+[p[1] for p in pts])-40
    xmax=max([w*1.25]+[p[0] for p in pts])+40
    ymax=max([h*1.25]+[p[1] for p in pts])+40
    vb=(xmin,ymin,xmax-xmin,ymax-ymin)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb[0]} {vb[1]} {vb[2]} {vb[3]}">',
         f'<rect x="{vb[0]}" y="{vb[1]}" width="{vb[2]}" height="{vb[3]}" fill="#fcfbf7"/>',
         f'<title>{escape(title)}</title>']
    if data.get('square'):
        if not (data.get('show_reference') and data.get('reference_curves')):
            out.append(f'<path d="M0,0 H{w} V{h} H0 Z" fill="none" stroke="#555" stroke-width="1.5"/>')
        family='Segoe UI, Arial, sans-serif' if data.get('show_reference') else 'Arial'
        for x,y,label in [(w/2,-12,'ℓ₁⁻ →'),(w/2,h+24,'ℓ₁⁺ →'),(-28,h/2,'m₁ ↑'),(w+28,h/2,'m₁ ↑')]:
            out.append(f'<text x="{x}" y="{y}" text-anchor="middle" font-family="{family}" font-size="14">{label}</text>')
    for name,(c,r) in data['holes'].items():
        out.append(f'<circle cx="{c[0]}" cy="{c[1]}" r="{r}" fill="#fcfbf7" stroke="#333" stroke-width="1.5"/>')
    if data.get('show_reference'):
        roles=data.get('roles') or {}
        for support, pieces in data.get('reference_curves', {}).items():
            role=roles.get(support)
            if role=='alpha':
                color,width,extra='#ba433f','2.2',''
            elif role=='beta':
                color,width,extra='#255dba','2.2',''
            elif str(support).startswith('l'):
                color,width,extra='#39444f','2.0',''
            else:
                color,width,extra='#7a6a9c','1.4',' stroke-dasharray="6 4"'
            kind='display-cut' if str(support).startswith('l') else 'reference'
            for piece in pieces:
                path='M'+' L'.join(f'{x},{y}' for x,y in piece['points'])+(' Z' if piece.get('closed') else '')
                out.append(
                    f'<path data-support="{escape(str(support), quote=True)}" data-kind="{kind}" '
                    f'd="{path}" fill="none" stroke="{color}" stroke-width="{width}"{extra}/>')
    for name,curve in data['curves'].items():
        role=curve.get('role')
        if role=='alpha':
            color='#ba433f'
        elif role=='beta':
            color='#255dba'
        else:
            color='#ba433f' if name.startswith(('a','α')) else '#255dba'
        path='M'+' L'.join(f'{x},{y}' for x,y in curve['points'])+(' Z' if curve['closed'] else '')
        if data.get('show_reference') or role:
            support=escape(str(curve.get('support_id', name.split('#',1)[0])), quote=True)
            out.append(
                f'<path data-support="{support}" data-role="{escape(role or "", quote=True)}" '
                f'd="{path}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        else:
            out.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        if dots:
            out.extend(f'<circle cx="{x}" cy="{y}" r="1.8" fill="{color}"/>' for x,y in curve['points'])
    for name,(c,r) in data['holes'].items():
        default = ('ℓ'+name[3:]+('⁺' if name.startswith('top') else '⁻')
                   if name.startswith(('top','bot')) and name[3:].isdigit() else name)
        label=escape(data.get('labels',{}).get(name,default))
        family='Segoe UI, Arial, sans-serif' if data.get('show_reference') else 'Arial'
        out.append(f'<text x="{c[0]}" y="{c[1]+4}" text-anchor="middle" font-family="{family}" font-size="12" fill="#253247">{label}</text>')
    out.append('</svg>')
    return '\n'.join(out)


def write_svg(data,path,**kwargs):
    Path(path).write_text(svg(data,**kwargs),encoding='utf-8')


def export_png(svg_path,png_path,*,inkscape='inkscape',width=1500):
    """Use installed Inkscape for SVG rasterization; no image synthesis."""
    subprocess.run([inkscape,str(svg_path),'--export-type=png',f'--export-filename={png_path}',
                    f'--export-width={width}'],check=True,capture_output=True)
