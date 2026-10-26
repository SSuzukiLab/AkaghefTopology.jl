"""Measure non-crossing branch clearance on the actual displayed polylines.

At an intended common crossing, only a small explicitly reported disk is exempt.
Consecutive polyline edges are merged into branches so mesh sampling density does
not manufacture 'near approaches'. Units match the SVG (600 * layout units).
"""
import json
import sys
from pathlib import Path
import numpy as np
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points


def branches(d):
    positions={int(v):np.asarray(p)*600 for v,p in d['positions'].items()}
    crossings={x['vertex']:x['id'] for x in d['crossings']}
    out=[]
    for curve in d['curves']:
        nxt={};prev={}
        for s in curve['segments']:
            for u,v in zip(s['vertices'],s['vertices'][1:]):
                assert u not in nxt and v not in prev
                nxt[u]=v;prev[v]=u
        stops=set(crossings)|set(nxt).symmetric_difference(prev)
        unused=set(nxt)
        semantic={s['vertices'][0] for s in curve['segments']}
        while unused:
            start=min(unused&stops) if unused&stops else min(unused&semantic) if unused&semantic else min(unused)
            path=[start];u=start
            while u in unused:
                unused.remove(u);u=nxt[u];path.append(u)
                if u in stops or u==start:break
            out.append(dict(curve=curve['id'],vertices=path,
                ends={v for v in (path[0],path[-1]) if v in crossings},
                line=LineString([positions[v] for v in path])))
    return positions,crossings,out


def measure(d,threshold=12.,crossing_radius=20.):
    p,xs,bs=branches(d);reports=[]
    for i,a in enumerate(bs):
        for j in range(i):
            b=bs[j];aa,bb=a['line'],b['line']
            for v in a['ends']&b['ends']:
                disk=Point(p[v]).buffer(crossing_radius,quad_segs=32)
                aa=aa.difference(disk);bb=bb.difference(disk)
            if aa.is_empty or bb.is_empty:continue
            gap=aa.distance(bb)
            if gap>=threshold:continue
            u,v=nearest_points(aa,bb)
            reports.append(dict(gap=float(gap),branches=[i,j],curves=[a['curve'],b['curve']],
                nearby_crossings=[[xs[v] for v in sorted(a['ends'])],[xs[v] for v in sorted(b['ends'])]],
                points=[list(u.coords[0]),list(v.coords[0])]))
    reports.sort(key=lambda r:r['gap'])
    return dict(units='SVG units',target=threshold,crossing_exemption_radius=crossing_radius,
                branch_count=len(bs),violations=len(reports),pairs=reports)


if __name__=='__main__':
    d=json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
    print(json.dumps(measure(d),ensure_ascii=False,indent=2))
