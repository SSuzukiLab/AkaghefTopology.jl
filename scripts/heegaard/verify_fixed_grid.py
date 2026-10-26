"""Check the displayed fixed-grid paths, independently of the layout objective."""
import json,math,sys
from pathlib import Path
import numpy as np
from shapely.geometry import LineString,Point,Polygon
from shapely.strtree import STRtree


def verify(d,source):
    p={int(v):np.round(np.array(x)*600,7) for v,x in d['positions'].items()}
    groups={}
    for b in d['boundaries']:
        center=p[b['cap']];groups.setdefault(b['cut'],[]).append(center)
        assert b['radius']>0
    assert len(d['boundaries'])==2*d['genus'] and len(groups)==d['genus']
    radii=[b['radius'] for b in d['boundaries']];assert not radii or max(radii)-min(radii)<1e-10
    # Source: layout_marked_document.coarse_view and marked_display.jl source
    # binding. Metadata is excluded from the mathematical digest, so infinity
    # and the exact D source partition need independent equality checks.
    for key in ('display_request','infinity_source_triangle','disk_region_marks','disk_pieces'):
        assert d.get(key)==source.get(key),('Source display selection changed',key)
    for i,label in enumerate(d['cuts']):
        pair=sorted(groups[label],key=lambda c:c[1])
        assert np.max(abs(np.array(pair)-np.array([[(i-(d['genus']-1)/2)*390,-250],[(i-(d['genus']-1)/2)*390,250]])))<1e-6
    assert d['crossings']==source['crossings'],'Crossing marks changed'
    for c,old in zip(d['curves'],source['curves']):
        assert c['id']==old['id']
        assert [s['source'] for s in c['segments']]==[s['source'] for s in old['segments']]
    for x in d['crossings']:
        rays={};v=x['vertex']
        for c in d['curves']:
            for s in c['segments']:
                vs=s['vertices']
                if vs[0]==v:rays[c['id']+'+']=p[vs[1]]-p[v]
                if vs[-1]==v:rays[c['id']+'-']=p[vs[-2]]-p[v]
        assert len(rays)==4,(x['id'],rays)
        order=sorted(rays,key=lambda k:math.atan2(rays[k][1],rays[k][0]))
        a,b=x['curves'];k=order.index(a+'+');order=order[k:]+order[:k]
        expected=[a+'+',b+('+' if x['sign']==1 else '-'),a+'-',b+('-' if x['sign']==1 else '+')]
        assert order==expected,(x['id'],order,expected)
    es=[(c['id'],u,v) for c in d['curves'] for s in c['segments'] for u,v in zip(s['vertices'],s['vertices'][1:])]
    ls=[LineString([p[u],p[v]]) for _,u,v in es];tree=STRtree(ls)
    # Source: prepare_marked_display's alpha-union-beta region semantics;
    # displayed auxiliary reference curves may cross the interior of D.
    mathematical={c['id'] for c in d['curves'] if c['role'] in ('alpha','beta')}
    mathematical_lines=[line for (cid,_,_),line in zip(es,ls) if cid in mathematical]
    for i,line in enumerate(ls):
        for j in tree.query(line,predicate='intersects'):
            if j<=i:continue
            common=set(es[i][1:])&set(es[j][1:]);hit=line.intersection(ls[j])
            assert common and hit.geom_type=='Point' and any(hit.distance(Point(p[v]))<1e-5 for v in common),('Extra intersection',i,int(j),es[i],es[j],hit.wkt)
        for b in d['boundaries']:
            assert line.distance(Point(p[b['cap']]))>=600*b['radius']-1e-5,('Hole entry',i,b['cut'])
    area=None
    if d.get('disk'):
        polygon=Polygon([p[v] for v in d['disk']]);assert polygon.is_valid and polygon.area>1e-8
        assert polygon.exterior.is_ccw
        for b in d['boundaries']:assert polygon.distance(Point(p[b['cap']]))>=600*b['radius']-1e-5
        interior=polygon.buffer(-1e-4)
        assert not any(interior.intersects(line) for line in mathematical_lines),'Mathematical curve in D interior'
        area=polygon.area
    from draw_marked_cut import disk_piece_cycles
    from shapely.geometry import box
    piece_areas=[]
    world=box(min(q[0] for q in p.values())-10,min(q[1] for q in p.values())-10,
        max(q[0] for q in p.values())+10,max(q[1] for q in p.values())+10)
    fills=[]
    for piece,cycles in disk_piece_cycles(d):
        fill=world if piece['contains_infinity'] else Polygon()
        signed=sum(float(np.sum(c[:,0]*np.roll(c[:,1],-1)-c[:,1]*np.roll(c[:,0],-1)))/2 for c in cycles)
        assert signed<0 if piece['contains_infinity'] else signed>0,'D piece orientation/infinity mismatch'
        for cycle in cycles:
            polygon=Polygon(cycle)
            assert polygon.is_valid and polygon.area>1e-8,'Invalid D piece boundary'
            fill=fill.symmetric_difference(polygon)
        for b in d['boundaries']:
            fill=fill.difference(Point(p[b['cap']]).buffer(600*b['radius'],quad_segs=1024))
        assert not fill.is_empty and fill.is_valid,'Empty or invalid D fill'
        assert not any(fill.buffer(-.01).intersects(line) for line in mathematical_lines),'Mathematical curve in D piece interior'
        assert not any(fill.intersection(old).area>.01 for old in fills),'Overlapping D pieces'
        fills.append(fill)
        # The infinity piece is unbounded in this planar chart; the temporary
        # verification window supplies no intrinsic finite area for that piece.
        piece_areas.append(None if piece['contains_infinity'] else fill.area)
    return dict(fixed_grid=f'2x{d["genus"]}',equal_radius=True,source_edges_preserved=True,signed_crossings=len(d['crossings']),extra_intersections=0,hole_entries=0,joining_disk_area=area,joining_disk_piece_areas=piece_areas)


if __name__=='__main__':
    path=Path(sys.argv[1]);data=json.loads(path.read_text(encoding='utf-8'))
    source=json.loads(path.with_name('KodaMarkedCut.json').read_text(encoding='utf-8'))
    print(json.dumps(verify(data,source),ensure_ascii=False))
