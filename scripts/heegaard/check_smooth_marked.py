"""Accept a smooth display only after checks independent of its objective."""
import json
import hashlib
from pathlib import Path
import sys

import numpy as np
from shapely.geometry import LineString, Point, MultiPoint
from shapely.strtree import STRtree

from measure_clearance import measure
from publish_marked_display import bind
from draw_marked_cut import render
from verify_fixed_grid import verify
from measure_clearance import branches


def split(q):
    a=(q[:-1]+q[1:])/2;b=(a[:-1]+a[1:])/2;c=(b[0]+b[1])/2
    return np.array([q[0],a[0],b[0],c]),np.array([c,b[1],a[2],q[3]])


def leaves(q, tolerance=.003):
    error=max(np.linalg.norm(q[1]-(2*q[0]+q[3])/3),np.linalg.norm(q[2]-(q[0]+2*q[3])/3))
    if error<=tolerance:
        yield q
    else:
        for half in split(q):yield from leaves(half,tolerance)


def check(d, old, prepared=None, base=None, *, minimum_separation=12., puncture_clearance=0., crossing_radius=20.):
    result=bind(d,prepared,base) if prepared is not None else dict(d,geometry_check=verify(d,old))
    for x in old['crossings']:
        assert np.linalg.norm(np.array(d['positions'][str(x['vertex'])])-old['positions'][str(x['vertex'])])<1e-12
    for a,b in zip(d['boundaries'],old['boundaries']):
        assert a==b
        assert d['positions'][str(a['cap'])]==old['positions'][str(b['cap'])]
        for port in a['ports']:
            v=str(port['vertex'])
            assert np.linalg.norm(np.array(d['positions'][v])-old['positions'][v])<1e-12
    hulls=[];qs=[];ids=[];rays={};max_join_error=0
    for curve in d['curves']:
        for s in curve['segments']:
            pieces=[np.array(q) for q in s['beziers']]
            vs=s['vertices']
            flattened=[pieces[0][0]]+[leaf[-1] for q in pieces for leaf in leaves(q)]
            actual=np.array([d['positions'][str(v)] for v in vs])*600
            assert np.shape(flattened)==np.shape(actual) and np.max(np.abs(flattened-actual))<1e-7,'Spline/sample mismatch'
            for end,direction in ((vs[0],pieces[0][1]-pieces[0][0]),(vs[-1],pieces[-1][-2]-pieces[-1][-1])):
                rays.setdefault((curve['id'],end),[]).append(direction/np.linalg.norm(direction))
            for a,b in zip(pieces,pieces[1:]):
                assert np.linalg.norm(a[-1]-b[0])<1e-7
                va=a[-1]-a[-2];vb=b[1]-b[0]
                max_join_error=max(max_join_error,np.linalg.norm(va/np.linalg.norm(va)-vb/np.linalg.norm(vb)))
            for q in pieces:
                for leaf in leaves(q):
                    assert np.all(np.diff(leaf,axis=0)@(leaf[-1]-leaf[0])>0),'Non-regular spline piece'
                    hulls.append(MultiPoint(leaf).convex_hull);qs.append(leaf);ids.append(curve['id'])
    # Stronger than sampling: a cubic is contained in its control hull.
    holes=[(np.array(d['positions'][str(b['cap'])])*600,b['radius']*600) for b in d['boundaries']]
    min_hole=min((h.distance(Point(c))-r for h in hulls for c,r in holes),default=float('inf'))
    assert min_hole>=-1e-6,('Bezier control hull enters hole',min_hole)
    tree=STRtree(hulls)
    for i,hull in enumerate(hulls):
        for j in tree.query(hull,predicate='intersects'):
            if j<=i:continue
            q,r=qs[i],qs[j]
            common=[p for p in (q[0],q[-1]) for s in (r[0],r[-1]) if np.linalg.norm(p-s)<1e-7]
            hit=hull.intersection(hulls[j])
            assert common and all(hit.hausdorff_distance(Point(p))<1e-5 for p in common),('Bezier hull collision',i,int(j),ids[i],ids[j],hit.wkt)
    for vectors in rays.values():
        if len(vectors)==2:max_join_error=max(max_join_error,np.linalg.norm(vectors[0]+vectors[1]))
    for hole in d['boundaries']:
        center=np.array(d['positions'][str(hole['cap'])])
        for port in hole['ports']:
            vertex=port['vertex'];radial=np.array(d['positions'][str(vertex)])-center;radial/=np.linalg.norm(radial)
            direction=next(vs[0] for (curve,v),vs in rays.items() if v==vertex)
            assert np.linalg.norm(direction-radial)<1e-7,'Port tangent is not normal to the hole'
    assert max_join_error<1e-7,('Non-smooth curve join',max_join_error)
    if d.get('disk'):
        disk=[np.array(d['disk_beziers'][0][0])]+[leaf[-1] for q in d['disk_beziers'] for leaf in leaves(np.array(q))]
        actual=np.array([d['positions'][str(v)] for v in d['disk']+[d['disk'][0]]])*600
        assert np.shape(disk)==np.shape(actual) and np.max(np.abs(disk-actual))<1e-7,'D spline/sample mismatch'
    def disk_word(view):
        lookup={}
        for c in view['curves']:
            for s in c['segments']:
                for u,v in zip(s['vertices'],s['vertices'][1:]):
                    lookup[u,v]=(c['id'],tuple(s['source']),1);lookup[v,u]=(c['id'],tuple(s['source']),-1)
        word=[];vs=view.get('disk',[])
        for edge in zip(vs,vs[1:]+vs[:1]):
            item=lookup[edge]
            if not word or word[-1]!=item:word.append(item)
        if len(word)>1 and word[0]==word[-1]:word.pop()
        return min((tuple(word[k:]+word[:k]) for k in range(len(word))),default=())
    assert disk_word(d)==disk_word(old),'D source boundary changed'
    positions,_,branch_list=branches(d);hole_gaps=[]
    for branch in branch_list:
        for hole in d['boundaries']:
            line=branch['line']
            for port in hole['ports']:
                if port['vertex'] in (branch['vertices'][0],branch['vertices'][-1]):
                    line=line.difference(Point(positions[port['vertex']]).buffer(19.99,quad_segs=64))
            if not line.is_empty:hole_gaps.append(line.distance(Point(positions[hole['cap']]))-600*hole['radius'])
    nonport_margin=min(hole_gaps,default=float('inf'))
    assert nonport_margin>=puncture_clearance+.006-1e-6 if puncture_clearance else nonport_margin>=-1e-6,('Non-port hole clearance',nonport_margin,puncture_clearance)
    if d.get('boundary_vertex_pairs'):
        portmap={p['source_vertex']:(b,p) for b in d['boundaries'] for p in b['ports']}
        for u,v in d['boundary_vertex_pairs']:
            if u not in portmap and v not in portmap:continue
            assert u in portmap and v in portmap,'Missing paired port'
            a,pa=portmap[u];b,pb=portmap[v]
            va=positions[pa['vertex']]-positions[a['cap']];vb=positions[pb['vertex']]-positions[b['cap']]
            assert np.linalg.norm(va-vb*np.array([1,-1]))<1e-6,'Paired port angle mismatch'
    clearance=measure(d,threshold=100)
    minimum=clearance['pairs'][0]['gap'] if clearance['pairs'] else None
    # Exempt at radius 20.01 rather than 20 to cover flatness at the cut circle.
    certified=measure(d,threshold=minimum_separation+.02,crossing_radius=crossing_radius-.01)
    assert not certified['violations'],certified
    result['clearance_check']=measure(d)
    result['smooth_check']=dict(crossing_coordinates_unchanged=len(d['crossings']),source_disk_boundary_unchanged=True,
        holes_and_ports_unchanged=True,bezier_control_hulls_checked=len(hulls),
        minimum_hole_clearance=min_hole,minimum_nonport_hole_clearance=nonport_margin,minimum_branch_spacing_svg=minimum,
        polyline_flatness_bound_svg=.003,maximum_tangent_join_error=max_join_error,
        no_extra_bezier_intersections=True,regular_spline_pieces=True,normal_port_tangents=True)
    result['geometry_provenance']='global cubic B-spline fairing; source intervals and D retained; Bezier control hulls checked'
    return result


if __name__=='__main__':
    root=Path(sys.argv[1]);read=lambda n:json.loads((root/n).read_text(encoding='utf-8'))
    result=check(read('KodaSmoothCandidate.json'),read('KodaSeparated.json'),read('KodaPrepared.json'),read('KodaMarkedCut.json'))
    isotopy=read('KodaSmoothIsotopy.json')
    for key,name in (('source_sha256','KodaSeparated.json'),('candidate_sha256','KodaSmoothCandidate.json')):
        assert isotopy[key]==hashlib.sha256((root/name).read_bytes()).hexdigest(),'Stale isotopy check'
    result['fairing_isotopy_check']=isotopy
    (root/'KodaSmooth.json').write_text(json.dumps(result),encoding='utf-8')
    (root/'KodaSmooth.svg').write_text(render(result),encoding='utf-8')
    print(json.dumps(result['smooth_check'],indent=2))
