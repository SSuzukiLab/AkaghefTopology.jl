"""Continuous-time collision check for the old-to-faired PL interpolation.

All orientation events of moving line segments are quadratic in time. Test
their roots and the intervening intervals, not just animation frames. This is
a floating-point geometric check, not a formal/exact-arithmetic theorem.
"""
import json
import hashlib
from pathlib import Path
import sys
import numpy as np
from shapely.geometry import MultiPoint, Point, LineString
from shapely.strtree import STRtree
from measure_clearance import branches


def split(q):
    a=(q[:-1]+q[1:])/2;b=(a[:-1]+a[1:])/2;c=(b[0]+b[1])/2
    return np.array([q[0],a[0],b[0],c]),np.array([c,b[1],a[2],q[3]])


def sample(q,lo,hi,out):
    err=max(np.linalg.norm(q[1]-(2*q[0]+q[3])/3),np.linalg.norm(q[2]-(q[0]+2*q[3])/3))
    if err<=.003:out.append((hi,q[3]));return
    a,b=split(q);mid=(lo+hi)/2;sample(a,lo,mid,out);sample(b,mid,hi,out)


def _parameter_reference_binding(geometry,reference):
    """Record original parameter labels separately from corrected PL coordinates."""
    pos,_,bs=branches(reference);grids=[]
    for branch in bs:
        q=np.array([pos[v] for v in branch['vertices']])
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(q,axis=0),axis=1))];arc/=arc[-1]
        assert np.isfinite(arc).all() and arc[0]==0. and arc[-1]==1. and np.all(np.diff(arc)>0)
        grids.append(dict(curve=branch['curve'],vertices=branch['vertices'],parameters=arc.tolist()))
    encode=lambda value:json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return dict(original_parameter_grid=grids,
        original_parameter_grid_sha256=hashlib.sha256(encode(grids)).hexdigest(),
        original_parameter_positions_sha256=hashlib.sha256(encode(reference['positions'])).hexdigest(),
        corrected_polyline_positions_sha256=hashlib.sha256(encode(geometry['positions'])).hexdigest())


def correspondence(old,new, *, with_branch_signatures=False, _parameter_source=None, with_source_parameters=False):
    pos,xs,bs=branches(old);starts=[];ends=[];edges=[];ids=[];signatures=[];source_parameters=[];nextid=max(pos)+1
    # Source: build_exact_anchor_beziers changes coordinates, not source vertex order.
    # A retained strictly monotone source grid parameterizes that corrected PL.
    reference=None
    if _parameter_source is not None:
        for key in ('source_document','display_request'):
            assert old.get(key)==_parameter_source.get(key), ('Parameter source identity',key)
        for curve,prior in zip(old['curves'],_parameter_source['curves'],strict=True):
            for key in ('id','role','closed','word','seam_word'):
                assert curve.get(key)==prior.get(key), ('Parameter source curve',key)
            for segment,prior_segment in zip(curve['segments'],prior['segments'],strict=True):
                for key in ('source','edge','edge_sign','vertices'):
                    assert segment.get(key)==prior_segment.get(key), ('Parameter source segment',key)
        rp,rx,rb=branches(_parameter_source)
        assert xs==rx and len(bs)==len(rb)
        for branch,prior in zip(bs,rb,strict=True):
            assert branch['curve']==prior['curve'] and branch['vertices']==prior['vertices']
        for curve,prior in zip(new['curves'],old['curves'],strict=True):
            assert curve['id']==prior['id']
            for segment,prior_segment in zip(curve['segments'],prior['segments'],strict=True):
                for key in ('source','edge','edge_sign'):
                    assert segment.get(key)==prior_segment.get(key), ('Emitted segment source order',key)
                assert segment['vertices'][0]==prior_segment['vertices'][0] and segment['vertices'][-1]==prior_segment['vertices'][-1]
        reference=(rp,rb)
    segs={c['id']:c['segments'] for c in new['curves']}
    for branch_index,branch in enumerate(bs):
        q=np.array([pos[v] for v in branch['vertices']])
        parameter_q=q if reference is None else np.array([reference[0][v] for v in branch['vertices']])
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(parameter_q,axis=0),axis=1))];length=arc[-1];arc/=length
        if reference is not None:
            assert np.isfinite(q).all() and np.isfinite(arc).all()
            assert arc[0]==0. and arc[-1]==1. and np.all(np.diff(arc)>0)
        vparam=dict(zip(branch['vertices'],arc))
        closed=branch['vertices'][0]==branch['vertices'][-1]
        if closed:vparam[branch['vertices'][0]]=0.
        n=max(6,int(np.ceil(length/55))+3)
        knots=np.r_[np.zeros(4),np.linspace(0,1,n-2)[1:-1],np.ones(4)]
        segments=[s for s in segs[branch['curve']] if s['vertices'][0] in vparam and s['vertices'][-1] in vparam]
        intervals=[];samples=[]
        for s in segments:
            a,z=vparam[s['vertices'][0]],vparam[s['vertices'][-1]]
            if closed and s['vertices'][-1]==branch['vertices'][0]:z=1.
            if 'spline_parameter_intervals' in s:
                spans=s['spline_parameter_intervals']
                assert spans and len(spans)==len(s['beziers'])
                assert spans[0][0]==a and spans[-1][1]==z
                assert all(a<=lo<hi<=z for lo,hi in spans)
                assert all(x[1]==y[0] for x,y in zip(spans,spans[1:]))
            else:
                assert reference is None, 'Diagnostic parameter transport requires actual emitted spans'
                ts=sorted(set([a,z]+[float(t) for t in knots if a+1e-12<t<z-1e-12]))
                assert len(ts)-1==len(s['beziers'])
                spans=list(zip(ts,ts[1:]))
            for (lo,hi),bezier in zip(spans,s['beziers']):
                curve=np.array(bezier);intervals.append((lo,hi,curve));samples.append((lo,curve[0]));sample(curve,lo,hi,samples)
        intervals.sort(key=lambda p:p[0]);assert abs(intervals[0][0])<1e-12 and abs(intervals[-1][1]-1)<1e-12
        if reference is not None:
            assert intervals and intervals[0][0]==0. and intervals[-1][1]==1.
            assert all(x[1]==y[0] for x,y in zip(intervals,intervals[1:]))
        # Include every old PL corner, so t=0 is exactly the old drawing.
        for t in arc:
            lo,hi,bez=next((lo,hi,bez) for lo,hi,bez in intervals if lo-1e-14<=t<=hi+1e-14)
            u=np.clip((t-lo)/(hi-lo),0,1);p=(1-u)**3*bez[0]+3*u*(1-u)**2*bez[1]+3*u*u*(1-u)*bez[2]+u**3*bez[3]
            samples.append((t,p))
        samples.sort(key=lambda p:p[0]);unique=[]
        for t,p in samples:
            if not unique or abs(t-unique[-1][0])>1e-13:unique.append((t,p))
        vs=[]
        for k,(t,p) in enumerate(unique):
            vertex=branch['vertices'][0] if k==0 else branch['vertices'][-1] if k==len(unique)-1 else nextid
            if k not in (0,len(unique)-1):nextid+=1
            start=np.array([np.interp(t,arc,q[:,axis]) for axis in (0,1)])
            if k==0 or k==len(unique)-1:assert np.linalg.norm(start-p)<1e-7
            vs.append((vertex,start,p))
        for edge_index,(a,b) in enumerate(zip(vs,vs[1:])):
            if np.linalg.norm(a[1]-b[1])+np.linalg.norm(a[2]-b[2])<1e-7:continue
            starts.append([a[1],b[1]]);ends.append([a[2],b[2]]);ids.append((a[0],b[0]));edges.append(branch['curve'])
            signatures.append((branch['curve'],branch['vertices'][0],branch['vertices'][-1]))
            source_parameters.append(dict(interval=[float(unique[edge_index][0]),float(unique[edge_index+1][0])],
                source_points_svg=[a[1].tolist()]+q[np.searchsorted(arc,unique[edge_index][0],side='right'):np.searchsorted(arc,unique[edge_index+1][0],side='left')].tolist()+[b[1].tolist()]))
    result=(np.array(starts),np.array(ends),ids,edges)
    if with_source_parameters:return (*result,signatures,source_parameters)
    return (*result,signatures) if with_branch_signatures else result


def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def polynomial(a,da,b,db,c,dc):
    u,v=b-a,c-a;du,dv=db-da,dc-da
    return np.array([cross(du,dv),cross(du,v)+cross(u,dv),cross(u,v)])


def roots(poly):
    a,b,c=poly
    if abs(a)<1e-12:
        return [] if abs(b)<1e-12 else [-c/b]
    disc=b*b-4*a*c
    if disc<0:return []
    r=np.sqrt(disc)
    return [(-b-r)/(2*a),(-b+r)/(2*a)]


def _motion_collision_location(p,q):
    # Source: the SAME evaluated moving edges at the rejected event time.
    from shapely.ops import nearest_points
    a,b=nearest_points(LineString(p),LineString(q))
    return ((np.array(a.coords[0])+np.array(b.coords[0]))/2).tolist()


def check(old,new, *, _parameter_source=None):
    a,b,ids,names,signatures,source_parameters=correspondence(old,new,with_branch_signatures=True,_parameter_source=_parameter_source,with_source_parameters=True);delta=b-a
    hulls=[MultiPoint(np.vstack((u,v))).convex_hull for u,v in zip(a,b)]
    holes=[(np.array(old['positions'][str(h['cap'])])*600,h['radius']*600) for h in old['boundaries']]
    for i,h in enumerate(hulls):
        for c,radius in holes:
            assert h.distance(Point(c))>=radius-1e-6,('Swept hull meets hole',i,names[i],h.distance(Point(c)),dict(branch_signatures=[signatures[i]]))
    print('Swept hole checks passed',len(hulls),flush=True)
    tree=STRtree(hulls);tested=0
    for i,h in enumerate(hulls):
        for j in tree.query(h,predicate='intersects'):
            if j<=i:continue
            common=set(ids[i])&set(ids[j]);u,v=a[i];du,dv=delta[i];w,z=a[j];dw,dz=delta[j]
            polynomials=[polynomial(u,du,v,dv,w,dw),polynomial(u,du,v,dv,z,dz),polynomial(w,dw,z,dz,u,du),polynomial(w,dw,z,dz,v,dv)]
            events=sorted(set([0.,1.]+[float(t) for poly in polynomials for t in roots(poly) if 0<t<1]))
            times=events+[(x+y)/2 for x,y in zip(events,events[1:])]
            tested+=1
            for t in times:
                p=a[i]+t*delta[i];q=a[j]+t*delta[j]
                if common:
                    vertex=next(iter(common));k=ids[i].index(vertex);l=ids[j].index(vertex)
                    r=p[1-k]-p[k];s=q[1-l]-q[l]
                    assert not (abs(cross(r,s))<1e-8*np.linalg.norm(r)*np.linalg.norm(s) and np.dot(r,s)>0),('Incident edges collapse',i,int(j),t,dict(branch_signatures=[signatures[i],signatures[j]],location_svg=_motion_collision_location(p,q)))
                else:
                    orientations=np.array([np.polyval(poly,t) for poly in polynomials])
                    if orientations[0]*orientations[1]>1e-15 or orientations[2]*orientations[3]>1e-15:continue
                    assert LineString(p).distance(LineString(q))>1e-7,('Interpolation collision',i,int(j),names[i],names[j],t,dict(branch_signatures=[signatures[i],signatures[j]],source_parameters=[source_parameters[i],source_parameters[j]],location_svg=_motion_collision_location(p,q)))
        if i%2000==0:print('Continuous edge check',i,len(hulls),tested,flush=True)
    report=dict(method='quadratic orientation events for straight PL interpolation; convex swept hull hole checks',
                arithmetic='floating-point, tolerance 1e-7 SVG units',moving_edges=len(hulls),candidate_pairs=tested,
                crossings_and_ports_fixed=True,collisions=0,hole_entries=0)
    if _parameter_source is not None:report['parameter_reference_binding']=_parameter_reference_binding(old,_parameter_source)
    return report


if __name__=='__main__':
    root=Path(sys.argv[1]);read=lambda f:json.loads((root/f).read_text(encoding='utf-8'))
    result=check(read('KodaSeparated.json'),read('KodaSmoothCandidate.json'))
    result['source_sha256']=hashlib.sha256((root/'KodaSeparated.json').read_bytes()).hexdigest()
    result['candidate_sha256']=hashlib.sha256((root/'KodaSmoothCandidate.json').read_bytes()).hexdigest()
    (root/'KodaSmoothIsotopy.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
