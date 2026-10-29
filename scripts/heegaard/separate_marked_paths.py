"""Open accidental near approaches by compactly supported ambient deformations.

The radial map is f(r)=r*(1+a*(1-r/R)^3) inside R, identity outside.
For 0<a<=2, f'(r)>=1-a/4>=1/2, so it is an orientation-preserving
homeomorphism. Its support avoids every puncture disk. ALL curves, crossings
and the D boundary are transported together; no curve is independently rerouted.
Polyline sampling is subsequently checked for intersections and source-D equality.
"""
import copy
import json
from pathlib import Path
import numpy as np
from measure_clearance import measure


def expand(d,center,radius,amount=2.):
    assert 0<amount<=2 and radius>0
    result=copy.deepcopy(d)
    old={int(v):np.asarray(p,dtype=float)*600 for v,p in d['positions'].items()}
    center=np.asarray(center,dtype=float)
    for b in d['boundaries']:
        assert np.linalg.norm(center-old[b['cap']])-600*b['radius']>radius
    def mapped(p):
        delta=p-center;t=np.linalg.norm(delta)/radius
        return p if t>=1 else center+delta*(1+amount*(1-t)**3)
    points={v:mapped(p) for v,p in old.items()};nextid=max(points)+1;cache={}
    def edge(u,v):
        nonlocal nextid
        if (u,v) in cache:return cache[u,v]
        a,b=old[u],old[v];delta=b-a
        def sample(t0,t1,depth=0):
            p,q=mapped(a+delta*t0),mapped(a+delta*t1)
            # Also split at support crossings; a long chord must not skip a disk.
            tests=[.25,.5,.75]
            error=max(np.linalg.norm(mapped(a+delta*(t0+(t1-t0)*s))-(p+(q-p)*s)) for s in tests)
            if error<.025 or depth>=18:return [q]
            mid=(t0+t1)/2;return sample(t0,mid,depth+1)+sample(mid,t1,depth+1)
        aa=np.dot(delta,delta);bb=2*np.dot(a-center,delta);cc=np.dot(a-center,a-center)-radius*radius
        disc=bb*bb-4*aa*cc;ts=[0.,1.]
        if disc>0 and aa>0:
            ts.extend(t for t in ((-bb-np.sqrt(disc))/(2*aa),(-bb+np.sqrt(disc))/(2*aa)) if 0<t<1)
        samples=[];ts=sorted(set(ts))
        for t0,t1 in zip(ts,ts[1:]):samples.extend(sample(t0,t1))
        ids=[u]
        for p in samples[:-1]:points[nextid]=p;ids.append(nextid);nextid+=1
        ids.append(v);cache[u,v]=ids;cache[v,u]=ids[::-1];return ids
    for c in result['curves']:
        for s in c['segments']:
            vs=s['vertices'];s['vertices']=[w for u,v in zip(vs,vs[1:]) for w in edge(u,v)[:-1]]+[vs[-1]]
    vs=result['disk'];result['disk']=[w for u,v in zip(vs,vs[1:]+vs[:1]) for w in edge(u,v)[:-1]]
    result['positions']={str(v):(p/600).tolist() for v,p in points.items()}
    result.setdefault('clearance_isotopies',[]).append(dict(center=center.tolist(),radius=radius,amount=amount,
        units='SVG units',minimum_radial_derivative=1-amount/4,sampling_tolerance=.025))
    return result


def separate(d,target=12.):
    report=measure(d,target);before=copy.deepcopy(report)
    energy=lambda r:sum((target-p['gap'])**2 for p in r['pairs'])
    for iteration in range(60):
        if not report['pairs']:break
        worst=report['pairs'][0];center=np.mean(worst['points'],axis=0)
        hole_clearance=min(np.linalg.norm(center-np.asarray(d['positions'][str(b['cap'])])*600)-b['radius']*600 for b in d['boundaries'])
        radius=min(90.,hole_clearance*.90)
        assert radius>target/2, 'Insufficient fixed-hole room for this local clearance deformation'
        candidate=expand(d,center,radius)
        proposed=measure(candidate,target)
        if energy(proposed)>=energy(report)-1e-5:
            candidate=expand(d,center,min(140.,hole_clearance*.95),amount=1.)
            proposed=measure(candidate,target)
            if energy(proposed)>=energy(report)-1e-5:break
        d,report=candidate,proposed
        print('Clearance',iteration+1,'violations',report['violations'],'minimum',
            report['pairs'][0]['gap'] if report['pairs'] else '>='+str(target),flush=True)
    d['clearance_before']=before;d['clearance_check']=report
    d['layout_method']+='; compact radial separation of non-crossing branches'
    return d


if __name__=='__main__':
    import argparse
    from publish_marked_display import bind
    from draw_marked_cut import render
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('directory',type=Path);args=parser.parse_args()
    read=lambda n:json.loads((args.directory/n).read_text(encoding='utf-8'))
    result=separate(read('KodaDisplay.json'))
    result=bind(result,read('KodaPrepared.json'),read('KodaMarkedCut.json'))
    result['geometry_provenance']='saved realization plus verified compact radial clearance deformations'
    assert result['clearance_check']['violations']==0, 'Clearance target not reached; do not publish'
    (args.directory/'KodaSeparated.json').write_text(json.dumps(result),encoding='utf-8')
    (args.directory/'KodaSeparated.svg').write_text(render(result),encoding='utf-8')
    print(json.dumps(result['clearance_check']))
