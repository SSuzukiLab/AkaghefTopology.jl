"""Open crossing angles by supported angular maps, without choosing new routes.

Inside disjoint disks, interpolate the increasing circle map taking the four
crossing rays to quarter turns with the identity on the support boundary.
This defines an ambient isotopy. Approximate its curve images adaptively and
run verify_fixed_grid on the actual resulting polylines before exporting.
"""
import copy,math
import numpy as np


def clarify(data):
    d=copy.deepcopy(data);old={int(v):np.array(p) for v,p in d['positions'].items()}
    maps=[]
    for x in d['crossings']:
        vertex=x['vertex'];center=old[vertex];rays=[]
        for curve in d['curves']:
            for s in curve['segments']:
                vs=s['vertices']
                if vs[0]==vertex:rays.append(old[vs[1]]-center)
                if vs[-1]==vertex:rays.append(old[vs[-2]]-center)
        assert len(rays)==4
        theta=np.sort([math.atan2(p[1],p[0]) for p in rays])
        target=np.arange(4)*math.pi/2;target+=np.mean(theta-target)
        radius=min(.11,*(.38*np.linalg.norm(center-old[y['vertex']]) for y in d['crossings'] if y!=x),*(.7*(np.linalg.norm(center-old[b['cap']])-b['radius']) for b in d['boundaries']))
        assert radius>0
        maps.append((center,radius,theta,target))
    def mapped(p):
        for center,radius,theta,target in maps:
            v=p-center;r=np.linalg.norm(v)
            if r>=radius or r<1e-14:continue
            a=math.atan2(v[1],v[0]);a=theta[0]+(a-theta[0])%(2*math.pi)
            angle=np.interp(a,np.r_[theta,theta[0]+2*math.pi],np.r_[target,target[0]+2*math.pi])
            t=np.clip((r/radius-.45)/.55,0,1);weight=1-t*t*(3-2*t)
            a+=weight*(angle-a)
            return center+r*np.array([math.cos(a),math.sin(a)])
        return p.copy()
    positions={v:mapped(p) for v,p in old.items()};cache={};nextid=max(old)+1
    def edge(u,v):
        nonlocal nextid
        if (u,v) in cache:return cache[(u,v)]
        a,b=old[u],old[v]
        ts=[0.,1.]
        # Cut at disk boundaries so a straight outside chord cannot skip a map.
        delta=b-a;aa=np.dot(delta,delta)
        for center,radius,_,_ in maps:
            for r in (radius,.45*radius):
                bb=2*np.dot(a-center,delta);cc=np.dot(a-center,a-center)-r*r;disc=bb*bb-4*aa*cc
                if disc>0:ts.extend(t for t in ((-bb-math.sqrt(disc))/(2*aa),(-bb+math.sqrt(disc))/(2*aa)) if 0<t<1)
        def sample(t0,t1,depth=0):
            p,q=mapped(a+delta*t0),mapped(a+delta*t1)
            fractions=(.25,.5,.75)
            error=max(np.linalg.norm(mapped(a+delta*(t0+(t1-t0)*f))-(p+(q-p)*f)) for f in fractions)
            if error<.015/600 or depth>=18:return [q]
            mid=(t0+t1)/2
            return sample(t0,mid,depth+1)+sample(mid,t1,depth+1)
        samples=[];ts=sorted(set(ts))
        for t0,t1 in zip(ts,ts[1:]):samples.extend(sample(t0,t1))
        ids=[u]
        for p in samples[:-1]:positions[nextid]=p;ids.append(nextid);nextid+=1
        ids.append(v);cache[(u,v)]=ids;cache[(v,u)]=ids[::-1];return ids
    for curve in d['curves']:
        for s in curve['segments']:
            vs=s['vertices'];s['vertices']=[z for u,v in zip(vs,vs[1:]) for z in edge(u,v)[:-1]]+[vs[-1]]
    vs=d['disk'];d['disk']=[z for u,v in zip(vs,vs[1:]+vs[:1]) for z in edge(u,v)[:-1]]
    d['pre_angle_mesh']={'positions':d['positions'],'triangles':d.pop('triangles'),'minimum_double_area':d.pop('minimum_double_area')}
    d['positions']={str(v):p.tolist() for v,p in positions.items()}
    d['crossing_angle_maps']=[dict(center=c.tolist(),radius=r,angles=a.tolist(),targets=t.tolist()) for c,r,a,t in maps]
    d['layout_method']+='; supported crossing-angle isotopies (sampled paths checked separately)'
    return d


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    from draw_marked_cut import render
    from verify_fixed_grid import verify
    path=Path(sys.argv[1]);d=clarify(json.loads(path.read_text(encoding='utf-8')))
    d['geometry_check']=verify(d,json.loads(path.with_name('KodaMarkedCut.json').read_text(encoding='utf-8')))
    out=path.with_name('KodaGrid.json');out.write_text(json.dumps(d),encoding='utf-8')
    out.with_suffix('.svg').write_text(render(d),encoding='utf-8')
    print(d['geometry_check'])
