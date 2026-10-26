"""Common marked-cut layout fairing; no named-example or genus-specific routing.

Reuses fair_marked_paths, check_smooth_marked, and the existing SVG renderer.
The input is an explicit feasible embedded layout, not endpoints to be routed.
"""
import argparse,copy,json,math
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from draw_marked_cut import render
from fair_marked_paths import optimise,export


def pair_ports(data):
    """One orientation-reversing circle parameter for each paired cut boundary.

Transport the entire drawing through disjoint annular angular homeomorphisms,
not just its endpoints. Outside the collars, crossings and geometry are fixed.
"""
    d=copy.deepcopy(data);pos={int(v):np.array(q)*600 for v,q in d['positions'].items()}
    pairing={a:b for u,v in d['boundary_vertex_pairs'] for a,b in ((u,v),(v,u))}
    for b in d['boundaries']:
        for p in b['ports']:
            assert abs(np.linalg.norm(pos[p['vertex']]-pos[b['cap']])-600*b['radius'])<1e-6,'Port is not on its hole boundary'
    maps=[];targets={}
    for cut in d['cuts']:
        pair=sorted([b for b in d['boundaries'] if b['cut']==cut],key=lambda b:b['row'])
        assert len(pair)==2 and [b['row'] for b in pair]==[0,1]
        top,bottom=pair
        bysource={p['source_vertex']:p for p in top['ports']}
        order=[v for v in top['vertices'] if v in bysource]
        # Old Koda geometry's boundary vertices include cap-subdivision nodes;
        # retain the native source boundary order explicitly when available.
        order=[v for v in top.get('source_vertices',top['vertices']) if v in bysource]
        assert len(order)==len(top['ports']) and len(order)==len(bottom['ports'])
        if not order:continue
        bp={p['source_vertex']:p for p in bottom['ports']}
        n=len(order);phases=[]
        for j,v in enumerate(order):
            for b,p,sgn in ((top,bysource[v],1),(bottom,bp[pairing[v]],-1)):
                q=pos[p['vertex']]-pos[b['cap']]
                phases.append(sgn*math.atan2(q[1],q[0])+2*math.pi*j/n)
        phase=math.atan2(sum(math.sin(a) for a in phases),sum(math.cos(a) for a in phases))
        for j,v in enumerate(order):
            angle=phase-2*math.pi*j/n
            targets[bysource[v]['vertex']]=angle;targets[bp[pairing[v]]['vertex']]=-angle
        already_paired=all(np.linalg.norm((pos[bysource[v]['vertex']]-pos[top['cap']])-(pos[bp[pairing[v]]['vertex']]-pos[bottom['cap']])*[1,-1])<1e-7 for v in order)
        if already_paired:
            for b in pair:
                for p in b['ports']:
                    q=pos[p['vertex']]-pos[b['cap']];targets[p['vertex']]=math.atan2(q[1],q[0])
        for b in pair:
            center=pos[b['cap']];r=b['radius']*600
            angle=sorted((math.atan2(*(pos[p['vertex']]-center)[::-1]),p['vertex']) for p in b['ports'])
            src=np.array([a for a,v in angle]);dst=np.array([targets[v] for a,v in angle])
            for j in range(1,len(dst)):
                dst[j]+=2*math.pi*math.ceil((dst[j-1]-dst[j]+1e-10)/(2*math.pi))
            assert dst[-1]-dst[0]<2*math.pi-1e-8,'Boundary port cyclic order mismatch'
            dst+=round(np.mean(src-dst)/(2*math.pi))*2*math.pi
            f=PchipInterpolator(np.r_[src-2*math.pi,src,src+2*math.pi],np.r_[dst-2*math.pi,dst,dst+2*math.pi])
            cross_dist=min((np.linalg.norm(pos[x['vertex']]-center) for x in d['crossings']),default=r+150)
            hole_dist=min((np.linalg.norm(pos[c['cap']]-center)-c['radius']*600 for c in d['boundaries'] if c is not b),default=r+150)
            support=min(r+100,cross_dist-5,(r+hole_dist)/2-2)
            assert support>r+3,'No collar room for changing ports without moving a crossing'
            maps.append((center,r,support,src,f))
    def mapped(p):
        for center,r,support,src,f in maps:
            q=p-center;rho=np.linalg.norm(q)
            if rho<1e-10 or rho>=support:continue
            a=src[0]+(math.atan2(q[1],q[0])-src[0])%(2*math.pi)
            t=np.clip((rho-r)/(support-r),0,1);weight=1-t*t*(3-2*t)
            angle=a+weight*(float(f(a))-a)
            return center+rho*np.array([math.cos(angle),math.sin(angle)])
        return p.copy()
    out={v:mapped(p) for v,p in pos.items()};nextid=max(pos)+1;cache={}
    def edge(u,v):
        nonlocal nextid
        if (u,v) in cache:return cache[u,v]
        a,b=pos[u],pos[v];delta=b-a
        def sample(t0,t1,depth=0):
            p,q=mapped(a+t0*delta),mapped(a+t1*delta)
            error=max(np.linalg.norm(mapped(a+(t0+(t1-t0)*f)*delta)-((1-f)*p+f*q)) for f in (.25,.5,.75))
            if error<.002 and np.linalg.norm(delta)*(t1-t0)<8:return [q]
            assert depth<25,'Port transport refinement failed'
            mid=(t0+t1)/2;return sample(t0,mid,depth+1)+sample(mid,t1,depth+1)
        points=sample(0,1);vs=[u]
        for q in points[:-1]:out[nextid]=q;vs.append(nextid);nextid+=1
        vs.append(v);cache[u,v]=vs;cache[v,u]=vs[::-1];return vs
    for c in d['curves']:
        for s in c['segments']:
            vs=s['vertices'];s['vertices']=[w for u,v in zip(vs,vs[1:]) for w in edge(u,v)[:-1]]+[vs[-1]]
            s.pop('beziers',None)
    vs=d.get('disk',[]);d['disk']=[w for u,v in zip(vs,vs[1:]+vs[:1]) for w in edge(u,v)[:-1]]
    d.pop('disk_beziers',None);d['positions']={str(v):(p/600).tolist() for v,p in out.items()}
    d['port_convention']='paired reflection (x,y) -> (x,-y), cyclic order retained, normal exits'
    d['port_transport']=dict(method='monotone angular maps in disjoint circular collars',crossings_fixed=True)
    return d


def plot(data, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20., correct_ports=True, maxiter=1800):
    from check_smooth_marked import check
    if data.get('schema')!='marked-cut-view/1':raise ValueError('Explicit marked-cut-view/1 embedded layout required')
    if not all(np.isfinite(x) and x>0 for x in (minimum_separation,puncture_clearance,crossing_radius)):raise ValueError('Positive finite drawing clearances required')
    if data.get('disk_pieces') or data.get('disks'):raise ValueError('Split joining-disk pieces need explicit boundary-arc representation; this input is not supported yet')
    if not data.get('boundary_vertex_pairs'):raise ValueError('Source boundary vertex pairing is required; do not infer it from endpoints')
    seed=pair_ports(data) if correct_ports else copy.deepcopy(data)
    bs,c,_=optimise(seed,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,maxiter=maxiter)
    result=export(seed,bs,c)
    result=check(result,seed,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius)
    from check_fairing_isotopy import check as check_motion
    result['layout_isotopy_check']=check_motion(seed,result)
    result['drawing_policy']=dict(minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius,port_collar_radius=20)
    result['svg']=render(result)
    # Source: draw_marked_cut.orientation_glyph_layer; optional display layer
    # uses the same directed source routes. The accepted default SVG is intact.
    try:
        result['orientation_svg']=render(result,orientation_glyphs=True)
        result.pop('orientation_display_error',None)
    except ValueError as exc:
        result.pop('orientation_svg',None)
        result['orientation_display_error']=str(exc)
    return result,seed


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--maxiter',type=int,default=1800)
    p.add_argument('--prepared',type=Path);p.add_argument('--minimum-separation',type=float,default=12.)
    p.add_argument('--recheck',action='store_true',help='Recheck existing spline output against its saved seed; do not optimise again')
    p.add_argument('--puncture-clearance',type=float,default=12.);p.add_argument('--crossing-radius',type=float,default=20.)
    a=p.parse_args();d=json.loads(a.input.read_text(encoding='utf-8'))
    if a.recheck:
        from check_smooth_marked import check
        from check_fairing_isotopy import check as check_motion
        seed=json.loads(a.input.with_suffix('.seed.json').read_text(encoding='utf-8'))
        d=check(d,seed,minimum_separation=a.minimum_separation,puncture_clearance=a.puncture_clearance,crossing_radius=a.crossing_radius)
        d['layout_isotopy_check']=check_motion(seed,d);d['svg']=render(d)
        a.output.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8');a.output.with_suffix('.svg').write_text(d['svg'],encoding='utf-8')
        print(json.dumps(d['smooth_check']),flush=True);raise SystemExit(0)
    if a.prepared:
        prepared=json.loads(a.prepared.read_text(encoding='utf-8'))
        assert prepared['source_cut_corner_classes']==d['source_cut_corner_classes']
        d['boundary_vertex_pairs']=prepared['boundary_vertex_pairs']
        for b in d['boundaries']:
            b['source_vertices']=next(c['vertices'] for c in prepared['boundaries'] if c['cut']==b['cut'] and c['side']==b['side'])
    result,seed=plot(d,maxiter=a.maxiter,minimum_separation=a.minimum_separation,puncture_clearance=a.puncture_clearance,crossing_radius=a.crossing_radius)
    a.output.write_text(json.dumps(result,ensure_ascii=False),encoding='utf-8')
    a.output.with_suffix('.svg').write_text(result.pop('svg'),encoding='utf-8')
    a.output.with_suffix('.seed.json').write_text(json.dumps(seed,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result['smooth_check']),flush=True)
