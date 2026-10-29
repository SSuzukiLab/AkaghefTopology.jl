"""Fair the source-edge paths from marked_cut_layout.fixed_mesh.

The regularized area denominator is the same numerical search device as in
fixed_mesh; scipy.optimize.minimize supplies L-BFGS-B. The first pass shortens
the paths; the second also spreads crossings. Neither pass replaces topology
checks or claims a global optimum. Final paths remain piecewise linear.
"""
import copy
from collections import defaultdict
import numpy as np
from scipy.optimize import minimize


def smooth(data,readable=True):
    d=copy.deepcopy(data)
    ids=sorted(map(int,d['positions']));ix={v:i for i,v in enumerate(ids)}
    P=np.array([d['positions'][str(v)] for v in ids])
    T=np.array([[ix[v] for v in t] for t in d['triangles']])
    counts=defaultdict(int)
    for t in d['triangles']:
        for a,b in zip(t,t[1:]+t[:1]):counts[tuple(sorted((a,b)))]+=1
    pinned={v for e,n in counts.items() if n==1 for v in e}
    pinned.update(b['cap'] for b in d['boundaries'])
    pinned.update(p['vertex'] for b in d['boundaries'] for p in b['ports'])
    free=np.array([ix[v] for v in ids if v not in pinned])
    edges=[];triples=[];outgoing={}
    for curve in d['curves']:
        before={};after={}
        for s in curve['segments']:
            for a,b in zip(s['vertices'],s['vertices'][1:]):
                edges.append((ix[a],ix[b]));after[a]=b;before[b]=a
        triples.extend((ix[before[v]],ix[v],ix[after[v]]) for v in before.keys()&after.keys())
        outgoing[curve['id']]=after
    E=np.array(edges);V=np.array(triples)
    X=np.array([(ix[x['vertex']],ix[outgoing[x['curves'][0]][x['vertex']]],ix[outgoing[x['curves'][1]][x['vertex']]]) for x in d['crossings']])
    adjacent=defaultdict(set)
    for u,v in edges:adjacent[u].add(v);adjacent[v].add(u)
    avoid=[];visible=sorted(adjacent)
    for j,v in enumerate(visible):
        near={v}|adjacent[v]|{w for u in adjacent[v] for w in adjacent[u]}
        avoid.extend((v,u) for u in visible[:j] if u not in near)
    R=np.array(avoid)
    def objective(x,turn_weight,eps):
        Q=P.copy();Q[free]=x.reshape((-1,2));grad=np.zeros_like(Q)
        a,b,c=(Q[T[:,j]] for j in range(3));u=b-a;v=c-a
        area=u[:,0]*v[:,1]-u[:,1]*v[:,0]
        root=np.sqrt(area*area+eps*eps)
        chi=np.empty_like(area);positive=area>=0
        chi[positive]=(area[positive]+root[positive])/2
        chi[~positive]=eps*eps/(2*(root[~positive]-area[~positive]))
        derivative=chi/root
        S=np.sum((a-b)**2+(b-c)**2+(c-a)**2,axis=1)
        barrier=.0003
        value=barrier*np.sum(S/chi)
        for j,r,s in ((0,b,c),(1,c,a),(2,a,b)):
            g=(-barrier*S*derivative/(chi*chi))[:,None]*np.column_stack((r[:,1]-s[:,1],s[:,0]-r[:,0]))
            np.add.at(grad,T[:,j],g)
        for j,r,s,t in ((0,a,b,c),(1,b,c,a),(2,c,a,b)):
            np.add.at(grad,T[:,j],2*barrier*(2*r-s-t)/chi[:,None])
        # Modest mesh energy keeps unused scaffold vertices well conditioned.
        value+=.02*np.sum((a-b)**2+(b-c)**2+(c-a)**2)
        for j,r,s,t in ((0,a,b,c),(1,b,c,a),(2,c,a,b)):
            np.add.at(grad,T[:,j],.04*(2*r-s-t))
        delta=Q[E[:,1]]-Q[E[:,0]];length=np.linalg.norm(delta,axis=1)
        value+=length.sum()*8
        g=8*delta/length[:,None]
        np.add.at(grad,E[:,0],-g);np.add.at(grad,E[:,1],g)
        u=Q[V[:,1]]-Q[V[:,0]];v=Q[V[:,2]]-Q[V[:,1]]
        lu=np.linalg.norm(u,axis=1);lv=np.linalg.norm(v,axis=1)
        U=u/lu[:,None];W=v/lv[:,None];diff=U-W
        value+=turn_weight*np.sum(diff*diff)
        gu=2*turn_weight*(diff-U*np.sum(U*diff,axis=1)[:,None])/lu[:,None]
        gv=-2*turn_weight*(diff-W*np.sum(W*diff,axis=1)[:,None])/lv[:,None]
        np.add.at(grad,V[:,0],-gu);np.add.at(grad,V[:,1],gu-gv);np.add.at(grad,V[:,2],gv)
        u=Q[X[:,1]]-Q[X[:,0]];v=Q[X[:,2]]-Q[X[:,0]]
        lu=np.linalg.norm(u,axis=1);lv=np.linalg.norm(v,axis=1)
        U=u/lu[:,None];W=v/lv[:,None];dot=np.sum(U*W,axis=1)
        angle_weight=10 if readable else 0
        value+=angle_weight*np.sum(dot*dot)
        gu=2*angle_weight*dot[:,None]*(W-U*dot[:,None])/lu[:,None]
        gv=2*angle_weight*dot[:,None]*(U-W*dot[:,None])/lv[:,None]
        np.add.at(grad,X[:,0],-gu-gv);np.add.at(grad,X[:,1],gu);np.add.at(grad,X[:,2],gv)
        delta=Q[R[:,1]]-Q[R[:,0]];dist=np.linalg.norm(delta,axis=1)
        gap=.075;short=np.maximum(0,gap-dist)
        repulsion=2 if readable else 0
        value+=repulsion*np.sum((short/gap)**2)
        g=-2*repulsion*short[:,None]*delta/(gap*gap*np.maximum(dist[:,None],1e-12))
        np.add.at(grad,R[:,0],-g);np.add.at(grad,R[:,1],g)
        return float(value),grad[free].ravel()
    x=P[free].ravel()
    stages=(2,8,25) if readable else (.02,.1,.4)
    for weight in stages:
        eps=1e-8
        # Regularized search may cross invalid iterates; only verified output is saved.
        opt=minimize(lambda x:objective(x,weight,eps),x,jac=True,method='L-BFGS-B',options={'maxiter':3500,'maxls':60,'ftol':1e-12,'maxcor':30})
        x=opt.x;P[free]=x.reshape((-1,2))
        print('Curve fairing',weight,opt.fun,opt.nit,opt.message,flush=True)
    d['positions']={str(v):P[ix[v]].tolist() for v in ids}
    a,b,c=(P[T[:,j]] for j in range(3));u=b-a;v=c-a
    d['minimum_double_area']=float(min(u[:,0]*v[:,1]-u[:,1]*v[:,0]))
    assert d['minimum_double_area']>0
    d['layout_method']+='; positive-mesh curve fairing'
    d.setdefault('fairing_stages',[]).extend(stages)
    return d


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    from draw_marked_cut import render
    from verify_fixed_grid import verify
    path=Path(sys.argv[1]);d=smooth(json.loads(path.read_text(encoding='utf-8')))
    report=verify(d,json.loads(path.with_name('KodaMarkedCut.json').read_text(encoding='utf-8')))
    d['geometry_check']=report
    out=path.with_name('KodaFixedMeshSmooth.json')
    out.write_text(json.dumps(d),encoding='utf-8')
    out.with_suffix('.svg').write_text(render(d),encoding='utf-8')
    print(report)
