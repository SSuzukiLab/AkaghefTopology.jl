"""Deterministic placement, ported from the owner-approved mock3.py.

Source: research workspace .research-scratch/PJ260922/
261109-lickorish-drawing-mock/mock/mock3.py, geometry/arcs/chord_check/tutte/Phi.
Contract: docs/specs/heegaard-lickorish-data.md sections 6.2 and 6.3.
Chord arithmetic, face order, sampling and the chord case's 160-by-160
infinity-center search are preserved. Fallback-only subdivision is documented
in docs/heegaard/261109_tutte_construction.md. This module consumes already
checked LC data; see lickorish_draw.py.
The returned coordinates use y-up, exactly as the reference JSON.
"""
import itertools
import math
from lickorish_faces import solve_faces
from lickorish_empty import store_intervals, add_empty
from lickorish_render import svg as render_svg
from lickorish_knots import attach_crossing_knots


def place(data, *, show_k=False):
    """Return (placement document, polyline SVG, diagnostic face information)."""
    if data['genus']==1:
        from lickorish_square import place_square
        doc,svg,details=place_square(data)
        add_empty(data,doc,details)
        for key in ('intervals','cuts','holes','map_point'): details.pop(key,None)
        return doc,render_svg(doc),details
    if data['genus']==0:
        if any(data['words'].values()):
            raise ValueError('sphere crossing placement is not part of mock3.py')
        doc=dict(W=520.,H=520.,holes={},labels={},curves={})
        details=dict(methods={'outer':'empty circles'},arcs=[],intervals={},holes=[],map_point=lambda f,p:p)
        add_empty(data,doc,details)
        for key in ('intervals','holes','map_point'): details.pop(key,None)
        return doc,render_svg(doc),details
    D = data
    g = D['genus']
    LD = {k: list(v) for k, v in D['words'].items()}
    CR={int(k):tuple(v) for k,v in D['crossings'].items()}
    Kn=[f'l{i}' for i in range(1,g+1)]+[f'm{i}' for i in range(1,g+1)]+[f'n{i}' for i in range(1,g)]
    CU=[k for k in LD if k not in Kn and LD[k]]
    intervals={}
    # ---------- display geometry (y-up): rectangle W x H, columns at x_i, triangular punctures
    colgap=520.;W=colgap*(g-1);H=520.
    xcol={i:colgap*(i-1) for i in range(1,g+1)}
    def count_l(i):return sum(1 for x in LD[f'l{i}'] if CR[abs(x)][0].startswith('l') and not CR[abs(x)][1].startswith(('m','n')))
    A=max(130.,min(220.,40.*max(1,max(count_l(i) for i in range(1,g+1)))+120))   # triangle size (mock rule)
    O=(W/2,H/2)
    def lerp(p,q,t):return (p[0]+(q[0]-p[0])*t,p[1]+(q[1]-p[1])*t)
    mirror=lambda p:(p[0],H-p[1])
    # K vertex labels on each l_i: find labels shared with m_i, n_{i-1}, n_i
    def kv(l,other):return next(abs(x) for x in LD[l] if abs(x) in {abs(y) for y in LD[other]})
    # triangles (top copies) with named vertices; pieces of l_i between K vertices in CCW(y-up) order,
    # each piece -> (point a, point b, face)
    tri={};lpieces={}
    for i in range(1,g+1):
        x=xcol[i]
        if i==1:
            Av=(x,H-A);Bv=(x+A,H);C=(x,H);tri[i]=[C,Av,Bv]
            vm=kv('l1','m1');vn=kv('l1','n1')
            # CCW: Bv(n1) -> C -> Av(m1) [outer sides, Finf], Av(m1) -> Bv(n1) [chamfer, F1]
            lpieces[i]=[(vm,vn,[(Av,Bv)],'F1'),(vn,vm,[(Bv,C),(C,Av)],'Finf')]
        elif i==g:
            Av=(x,H-A);Bv=(x-A,H);C=(x,H);tri[i]=[C,Av,Bv]
            vm=kv(f'l{g}',f'm{g}');vn=kv(f'l{g}',f'n{g-1}')
            # CCW: C -> Bv(n) -> Av(m) [chamfer, F_{g-1}], Av(m) -> C -> Bv(n) [outer, Finf]
            lpieces[i]=[(vn,vm,[(Bv,Av)],f'F{g-1}'),(vm,vn,[(Av,C),(C,Bv)],'Finf')]
        else:
            BL=(x-A,H);BR=(x+A,H);Av=(x,H-A*0.9);tri[i]=[BL,BR,Av]
            vm=kv(f'l{i}',f'm{i}');vl=kv(f'l{i}',f'n{i-1}');vr=kv(f'l{i}',f'n{i}')
            # CCW: BR(n_i) -> BL(n_{i-1}) [top, Finf], BL -> Av(m_i) [left chamfer, F_{i-1}], Av -> BR [right chamfer, F_i]
            lpieces[i]=[(vr,vl,[(BR,BL)],'Finf'),(vl,vm,[(BL,Av)],f'F{i-1}'),(vm,vr,[(Av,BR)],f'F{i}')]
    # K edges m_i, n_i as straight pieces with traversal order and left/right faces (y-up)
    edges={};side={}
    for i in range(1,g+1):
        tp=[v for v in tri[i] if v[1]<H][0] if i in(1,g) else tri[i][2]   # m junction on the top triangle
        edges[f'm{i}']=[(tp,mirror(tp))]
        left='F%d'%i if i<g else 'Finf';right='F%d'%(i-1) if i>1 else 'Finf'
        side[f'm{i}']={+1:left,-1:right}
    for i in range(1,g):
        a=[v for v in tri[i] if v[1]==H and (v[0]>xcol[i] or i==1)]; a=max(a,key=lambda v:v[0])  # right-most top vertex of l_i triangle
        b=[v for v in tri[i+1] if v[1]==H]; b=min(b,key=lambda v:v[0])                              # left-most top vertex of l_{i+1}
        edges[f'n{i}']=[(a,b),(mirror(b),mirror(a))]
        side[f'n{i}']={+1:'Finf',-1:f'F{i}'}
    KV={abs(x) for l in Kn if l.startswith('l') for x in LD[l] if CR[abs(x)][1].startswith(('m','n')) and CR[abs(x)][0].startswith('l')}
    # ---------- points on K edges
    pos={};lface={}
    for e,pieces in edges.items():
        # cf:docs/heegaard/261109_cycle_knots.md K 辺の始点 — anchor each directed interval at its positive K crossing
        labs=[abs(x) for x in LD[e]];i0=next(i for i,x in enumerate(LD[e]) if x>0 and x in KV);labs=labs[i0:]+labs[:i0]
        groups=[];cur=[]
        for l in labs[1:]+[labs[0]]:
            if l in KV:groups.append(cur);cur=[]
            else:cur.append(l)
        if len(groups)!=len(pieces): raise ValueError((e,groups))
        for vi,((p,q),gr) in enumerate(zip(pieces,groups)):
            start=[l for l in labs if l in KV][vi]
            signed=[next(x for x in LD[e] if abs(x)==l) for l in [start]+gr]
            store_intervals(intervals,e,signed,[p,q])
            for j,l in enumerate(gr):pos[l]=lerp(p,q,(j+1)/(len(gr)+1))
    for i in range(1,g+1):
        w=[abs(x) for x in LD[f'l{i}']]
        for v0,v1,segs,face in lpieces[i]:
            a=w.index(v0);ww=w[a:]+w[:a];b=ww.index(v1);labs=ww[1:b]
            L=sum(math.dist(p,q) for p,q in segs);n=len(labs)
            signed=[next(x for x in LD[f'l{i}'] if abs(x)==l) for l in [v0]+labs]
            path=[segs[0][0]]+[q for p,q in segs]
            store_intervals(intervals,f'l{i}',signed,path)
            store_intervals(intervals,f'l{i}',signed,[mirror(p) for p in path])
            for j,l in enumerate(labs):
                t=(j+1)/(n+1)*L
                for p,q in segs:
                    d=math.dist(p,q)
                    if t<=d+1e-9:pos[l]=lerp(p,q,t/d);break
                    t-=d
                lface[l]=face
    # ---------- arcs
    labedge={}
    for e in Kn:
        for x in LD[e]:
            if abs(x) not in KV:labedge[abs(x)]=e
    def ksign(e,lab):return next(x for x in LD[e] if abs(x)==lab)>0
    def csign(c,lab):return next(x for x in LD[c] if abs(x)==lab)>0
    arcs=[]
    for c in CU:
        w=LD[c];n=len(w);kidx=[i for i,x in enumerate(w) if abs(x) in labedge]
        if not kidx: raise ValueError('curve without K crossings: '+c)
        for a,b in zip(kidx,kidx[1:]+[kidx[0]+n]):
            s,e=abs(w[a%n]),abs(w[b%n]);mids=[abs(w[i%n]) for i in range(a+1,b)]
            es=labedge[s];face=lface[s] if es.startswith('l') else side[es][+1 if ksign(es,s) else -1]
            ee=labedge[e];before=lface[e] if ee.startswith('l') else side[ee][-1 if ksign(ee,e) else +1]
            if before!=face: raise ValueError(('face mismatch',c,s,e,face,before))
            arc=dict(curve=c,start=s,end=e,mids=mids,face=face)
            arc['P']=(mirror(pos[s]) if not csign(c,s) else pos[s]) if es.startswith('l') else pos[s]
            arc['Q']=(pos[e] if not csign(c,e) else mirror(pos[e])) if ee.startswith('l') else pos[e]
            arcs.append(arc)
    # ---------- faces: boundary parametrisation
    faces=[f'F{i}' for i in range(1,g)]+['Finf']
    L=2*(W+H)
    def frac_of(p):
        x,y=p
        if abs(y)<1e-6:return x/L
        if abs(x-W)<1e-6:return (W+y)/L
        if abs(y-H)<1e-6:return (W+H+(W-x))/L
        return (2*W+H+(H-y))/L
    def outline_pt(f):
        t=(f%1)*L
        if t<W:return (t,0)
        t-=W
        if t<H:return (W,t)
        t-=H
        if t<W:return (W-t,H)
        t-=W;return (0,H-t)
    for a in arcs:
        if a['face']=='Finf':
            for key,P in (('mP',a['P']),('mQ',a['Q'])):
                th=2*math.pi*frac_of(P);a[key]=(math.cos(th),math.sin(th))
        else:a['mP'],a['mQ']=a['P'],a['Q']
    # cf:docs/heegaard/261109_tutte_construction.md 組合せ的な細分 — fixed convex boundary corners for the fallback only
    boundaries={'Finf':[(1.,0.),(0.,1.),(-1.,0.),(0.,-1.)]}
    for i in range(1,g):
        p=edges[f'm{i}'][0][0];q=edges[f'm{i+1}'][0][0]
        a,b=edges[f'n{i}'][0]
        boundaries[f'F{i}']=[p,mirror(p),mirror(a),mirror(b),mirror(q),q,b,a]
    meshes={}
    method=solve_faces(arcs,LD,CU,faces,boundaries,meshes)
    # ---------- map model -> plane; Finf through the radial map from c
    def dseg(p,a,b):
        ab=(b[0]-a[0],b[1]-a[1]);t=max(0,min(1,((p[0]-a[0])*ab[0]+(p[1]-a[1])*ab[1])/(ab[0]**2+ab[1]**2)))
        q=(a[0]+ab[0]*t,a[1]+ab[1]*t);return math.hypot(p[0]-q[0],p[1]-q[1])
    finf=[a for a in arcs if a['face']=='Finf'];fsegs=[(a['poly_model'][k],a['poly_model'][k+1]) for a in finf for k in range(len(a['poly_model'])-1)]
    best=None;N=160
    if method['Finf']=='tutte':
        # cf:docs/heegaard/261109_tutte_construction.md 外側の中心 — triangle-interior candidates, same clearance objective
        mesh=meshes['Finf'];points=mesh['points']
        candidates=sorted(tuple(sum(points[v][k] for v in t)/3 for k in (0,1)) for t in mesh['triangles'])
        for p in candidates:
            d=min([dseg(p,a,b) for a,b in fsegs]+[1-math.hypot(*p)])
            if best is None or d>best[0]:best=(d,p)
    else:
        for i in range(N):
            for j in range(N):
                p=(-1+2*(i+.5)/N,-1+2*(j+.5)/N)
                if math.hypot(*p)>=0.97:continue
                d=min([dseg(p,a,b) for a,b in fsegs]+[1-math.hypot(*p)])
                if best is None or d>best[0]:best=(d,p)
    delta,CC=best;print('c =',CC,'delta =',round(delta,3));kappa=1.5
    def Phi(p):
        v=(p[0]-CC[0],p[1]-CC[1]);nv=math.hypot(*v);u=(v[0]/nv,v[1]/nv)
        bq=CC[0]*u[0]+CC[1]*u[1];cq=CC[0]**2+CC[1]**2-1;t=-bq+math.sqrt(bq*bq-cq)
        q=(CC[0]+t*u[0],CC[1]+t*u[1]);s=nv/t
        th=math.atan2(q[1],q[0])%(2*math.pi);P=outline_pt(th/(2*math.pi))
        d=math.hypot(P[0]-O[0],P[1]-O[1]);uu=((P[0]-O[0])/d,(P[1]-O[1])/d)
        f=1+(kappa-1)*(1-s)/(1-delta) if s<1 else 1.
        return (O[0]+d*f*uu[0],O[1]+d*f*uu[1])
    for a in arcs:
        ch=a['poly_model']
        if a['face']=='Finf':
            pts=[]
            for p,q in zip(ch,ch[1:]):pts+=[Phi(lerp(p,q,k/300)) for k in range(300)]
            pts.append(Phi(ch[-1]));pts[0]=a['P'];pts[-1]=a['Q'];a['poly']=pts
        else:a['poly']=list(ch)
    # ---------- punctures: incircles, extension
    def incircle(t):
        a=math.dist(t[1],t[2]);b=math.dist(t[0],t[2]);cc=math.dist(t[0],t[1]);p=a+b+cc
        I=((a*t[0][0]+b*t[1][0]+cc*t[2][0])/p,(a*t[0][1]+b*t[1][1]+cc*t[2][1])/p)
        s=p/2;area=math.sqrt(s*(s-a)*(s-b)*(s-cc));return I,area/s
    holes={}
    for i in range(1,g+1):
        holes[f'top{i}']=tri[i];holes[f'bot{i}']=[mirror(p) for p in tri[i]]
    inc={k:incircle(t) for k,t in holes.items()}
    def extend(p):
        for k,t in holes.items():
            if min(dseg(p,t[i],t[(i+1)%3]) for i in range(3))<1e-6:
                I,r=inc[k];v=(p[0]-I[0],p[1]-I[1]);n=math.hypot(*v);return (I[0]+r*v[0]/n,I[1]+r*v[1]/n)
        raise ValueError(p)
    # ---------- SVG + JSON export
    def Y(p):return (p[0],H-p[1])
    margin=60;vb=(-(kappa-1)*W/2-margin,-(kappa-1)*H/2-margin,kappa*W+2*margin,kappa*H+2*margin)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb[0]} {vb[1]} {vb[2]} {vb[3]}" font-family="Arial,sans-serif">',
         f'<rect x="{vb[0]}" y="{vb[1]}" width="{vb[2]}" height="{vb[3]}" fill="#fcfbf7"/>',
         f'<text x="{vb[0]+20}" y="{vb[1]+35}" font-size="23" fill="#253247">Heegaard diagram · g = {g} (Lickorish Data, placement)</text>']
    color=lambda c:'#ba433f' if c.startswith(('a','alpha')) else '#255dba'
    export={'W':W,'H':H,'holes':{k:[list(I),r] for k,(I,r) in inc.items()},'labels':D.get('holes',{}),'curves':{}}
    def path_d(segs):
        # trial mode: plain polyline, no fillets
        pts=[p for sg in segs for p in sg]
        return 'M'+' L'.join(f'{Y(p)[0]:.1f},{Y(p)[1]:.1f}' for p in pts)
    
    owners={}
    # cf:docs/heegaard/261109_smoothing_gates.md LC gate — signed visits in each drawn run; coordinates unchanged
    export['crossing_structure']={}
    arc_of={id(a['poly']):a for a in arcs}
    for c in CU:
        segs=[a['poly'] for a in arcs if a['curve']==c]
        # split into runs at puncture passages (discontinuities between consecutive arcs)
        runs=[[segs[0]]]
        for sg in segs[1:]:
            if math.dist(runs[-1][-1][-1],sg[0])<1e-6:runs[-1].append(sg)
            else:runs.append([sg])
        if len(runs)>1 and math.dist(runs[-1][-1][-1],runs[0][0][0])<1e-6:runs[0]=runs.pop()+runs[0]
        closed=len(runs)==1 and math.dist(runs[0][0][0],runs[0][-1][-1])<1e-6
        for ri,run in enumerate(runs):
            s0,e0=run[0][0],run[-1][-1]
            segs2=run+[run[0][:2]] if closed else [[extend(s0),s0]]+run+[[e0,extend(e0)]]
            out.append(f'<path d="{path_d(segs2)}" fill="none" stroke="{color(c)}" stroke-width="2.5"/>')
            pts=[]
            for sg in (run if closed else [[extend(s0),s0]]+run+[[e0,extend(e0)]]):
                for p in sg:
                    if not pts or math.dist(pts[-1],p)>1e-6:pts.append(list(p))
            export['curves'][c if closed else f'{c}#{ri}']={'points':pts,'closed':closed}
            owners[c if closed else f'{c}#{ri}']=c
            export['crossing_structure'][c if closed else f'{c}#{ri}']=[
                next(x for x in LD[c] if abs(x)==m)
                for sg in run for m in arc_of[id(sg)]['mids']]
    labels=D.get('holes',{})
    for k,(I,r) in inc.items():
        out.append(f'<circle cx="{Y(I)[0]:.1f}" cy="{Y(I)[1]:.1f}" r="{r:.1f}" fill="#fcfbf7" stroke="#333" stroke-width="1.5"/><text x="{Y(I)[0]:.1f}" y="{Y(I)[1]+4:.1f}" font-size="13" text-anchor="middle" fill="#253247">{labels.get(k,k)}</text>')
    if show_k:
        for e,pieces in edges.items():
            for p,q in pieces:out.append(f'<path d="M{Y(p)[0]},{Y(p)[1]}L{Y(q)[0]},{Y(q)[1]}" stroke="#7a6a9c" stroke-width="1" stroke-dasharray="6 4"/>')
        for k,t in holes.items():out.append('<path d="M'+' L'.join(f'{Y(p)[0]},{Y(p)[1]}' for p in t)+' Z" fill="none" stroke="#7a6a9c" stroke-width="1" stroke-dasharray="6 4"/>')
    out.append('</svg>')
    attach_crossing_knots(export,arcs,owners,lambda face,p:Phi(p) if face=='Finf' else p)
    details=dict(methods=method,center=CC,delta=delta,arcs=arcs,tutte=meshes,intervals=intervals,
                 holes=list(holes.values()),map_point=lambda face,p:Phi(p) if face=='Finf' else p)
    add_empty(data,export,details)
    for key in ('intervals','holes','map_point'): details.pop(key,None)
    return export, render_svg(export), details
