"""Genus-one cut square (specification section 6.2).

Uses the common chord/subdivision solver and the specified square boundary.
Both copies of each cut share the same ordered ports.
"""
from lc_faces import solve_faces
from lc_knots import attach_crossing_knots


def place_square(data):
    words=data['words']; size=520.
    ports={}; intervals={}; cuts=[]
    for name in ('l1','m1'):
        w=words[name]; k=next(i for i,x in enumerate(w) if abs(x)==1)
        w=w[k:]+w[:k]
        # l points right, m downwards: det(tangent_m, tangent_l)>0 at K.
        copies=[]
        for side in (0.,size):
            pts=[(size*j/len(w),side) if name=='l1' else
                 (side,size-size*j/len(w)) for j in range(len(w)+1)]
            copies.append(pts); cuts.append(pts)
            for j,x in enumerate(w): intervals.setdefault((name,x),[]).append([pts[j],pts[j+1]])
        for j,x in enumerate(w[1:],1): ports[abs(x)]=(name,copies[0][j],copies[1][j])
    curves=[c for c in words if c not in ('l1','m1') and words[c]]
    arcs=[]
    def endpoint(x,outgoing):
        _,low,high=ports[abs(x)]
        return high if (x>0)==outgoing else low
    owners={}
    for c in curves:
        w=words[c]; n=len(w); indices=[i for i,x in enumerate(w) if abs(x) in ports]
        if not indices: raise ValueError(f'F1: curve without K crossings: {c}; word={w}')
        for a,b in zip(indices,indices[1:]+[indices[0]+n]):
            p,q=endpoint(w[a],True),endpoint(w[b%n],False)
            arcs.append(dict(curve=c,start=abs(w[a]),end=abs(w[b%n]),
                             mids=[abs(w[i%n]) for i in range(a+1,b)],face='F1',
                             P=p,Q=q,mP=p,mQ=q))
    meshes={}
    methods=solve_faces(arcs,words,curves,['F1'],
                        {'F1':[(0.,0.),(size,0.),(size,size),(0.,size)]},meshes)
    # cf:docs/heegaard/261109_smoothing_gates.md LC gate — signed visits of each square arc
    doc=dict(W=size,H=size,holes={},labels={},curves={},square=True,crossing_structure={})
    for c in curves:
        for j,a in enumerate(a for a in arcs if a['curve']==c):
            a['poly']=a['poly_model']
            doc['curves'][f'{c}#{j}']=dict(points=[list(p) for p in a['poly']],closed=False)
            owners[f'{c}#{j}']=c
            doc['crossing_structure'][f'{c}#{j}']=[next(x for x in words[c] if abs(x)==m) for m in a['mids']]
    attach_crossing_knots(doc,arcs,owners,lambda face,p:p)
    return doc,'',dict(methods=methods,arcs=arcs,tutte=meshes,intervals=intervals,cuts=cuts,
                       holes=[],map_point=lambda face,p:p)
