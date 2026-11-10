"""Keep mock3.py's chord arithmetic; subdivide only the fallback faces.
cf:docs/heegaard/261109_tutte_construction.md 弦分岐の出典 — unchanged mock chord arithmetic
The Tutte construction and its source comments live in lc_tutte.py.
"""
import itertools


def solve_faces(arcs, LC, CU, faces, boundaries, meshes):
    def csign(c,lab):return next(x for x in LC[c] if abs(x)==lab)>0
    def seg_int(a,b,c,d):
        r=(b[0]-a[0],b[1]-a[1]);s=(d[0]-c[0],d[1]-c[1]);den=r[0]*s[1]-r[1]*s[0]
        if abs(den)<1e-12:return None
        q=(c[0]-a[0],c[1]-a[1]);u=(q[0]*s[1]-q[1]*s[0])/den;v=(q[0]*r[1]-q[1]*r[0])/den
        return (u,v) if 0<u<1 and 0<v<1 else None
    def other_curve(a,m):return [c for c in CU if m in [abs(x) for x in LC[c]] and c!=a['curve']][0]
    def chord_check(face):
        fa=[a for a in arcs if a['face']==face];found={}
        for a,b in itertools.combinations(fa,2):
            r=seg_int(a['mP'],a['mQ'],b['mP'],b['mQ'])
            if r:
                u,v=r;found.setdefault(id(a),[]).append((u,b));found.setdefault(id(b),[]).append((v,a))
        ok=True
        for a in fa:
            seq=[x[1]['curve'] for x in sorted(found.get(id(a),[]),key=lambda x:x[0])]
            exp=[other_curve(a,m) for m in a['mids']]
            if seq!=exp:ok=False;print('  chord mismatch',face,a['curve'],a['start'],a['end'],seq,exp)
        return ok
    method={}
    for face in faces:
        if chord_check(face):
            method[face]='chords'
            for a in arcs:
                if a['face']==face:a['poly_model']=[a['mP'],a['mQ']]
        else:
            from lc_tutte import embed_face
            fa=[a for a in arcs if a['face']==face]
            meshes[face]=embed_face(fa,LC,boundaries[face],exterior=face=='Finf')
            method[face]='tutte'
    print('methods:',method)
    # Event locations are bookkeeping for section 3.1 empty-component faces.
    # Tutte keeps the original events separately from subdivision vertices.
    for a in arcs:
        if method[a['face']]!='tutte':
            hits=[]
            for b in arcs:
                if b is a or b['face']!=a['face']: continue
                hit=seg_int(a['mP'],a['mQ'],b['mP'],b['mQ'])
                if hit: hits.append(hit[0])
            p,q=a['mP'],a['mQ']
            a['events_model']=[p]+[(p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])) for t in sorted(hits)]+[q]
    return method
