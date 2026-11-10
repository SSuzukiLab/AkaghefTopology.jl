"""Output defect checks, deliberately outside the default drawing path.

Source: specification sections 1.2 and 6.3; mock3.py's planarized face graph.
In particular, a Tutte solve must not collapse parallel edges to one segment.
"""
import itertools
import math


def audit(details,words):
    failures=[]; events={}; segments=[]
    for a in details['arcs']:
        pts=a['events_model']; mids=a['mids']
        if len(pts)!=len(mids)+2:
            failures.append(f"{a['face']}: {a['curve']} event count");continue
        for i,m in enumerate(mids,1):
            local=a['poly_model'] if 'event_indices' in a else pts
            k=a['event_indices'][i] if 'event_indices' in a else i
            tangent=(local[k+1][0]-local[k-1][0],local[k+1][1]-local[k-1][1])
            sign=next(x for x in words[a['curve']] if abs(x)==m)
            events.setdefault(m,[]).append((pts[i],tangent,sign,a['face']))
        actual=a['poly_model'] if 'event_indices' in a else pts
        for p,q in zip(actual,actual[1:]):
            segments.append((p,q,a['face'],a['curve']))
    cross=lambda u,v:u[0]*v[1]-u[1]*v[0]
    sub=lambda u,v:(u[0]-v[0],u[1]-v[1])
    for m,visits in events.items():
        if len(visits)!=2:
            failures.append(f'crossing {m}: {len(visits)} visits');continue
        a,b=sorted(visits,key=lambda v:-v[2])
        # Phi sends the disk interior to the outside annulus; radial distance
        # decreases, so it reverses orientation in Finf (mock3.py::Phi).
        orientation=-1 if a[3]=='Finf' else 1
        if math.dist(a[0],b[0])>1e-6 or orientation*cross(a[1],b[1])<=1e-9:
            failures.append(f'{a[3]}: crossing {m} position/sign/degeneracy')
    for (a,b,f,c),(p,q,g,d) in itertools.combinations(segments,2):
        if f!=g: continue
        u,v=sub(b,a),sub(q,p); length=math.hypot(*u)
        if length<1e-8: failures.append(f'{f}: zero-length segment on {c}');continue
        if abs(cross(u,v))<1e-8*max(1,length*math.hypot(*v)) and abs(cross(u,sub(p,a)))<1e-8*max(1,length):
            ts=[(sub(z,a)[0]*u[0]+sub(z,a)[1]*u[1])/length**2 for z in (p,q)]
            if min(1,max(ts))-max(0,min(ts))>1e-7:
                failures.append(f'{f}: overlapping segments {c}/{d}')
    return list(dict.fromkeys(failures))
