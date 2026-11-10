"""Empty-word circles in named rotation faces (specification section 3.1).

An oriented edge midpoint has its named face on its left. A disk smaller
than the distance to the rest of the embedded graph gives a local placement
inside that face. This is geometric clearance, not mesh/path optimization.
"""
import math
import lc_input as lc


def distance(p,a,b):
    v=(b[0]-a[0],b[1]-a[1]); den=v[0]**2+v[1]**2
    t=max(0.,min(1.,((p[0]-a[0])*v[0]+(p[1]-a[1])*v[1])/den)) if den else 0.
    return math.dist(p,(a[0]+t*v[0],a[1]+t*v[1]))


def store_intervals(intervals,name,labels,path):
    """Split an oriented K path at evenly spaced labels, keeping corners."""
    lengths=[math.dist(a,b) for a,b in zip(path,path[1:])]; total=sum(lengths)
    stations=[0.]
    for d in lengths: stations.append(stations[-1]+d)
    def point(t):
        for k,(s,d) in enumerate(zip(stations,lengths)):
            if t<=s+d+1e-9:
                u=(t-s)/d; a,b=path[k:k+2]
                return (a[0]+u*(b[0]-a[0]),a[1]+u*(b[1]-a[1]))
        return path[-1]
    for j,label in enumerate(labels):
        a,b=total*j/len(labels),total*(j+1)/len(labels)
        points=[point(a)]+[p for p,s in zip(path[1:-1],stations[1:-1]) if a<s<b]+[point(b)]
        intervals.setdefault((name,label),[]).append(points)


def add_empty(data,doc,details):
    words=data['words']; empty=[c for c,w in words.items() if not w]
    if not empty: return
    intervals=details['intervals']; map_point=details['map_point']
    for a in details['arcs']:
        labels=[a['start']]+a['mids']; events=a['events_model']
        for j,(label,p,q) in enumerate(zip(labels,events,events[1:])):
            signed=next(x for x in words[a['curve']] if abs(x)==label)
            count=300 if a['face']=='Finf' else 1
            # Subdivision vertices belong to the same original word interval.
            model=a.get('event_paths_model',[])
            if model:
                path=[]
                for p,q in zip(model[j],model[j][1:]):
                    path.extend(map_point(a['face'],(p[0]+(q[0]-p[0])*k/count,p[1]+(q[1]-p[1])*k/count)) for k in range(count))
                path.append(map_point(a['face'],q))
            else:
                path=[map_point(a['face'],(p[0]+(q[0]-p[0])*k/count,p[1]+(q[1]-p[1])*k/count)) for k in range(count+1)]
            intervals.setdefault((a['curve'],signed),[]).append(path)
    segments=[]
    for paths in intervals.values():
        for path in paths: segments.extend(zip(path,path[1:]))
    for curve in doc['curves'].values():
        pts=curve['points'];segments.extend(zip(pts,pts[1:]+([pts[0]] if curve['closed'] else [])))
    faces,_=lc.face_names({c:w for c,w in words.items() if w})
    occ=lc.rotation_occ({c:w for c,w in words.items() if w})
    circles={}; circle_disks={}; inside={}; required=set()
    def allowed(center):
        return all((math.dist(center,p)<r)==(name in required)
                   for name,(p,r) in circle_disks.items())
    def free_disk(paths):
        candidates=[]
        for path in paths:
            for a,b in zip(path,path[1:]):
                length=math.dist(a,b)
                if length<1e-7: continue
                for t in (.5,.38196601125,.61803398875):
                    p=(a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1]))
                    distances=[distance(p,c,d) for c,d in segments]
                    clearance=min([length*t,length*(1-t)]+[d for d in distances if d>1e-7])
                    offset=clearance/2
                    center=(p[0]-(b[1]-a[1])*offset/length,p[1]+(b[0]-a[0])*offset/length)
                    if not allowed(center): continue
                    if doc.get('square') and not (0<center[0]<doc['W'] and 0<center[1]<doc['H']): continue
                    if any(math.dist(center,c)<=r for c,r in doc['holes'].values()): continue
                    # Triangle cut faces exclude their interiors, even though
                    # the final hole is displayed as its smaller incircle.
                    inside=False
                    for tri in details.get('holes',[]):
                        cross=[(v[0]-u[0])*(center[1]-u[1])-(v[1]-u[1])*(center[0]-u[0]) for u,v in zip(tri,tri[1:]+tri[:1])]
                        if min(cross)>=0 or max(cross)<=0: inside=True;break
                    if inside: continue
                    radius=min([offset*.8]+[distance(center,c,d)*.8 for c,d in segments])
                    if radius>1e-7: candidates.append((radius,center))
        if not candidates: raise ValueError('no nondegenerate local disk in named face')
        return max(candidates,key=lambda x:x[0])[::-1]
    def face_paths(fid):
        paths=[]
        for label,sign,way in faces[fid]:
            c,i=occ[label][sign]
            w=words[c]
            if way=='out': paths.extend(intervals.get((c,w[i]),[]))
            else: paths.extend([list(reversed(p)) for p in intervals.get((c,w[(i-1)%len(w)]),[])])
        return paths
    pending=list(empty)
    while pending:
        advanced=False
        for c in pending[:]:
            ref=data['places'][c]
            parent=ref[:-1] if ref.endswith(('+','-')) else None
            if parent in empty and parent not in circles: continue
            required=set(inside.get(parent,set()))
            if parent in circles and ref.endswith('+'): required.add(parent)
            if ref=='outer':
                if circles: center,radius=free_disk([list(reversed(p)) for p in circles.values()])
                else: center,radius=(260.,260.),160.
            elif ref in faces: center,radius=free_disk(face_paths(ref))
            else:
                parent,side=ref[:-1],ref[-1]
                if parent in circles:
                    path=circles[parent]
                    if side=='+':
                        # Pack direct children inside the parent disk, without
                        # making a later sibling's size depend on polygon edges.
                        siblings=[name for name in empty if data['places'][name]==ref]
                        origin,extent=circle_disks[parent]; n=len(siblings)
                        angle=2*math.pi*siblings.index(c)/n
                        center=origin if n==1 else (origin[0]+extent*.5*math.cos(angle),origin[1]+extent*.5*math.sin(angle))
                        radius=min(extent/3 if n==1 else extent*math.sin(math.pi/n)/3,
                                   min(distance(center,a,b) for a,b in segments)*.8)
                        if radius<1e-7 or not allowed(center): center,radius=free_disk([path])
                    else: center,radius=free_disk([list(reversed(path))])
                elif parent in words and words[parent]:
                    paths=[p for (name,label),parts in intervals.items() if name==parent for p in parts]
                    center,radius=free_disk(paths if side=='+' else [list(reversed(p)) for p in paths])
                else: continue
            pts=[(center[0]+radius*math.cos(2*math.pi*k/64),center[1]+radius*math.sin(2*math.pi*k/64)) for k in range(64)]
            doc['curves'][c]=dict(points=[list(p) for p in pts],closed=True)
            circle_disks[c]=(center,radius);inside[c]=set(required)
            circles[c]=pts+[pts[0]];segments.extend(zip(circles[c],circles[c][1:]))
            pending.remove(c);advanced=True
        if not advanced: raise ValueError(f'unsupported placement references: {pending}')
