"""Inspect emitted polylines independently of the mechanical gate.

cf:docs/specs/lickorish-code.md §6.3(5) — unshared crossings and hole intrusion diagnostics, outside the default path
"""
import numpy as np

def audit_document(document):
    segments=[];owners=[]
    for name,curve in document['curves'].items():
        points=list(curve['points'])
        if curve['closed'] and np.linalg.norm(np.subtract(points[0],points[-1]))<1e-9:points.pop()
        for a,b in zip(points,points[1:]+([points[0]] if curve['closed'] else [])):
            segments.append((a,b));owners.append(name)
    if not segments:return dict(unshared_intersections=[],hole_intrusions=[],boundary_exits=[])
    seg=np.asarray(segments);a=seg[:,0];r=seg[:,1]-a
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    intersections=[]
    for i in range(len(a)):
        q=a[i+1:]-a[i];s=r[i+1:];den=cross(r[i],s)
        with np.errstate(divide='ignore',invalid='ignore'):
            u=cross(q,s)/den;v=cross(q,r[i])/den
        hits=np.nonzero((np.abs(den)>1e-12)&(u>1e-8)&(u<1-1e-8)&(v>1e-8)&(v<1-1e-8))[0]
        for k in hits:
            j=i+1+int(k)
            intersections.append(dict(curves=[owners[i],owners[j]],point=(a[i]+u[k]*r[i]).tolist()))
    intrusions=[]
    for name,(center,radius) in document['holes'].items():
        den=np.sum(r*r,axis=1)
        t=np.clip(np.sum((np.asarray(center)-a)*r,axis=1)/np.maximum(den,1e-30),0,1)
        distances=np.linalg.norm(a+t[:,None]*r-center,axis=1)
        for i in np.nonzero(distances<radius-1e-6)[0]:
            intrusions.append(dict(curve=owners[int(i)],hole=name,depth=float(radius-distances[i])))
    exits=[]
    if document.get('square'):
        for name,curve in document['curves'].items():
            for i,(x,y) in enumerate(curve['points']):
                depth=max(-x,x-document['W'],-y,y-document['H'])
                if depth>1e-6:exits.append(dict(curve=name,index=i,point=[x,y],depth=depth))
    return dict(unshared_intersections=intersections,hole_intrusions=intrusions,boundary_exits=exits)
