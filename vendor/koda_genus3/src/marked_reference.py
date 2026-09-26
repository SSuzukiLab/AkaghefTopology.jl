"""Executable reference specification for marked PL surface diagrams.
Python >=3.11, standard library only. The installation target is Julia.
A digest identifies a representation, NEVER an ambient-isotopy class.
"""
from __future__ import annotations
from collections import defaultdict, Counter
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import argparse, copy, hashlib, html, json, math, re, tomllib

class DiagramError(ValueError):
    def __init__(self, code: str, detail: str):
        self.code, self.detail = code, detail
        super().__init__(f'{code}: {detail}')

def need(ok, code, detail):
    if not ok: raise DiagramError(code, detail)

def edge(u, v): return tuple(sorted((u, v)))
def directed_triangle(t): return [(t[i],t[(i+1)%3]) for i in range(3)]
def components(nodes, adj):
    unseen=set(nodes); out=[]
    while unseen:
        seed=min(unseen); unseen.remove(seed); comp={seed}; stack=[seed]
        while stack:
            v=stack.pop()
            for w in adj.get(v, []):
                if w in unseen: unseen.remove(w); comp.add(w); stack.append(w)
        out.append(comp)
    return out

class DSU:
    def __init__(self, nodes): self.p={v:v for v in nodes}
    def find(self, x):
        p=self.p[x]
        if p!=x: self.p[x]=self.find(p)
        return self.p[x]
    def union(self, a,b):
        a,b=self.find(a),self.find(b)
        if a!=b: self.p[max(a,b)]=min(a,b)

@dataclass(frozen=True)
class Dart:
    face:str
    u:str
    v:str

@dataclass
class State:
    doc:dict
    patches:dict
    boundaries:dict
    incidence:dict
    triangles:dict
    seam_partner:dict
    vertex_dsu:DSU
    curves:dict
    used:dict
    crossing_at:dict
    surface_report:dict
    fingerprint:str

ID=r'[^\W\d][\w.-]*'
BLOCK=re.compile(rf'({ID})\{{\s*({ID}(?:\s+{ID})+)\s*\}}',re.UNICODE)
def parse_word(s:str)->list[Dart]:
    """C{a b x c} | D{q r}: named paths in prescribed triangulated disks."""
    need(isinstance(s,str) and s.strip(), 'WORD_EMPTY', repr(s))
    out=[]
    for block in s.split('|'):
        m=BLOCK.fullmatch(block.strip())
        need(m is not None,'WORD_SYNTAX',block)
        f=m.group(1); vs=m.group(2).split()
        out.extend(Dart(f,u,v) for u,v in zip(vs,vs[1:]))
    return out

def print_word(ds:list[Dart])->str:
    blocks=[]
    for d in ds:
        if blocks and blocks[-1][0]==d.face and blocks[-1][1][-1]==d.u:
            blocks[-1][1].append(d.v)
        else: blocks.append([d.face,[d.u,d.v]])
    return ' | '.join(f'{f}{{{" ".join(vs)}}}' for f,vs in blocks)

def patch_info(p):
    fid=p['face']; order=p['side_order']; sides=p['sides']
    need(len(order)>=3 and len(set(order))==len(order),'PATCH_SIDES',fid)
    need(set(order)==set(sides),'PATCH_SIDES',fid)
    boundary=[]
    for i,s in enumerate(order):
        vs=sides[s]
        need(len(vs)>=2 and len(vs)==len(set(vs)),'SIDE_PATH',f'{fid}/{s}')
        need(vs[-1]==sides[order[(i+1)%len(order)]][0],'SIDE_CORNERS',f'{fid}/{s}')
        boundary.extend(vs[:-1])
    need(len(set(boundary))==len(boundary),'BOUNDARY_NOT_SIMPLE',fid)
    inc=defaultdict(list); tris={}; vertices=set()
    for i,t in enumerate(p['triangles']):
        t=tuple(t); need(len(t)==3 and len(set(t))==3,'TRIANGLE_INVALID',f'{fid}/{i}')
        need(tuple(sorted(t)) not in tris,'TRIANGLE_DUPLICATE',fid)
        tris[tuple(sorted(t))]=i
        vertices.update(t)
        for u,v in directed_triangle(t): inc[edge(u,v)].append((i,u,v))
    need(vertices>=set(boundary),'MISSING_BOUNDARY_VERTEX',fid)
    bd={edge(u,v):(u,v) for u,v in zip(boundary,boundary[1:]+boundary[:1])}
    for e,occ in inc.items():
        if e in bd:
            need(len(occ)==1 and occ[0][1:]==bd[e],'BOUNDARY_ORIENTATION',f'{fid}/{e}')
        else:
            need(len(occ)==2 and occ[0][1:]==tuple(reversed(occ[1][1:])),
                 'EDGE_MANIFOLD',f'{fid}/{e}')
    need(set(bd)<=set(inc),'BOUNDARY_EDGE_MISSING',fid)
    adj=defaultdict(set)
    for e,occ in inc.items():
        if len(occ)==2:
            a,b=occ[0][0],occ[1][0]; adj[a].add(b);adj[b].add(a)
    need(len(components(range(len(p['triangles'])),adj))==1,'PATCH_DISCONNECTED',fid)
    need(len(vertices)-len(inc)+len(p['triangles'])==1,'PATCH_NOT_DISK',fid)
    # Link tests rule out pinched pseudomanifolds despite a plausible Euler number.
    for v in vertices:
        link=defaultdict(list)
        for t in p['triangles']:
            if v in t:
                a,b=[w for w in t if w!=v];link[a].append(b);link[b].append(a)
        need(len(components(link,link))==1,'VERTEX_LINK_DISCONNECTED',f'{fid}/{v}')
        degrees=sorted(len(x) for x in link.values())
        if v in boundary:
            need(degrees.count(1)==2 and all(x in (1,2) for x in degrees),
                 'BOUNDARY_LINK',f'{fid}/{v}')
        else: need(all(x==2 for x in degrees),'INTERIOR_LINK',f'{fid}/{v}')
    return boundary,dict(inc),vertices

def check_keys(x, required, optional=()):
    need(isinstance(x,dict),'SCHEMA_OBJECT',repr(x))
    need(set(required)<=set(x) and set(x)<=set(required)|set(optional),'SCHEMA_KEYS',f'keys={list(x)} required={required}')

def check_shape(doc):
    check_keys(doc,('version','surface','patches','curves','crossings'),('metadata','region_marks'))
    check_keys(doc['surface'],('marking_id','genus','gluings','boundary'))
    for p in doc['patches']:
        check_keys(p,('face','side_order','sides','triangles'),('side_parameters',))
    for c in doc['curves']:check_keys(c,('id','role','closed','word'))
    for x in doc['crossings']:check_keys(x,('id','at','curves','sign'))
    for b in doc['surface']['boundary']:check_keys(b,('side','label'))
    for r in doc.get('region_marks',[]):
        check_keys(r,('id','seed','claim'))
        z=r['seed'];need(isinstance(z,list) and len(z)==2 and isinstance(z[0],str) and type(z[1]) is int and z[1]>=0,'REGION_SEED',str(z))

def check_state(st):
    need(representation_digest(st.doc)==st.fingerprint,'MATHEMATICAL_STATE_MUTATED','revalidate the edited document')

def _validate(doc:dict)->State:
    check_shape(doc);doc=copy.deepcopy(doc)
    need(type(doc.get('version')) is int and doc.get('version')==1,'VERSION','expected version=1')
    surf=doc['surface']; need(bool(surf.get('marking_id')),'MARKING_REQUIRED','surface.marking_id')
    patches={p['face']:p for p in doc['patches']}
    need(len(patches)==len(doc['patches']),'DUPLICATE_FACE','patches')
    need(bool(patches),'SURFACE_EMPTY','patches')
    need(len(patches)<=64,'RESOURCE_LIMIT','>64 patches')
    bds={}; incs={}; nodes=[]; tri={}; sidepaths={}; corners=set()
    for f,p in patches.items():
        bd,inc,vs=patch_info(p); bds[f]=bd;incs[f]=inc
        need(len(vs)<=256,'RESOURCE_LIMIT',f'{f}: >256 vertices in reference exact kernel')
        supplied=p.get('side_parameters',{})
        need(set(supplied)<=set(p['sides']),'SIDE_PARAMETERS',f)
        for sid,path in p['sides'].items():
            ts=[qnum(x) for x in supplied[sid]] if sid in supplied else [Q(k,len(path)-1) for k in range(len(path))]
            need(len(ts)==len(path) and ts[0]==0 and ts[-1]==1 and all(a<b for a,b in zip(ts,ts[1:])), 'SIDE_PARAMETERS',f'{f}/{sid}')
        nodes.extend((f,v) for v in vs)
        tri.update({(f,i):tuple(t) for i,t in enumerate(p['triangles'])})
        for sid in p['side_order']:
            sidepaths[(f,sid)]=p['sides'][sid]
            corners.add((f,p['sides'][sid][0]))
    dsu=DSU(nodes); partner={}; faceadj=defaultdict(set); paired=set()
    # Local fine boundary edge -> reversed local edge in the glued side.
    seam={}; corneradj=defaultdict(list)
    for a,b in surf.get('gluings',[]):
        a,b=tuple(a),tuple(b)
        need(a in sidepaths and b in sidepaths and a!=b,'GLUING_SIDE',f'{a},{b}')
        need(a not in paired and b not in paired,'GLUING_DUPLICATE',f'{a},{b}')
        paired.update((a,b));partner[a]=b;partner[b]=a
        va,vb=sidepaths[a],sidepaths[b]
        need(len(va)==len(vb),'GATE_COUNT',f'{a},{b}')
        pa=side_parameters(patches[a[0]],a[1]);pb=side_parameters(patches[b[0]],b[1])
        need(pa==[1-t for t in reversed(pb)],'GATE_PARAMETERS',f'{a},{b}')
        for u,v in zip(va,reversed(vb)): dsu.union((a[0],u),(b[0],v))
        for ia in range(len(va)-1):
            ka=(a[0],edge(va[ia],va[ia+1]));ib=len(vb)-2-ia
            kb=(b[0],edge(vb[ib],vb[ib+1])); seam[ka]=kb;seam[kb]=ka
        # Corner sectors are link edges; these joins connect sector ends.
        for u,v in [(va[0],vb[-1]),(va[-1],vb[0])]:
            ca,cb=(a[0],u),(b[0],v);corneradj[ca].append(cb);corneradj[cb].append(ca)
        faceadj[a[0]].add(b[0]);faceadj[b[0]].add(a[0])
    need(len(components(patches,faceadj))==1,'SURFACE_DISCONNECTED','base cells')
    unpaired=set(sidepaths)-paired
    labels={tuple(x['side']):x['label'] for x in surf.get('boundary',[])}
    need(set(labels)==unpaired and len(labels)==len(surf.get('boundary',[])),
         'BOUNDARY_LABELS','label each unglued side exactly once')
    # Every quotient corner link must be a circle or an interval, never a wedge.
    groups=defaultdict(set)
    for c in corners: groups[dsu.find(c)].add(c)
    for root,cs in groups.items():
        need(len(components(cs,corneradj))==1,'SURFACE_VERTEX_LINK',str(root))
        degrees=[len(corneradj[c]) for c in cs]
        defect=sum(2-x for x in degrees)
        need(all(x<=2 for x in degrees) and defect in (0,2),'SURFACE_VERTEX_LINK',str(root))
    badj=defaultdict(list);bv=set();side_ends={}
    for s in unpaired:
        a,b=dsu.find((s[0],sidepaths[s][0])),dsu.find((s[0],sidepaths[s][-1]))
        badj[a].append(b);badj[b].append(a);bv.update((a,b));side_ends[s]=(a,b)
    need(all(len(badj[v])==2 for v in bv),'SURFACE_BOUNDARY_LINK','boundary')
    bc=components(bv,badj); usedlabels=[]
    for c in bc:
        names={labels[s] for s,ends in side_ends.items() if ends[0] in c}
        need(len(names)==1,'BOUNDARY_COMPONENT_LABEL','one label per component')
        usedlabels.extend(names)
    need(len(set(usedlabels))==len(usedlabels),'BOUNDARY_LABEL_REUSED','distinct components')
    V=len(groups); E=len(sidepaths)-len(paired)//2;F=len(patches);chi=V-E+F
    g=surf['genus'];need(isinstance(g,int) and not isinstance(g,bool) and g>=0 and chi==2-2*g-len(bc),
                       'SURFACE_GENUS',f'chi={chi}, g={g}, b={len(bc)}')
    report={'genus':g,'boundary_components':len(bc),'chi':chi,'faces':F}
    curves={}; used=defaultdict(list); visits=defaultdict(dict)
    def globalv(f,v): return dsu.find((f,v))
    physical=set()
    for s in unpaired: physical.update((s[0],v) for v in sidepaths[s][1:-1])
    physical_global={dsu.find(v) for v in physical}
    for c in doc.get('curves',[]):
        cid=c['id'];need(type(c['closed']) is bool,'CURVE_CLOSED','closed must be boolean')
        need(cid not in curves,'CURVE_DUPLICATE',cid)
        need(c['role'] in ('alpha','beta','aux'),'CURVE_ROLE',cid)
        ds=parse_word(c['word']); curves[cid]=(c,ds)
        localvertices=[]
        for d in ds:
            need(d.face in incs and edge(d.u,d.v) in incs[d.face],
                 'ROUTE_EDGE_MISSING',str(d))
            need(len(incs[d.face][edge(d.u,d.v)])==2,'ROUTE_ON_MARKING',str(d))
            need((d.face,d.u) not in corners and (d.face,d.v) not in corners,
                 'ROUTE_AT_CORNER',str(d))
            key=(d.face,edge(d.u,d.v));used[key].append(cid)
            localvertices.extend([(d.face,d.u),(d.face,d.v)])
        # A path may leave a disk only via its declared paired gate.
        pairs=list(zip(ds,ds[1:]))+([(ds[-1],ds[0])] if c['closed'] else [])
        for a,b in pairs:
            need(globalv(a.face,a.v)==globalv(b.face,b.u),'ROUTE_DISCONTINUITY',cid)
            if (a.face,a.v)!=(b.face,b.u):
                need(a.v in bds[a.face] and b.u in bds[b.face],'ILLEGAL_CHART_JUMP',cid)
            else:
                need(a.v not in bds[a.face],'GATE_TURNBACK',f'{cid}/{a.v}')
        if not c['closed']:
            need(globalv(ds[0].face,ds[0].u) in physical_global and
                 globalv(ds[-1].face,ds[-1].v) in physical_global,'ARC_ENDPOINT',cid)
        ordered=[globalv(d.face,d.u) for d in ds]
        if not c['closed']: ordered.append(globalv(ds[-1].face,ds[-1].v))
        need(len(ordered)==len(set(ordered)),'SELF_INTERSECTION',cid)
        # A crossing has two oriented neighbors, both in one disk.
        for a,b in pairs:
            root=globalv(a.face,a.v)
            need(cid not in visits[root],'CURVE_REPEATED_VERTEX',cid)
            visits[root][cid]=(a,b)
        if not c['closed']:
            visits[ordered[0]][cid]=(None,ds[0]);visits[ordered[-1]][cid]=(ds[-1],None)
    need(all(len(v)==1 for v in used.values()),'SHARED_ROUTE_EDGE','curves overlap along an edge')
    declared={}; crossingids=set()
    for x in doc.get('crossings',[]):
        fid,v=x['at']; key=(fid,v)
        need(key in dsu.p,'CROSSING_VERTEX',str(key))
        need(v not in bds[fid],'CROSSING_ON_MARKING',str(key))
        need(x['id'] not in crossingids and key not in declared,'CROSSING_DUPLICATE',str(key))
        crossingids.add(x['id']);declared[key]=x
    seen=set()
    for root,at in visits.items():
        if len(at)==1:
            need(root not in declared,'FALSE_CROSSING',str(root));continue
        need(len(at)==2,'TRIPLE_POINT',str(root))
        # Interior vertices are not glued, hence root is their local key.
        need(root in declared,'UNDECLARED_INTERSECTION',str(root));x=declared[root];seen.add(root)
        first,second=x['curves'];need(set((first,second))==set(at),'CROSSING_CURVES',x['id'])
        r1,r2=curves[first][0]['role'],curves[second][0]['role']
        need(not(r1==r2 and r1 in ('alpha','beta')),'SAME_FAMILY_INTERSECTION',x['id'])
        a0,a1=at[first];b0,b1=at[second]
        need(all(d is not None for d in (a0,a1,b0,b1)),'ENDPOINT_INTERSECTION',x['id'])
        links={}
        for t in patches[root[0]]['triangles']:
            if root[1] in t:
                i=t.index(root[1]);links[t[(i+1)%3]]=t[(i+2)%3]
        ao,ai,bo,bi=a1.v,a0.u,b1.v,b0.u
        order=[];v=ao
        while not order or v!=ao:
            order.append(v);v=links[v]
            need(len(order)<=len(links),'CROSSING_LINK',x['id'])
        reduced=[v for v in order if v in (ao,ai,bo,bi)]
        need(reduced in ([ao,bo,ai,bi],[ao,bi,ai,bo]),'NONTRANSVERSE_INTERSECTION',x['id'])
        sign=1 if reduced==[ao,bo,ai,bi] else -1
        need(type(x['sign']) is int and x['sign']==sign,'CROSSING_SIGN',f'{x["id"]}: expected {sign}')
    need(seen==set(declared),'FALSE_CROSSING','unused crossing declaration')
    st=State(doc,patches,bds,incs,tri,seam,dsu,curves,dict(used),declared,report,representation_digest(doc))
    ids=set()
    for mark in doc.get('region_marks',[]):
        need(mark['id'] not in ids,'REGION_MARK_DUPLICATE',mark['id']);ids.add(mark['id'])
        need(mark['claim'] in ('region','embedded_closed_disk','koda_joining'),'REGION_CLAIM',mark['claim'])
        audit=check_joining_disk(st,mark['seed'])
        need(mark['claim']=='region' or audit[mark['claim']],'REGION_CLAIM_FAILED',mark['id'])
    return st

def validate(doc:dict)->State:
    """Public fail-closed boundary; unexpected structural inputs are schema errors."""
    try:return _validate(doc)
    except DiagramError:raise
    except (KeyError,TypeError,IndexError,ValueError,AttributeError) as e:
        raise DiagramError('SCHEMA_SHAPE',str(e)) from e

def side_parameters(p,sid):
    xs=p.get('side_parameters',{}).get(sid)
    n=len(p['sides'][sid])-1
    return [qnum(x) for x in xs] if xs is not None else [Q(k,n) for k in range(n+1)]

def endpoint_anchors(st:State,cid:str):
    check_state(st)
    c,ds=st.curves[cid]
    if c['closed']:return []
    out=[]
    for f,v in [(ds[0].face,ds[0].u),(ds[-1].face,ds[-1].v)]:
        p=st.patches[f];sid=next(s for s in p['side_order'] if v in p['sides'][s][1:-1])
        i=p['sides'][sid].index(v)
        out.append({'side':[f,sid], 'parameter':str(side_parameters(p,sid)[i])})
    return out

def qnum(x):
    need(isinstance(x,(str,int)) and not isinstance(x,bool),'RATIONAL_LITERAL',repr(x))
    need(re.fullmatch(r'-?[0-9]+(?:/[1-9][0-9]*)?',str(x)) is not None,'RATIONAL_LITERAL',repr(x))
    return Q(x)
def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
def onsegment(a,b,p):
    return orient(a,b,p)==0 and min(a[0],b[0])<=p[0]<=max(a[0],b[0]) and min(a[1],b[1])<=p[1]<=max(a[1],b[1])
def segments_meet(a,b,c,d):
    s,t,u,v=orient(a,b,c),orient(a,b,d),orient(c,d,a),orient(c,d,b)
    return (s*t<0 and u*v<0) or any((z==0 and onsegment(p,q,r)) for z,p,q,r in
                                 [(s,a,b,c),(t,a,b,d),(u,c,d,a),(v,c,d,b)])

def circle_boundary(n):
    """n strictly convex, rational unit-circle points in CCW order."""
    out=[]
    for i in range(n):
        quadrant=(4*i)//n;t=Q((4*i)%n,n)
        x,y=(1-t*t)/(1+t*t),2*t/(1+t*t)
        for _ in range(quadrant):x,y=-y,x
        out.append((x,y))
    return out

def solve_exact(A,B):
    n=len(A)
    if n==0:return []
    M=[list(map(Q,A[i]))+list(map(Q,B[i])) for i in range(n)]
    for i in range(n):
        k=next((k for k in range(i,n) if M[k][i]),None)
        need(k is not None,'LAYOUT_SINGULAR','Dirichlet matrix')
        M[i],M[k]=M[k],M[i];pivot=M[i][i];M[i]=[z/pivot for z in M[i]]
        for j in range(n):
            if i!=j and M[j][i]:
                a=M[j][i];M[j]=[x-a*y for x,y in zip(M[j],M[i])]
    return [r[n:] for r in M]

def certify_layout(st:State, positions:dict):
    check_state(st)
    need(set(positions)==set(st.patches),'LAYOUT_FACES','exact same face set required')
    for f,p in st.patches.items():
        pts=positions[f]
        need(all(isinstance(p,tuple) and len(p)==2 and all(isinstance(x,Q) for x in p) for p in pts.values()),'LAYOUT_RATIONAL',f)
        allv={v for t in p['triangles'] for v in t}
        need(set(pts)==allv,'LAYOUT_VERTICES',f)
        need(len(set(pts.values()))==len(pts),'COINCIDENT_VERTICES',f)
        for i,t in enumerate(p['triangles']):
            need(orient(*(pts[v] for v in t))>0,'FLIPPED_TRIANGLE',f'{f}/{i}')
        es=list(st.incidence[f])
        for i,(a,b) in enumerate(es):
            for c,d in es[i+1:]:
                common=set((a,b))&set((c,d))
                if not common:
                    need(not segments_meet(pts[a],pts[b],pts[c],pts[d]),'MESH_EDGE_CROSSING',f'{f}/{a},{b}/{c},{d}')
                else:
                    z=next(iter(common));u=next(v for v in (a,b) if v!=z);v=next(v for v in (c,d) if v!=z)
                    need(not onsegment(pts[z],pts[u],pts[v]) and not onsegment(pts[z],pts[v],pts[u]),
                         'OVERLAPPING_EDGES',f'{f}/{z}')
    return True

def realize(st:State)->dict:
    check_state(st)
    out={}
    for f,p in st.patches.items():
        bd=st.boundaries[f];pts=dict(zip(bd,circle_boundary(len(bd))))
        adj=defaultdict(set)
        for u,v in st.incidence[f]:adj[u].add(v);adj[v].add(u)
        interior=sorted(set(adj)-set(bd));index={v:i for i,v in enumerate(interior)}
        n=len(interior);A=[[Q(0) for _ in range(n)] for _ in range(n)];B=[[Q(0),Q(0)] for _ in range(n)]
        for v,i in index.items():
            A[i][i]=len(adj[v])
            for w in adj[v]:
                if w in index:A[i][index[w]]-=1
                else:
                    for j in (0,1):B[i][j]+=pts[w][j]
        for v,row in zip(interior,solve_exact(A,B)):pts[v]=tuple(row)
        out[f]=pts
    certify_layout(st,out)
    return out

def representation_digest(doc):
    # Exact document hash, not isomorphism- or isotopy-canonicalization.
    semantic={k:v for k,v in doc.items() if k not in ('metadata',)}
    return hashlib.sha256(json.dumps(semantic,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()

def marking_digest(doc):
    obj={'surface':doc['surface'], 'faces':[{ 'id':p['face'],'side_order':p['side_order']} for p in doc['patches']]}
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()

def seam_word(st:State,cid:str):
    check_state(st)
    c,ds=st.curves[cid];pair={}
    for a,b in st.doc['surface'].get('gluings',[]):pair[tuple(a)]=tuple(b);pair[tuple(b)]=tuple(a)
    token=[];pairs=list(zip(ds,ds[1:]))+([(ds[-1],ds[0])] if c['closed'] else [])
    for a,b in pairs:
        if (a.face,a.v)==(b.face,b.u):continue
        p=st.patches[a.face]
        side=next(s for s in p['side_order'] if a.v in p['sides'][s][1:-1])
        s=(a.face,side);t=pair[s];base=min(s,t);sign=1 if s==base else -1
        vs=st.patches[base[0]]['sides'][base[1]]
        # Base-side rank is independent of which copy was traversed.
        k=p['sides'][side].index(a.v)
        if s!=base:k=len(vs)-1-k
        token.append({'side':list(base),'direction':sign,'gate_rank':k})
    return token

def intersection_cycles(st:State):
    check_state(st)
    ans={}
    for cid,(c,ds) in st.curves.items():
        seq=[]
        for d in ds:
            x=st.crossing_at.get((d.face,d.u))
            if x:seq.append({'id':x['id'],'sign':x['sign'],'curves':x['curves']})
        ans[cid]=seq
    return ans

def legacy_projection(st:State):
    """Coordinates-free compatibility projection; never an inverse migration."""
    check_state(st)
    need(all(c['role'] in ('alpha','beta') and c['closed'] for c,_ in st.curves.values()),'LEGACY_NOT_HEEGAARD','closed alpha/beta curves only')
    events=[]
    for x in st.doc['crossings']:
        a,b=x['curves'];sign=x['sign']
        if st.curves[a][0]['role']=='beta':a,b=b,a;sign=-sign
        need(st.curves[a][0]['role']=='alpha' and st.curves[b][0]['role']=='beta','LEGACY_NOT_HEEGAARD',x['id'])
        events.append({'id':x['id'],'alpha':a,'beta':b,'sign':sign})
    cyc=intersection_cycles(st)
    return {'genus':st.surface_report['genus'],
      'alpha':[{'id':c['id'],'cycle':[x['id'] for x in cyc[c['id']]]} for c in st.doc['curves'] if c['role']=='alpha'],
      'beta':[{'id':c['id'],'cycle':[x['id'] for x in cyc[c['id']]]} for c in st.doc['curves'] if c['role']=='beta'],
      'intersections':events,'source_representation_sha256':st.fingerprint,'marking_sha256':marking_digest(st.doc)}

def fine_edge_occurrences(st):
    """Pair triangle-edge occurrences, including seams. Orientation is reversed."""
    pairs={};one={}
    for f,inc in st.incidence.items():
        for e,occ in inc.items():
            ids=[(f,o[0],o[1],o[2]) for o in occ]
            if len(ids)==2:pairs[ids[0]]=ids[1];pairs[ids[1]]=ids[0]
            else:one[(f,e)]=ids[0]
    for a,b in st.seam_partner.items():pairs[one[a]]=one[b]
    return pairs

def regions(st:State, barriers=None):
    """Normalized cut-complement, NOT the potentially pinched image closure."""
    check_state(st)
    blocked=set(st.used) if barriers is None else set(barriers)
    epairs=fine_edge_occurrences(st);adj=defaultdict(set)
    def cut(o):return (o[0],edge(o[2],o[3])) in blocked
    for a,b in epairs.items():
        if not cut(a) and not cut(b):adj[a[:2]].add(b[:2])
    comps=components(st.triangles,adj);out=[]
    for comp in comps:
        cn=[(f,i,j) for f,i in comp for j in range(3)]; vd=DSU(cn)
        eo=[(f,i,u,v) for f,i in comp for u,v in directed_triangle(st.triangles[(f,i)])]
        ed=DSU(eo)
        for a in eo:
            b=epairs.get(a)
            if b is not None and b[:2] in comp and not cut(a) and not cut(b):
                ed.union(a,b)
                ta,tb=st.triangles[a[:2]],st.triangles[b[:2]]
                vd.union((a[0],a[1],ta.index(a[2])),(b[0],b[1],tb.index(b[3])))
                vd.union((a[0],a[1],ta.index(a[3])),(b[0],b[1],tb.index(b[2])))
        bd=[a for a in eo if a not in epairs or cut(a) or cut(epairs[a])]
        bnext={};bverts=set()
        for a in bd:
            t=st.triangles[a[:2]];u=vd.find((a[0],a[1],t.index(a[2])));v=vd.find((a[0],a[1],t.index(a[3])))
            need(u not in bnext,'REGION_BOUNDARY_BRANCH','normalized cut complex')
            bnext[u]=(v,a);bverts.update((u,v))
        cycles=[];unseen=set(bnext)
        while unseen:
            seed=min(unseen);z=seed;cycle=[]
            while True:
                need(z in unseen,'REGION_BOUNDARY_CYCLE','normalized cut complex')
                unseen.remove(z);z,a=bnext[z];cycle.append(a)
                if z==seed:break
            cycles.append(cycle)
        chi=len({vd.find(x) for x in cn})-len({ed.find(x) for x in eo})+len(comp)
        b=len(cycles);g=(2-b-chi)//2
        need(2-2*g-b==chi and g>=0,'REGION_TOPOLOGY','normalized cut complex')
        image_boundary_simple=False
        if b==1:
            ids=[st.vertex_dsu.find((a[0],a[2])) for a in cycles[0]]
            image_boundary_simple=len(ids)==len(set(ids))
        out.append({'seed':list(min(comp)),'triangles':[list(x) for x in sorted(comp)],'chi':chi,'genus':g,
                    'boundary_components':b,'normalized_disk':chi==1 and b==1,
                    'embedded_closed_disk':chi==1 and b==1 and image_boundary_simple,
                    'boundary_cycles':[[list(a) for a in cyc] for cyc in cycles]})
    return out

def check_joining_disk(st:State, seed:list, orientations=None):
    """Check geometric disk and Koda's boundary-outward-positive convention.
    orientations is a non-mutating audit of an explicit proposed curve reversal.
    It is never silently applied to the stored crossings or exported core.
    """
    orientations=orientations or {}
    need(set(orientations)<=set(st.curves) and all(type(x) is int and x in (-1,1) for x in orientations.values()),'ORIENTATION_AUDIT','expected known curve -> +/-1')
    r=next((r for r in regions(st) if seed in r['triangles']),None)
    need(r is not None,'REGION_SEED',str(seed))
    result={'seed':seed,'embedded_closed_disk':r['embedded_closed_disk'],'arcs_all_nonempty_connected':False,
            'boundary_outward_positive':False,'koda_joining':False}
    if not r['embedded_closed_disk']:return result
    labels=[];against=[]
    for f,i,u,v in r['boundary_cycles'][0]:
        owners=st.used.get((f,edge(u,v)),[])
        if len(owners)!=1:return result
        cid=owners[0];labels.append(cid)
        d=next(d for d in st.curves[cid][1] if d.face==f and edge(d.u,d.v)==edge(u,v))
        # Region boundary traversed CCW (disk on LEFT); curve must point the OTHER way.
        against.append((1 if (d.u,d.v)==(u,v) else -1)*orientations.get(cid,1)==-1)
    required={cid for cid,(c,_) in st.curves.items() if c['role'] in ('alpha','beta')}
    runs=[]
    for x in labels:
        if not runs or x!=runs[-1]:runs.append(x)
    if len(runs)>1 and runs[0]==runs[-1]:runs.pop()
    result['arcs_all_nonempty_connected']=set(runs)==required and len(runs)==len(required)
    result['boundary_outward_positive']=all(against)
    result['koda_joining']=result['arcs_all_nonempty_connected'] and all(against)
    return result

def family_cut_report(st:State, role:str):
    ids={cid for cid,(c,_) in st.curves.items() if c['role']==role}
    rs=regions(st,{e for e,cs in st.used.items() if cs[0] in ids})
    g=st.surface_report['genus'];closed=st.surface_report['boundary_components']==0
    return {'role':role,'curves':len(ids),'components':len(rs),
            'cut_system':closed and len(ids)==g and len(rs)==1 and rs[0]['genus']==0 and rs[0]['boundary_components']==2*g,
            'topology':[{'genus':r['genus'],'boundary_components':r['boundary_components']} for r in rs]}

def move_vertex(st:State,pos:dict,face:str,vertex:str,target):
    """Certified one-vertex linear isotopy, fixed boundary, unchanged core.
    Incident signed areas are affine in time; positive endpoints certify all t.
    """
    certify_layout(st,pos)
    need(face in st.patches,'MOVE_FACE',face)
    need(vertex not in st.boundaries[face],'MOVE_BOUNDARY','boundary motion requires a separate marked transport')
    need(vertex in pos[face],'MOVE_VERTEX',vertex)
    need(len(target)==2,'MOVE_TARGET','two exact coordinates')
    out=copy.deepcopy(pos);out[face][vertex]=tuple(map(qnum,target));certify_layout(st,out)
    return out,{'kind':'interior_vertex_linear','face':face,'vertex':vertex,'old':[str(x) for x in pos[face][vertex]],
               'new':[str(x) for x in out[face][vertex]],'same_representation':representation_digest(st.doc)}

def rounded(q:Q,digits:int)->Q:
    n=10**digits;return Q((q.numerator*n*2+q.denominator)//(2*q.denominator),n)
def dec(q:Q,digits:int)->str:
    q=rounded(q,digits);sgn='-' if q<0 else '';n=abs(q.numerator)*(10**digits)//q.denominator
    return f'{sgn}{n//(10**digits)}.{n%(10**digits):0{digits}d}'

def write_svg(st:State,pos:dict,path,mesh=True):
    """Round-trip certify decimal centerlines before serializing SVG."""
    certify_layout(st,pos);digits=None;rpos=None
    for d in (6,12,24,48,96,192,384,768):
        candidate={f:{v:tuple(rounded(x,d) for x in p) for v,p in ps.items()} for f,ps in pos.items()}
        try:certify_layout(st,candidate);digits=d;rpos=candidate;break
        except DiagramError:pass
    need(digits is not None,'SVG_PRECISION_LIMIT','no certified decimal realization')
    n=len(st.patches);cols=min(n,2);rows=(n+cols-1)//cols;w=520*cols;h=500*rows
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img">',
         '<style>.mesh{stroke:#bcc5d3;stroke-width:0.004;fill:none}.curve{fill:none;stroke-width:0.014;stroke-linejoin:round}.alpha{stroke:#b63837}.beta{stroke:#275db2}.aux{stroke:#7255a3}.label{font:0.065px monospace;fill:#152435}.gate{fill:white;stroke:#566479;stroke-width:0.006}</style>',
         f'<metadata>representation_sha256={representation_digest(st.doc)};decimal_digits={digits};centerlines_certified=true</metadata>',
         '<rect width="100%" height="100%" fill="white"/>']
    for ix,(f,p) in enumerate(st.patches.items()):
        ox=260+520*(ix%cols);oy=250+500*(ix//cols)
        out.append(f'<text x="{ox-225}" y="{oy-222}" font-size="18" font-family="sans-serif">{html.escape(f)} · marked disk</text>')
        out.append(f'<g transform="translate({ox} {oy}) scale(195 -195)">')
        pts=rpos[f]
        def point(v):return ','.join(dec(x,digits) for x in pts[v])
        out.append(f'<polygon points="{" ".join(point(v) for v in st.boundaries[f])}" fill="#f6f8fc" stroke="#53647d" stroke-width="0.01"/>')
        for mark in st.doc.get('region_marks',[]):
            region=next(r for r in regions(st) if mark['seed'] in r['triangles'])
            for ff,i in region['triangles']:
                if ff==f:
                    poly=' '.join(point(v) for v in p['triangles'][i])
                    out.append(f'<polygon class="region" data-region="{html.escape(mark["id"])}" points="{poly}" fill="#f1c94c" fill-opacity="0.55"><title>{html.escape(mark["id"])}: {mark["claim"]}</title></polygon>')
        if mesh:
            for u,v in st.incidence[f]:out.append(f'<path class="mesh" d="M {point(u)} L {point(v)}"/>')
        for cid,(c,ds) in st.curves.items():
            for d in ds:
                if d.face==f:out.append(f'<path class="curve {c["role"]}" data-curve="{html.escape(cid)}" d="M {point(d.u)} L {point(d.v)}"><title>{html.escape(cid)}</title></path>')
        for v in st.boundaries[f]:
            x,y=pts[v];out.append(f'<circle class="gate" cx="{dec(x,digits)}" cy="{dec(y,digits)}" r="0.012"/>')
        labelv=set(st.boundaries[f])|{v for ff,v in st.crossing_at if ff==f}
        for v in sorted(labelv):
            x,y=pts[v];lab=st.crossing_at.get((f,v),{}).get('id',v)
            out.append(f'<g transform="translate({dec(x,digits)} {dec(y,digits)}) scale(1 -1)"><text class="label" x="0.025" y="-0.025">{html.escape(lab)}</text></g>')
        for (ff,v),cr in st.crossing_at.items():
            if ff==f:
                x,y=pts[v];out.append(f'<circle cx="{dec(x,digits)}" cy="{dec(y,digits)}" r="0.018" fill="#172333"/>')
        out.append('</g>')
    out.append('</svg>');Path(path).write_text('\n'.join(out),encoding='utf-8')
    return {'digits':digits,'certified_centerlines':True,'representation_sha256':representation_digest(st.doc)}

def read_document(path):
    path=Path(path)
    return tomllib.loads(path.read_text()) if path.suffix=='.toml' else json.loads(path.read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('document');ap.add_argument('--svg');ap.add_argument('--report');args=ap.parse_args()
    try:
        d=read_document(args.document);s=validate(d);p=realize(s)
        result={'surface':s.surface_report,'representation_sha256':representation_digest(d),'marking_sha256':marking_digest(d),
                'curves':{c:{'seam_word':seam_word(s,c),'endpoint_anchors':endpoint_anchors(s,c),'crossings':intersection_cycles(s)[c]} for c in s.curves},
                'regions':regions(s),'alpha_cut':family_cut_report(s,'alpha'),'beta_cut':family_cut_report(s,'beta'),
                'rational_positions':{f:{v:[str(x) for x in q] for v,q in ps.items()} for f,ps in p.items()}}
        if args.svg:result['svg']=write_svg(s,p,args.svg)
        text=json.dumps(result,indent=2,ensure_ascii=False)
        if args.report:Path(args.report).write_text(text,encoding='utf-8')
        else:print(text)
    except (DiagramError,KeyError,TypeError,ValueError) as e:
        raise SystemExit(str(e))
if __name__=='__main__':main()
