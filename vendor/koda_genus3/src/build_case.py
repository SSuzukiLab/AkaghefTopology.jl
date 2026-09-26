"""Reconstruct the supplied drawing, derive its o-data, and validate a marked mesh.

Input: trace_input.py only.  The user's proposed o-data is NOT read here.
Python 3.10+, standard library.  All sign decisions are exact integer arithmetic.
"""
from __future__ import annotations
import argparse, copy, hashlib, json
from collections import defaultdict, Counter
from functools import cmp_to_key
from pathlib import Path
import trace_input as T
import marked_reference as R

ROOT = Path(__file__).resolve().parents[1]

def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')

def point(v): return T.POINTS[v] if isinstance(v,str) else tuple(v)
def det(a,b): return a[0]*b[1]-a[1]*b[0]
def sub(a,b): return a[0]-b[0],a[1]-b[1]
def orient(a,b,c): return det(sub(b,a),sub(c,a))

def angular_cmp(a,b):
    # a,b = (dart, oriented chart vector); CCW order starting at +x.
    def half(v):return 0 if (v[1]>0 or (v[1]==0 and v[0]>0)) else 1
    av,bv=a[1],b[1]
    if half(av)!=half(bv):return -1 if half(av)<half(bv) else 1
    d=det(av,bv)
    if d==0: raise ValueError('nontransverse / coincident rays')
    return -1 if d>0 else 1

def input_data():
    """Derive cyclic orders and signed rotation from the explicit image traces."""
    cycles={}; event_sheets=defaultdict(set); rays=defaultdict(dict); transitions={}
    edge_sheets={}; trace_segments=[]
    for cid,pieces in T.PIECES.items():
        seq=[]; flat=[]; fsheets=[]
        for k,(sheet,vs) in enumerate(pieces):
            assert vs[-1]==pieces[(k+1)%len(pieces)][1][0],(cid,k)
            for a,b in zip(vs,vs[1:]):
                assert point(a)!=point(b),(cid,a,b)
                trace_segments.append((cid,sheet,a,b))
            flat+=vs[:-1];fsheets+=[sheet]*(len(vs)-1)
        events=[]
        for i,v in enumerate(flat):
            if not isinstance(v,str):continue
            seq.append(v);event_sheets[v].add(fsheets[i]);events.append(i)
            before=fsheets[(i-1)%len(flat)];after=fsheets[i]
            if before!=after:transitions[v]=(before,after,cid)
            else:
                # Front chart orientation dx wedge d(-y), back its negative.
                flip=-1 if after=='F' else 1
                for direction,j in ((-1,i-1),(1,(i+1)%len(flat))):
                    dx,dy=sub(point(flat[j]),point(v))
                    rays[v][(cid,v,direction)]=(dx,flip*dy)
        for i,j in zip(events,events[1:]+[events[0]+len(flat)]):
            sheetseq={fsheets[k%len(flat)] for k in range(i,j)}
            assert len(sheetseq)==1,(cid,flat[i],sheetseq)
            edge_sheets[(cid,flat[i],1)]=next(iter(sheetseq))
            edge_sheets[(cid,flat[j%len(flat)],-1)]=next(iter(sheetseq))
        assert len(seq)==len(set(seq)),cid
        cycles[cid]=seq
    pairs={}; signs={}; source_signs={}
    for v,rr in rays.items():
        assert len(rr)==4,(v,rr)
        alpha=next(c for c,_,d in rr if c.startswith('alpha'))
        beta=next(c for c,_,d in rr if c.startswith('beta'))
        ordered=[d for d,r in sorted(rr.items(),key=cmp_to_key(angular_cmp))]
        i=ordered.index((alpha,v,1)); ordered=ordered[i:]+ordered[:i]
        assert ordered[2]==(alpha,v,-1),(v,ordered)
        assert ordered[1][0]==beta and ordered[3][0]==beta,(v,ordered)
        pairs[v]=(alpha,beta);signs[v]=ordered[1][2]
        source_signs[v]={'sheet':next(iter(event_sheets[v])), 'pixel':T.POINTS[v],
                         'oriented_chart_rays':[{'dart':list(d),'ray':r} for d,r in sorted(rr.items(),key=cmp_to_key(angular_cmp))]}
    for oid,vs in T.OUTLINES.items():
        cycles[oid]=list(vs)
        for v in vs:
            before,after,cid=transitions[v]
            assert (before,after) in [('F','B'),('B','F')]
            pairs[v]=(cid,oid);signs[v]=1 if after=='B' else -1
            # Silhouette direction has front on left and back on right.
            edge_sheets[(oid,v,1)]='F';edge_sheets[(oid,v,-1)]='B'
    assert set(pairs)==set(T.POINTS)
    return cycles,pairs,signs,edge_sheets,source_signs,trace_segments


def ribbon(cycles,pairs,signs):
    rotations={}
    for v,(a,b) in pairs.items():
        e=signs[v]
        rotations[v]=[(a,v,1),(b,v,e),(a,v,-1),(b,v,-e)]
    sigma={};twin={}
    for rot in rotations.values():
        sigma.update({d:rot[(i+1)%4] for i,d in enumerate(rot)})
    for c,vs in cycles.items():
        for x,y in zip(vs,vs[1:]+vs[:1]):
            twin[c,x,1]=(c,y,-1);twin[c,y,-1]=(c,x,1)
    assert set(sigma)==set(twin)
    sigmainv={v:k for k,v in sigma.items()}
    nxt={d:sigmainv[twin[d]] for d in twin}
    seen=set();faces=[]
    for d in twin:
        if d in seen:continue
        f=[];z=d
        while z not in seen:
            f.append(z);seen.add(z);z=nxt[z]
        assert z==d
        faces.append(f)
    return rotations,twin,faces


def make_document(cycles,pairs,signs,rotations,twin,faces,edge_sheets):
    """Thicken each vertex to a disk and cap each face with a disk.
    This is a complete finite cellulation, not just a signed intersection code.
    Auxiliary silhouette curves retain the front/back marking of the picture.
    """
    patches=[];gluings=[];corner_side={};port={}
    for v,rot in rotations.items():
        sides={};order=[];bd=[]
        for i,d in enumerate(rot):
            sides[f'a{i}']=[f'b{2*i}',f'p{i}',f'b{2*i+1}']
            sides[f'g{i}']=[f'b{2*i+1}',f'b{(2*i+2)%8}']
            order += [f'a{i}',f'g{i}']
            bd += [f'b{2*i}',f'p{i}',f'b{2*i+1}']
            port[d]=('V_'+v,f'a{i}',f'p{i}')
            corner_side[d]=('V_'+v,f'g{i}')
        patches.append({'face':'V_'+v,'side_order':order,'sides':sides,
                        'triangles':[['z',bd[i],bd[(i+1)%len(bd)]] for i in range(len(bd))]})
    done=set()
    for d,e in twin.items():
        if d in done:continue
        done.update((d,e));gluings.append([list(port[d][:2]),list(port[e][:2])])
    face_sheets=[]
    for j,face in enumerate(faces):
        sheets={edge_sheets[d] for d in face}
        assert len(sheets)==1,('inconsistent front/back face',j,sheets,face)
        face_sheets.append(next(iter(sheets)))
        n=len(face);sid=[f'e{k}' for k in range(n)]
        patches.append({'face':f'F_{j}','side_order':sid,
           'sides':{f'e{k}':[f't{k}',f't{(k+1)%n}'] for k in range(n)},
           'triangles':[['z',f't{k}',f't{(k+1)%n}'] for k in range(n)]})
        for k,d in enumerate(face):
            gluings.append([[f'F_{j}',f'e{k}'],list(corner_side[d])])
    crvs=[]
    for c,vs in cycles.items():
        blocks=[]
        for v in vs:
            pin=port[c,v,-1][2];pout=port[c,v,1][2]
            blocks.append(f'V_{v}'+'{'+pin+' z '+pout+'}')
        crvs.append({'id':c,'role':'alpha' if c.startswith('alpha') else 'beta' if c.startswith('beta') else 'aux',
                     'closed':True,'word':' | '.join(blocks)})
    dset={(c,b,-1) for c,a,b in T.D_RIGHT}
    di=next(i for i,f in enumerate(faces) if set(f)==dset)
    doc={'version':1,'surface':{'marking_id':'supplied-koda-fig6-front-back-v1','genus':3,'gluings':gluings,'boundary':[]},
         'patches':patches,'curves':crvs,
         'crossings':[{'id':v,'at':['V_'+v,'z'],'curves':list(pairs[v]),'sign':signs[v]} for v in pairs],
         'region_marks':[{'id':'D','seed':[f'F_{di}',0],'claim':'koda_joining'}],
         'metadata':{'source':'data/source.png; hand-traced, front/back distinguished',
                     'source_sha256':hashlib.sha256((ROOT/'data/source.png').read_bytes()).hexdigest(),
                     'provenance':'manual image transcription, NOT reconstruction from expected o-data',
                     'outline_roles':{'O':'outer silhouette','L':'left aperture','M':'middle aperture','R':'right aperture'},
                     'face_sheets':{f'F_{i}':s for i,s in enumerate(face_sheets)},
                     'orientation':'front dx wedge d(-y), back opposite; all alpha/beta have D on RIGHT'}}
    return doc


def koda(core):
    """Delete D arcs, follow all complementary directed edges through corners.
    For each retained vertex: alpha visit = negative, beta visit = positive.
    Local normal-o-graph type l=-1/r=+1 equals minus eps_H in this convention.
    """
    cyc={q['id']:q['cycle'] for fam in ('alpha','beta') for q in core[fam]}
    xs={x['id']:x for x in core['intersections']}
    deleted={(c,a,b) for c,a,b in T.D_RIGHT}
    corners={a for c,a,b in T.D_RIGHT}
    edges={}; succ=defaultdict(list)
    for c,vs in cyc.items():
        for a,b in zip(vs,vs[1:]+vs[:1]):
            if (c,a,b) in deleted:continue
            edges[(c,a)]=(c,a,b);succ[a].append((c,a))
    assert all(len(succ[c])==1 for c in corners)
    start=('alpha1','c1');z=start;walk=[];word=[];visit_edges=[];seen=set();blocks=[]
    while z not in seen:
        seen.add(z);c,a,b=edges[z];walk.append((c,a,b))
        if a not in corners:
            n=int(a[1:]);word.append(-n if c.startswith('alpha') else n);visit_edges.append(len(walk)-1)
        z=succ[b][0] if b in corners else (c,b)
    assert z==start and len(seen)==len(edges),(z,len(seen),len(edges))
    for c,a,b in walk:
        if not blocks or blocks[-1]['curve']!=c:blocks.append({'curve':c,'path':[a,b]})
        else:blocks[-1]['path'].append(b)
    # start chosen before first retained vertex; visit order already starts at x1.
    types=[-xs[f'x{i}']['sign'] for i in range(1,6)]
    oedges=[]
    for j,starti in enumerate(visit_edges):
        endi=visit_edges[(j+1)%len(visit_edges)]
        inds=[];k=starti
        while True:
            inds.append(k);k=(k+1)%len(walk)
            if k==endi:break
        oedges.append({'id':f'e{j+1}','from_visit':word[j],'to_visit':word[(j+1)%len(word)],
                        'heegaard_edges':[list(walk[k]) for k in inds]})
    return {'convention':{'visit_sign':{'alpha':-1,'beta':1},'vertex_type':{'l':-1,'r':1},
                           'vertex_type_from_heegaard':'sigma_o = -epsilon(alpha,beta)'},
            'removed_corner_ids':sorted(corners),'removed_D_edges':[list(x) for x in T.D_RIGHT],
            'blocks':blocks,'o_data':[[word],types],
            'vertices':[{'id':i,'heegaard_crossing':f'x{i}','type_sign':types[i-1],
                         'type':'l' if types[i-1]<0 else 'r'} for i in range(1,6)],
            'edges':oedges}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--render',action='store_true');args=ap.parse_args()
    cycles,pairs,signs,edge_sheets,source_signs,trace_segments=input_data()
    ro,twin,faces=ribbon(cycles,pairs,signs)
    doc=make_document(cycles,pairs,signs,ro,twin,faces,edge_sheets)
    full=R.validate(doc)
    mathdoc=copy.deepcopy(doc)
    mathdoc['curves']=[c for c in doc['curves'] if c['role']!='aux']
    mathdoc['crossings']=[x for x in doc['crossings'] if x['id'].startswith(('c','x'))]
    st=R.validate(mathdoc)
    core=R.legacy_projection(st)
    K=koda(core)
    rs=R.regions(st)
    outlineids=set(T.OUTLINES)
    sheets=R.regions(full,{e for e,cs in full.used.items() if cs[0] in outlineids})
    report={'surface':st.surface_report,'complete_alpha':R.family_cut_report(st,'alpha'),
            'complete_beta':R.family_cut_report(st,'beta'),
            'joining_disk':R.check_joining_disk(st,mathdoc['region_marks'][0]['seed']),
            'crossings':len(core['intersections']),'retained_o_vertices':len(K['vertices']),
            'ribbon_counts':{'V':25,'E':len(twin)//2,'F':len(faces),'chi':25-len(twin)//2+len(faces)},
            'complement_regions':[{'id':i,**{k:r[k] for k in ('seed','chi','genus','boundary_components','normalized_disk','embedded_closed_disk')}} for i,r in enumerate(rs)],
            'outline_cut_components':[{k:r[k] for k in ('chi','genus','boundary_components')} for r in sheets],
            'source_sign_witnesses':source_signs}
    assert report['complete_alpha']['cut_system'] and report['complete_beta']['cut_system']
    assert report['joining_disk']['koda_joining']
    assert len(rs)==7 and all(r['normalized_disk'] for r in rs)
    assert len(sheets)==2 and all(r['genus']==0 and r['boundary_components']==4 for r in sheets)
    dump(ROOT/'data/figure6.fullmarked.json',doc)
    dump(ROOT/'data/figure6.heegaard.marked.json',mathdoc)
    dump(ROOT/'data/figure6.core.json',core)
    dump(ROOT/'data/figure6.koda.json',K)
    dump(ROOT/'data/figure6.odata.json',K['o_data'])
    dump(ROOT/'artifacts/validation.json',report)
    dump(ROOT/'data/image_transcription.json',{'points':T.POINTS,'pieces':T.PIECES,'outlines':T.OUTLINES,
         'D_right':T.D_RIGHT,'augmented_cycles':cycles,'signs':signs,
         'ribbon_faces':[[list(d) for d in f] for f in faces]})
    if args.render:
        pos=R.realize(st)
        R.write_svg(st,pos,ROOT/'artifacts/certified_atlas.svg',mesh=True)
        dump(ROOT/'artifacts/certified_atlas.coordinates.json',{f:{v:[str(z) for z in p] for v,p in ps.items()} for f,ps in pos.items()})
        report['exact_layout_certified']=True;dump(ROOT/'artifacts/validation.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sign_witnesses','complement_regions')},ensure_ascii=False,indent=2))
    print('Computed o-data:',K['o_data'])

if __name__=='__main__':main()
