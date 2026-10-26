"""Realise the stored Fig.12.3 reference, with l_i as the explicit display cuts.

This is the standard marked reference, not a reconstructed Koda marking.
Connectivity is read from the original disk words and the native cut quotient.
"""
import argparse,json,sys,hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
from draw_marked_cut import render


def build(source, prepared, kernel):
    sys.path.insert(0,str(kernel));import marked_reference as m
    doc=m.read_document(source);st=m.validate(doc)
    g=prepared['genus'];cuts=set(prepared['cut_curves'])
    classes=prepared['source_cut_corner_classes']
    def vertex(f,i,v):return classes[f'{f},{i},{st.triangles[f,i].index(v)}']
    curves=[];positions={};boundaries=[];port_labels=[]
    for c in doc['curves']:
        if c['id'] in cuts:continue
        segments=[]
        for dart in st.curves[c['id']][1]:
            f,u,v=dart.face,dart.u,dart.v;i=st.incidence[f][m.edge(u,v)][0][0]
            segments.append(dict(source=[f,u,v],vertices=[vertex(f,i,u),vertex(f,i,v)]))
        curves.append(dict(id=c['id'],role='aux',color='#255dba' if c['id'].startswith('m') else '#a24743',segments=segments))
    crossings=defaultdict(set)
    for x in doc['crossings']:
        f,v=x['at']
        for (face,i),triangle in st.triangles.items():
            if face==f and v in triangle:crossings[x['id']].add(vertex(face,i,v))
    nextid=max(classes.values())+1
    for b in prepared['boundaries']:
        col=b['column'];i=col+1;row=b['row'];center=np.array([390*(col-(g-1)/2),250*(1-2*row)])
        cap=nextid;nextid+=1;positions[str(cap)]=(center/600).tolist();ports=[]
        for label,vertices in crossings.items():
            hits=set(b['vertices'])&vertices
            if not hits:continue
            assert len(hits)==1;v=hits.pop()
            angle=-np.pi/2 if label.startswith('lm_') else 0 if label==f'ln_{i}_{i}' else np.pi
            if row==1:angle=-angle
            positions[str(v)]=((center+55*np.array([np.cos(angle),np.sin(angle)]))/600).tolist()
            ports.append(dict(vertex=v,source_vertex=v,crossing=label,angle=float(angle)))
            parts=label.split('_');human=f'l{parts[1]}∩m{parts[1]}' if parts[0]=='lm' else f'l{parts[1]}∩n{parts[2]}'
            port_labels.append(dict(vertex=v,label=human))
        boundaries.append(dict(**b,cap=cap,radius=55/600,ports=ports,label=f'l{i}'+('⁺' if row==0 else '⁻')))
    # Follow each cut-open source curve, never connect unpaired endpoints.
    for curve in curves:
        nxt={u:v for s in curve['segments'] for u,v in zip(s['vertices'],s['vertices'][1:])}
        prev={v:u for u,v in nxt.items()};starts=set(nxt)-set(prev);paths=[]
        for start in sorted(starts):
            path=[start]
            while path[-1] in nxt:path.append(nxt[path[-1]])
            paths.append(path)
            a,z=np.array(positions[str(path[0])]),np.array(positions[str(path[-1])])
            # Reference m arcs join matched rows; n arcs join adjacent columns.
            assert abs(a[0]-z[0])<1e-10 if curve['id'].startswith('m') else abs(a[1]-z[1])<1e-10
            for t,v in zip(np.linspace(0,1,len(path)),path):positions[str(v)]=((1-t)*a+t*z).tolist()
        assert sum(len(p)-1 for p in paths)==len(nxt)
        middle=np.mean([positions[str(v)] for v in paths[0]],axis=0)
        curve['label_at']=(middle+np.array([.035,0]) if curve['id'].startswith('m') else middle+np.array([0,.035])).tolist()
    for pair in prepared['boundary_vertex_pairs']:
        if all(str(v) in positions for v in pair):
            a=next(b for b in boundaries if pair[0] in b['vertices']);b=next(b for b in boundaries if pair[1] in b['vertices'])
            p=np.array(positions[str(pair[0])])-positions[str(a['cap'])];q=np.array(positions[str(pair[1])])-positions[str(b['cap'])]
            assert np.linalg.norm(p-q*np.array([1,-1]))<1e-10,'Pairing is not reflected across the row gap'
    return dict(schema='marked-cut-view/1',title=f'Lickorish · g = {g} · l: paired boundaries · m: blue · n: red',genus=g,
        source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cuts=prepared['cut_curves'],
        source_figure=doc['metadata']['source_figure'],source_status=doc['metadata']['source_status'],
        grid=prepared['grid'],positions=positions,curves=curves,boundaries=boundaries,crossings=[],disk=[],
        port_labels=port_labels,boundary_vertex_pairs=prepared['boundary_vertex_pairs'],
        source_cut_corner_classes=classes,reference_only=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--genus',type=int,default=3)
    p.add_argument('--kernel',type=Path,required=True);p.add_argument('--register',action='store_true');a=p.parse_args();name=f'LickorishG{a.genus}'
    d=build(a.directory/(name+'.marked.json'),json.loads((a.directory/(name+'.display.json')).read_text()),a.kernel)
    (a.directory/(name+'.view.json')).write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8')
    (a.directory/(name+'.svg')).write_text(render(d),encoding='utf-8')
    if a.register:
        assert (a.directory/(name+'Common.json')).exists(),'Generate the common plot output before registering it'
        path=a.directory/'manifest.json';items=json.loads(path.read_text(encoding='utf-8'))
        entry=dict(id=name+'Common',title=f'Lickorish基準曲線 · 種数{a.genus}',genus=a.genus)
        items=[x for x in items if x['id']!=entry['id']]+[entry]
        path.write_text(json.dumps(items,ensure_ascii=False),encoding='utf-8')
    print(name,'paired boundary positions and source-edge paths checked',flush=True)
