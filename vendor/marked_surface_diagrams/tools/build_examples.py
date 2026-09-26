from pathlib import Path
import copy,json,sys
sys.path.insert(0,str(Path(__file__).parent))
from reference import validate,realize,regions,intersection_cycles,seam_word,check_joining_disk
ROOT=Path(__file__).resolve().parents[1]

def grid(n=4, surface='torus', flip=False):
    v=lambda x,y:f'v{x}_{y}'
    tris=[]
    for y in range(n):
        for x in range(n):
            a,b,c,d=v(x,y),v(x+1,y),v(x+1,y+1),v(x,y+1)
            tris += [[a,b,d],[b,c,d]] if flip else [[a,b,c],[a,c,d]]
    p={'face':'P','side_order':['south','east','north','west'],'sides':{
       'south':[v(x,0) for x in range(n+1)],'east':[v(n,y) for y in range(n+1)],
       'north':[v(x,n) for x in range(n,-1,-1)],'west':[v(0,y) for y in range(n,-1,-1)]},'triangles':tris}
    gluings=[];bd=[]
    if surface in ('torus','annulus'):gluings.append([['P','east'],['P','west']])
    if surface=='torus':gluings.append([['P','north'],['P','south']])
    elif surface=='annulus':bd=[{'side':['P','south'],'label':'inner'},{'side':['P','north'],'label':'outer'}]
    else:bd=[{'side':['P',s],'label':'outer'} for s in p['side_order']]
    return {'version':1,'surface':{'marking_id':surface+'-square-v1','genus':int(surface=='torus'),'gluings':gluings,'boundary':bd},
            'patches':[p],'curves':[],'crossings':[]},v

def word(v,blocks):return ' | '.join('P{'+' '.join(v(x,y) for x,y in b)+'}' for b in blocks)
def curve(d,cid,role,closed,w):d['curves'].append({'id':cid,'role':role,'closed':closed,'word':w})
def line(a,b):
    dx=b[0]-a[0];dy=b[1]-a[1];n=max(abs(dx),abs(dy))
    assert n and dx%n==0 and dy%n==0
    return [(a[0]+i*dx//n,a[1]+i*dy//n) for i in range(n+1)]
def polyline(pts):
    out=[]
    for a,b in zip(pts,pts[1:]):out.extend(line(a,b)[:-1])
    return out+[pts[-1]]

def write(name,d):
    validate(d)
    (ROOT/'examples'/f'{name}.json').write_text(json.dumps(d,indent=2),encoding='utf8')
    # A deterministic TOML serialization using only this schema's fields.
    q=lambda x:json.dumps(x,ensure_ascii=False)
    lines=['version = 1','', '[surface]',f'marking_id = {q(d["surface"]["marking_id"])}',f'genus = {d["surface"]["genus"]}',
           'gluings = '+q(d['surface']['gluings'])]
    if not d['surface']['boundary']:lines.append('boundary = []')
    for b in d['surface']['boundary']:lines+=['[[surface.boundary]]','side = '+q(b['side']),'label = '+q(b['label'])]
    for p in d['patches']:
        lines+=['','[[patches]]','face = '+q(p['face']),'side_order = '+q(p['side_order']),
                'triangles = [']+['  '+q(t)+',' for t in p['triangles']]+[']','[patches.sides]']
        lines.extend(q(k)+' = '+q(vs) for k,vs in p['sides'].items())
    for c in d['curves']:
        lines+=['','[[curves]]','id = '+q(c['id']),'role = '+q(c['role']),'closed = '+str(c['closed']).lower(),'word = '+q(c['word'])]
    for x in d['crossings']:
        lines+=['','[[crossings]]','id = '+q(x['id']),'at = '+q(x['at']),'curves = '+q(x['curves']),'sign = '+str(x['sign'])]
    for r in d.get('region_marks',[]):
        lines+=['','[[region_marks]]','id = '+q(r['id']),'seed = '+q(r['seed']),'claim = '+q(r['claim'])]
    # Empty arrays have to be root values, not entries in the last table.
    if not d['curves']:lines.insert(1,'curves = []')
    if not d['crossings']:lines.insert(1,'crossings = []')
    (ROOT/'examples'/f'{name}.toml').write_text('\n'.join(lines)+'\n',encoding='utf8')

for k in (-1,0,1):
    d,v=grid(4,'annulus',flip=k==-1)
    if k==0:blocks=[line((1,0),(1,4))]
    elif k==1:blocks=[line((1,0),(4,3)),line((0,3),(1,4))]
    else:blocks=[line((1,0),(0,1)),line((4,1),(1,4))]
    curve(d,'gamma','aux',False,word(v,blocks));write(f'annulus_winding_{k:+d}'.replace('+','p').replace('-','m'),d)

d,v=grid();curve(d,'alpha1','alpha',True,word(v,[line((0,2),(4,2))]));curve(d,'beta1','beta',True,word(v,[line((2,0),(2,4))]));d['crossings']=[{'id':'x','at':['P',v(2,2)],'curves':['alpha1','beta1'],'sign':1}];write('torus_mer_lon',d)

d,v=grid();curve(d,'alpha1','alpha',True,word(v,[line((0,1),(4,1))]));curve(d,'beta1','beta',True,word(v,[line((0,3),(4,3))]));write('torus_parallel_zero_crossings',d)

d,v=grid(8);curve(d,'alpha1','alpha',True,word(v,[line((0,4),(8,4))]));curve(d,'beta1','beta',True,word(v,[polyline([(2,0),(2,5),(4,5),(4,3),(6,3),(6,7),(2,7),(2,8)])]));d['crossings']=[{'id':f'x{i}','at':['P',v(x,4)],'curves':['alpha1','beta1'],'sign':sign} for i,x,sign in [(1,2,1),(2,4,-1),(3,6,1)]];d['region_marks']=[{'id':'D0','seed':['P',68],'claim':'embedded_closed_disk'}];write('torus_three_crossings',d)

d,v=grid(6,'disk');curve(d,'outer','alpha',True,word(v,[polyline([(1,1),(5,1),(5,5),(1,5),(1,1)])]));curve(d,'inner','beta',True,word(v,[polyline([(2,2),(4,2),(4,4),(2,4),(2,2)])]));write('nested_zero_crossings',d)

for g in (1,2,3,4):
    n=4*g;sides=[f's{i}' for i in range(n)];p={'face':'P','side_order':sides,'sides':{s:[f'c{i}',f'c{(i+1)%n}'] for i,s in enumerate(sides)},'triangles':[[f'c{i}',f'c{(i+1)%n}','center'] for i in range(n)]}
    gl=[]
    for k in range(g):gl.extend([[['P',f's{4*k}'],['P',f's{4*k+2}']],[['P',f's{4*k+1}'],['P',f's{4*k+3}']]])
    d={'version':1,'surface':{'marking_id':f'polygon-genus-{g}-v1','genus':g,'gluings':gl,'boundary':[]},'patches':[p],'curves':[],'crossings':[]};write(f'genus_{g}_empty',d)

# Two triangles, oppositely glued along every edge, give an oriented sphere.
p=lambda f:{'face':f,'side_order':['a','b','c'],'sides':{'a':['u','v'],'b':['v','w'],'c':['w','u']},'triangles':[['u','v','w']]}
d={'version':1,'surface':{'marking_id':'sphere-two-triangles-v1','genus':0,'gluings':[[['P',s],['Q',t]] for s,t in [('a','a'),('b','c'),('c','b')]],'boundary':[]},'patches':[p('P'),p('Q')],'curves':[],'crossings':[]};write('sphere_empty',d)
print('Generated',len(list((ROOT/'examples').glob('*.json'))),'fixtures')
