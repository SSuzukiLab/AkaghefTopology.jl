"""Run: python src/verify_case.py.  Standard library only.
Tests actual image-derived data against the proposed answer, without using it
as an input to the reconstruction. A successful verification intentionally
reports target_equal=false.
"""
import copy, json, math, itertools
from pathlib import Path
import build_case as B
import marked_reference as R
import trace_input as T
ROOT=B.ROOT

def read(n):return json.loads((ROOT/n).read_text(encoding='utf-8'))
def sign(v):return 1 if v>0 else -1

def max_run(word):
    if len({sign(v) for v in word})==1:return len(word)
    best=0
    for i in range(len(word)):
        k=1
        while k<len(word) and sign(word[(i+k)%len(word)])==sign(word[i]):k+=1
        best=max(best,k)
    return best

def signed_runs(word):
    i=next(i for i in range(len(word)) if sign(word[i])!=sign(word[i-1]))
    rotated=word[i:]+word[:i];runs={1:[],-1:[]};prev=sign(rotated[0]);n=0
    for v in rotated:
        s=sign(v)
        if s!=prev:runs[prev].append(n);n=0;prev=s
        n+=1
    runs[prev].append(n)
    return {'positive':sorted(runs[1],reverse=True),'negative':sorted(runs[-1],reverse=True)}

def same_up_to_cyclic_relabel(a,b,reverse=False,global_visit_reverse=False):
    wa,ta=a[0][0],a[1];wb,tb=b[0][0],b[1]
    variants=[wa]
    if reverse:variants.append(wa[::-1])
    for raw in variants:
        for sgn in ([1,-1] if global_visit_reverse else [1]):
            w=[sgn*x for x in raw]
            for k in range(len(w)):
                q=w[k:]+w[:k];rename={};inv={};ok=True
                for x,y in zip(q,wb):
                    if sign(x)!=sign(y):ok=False;break
                    u,v=abs(x),abs(y)
                    if u in rename and rename[u]!=v:ok=False;break
                    if v in inv and inv[v]!=u:ok=False;break
                    if ta[u-1]!=tb[v-1]:ok=False;break
                    rename[u]=v;inv[v]=u
                if ok:return {'shift':k,'global_visit_sign':sgn,'relabel':rename}
    return None

def pixel_intersection_audit():
    _,_,_,_,witness,segs=B.input_data()
    actual=set();projected_ignored=0
    for a,b in itertools.combinations(segs,2):
        c1,s1,u,v=a;c2,s2,w,z=b
        p,q,r,s=map(B.point,(u,v,w,z))
        meet=R.segments_meet(p,q,r,s)
        if not meet:continue
        common=set((p,q))&set((r,s))
        if s1!=s2:
            if not common:projected_ignored+=1
            continue
        if not common:raise AssertionError(('unexpected same-sheet intersection',a,b))
        assert len(common)==1,(a,b)
        pt=next(iter(common))
        if c1==c2:continue # adjacency; separate topology validator checks simplicity
        names={x for x in (u,v,w,z) if isinstance(x,str) and T.POINTS[x]==pt}
        assert len(names)==1,(a,b)
        x=next(iter(names));assert x in witness,(a,b)
        actual.add(x)
    assert actual==set(witness)
    return {'same_sheet_crossing_ids':sorted(actual),
            'front_crossings':sorted(x for x in actual if witness[x]['sheet']=='F'),
            'back_crossings':sorted(x for x in actual if witness[x]['sheet']=='B'),
            'unplanned_same_sheet_intersections':0,
            'different_sheet_segment_intersections_ignored':projected_ignored,
            'scope':'exact check of the recorded PL tracing witness, not automatic image recognition'}

def main():
    passed=[]
    def check(name,condition):
        assert condition,name;passed.append(name)
    doc=read('data/figure6.heegaard.marked.json');full=read('data/figure6.fullmarked.json')
    st=R.validate(doc);sf=R.validate(full)
    check('closed oriented genus 3',st.surface_report['genus']==3 and st.surface_report['chi']==-4 and st.surface_report['boundary_components']==0)
    check('complete alpha meridian system',R.family_cut_report(st,'alpha')['cut_system'])
    check('complete beta meridian system',R.family_cut_report(st,'beta')['cut_system'])
    check('joining disk geometric and oriented',R.check_joining_disk(st,doc['region_marks'][0]['seed'])['koda_joining'])
    check('eleven Heegaard crossings',len(doc['crossings'])==11)
    check('fourteen extra front/back gates',len(full['crossings'])-len(doc['crossings'])==14)
    check('same surface marking with/without outlines',R.marking_digest(doc)==R.marking_digest(full))
    check('seven disk complementary regions',len(R.regions(st))==7 and all(r['normalized_disk'] for r in R.regions(st)))
    core=R.legacy_projection(st)
    check('legacy core independently derived',core==read('data/figure6.core.json'))
    k=B.koda(core)
    check('exact repeat of Koda conversion',k==read('data/figure6.koda.json'))
    check('five retained vertices and ten directed edges',len(k['vertices'])==5 and len(k['edges'])==10)
    check('five positive and five negative visits',sorted(k['o_data'][0][0])==list(range(-5,0))+list(range(1,6)))
    check('six complementary alpha/beta blocks',len(k['blocks'])==6)
    check('correct local vertex types',k['o_data'][1]==[-1,1,1,1,-1])
    # Printed Figure 6 E-data is a separate external regression, not reconstruction input.
    paper=[[[-1,-2,3,1,4,-5,-4,5,-3,2]],[-1,1,1,1,-1]]
    check('matches Koda printed Figure 6 E-data',k['o_data']==paper)
    target=read('data/target_odata.json')
    check('does not equal the user target',k['o_data']!=target)
    check('not cyclic-relabel equivalent',same_up_to_cyclic_relabel(k['o_data'],target) is None)
    check('not equivalent even after extra reversals',same_up_to_cyclic_relabel(k['o_data'],target,True,True) is None)
    check('renaming preserves maximum signed run',max_run(k['o_data'][0][0])==3 and max_run(target[0][0])==2)
    # Verify all local integer orientation witnesses, including the back-sheet sign.
    _,pairs,signs,_,witness,_=B.input_data()
    check('all signs derived from chart rays',[signs[x['id']] for x in doc['crossings']]==[x['sign'] for x in doc['crossings']])
    check('x2 is the rear crossing',witness['x2']['sheet']=='B' and signs['x2']==-1)
    px=pixel_intersection_audit();check('same-sheet trace has no unplanned crossings',px['unplanned_same_sheet_intersections']==0)
    # Deliberately corrupt the back sign: kernel must reject, not repair silently.
    bad=copy.deepcopy(doc)
    next(x for x in bad['crossings'] if x['id']=='x2')['sign']*= -1
    try:R.validate(bad)
    except R.DiagramError:check('wrong rear sign rejected',True)
    else:raise AssertionError('wrong rear sign accepted')
    # Deliberately omit a genuine intersection: no undeclared crossing accepted.
    bad=copy.deepcopy(doc);bad['crossings']=[x for x in bad['crossings'] if x['id']!='x4']
    try:R.validate(bad)
    except R.DiagramError:check('omitted x4 rejected',True)
    else:raise AssertionError('undeclared x4 accepted')
    coords=read('artifacts/certified_atlas.coordinates.json')
    pos={f:{v:tuple(R.qnum(z) for z in p) for v,p in ps.items()} for f,ps in coords.items()}
    check('exact rational atlas rechecked',R.certify_layout(st,pos))
    diff=[{'position_1_based':i+1,'computed':a,'target':b} for i,(a,b) in enumerate(zip(k['o_data'][0][0],target[0][0])) if a!=b]
    comparison={'computed':k['o_data'],'target':target,'literal_equal':False,'cyclic_relabel_equal':False,
         'even_reversal_and_global_visit_flip_equal':False,'vertex_type_lists_equal':True,
         'differences':diff,'computed_run_lengths':signed_runs(k['o_data'][0][0]),
         'target_run_lengths':signed_runs(target[0][0]),
         'scope':'equality of the oriented typed o-graph code, NOT a 3-manifold inequivalence claim'}
    B.dump(ROOT/'artifacts/comparison.json',comparison)
    B.dump(ROOT/'artifacts/pixel_trace_audit.json',px)
    B.dump(ROOT/'artifacts/tests.json',{'passed':len(passed),'failures':0,'checks':passed,
               'python_executed':True,'julia_executed':False,'existing_application_integrated':False})
    print(json.dumps({'passed':len(passed),'failures':0,'computed_o_data':k['o_data'],'target_equal':False,
       'cyclic_relabel_equal':False,'trace_audit':px},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
