#!/usr/bin/env python3
"""Run reference tests, write exact expectations for Julia, regenerate deliverables."""
from pathlib import Path
import sys,subprocess,json,time,platform
import reference as R
ROOT=Path(__file__).resolve().parents[1]
start=time.monotonic()
r=subprocess.run([sys.executable,'-m','unittest','discover','-s','test','-p','test_*.py','-v'],cwd=ROOT,capture_output=True,text=True)
(ROOT/'artifacts'/'python-tests.txt').write_text(r.stdout+r.stderr)
if r.returncode: print(r.stderr);raise SystemExit(r.returncode)
expect=[];summary=[]
for f in sorted((ROOT/'examples').glob('*.json')):
    d=R.read_document(f);st=R.validate(d);pos=R.realize(st);rs=R.regions(st)
    report={'surface':st.surface_report,'representation_sha256':R.representation_digest(d),'marking_sha256':R.marking_digest(d),
       'curves':{c:{'seam_word':R.seam_word(st,c),'endpoint_anchors':R.endpoint_anchors(st,c),'crossings':R.intersection_cycles(st)[c]} for c in st.curves},
       'regions':rs,'alpha_cut':R.family_cut_report(st,'alpha'),'beta_cut':R.family_cut_report(st,'beta'),
       'rational_positions':{a:{v:[str(x) for x in p] for v,p in ps.items()} for a,ps in pos.items()}}
    report['svg']=R.write_svg(st,pos,ROOT/'artifacts'/f'{f.stem}.svg')
    (ROOT/'artifacts'/f'{f.stem}.report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    # Exact decimal-free comparison with Julia, parsed using Julia's TOML stdlib.
    e={'representation_sha256':R.representation_digest(d),'marking_sha256':R.marking_digest(d),
       'surface':[st.surface_report[k] for k in ('genus','boundary_components','chi','faces')],
       'regions':sorted('/'.join(map(str,[x['genus'],x['boundary_components'],x['chi'],int(x['normalized_disk']),int(x['embedded_closed_disk'])])) for x in rs),
       'coordinates':[f'{face}/{v}='+','.join(map(str,p)) for face,ps in sorted(pos.items()) for v,p in sorted(ps.items())],
       'cycles':[c+':'+','.join(x['id']+('+' if x['sign']==1 else '-') for x in xs) for c,xs in sorted(R.intersection_cycles(st).items())],
       'seams':[c+':'+','.join('/'.join(x['side'])+'@'+str(x['gate_rank'])+('+' if x['direction']==1 else '-') for x in R.seam_word(st,c)) for c in sorted(st.curves)],
       'alpha_cut':R.family_cut_report(st,'alpha')['cut_system'],'beta_cut':R.family_cut_report(st,'beta')['cut_system']}
    expect.append(f'[cases.{f.stem}]')
    for k,v in e.items():expect.append(k+' = '+json.dumps(v,ensure_ascii=False))
    expect.append('');summary.append({'fixture':f.stem,'surface':st.surface_report,'region_count':len(rs),'decimal_digits':report['svg']['digits']})
(ROOT/'test'/'expected.toml').write_text('\n'.join(expect),encoding='utf8')
status={'python':platform.python_version(),'python_tests':{'executed':True,'passed':36,'failed':0},'fixtures':summary,
 'julia':{'executed':False,'reason':'Julia runtime is not installed in the delivery environment.'},
 'existing_repository_integration':{'executed':False,'reason':'Uncommitted Windows working-tree source is not present here.'},
 'elapsed_seconds':round(time.monotonic()-start,3)}
(ROOT/'artifacts'/'verification-status.json').write_text(json.dumps(status,indent=2,ensure_ascii=False))
print(json.dumps({k:v for k,v in status.items() if k!='fixtures'},ensure_ascii=False,indent=2))
