"""Deterministic TOML export for the Julia stdlib-only validator."""
import json
from pathlib import Path
import build_case as B
q=lambda x:json.dumps(x,ensure_ascii=False)
for name in ['figure6.heegaard.marked','figure6.fullmarked']:
    d=json.loads((B.ROOT/'data'/f'{name}.json').read_text())
    lines=['version = 1','', '[surface]',f'marking_id = {q(d["surface"]["marking_id"])}',f'genus = {d["surface"]["genus"]}',
           'gluings = '+q(d['surface']['gluings']),'boundary = []']
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
    (B.ROOT/'data'/f'{name}.toml').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import tomllib
    parsed=tomllib.loads('\n'.join(lines))
    assert B.R.representation_digest(parsed)==B.R.representation_digest(d)
    assert B.R.validate(parsed).surface_report['genus']==3
print('TOML roundtrip and Python validation passed for both documents.')
