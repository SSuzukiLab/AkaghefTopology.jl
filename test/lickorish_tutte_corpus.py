"""Replay the previously qualified 99 LC inputs; no generated input fixtures.

cf:docs/heegaard/261109_tutte_construction.md 検証の境界 — corpus and construction diagnostics outside the drawing path
"""
import contextlib
import io
import json
from pathlib import Path
import sys

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts/heegaard'))
from lickorish_draw import read_ld
from lickorish_placement import place
from lickorish_render import write_svg
from lickorish_tutte import audit_mesh
from lickorish_placement_geometry import audit


def run(research, output):
    baseline=research/'.research-scratch/PJ260922/ld-draw-astra'
    records=json.loads((baseline/'261109_corpus-inputs.json').read_text())
    previous={r['path']:r for r in json.loads((baseline/'261109_corpus-placements.json').read_text())}
    output.mkdir(parents=True,exist_ok=True)
    result=[]
    for record in records:
        if record['input_status']!='pass': continue
        source=record['path']; data=read_ld(research/source)
        name=source.split('data/',1)[1].replace('/','__')[:-3]
        row=dict(path=source,genus=record['genus'],vertices=record['vertices'])
        log=io.StringIO()
        try:
            with contextlib.redirect_stdout(log): doc,svg,details=place(data)
            row['methods']=details['methods']
            row['construction_problems']={f:audit_mesh(m) for f,m in details.get('tutte',{}).items() if audit_mesh(m)}
            row['geometry_problems']=audit(details,data['words'])
            row['status']='diagnostic' if row['construction_problems'] or row['geometry_problems'] else 'pass'
            if details.get('tutte'):
                with contextlib.redirect_stdout(log): replay,_,_=place(data)
                row['deterministic']=doc==replay
                (output/(name+'.details.json')).write_text(json.dumps(details),encoding='utf-8')
            else:
                old=previous[source]
                row['unchanged']=doc==json.loads((baseline/old['artifact']).with_suffix('.json').read_text())
            row['artifact']=name+'.svg'
            (output/(name+'.json')).write_text(json.dumps(doc),encoding='utf-8')
            write_svg(doc,output/row['artifact'],title='LC · '+source.split('data/',1)[1][:-3])
            if row['status']!='pass': row['words']=data['words']
        except Exception as error:
            row.update(status='error',error=repr(error),words=data['words'])
        row['log']=log.getvalue();result.append(row)
        (output/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
        print(source,row['status'],row.get('construction_problems',''),flush=True)
    assert len(result)==99
    assert len({r.get('artifact') for r in result if 'artifact' in r})==sum('artifact' in r for r in result)
    assert all(r.get('unchanged',True) for r in result), 'a chord-only drawing changed'
    assert all(r.get('deterministic',True) for r in result), 'non-deterministic placement'
    return result


if __name__=='__main__':
    result=run(Path(sys.argv[1]),Path(sys.argv[2]))
    raise SystemExit(1 if any(r['status']!='pass' for r in result) else 0)
