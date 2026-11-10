"""Replay existing LC data and rotate K words without creating new inputs.

cf:docs/heegaard/261109_cycle_knots.md 確認 — coordinate preservation and cyclic-word invariance
"""
import contextlib,copy,io,json,sys
from pathlib import Path
REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'scripts/heegaard'))
from lickorish_draw import read_ld
from lickorish_placement import place
from lickorish_render import write_svg
from lickorish_tutte import audit_mesh
from lickorish_placement_geometry import audit

def coordinates(d):
    return {k:d[k] for k in ('W','H','holes','curves')}

def run(research,output):
    prior=research/'.research-scratch/PJ260922/ld-draw-astra'
    records=json.loads((prior/'261109_tutte-corpus/results.json').read_text())
    output.mkdir(parents=True,exist_ok=True);result=[]
    for old in records:
        data=read_ld(research/old['path']);log=io.StringIO()
        row={k:old[k] for k in ('path','genus','vertices')}
        row['previous_status']=old['status'];row['previous_tutte']='tutte' in old['methods'].values()
        try:
            with contextlib.redirect_stdout(log):doc,_,details=place(data)
            problems=audit(details,data['words'])
            for f,m in details.get('tutte',{}).items():problems.extend(f+': '+p for p in audit_mesh(m))
            # The same cyclic words, with only their textual starts moved.
            rotated=copy.deepcopy(data)
            g=data['genus'];names=[f'l{i}' for i in range(1,g+1)]+[f'm{i}' for i in range(1,g+1)]+[f'n{i}' for i in range(1,g)]
            for name in names:
                w=rotated['words'][name];rotated['words'][name]=w[-1:]+w[:-1]
            with contextlib.redirect_stdout(log):again,_,_=place(rotated)
            row['cyclic_identical']=doc==again
            expected=json.loads((prior/'261109_tutte-corpus'/old['artifact']).with_suffix('.json').read_text())
            row['coordinates_unchanged']=coordinates(doc)==coordinates(expected)
            row['problems']=problems;row['methods']=details['methods']
            row['status']='pass' if not problems and row['cyclic_identical'] else 'fail'
            row['artifact']=old['artifact'];row['declared_crossings']=len(doc.get('crossing_knots',[]))
            (output/old['artifact']).with_suffix('.json').write_text(json.dumps(doc))
            write_svg(doc,output/old['artifact'],title='LC · '+old['path'].split('data/',1)[1][:-3])
            if row['previous_tutte']:
                (output/old['artifact']).with_suffix('.details.json').write_text(json.dumps(details))
        except Exception as error:row.update(status='error',error=repr(error))
        row['log']=log.getvalue();result.append(row)
        (output/'results.json').write_text(json.dumps(result,indent=2))
        print(row['path'],row['status'],'unchanged',row.get('coordinates_unchanged'),flush=True)
    assert len(result)==99
    assert all(r.get('coordinates_unchanged') for r in result if r['previous_status']=='pass')
    assert all(r['status']=='pass' for r in result)
    return result

if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))
