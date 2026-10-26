"""Diagnostic closed-component/D regression and malformed-port negative control.

The added contractible circle is a geometry regression, not a manifold example
or acceptance evidence for arbitrary Heegaard diagrams.
"""
import copy,json
from pathlib import Path
import numpy as np
from marked_display_pipeline import plot,pair_ports

root=Path(__file__).resolve().parents[2]/'artifacts/heegaard/full_curves'
d=json.loads((root/'LickorishG1.view.json').read_text(encoding='utf-8'))
first=max(map(int,d['positions']))+1;vs=list(range(first,first+16))
for i,v in enumerate(vs):
    p=np.array([160.,0.])+45*np.array([np.cos(2*np.pi*i/16),np.sin(2*np.pi*i/16)])
    d['positions'][str(v)]=(p/600).tolist()
d['curves'].append(dict(id='closed component',role='alpha',segments=[dict(source=['diagnostic',str(u),str(v)],vertices=[u,v]) for u,v in zip(vs,vs[1:]+vs[:1])]))
d['disk']=vs
result,_=plot(d,maxiter=900)
assert result['disk'] and result['smooth_check']['source_disk_boundary_unchanged']
assert result['layout_isotopy_check']['collisions']==0
bad=copy.deepcopy(d);v=bad['boundaries'][0]['ports'][0]['vertex'];bad['positions'][str(v)][0]+=.01
try:pair_ports(bad)
except AssertionError as e:assert 'not on its hole boundary' in str(e)
else:raise AssertionError('Malformed port accepted')
print('Diagnostic: closed component, optional D, malformed port rejection passed.',flush=True)
