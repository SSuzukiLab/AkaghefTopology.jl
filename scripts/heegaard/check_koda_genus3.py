"""Thin end-to-end verification of the actual app import/draw/convert workflow."""
import copy
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/heegaard/cellulated'
BASE='http://localhost:8741'
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1500,'height':1200})
    errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(BASE+'/surface.html?cellulated=KodaGenus3')
    page.wait_for_load_state('networkidle')
    assert 'genus 3' in page.locator('#subtitle').inner_text()
    assert '11 signed intersections' in page.locator('#subtitle').inner_text()
    assert 'Koda joining = true' in page.locator('#joining-description').inner_text()
    assert page.locator('#drawing [data-point="x2"]').count()==1
    assert page.locator('#drawing [data-sheet="B"]').count()==7
    page.locator('#legend button[data-curve="alpha1"]').click()
    assert page.locator('#drawing [data-curve="beta1"].dim').count()>0
    page.locator('#clear-selection').click()
    page.locator('#toggle-disk').click()
    assert page.locator('#drawing .region').first.is_hidden()
    page.locator('#toggle-disk').click()
    page.locator('#cell-atlas').click()
    assert page.locator('#drawing .mesh').count()>0
    page.locator('#cell-projection').click()
    page.locator('#drawing').screenshot(path=str(OUT/'KodaGenus3.projection.png'))
    # Download the exact marked input and upload using the existing app dialog.
    with page.expect_download() as d:page.locator('#cell-project').click()
    project=OUT/'KodaGenus3.cellulated-project.json';d.value.save_as(project)
    page.locator('#open-import').click();page.locator('#core-file').set_input_files(project)
    page.locator('#render').click()
    page.wait_for_function("!document.querySelector('#import-dialog').open",timeout=240000)
    saved=json.loads(project.read_text(encoding='utf-8'))
    assert page.evaluate('current.diagram')==saved['diagram']
    assert page.evaluate('current.layout')==saved['layout']
    assert page.locator('#drawing [data-point="x2"]').count()==1
    print('import/draw/save/reload passed',flush=True)
    # Use the app's real conversion button, not a canned result from the bundle.
    page.locator('#convert-koda').click()
    page.wait_for_function("document.querySelector('#operation-status').textContent==='Koda変換完了' || !document.querySelector('#convert-koda').disabled",timeout=300000)
    assert page.locator('#operation-status').inner_text()=='Koda変換完了',page.locator('#operation-status').inner_text()
    result=page.evaluate('current.koda')
    assert result['closedness']=={'closed':True,'criteria':[True,True,True]}
    assert result['ograph']['orientation']==[-1,1,1,1,-1]
    word=result['ograph']['gauss'][0]
    source=json.loads((ROOT/'vendor/koda_genus3/data/figure6.odata.json').read_text())
    reverse=source[0][0][::-1]
    assert any(word==reverse[i:]+reverse[:i] for i in range(len(word)))
    assert len(result['complex']['edge_incidence'])==10 and len(result['complex']['regions'])==6
    assert len(result['cellulated_projection']['region_map'])==7
    page.locator('#koda-vertices button').nth(1).click()
    assert 'x2' in page.locator('#status').inner_text()
    (OUT/'KodaGenus3.koda.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'KodaGenus3.repository.odata.json').write_text(json.dumps([result['ograph']['gauss'],result['ograph']['orientation']]),encoding='utf-8')
    print('actual Koda conversion/closedness/provenance passed',flush=True)
    bad=copy.deepcopy(saved['diagram'])
    next(x for x in bad['crossings'] if x['id']=='x2')['sign']*=-1
    r=page.request.post(BASE+'/api/render',data={'diagram':bad},timeout=240000)
    assert r.status==422,r.text()
    false_disk=json.loads((OUT/'torus_three_crossings_D.json').read_text())['diagram']
    r=page.request.post(BASE+'/api/convert',data={'diagram':false_disk},timeout=240000)
    assert r.status==422 and 'stored-orientation' in r.json()['error']
    page.locator('#examples button[data-id="L2_1"]').click()
    page.wait_for_function("current.id==='L2_1'")
    assert not page.locator('body').evaluate("e=>e.classList.contains('cellulated')")
    page.locator('#cellulated-nav button[data-id="KodaGenus3"]').click()
    page.wait_for_function("current.id==='KodaGenus3'")
    assert not errors,errors
    (OUT/'KodaGenus3.app-verification.json').write_text(json.dumps({
        'genus':3,'intersections':11,'D_corners':6,'o_vertices':5,'o_edges':10,'o_regions':6,
        'closedness':result['closedness'],'source_convention':'reverse complete traversal; unchanged vertex signs',
        'repository_odata':[result['ograph']['gauss'],result['ograph']['orientation']],
        'source_representation_sha256':result['cellulated_projection']['source_representation_sha256'],
        'ui':['front/back projection','curve selection','D toggle','certified atlas switch','save/upload roundtrip','live Koda button','vertex provenance','legacy switch'],
        'negative_controls':['wrong backside sign','geometric but not oriented joining disk'],
        'page_errors':errors},ensure_ascii=False,indent=2),encoding='utf-8')
    page.screenshot(path=str(OUT/'KodaGenus3.app.png'),full_page=True)
    browser.close()
    print('Koda genus3 app checks passed',flush=True)
