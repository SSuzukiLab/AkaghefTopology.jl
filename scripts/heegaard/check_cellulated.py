"""Actual common app API + browser roundtrip; generated fixtures are not proof oracles."""
import copy
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/heegaard/cellulated'
BASE = 'http://localhost:8741'

with sync_playwright() as p:
    browser = p.chromium.launch(channel='msedge', headless=True)
    page = browser.new_page(viewport={'width':1440, 'height':1100})
    errors = []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(BASE+'/surface.html?cellulated=annulus_winding_p1')
    page.wait_for_load_state('networkidle')
    assert page.locator('#cell-words').inner_text().startswith('gamma [arc]')
    assert page.locator('#cellulated-nav button').count() == len(json.loads((OUT/'manifest.json').read_text(encoding='utf-8')))
    def post(body, status=200):
        response = page.request.post(BASE+'/api/render', data=body, timeout=180000)
        assert response.status == status, response.text()[-3500:]
        return response.json()
    checked=[]
    for item in json.loads((OUT/'manifest.json').read_text(encoding='utf-8')):
        saved=json.loads((OUT/(item['id']+'.json')).read_text(encoding='utf-8'))
        project={'schema':'akaghef.heegaard.cellulated-project/1','diagram':saved['diagram'],'layout':saved['layout']}
        rendered=post({'diagram':project})
        assert rendered['diagram']==saved['diagram']
        assert rendered['layout']==saved['layout']
        assert rendered['report']==saved['report']
        checked.append(item['id'])
        print('roundtrip',item['id'],flush=True)
    small=json.loads((OUT/'annulus_winding_p1.json').read_text(encoding='utf-8'))
    plain=post({'diagram':small['diagram']})
    assert plain['layout']==small['layout']  # no geometry input, no route inference
    wrong=copy.deepcopy(small['layout']);wrong['representation_sha256']='0'*64
    assert 'CELL_LAYOUT_HASH_MISMATCH' in post({'diagram':small['diagram'],'layout':wrong},422)['error']
    wrong=copy.deepcopy(small['diagram']);wrong['curves'][0]['word']='P{v1_0 absent}'
    assert 'ROUTE_EDGE_MISSING' in post({'diagram':wrong},422)['error']
    wrong=copy.deepcopy(small['layout']);wrong['positions']['P']['v1_1']=wrong['positions']['P']['v2_2']
    assert 'COINCIDENT_VERTICES' in post({'diagram':small['diagram'],'layout':wrong},422)['error']
    marked=json.loads((OUT/'torus_three_crossings_D.json').read_text(encoding='utf-8'))
    wrong=copy.deepcopy(marked['diagram']);wrong['region_marks'][0]['claim']='koda_joining'
    assert 'REGION_CLAIM_FAILED' in post({'diagram':wrong},422)['error']
    response=page.request.post(BASE+'/api/convert',data={'diagram':marked['diagram']},timeout=180000)
    assert response.status==422 and 'stored-orientation' in response.json()['error']
    # Actual UI save, upload and render through the common dialog.
    with page.expect_download() as d:
        page.locator('#cell-project').click()
    project_path=OUT/'browser-roundtrip.cellulated-project.json';d.value.save_as(project_path)
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files(project_path)
    page.locator('#render').click()
    page.wait_for_function("!document.querySelector('#import-dialog').open",timeout=180000)
    page.locator('#legend button').first.click()
    assert page.locator('#legend button').first.get_attribute('aria-pressed')=='true'
    page.locator('#toggle-mesh').click()
    assert page.locator('#drawing').evaluate("e=>e.classList.contains('hide-mesh')")
    # Empty-intersection, disconnected and nested objects use exactly the same UI.
    for name in ['torus_parallel_zero_crossings','nested_zero_crossings','genus_2_empty','torus_three_crossings_D']:
        page.locator(f'#cellulated-nav button[data-id="{name}"]').click()
        page.wait_for_function('(id)=>current.id===id',arg=name)
        assert page.locator('#drawing svg').count()==1
    assert 'Koda joining = false' in page.locator('#joining-description').inner_text()
    page.screenshot(path=str(OUT/'cellulated-D.png'),full_page=True)
    page.locator('#examples button[data-id="L2_1"]').click()
    page.wait_for_function("current.id==='L2_1'")
    assert not page.locator('body').evaluate("e=>e.classList.contains('cellulated')")
    assert page.locator('#drawing svg').count()>0
    page.locator('#cellulated-nav button[data-id="annulus_winding_p1"]').click()
    page.wait_for_function("current.id==='annulus_winding_p1'")
    page.screenshot(path=str(OUT/'cellulated-winding.png'),full_page=True)
    assert not errors,errors
    (OUT/'app-verification.json').write_text(json.dumps({'roundtrip':checked,'negative_controls':['hash','word','overlap','Koda claim','invalid joining orientation'],'ui':['download/upload','curve selection','mesh toggle','D condition distinction','legacy L2_1 switch'],'page_errors':errors},indent=2),encoding='utf-8')
    browser.close()
