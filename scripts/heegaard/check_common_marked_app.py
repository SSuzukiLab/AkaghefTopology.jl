"""Exercise the existing application's marked-display import and plot endpoint."""
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/heegaard/full_curves'
with sync_playwright() as pw:
    browser=pw.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1450,'height':1100});errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto('http://localhost:8741/surface_LickorishG4Common.html');page.wait_for_load_state('networkidle')
    assert page.locator('[data-puncture]').count()==8
    assert page.locator('[data-curve]').count()==7
    assert page.locator('#toggle-disk').is_hidden()
    assert '図12.3' in page.locator('#marked-source').text_content()
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files(OUT/'LickorishG2.view.json')
    with page.expect_response(lambda r:r.url.endswith('/api/render'),timeout=180000) as pending:
        page.locator('#render').click()
    response=pending.value;assert response.ok,response.text();data=response.json()
    page.wait_for_function("document.querySelectorAll('[data-puncture]').length===4")
    assert page.locator('[data-curve]').count()==3
    assert data['layout_isotopy_check']['collisions']==0
    assert data['smooth_check']['minimum_nonport_hole_clearance']>=12
    assert data['genus']==2
    assert data['source_figure']['cite_key']=='AnonND60'
    assert data['source_figure']['asset'].endswith('/p0251-0001.jpeg')
    with page.expect_download() as download:page.locator('#download-svg').click()
    assert Path(download.value.path()).read_text(encoding='utf-8')==data['svg']
    with page.expect_response(lambda r:r.url.endswith('/api/render'),timeout=180000) as pending:
        page.locator('#regenerate').click()
    second=pending.value;assert second.ok,second.text()
    assert second.json()['genus']==2
    assert second.json()['source_figure']==data['source_figure']
    # Reject a missing pairing before any substitute routing is attempted.
    bad=json.loads((OUT/'LickorishG1.view.json').read_text(encoding='utf-8'));bad.pop('boundary_vertex_pairs')
    rejected=page.request.post('http://localhost:8741/api/render',data={'diagram':bad},timeout=180000)
    assert rejected.status==422
    assert 'boundary vertex pairing is required' in rejected.json()['error']
    assert not errors,errors
    page.locator('svg').screenshot(path=str(OUT/'LickorishG2Common-app.png'))
    (OUT/'common-marked-app-check.json').write_text(json.dumps(dict(
        entry='existing surface.html and /api/render',genus4_display=True,genus2_import=True,
        genus2_regeneration=True,svg_download_exact=True,missing_pairing_rejected=True,
        library_source_displayed=True,library_source_roundtrip=True,page_errors=errors),indent=2),encoding='utf-8')
    print('Existing app: genus 4 display; genus 2 import, regeneration, SVG download; missing pairing rejected.',flush=True)
    browser.close()
