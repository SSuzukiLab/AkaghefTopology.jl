"""App entry points: independent Koda panel, fixture-free regeneration and D preservation."""
from pathlib import Path
import json
import sys
from playwright.sync_api import sync_playwright, expect

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/heegaard/full_curves'
base=sys.argv[1] if len(sys.argv)>1 else 'http://localhost:8741'
errors=[]
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1450,'height':1100})
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(base+'/surface_AbaloneAuto.html')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#title')).to_contain_text('自動配置')
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(3)
    page.locator('#convert-koda').click()
    expect(page.locator('#koda-summary')).to_contain_text('1 頂点 · 2 辺 · 2 領域 · 1 blocks')
    with page.expect_download() as download:
        page.locator('#koda-section summary').click()
        page.locator('#download-koda').click()
    koda=json.loads(Path(download.value.path()).read_text(encoding='utf-8'))
    assert koda['ograph']['gauss']==[[1,-1]]
    assert koda['complex']['checks']['oriented_branching']
    # Real user button -> API, with no supplied layout.
    with page.expect_response('**/api/render',timeout=900_000) as response:
        page.locator('#regenerate').click()
    assert response.value.ok,response.value.text()
    result=response.value.json()
    assert result['layout']['backend_trace']['fixture_used'] is False
    assert result['diagram']['joining_disk']['boundary_darts']==[1,8]
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(3)
    page.screenshot(path=str(OUT/'AbaloneAuto-page.png'),full_page=True)
    # Invalid disk is rejected, not changed or guessed by the conversion endpoint.
    bad=json.loads(json.dumps(result['diagram']))
    bad['joining_disk']['boundary_darts']=[1,7]
    rejection=page.request.post(base+'/api/convert',data={'diagram':bad},timeout=900_000)
    assert rejection.status==422
    assert 'boundary' in rejection.json()['error']
    assert not errors,errors
    browser.close()
print('PASS: fixture-free app regeneration, unchanged D, Koda panel/download, invalid D rejection')
