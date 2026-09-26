"""The actual Poincare app regeneration entry, without any layout in the request."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright, expect
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/heegaard/full_curves'
errors=[]
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1450,'height':1200})
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto('http://localhost:8741/surface_PoincareD.html')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#title')).to_have_text('Poincaré + D₀')
    page.locator('#convert-koda').click()
    expect(page.locator('#koda-summary')).to_contain_text('8 頂点 · 16 辺 · 9 領域 · 2 blocks')
    page.get_by_role('button',name='v5 ← x7 (+)',exact=True).click()
    expect(page.locator('#status')).to_contain_text('交点 x7+')
    page.locator('#koda-section summary').click()
    with page.expect_download() as downloaded:
        page.locator('#download-koda').click()
    koda=json.loads(Path(downloaded.value.path()).read_text(encoding='utf-8'))
    assert len(koda['complex']['regions'])==9
    assert koda['ograph']['metadata']['o_block_count']==2
    with page.expect_response('**/api/render',timeout=900_000) as response:
        page.locator('#regenerate').click()
    assert response.value.ok,response.value.text()
    request=response.value.request.post_data_json
    assert 'layout' not in request
    result=response.value.json()
    assert result['diagram']==request['diagram']
    trace=result['layout']['backend_trace']
    assert trace['fixture_used'] is False
    assert trace['post_router']=='terminal-planar-realization'
    assert 'exterior_compaction' in trace
    assert trace['plot']['profile']=='fixed-100'
    points=[p for r in result['layout']['routes'] for p in r['points']]
    assert all(max(p[k] for p in points)-min(p[k] for p in points)<8 for k in (0,1))
    assert result['chart']['joining_disk']['disk']['boundary_darts']==[31,24,43,8]
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(12)
    expect(page.locator('#drawing circle[data-curve-family="alpha"]')).to_have_count(2)
    expect(page.locator('#drawing path[data-curve-family="beta"]')).to_have_count(24)
    page.locator('#joining-boundary').get_by_role('button',name='x6',exact=True).click()
    expect(page.locator('#status')).to_contain_text('交点 x6+')
    page.locator('#toggle-disk').click()
    assert page.locator('#drawing [data-joining-disk]').evaluate_all(
        '(nodes)=>nodes.every(n=>getComputedStyle(n).display==="none")')
    page.locator('#toggle-disk').click()
    page.locator('#clear-selection').click()
    page.locator('#drawing').screenshot(path=str(OUT/'PoincareAuto-chart.png'))
    page.screenshot(path=str(OUT/'PoincareAuto-page.png'),full_page=True)
    (OUT/'PoincareAuto-api.json').write_text(json.dumps(result,ensure_ascii=False),encoding='utf-8')
    assert not errors,errors
    browser.close()
print('PASS: Poincare Koda trace/download; fixture-free API regeneration; unchanged labelled core/D; both curve families, D toggle/selection; no JS errors')
