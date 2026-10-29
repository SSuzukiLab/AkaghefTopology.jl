"""Actual app route: 2g holes, D, view switching and fixture-free Koda redraw."""
import json
from pathlib import Path
import sys
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/heegaard/full_curves'
base = sys.argv[1] if len(sys.argv) > 1 else 'http://localhost:8741'
errors = []
with sync_playwright() as p:
    browser = p.chromium.launch(channel='msedge', headless=True)
    page = browser.new_page(viewport={'width': 1450, 'height': 1100})
    page.on('pageerror', lambda e: errors.append(str(e)))
    # The owner's old URL must land on the requested Koda 6-puncture example.
    page.goto(base + '/surface.html?cellulated=KodaGenus3')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#drawing svg[data-heegaard-view="alpha_cut"]')).to_have_count(1)
    expect(page.locator('#drawing circle[data-boundary-id]')).to_have_count(6)
    expect(page.locator('#drawing [data-chart="annulus"]')).to_have_count(0)
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_have_count(1)
    expect(page.locator('#drawing path[data-route-segment]')).to_have_count(11)
    page.locator('#toggle-disk').click()
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_be_hidden()
    page.locator('#toggle-disk').click()
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_be_visible()
    page.locator('#atlas').click()
    expect(page.locator('#drawing [data-chart="annulus"]')).to_have_count(3)
    page.locator('#cut').click()
    expect(page.locator('#drawing [data-chart="annulus"]')).to_have_count(0)
    with page.expect_download() as dl:
        page.locator('#download-svg').click()
    saved = Path(dl.value.path()).read_text(encoding='utf-8')
    assert 'data-heegaard-view="alpha_cut"' in saved
    assert 'data-disk-piece="base"' in saved
    print('PASS: Koda entry, 6 holes, D toggle, view switch and SVG download', flush=True)
    # Button -> API -> Julia, with no layout in the request.
    with page.expect_response('**/api/render', timeout=900_000) as response:
        page.locator('#regenerate').click()
    assert response.value.ok, response.value.text()
    result = response.value.json()
    original = json.loads((OUT/'KodaGenus3.punctured.hgd.json').read_text(encoding='utf-8'))
    assert result['diagram']['joining_disk'] == original['joining_disk']
    assert result['layout']['backend_trace']['fixture_used'] is False
    assert len(result['layout']['punctures']) == 6
    assert len(result['layout']['routes']) == 11
    assert result['view'] == 'alpha_cut'
    expect(page.locator('#drawing circle[data-boundary-id]')).to_have_count(6)
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_be_visible()
    (OUT/'KodaGenus3-alpha-cut-api.json').write_text(json.dumps(result), encoding='utf-8')
    page.locator('#drawing').screenshot(path=str(OUT/'KodaGenus3-alpha-cut.png'))
    page.screenshot(path=str(OUT/'KodaGenus3-alpha-cut-page.png'), full_page=True)
    print('PASS: fixture-free Koda regeneration, unchanged joining disk', flush=True)
    bad = json.loads(json.dumps(original))
    bad['joining_disk']['boundary_darts'][0] = 1
    rejected = page.request.post(base+'/api/render', data={'diagram': bad, 'view': 'alpha_cut'}, timeout=900_000)
    assert rejected.status == 422, rejected.text()
    assert 'boundary' in rejected.json()['error']
    for id, holes in [('AbaloneAuto', 2), ('PoincareDAuto', 4), ('S3', 2)]:
        page.goto(base+'/surface_'+id+'.html')
        page.wait_for_load_state('networkidle')
        expect(page.locator('#drawing circle[data-boundary-id]')).to_have_count(holes)
        expect(page.locator('#drawing [data-chart="annulus"]')).to_have_count(0)
    assert not errors, errors
    browser.close()
print('PASS: invalid D rejection; abalone, Poincare and ordinary diagram entry points')
