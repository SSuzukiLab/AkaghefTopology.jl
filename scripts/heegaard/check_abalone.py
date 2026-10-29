"""Check abalone through the existing surface app, including file re-import."""
from pathlib import Path
import json
import sys
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/heegaard/full_curves/abalone'
base = sys.argv[1] if len(sys.argv) > 1 else 'http://localhost:8741'
errors = []
with sync_playwright() as p:
    browser = p.chromium.launch(channel='msedge', headless=True)
    page = browser.new_page(viewport={'width': 1400, 'height': 1150}, device_scale_factor=1)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(base + '/abalone.html')
    page.wait_for_load_state('networkidle')
    page.wait_for_url(base + '/surface_Abalone.html')
    expect(page.locator('#title')).to_have_text('S³ · abalone + D')
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(3)
    expect(page.locator('#drawing circle[data-curve-family="alpha"]')).to_have_count(1)
    expect(page.locator('#drawing path[data-curve-family="beta"]')).to_have_count(6)
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_have_count(1)
    expect(page.locator('#drawing path[data-disk-piece="annulus"]')).to_have_count(1)
    assert '1 → α1 → 2 → β1⁻¹ → 1' in page.locator('#joining-boundary').inner_text()
    page.screenshot(path=str(OUT / 'abalone-app.png'), full_page=True)
    page.locator('#drawing').screenshot(path=str(OUT / 'abalone-chart.png'))
    original_svg = page.locator('#drawing').inner_html()

    page.locator('#toggle-disk').click()
    assert page.locator('#drawing [data-joining-disk]').evaluate_all(
        '(items)=>items.every(e=>getComputedStyle(e).display==="none")')
    page.locator('#toggle-disk').click()
    expect(page.locator('#drawing path[data-disk-piece="annulus"]')).to_be_visible()
    page.locator('#joining-boundary button').first.click()
    assert '交点 1+' in page.locator('#status').inner_text()
    page.locator('#legend button[data-family="alpha"]').click()
    assert 'α1' in page.locator('#status').inner_text()
    page.locator('#clear-selection').click()

    downloads = {}
    for selector, key in [('#download-core', 'diagram'), ('#download-feasible', 'layout')]:
        with page.expect_download() as download:
            page.locator(selector).click()
        downloads[key] = download.value.path()
    diagram = json.loads(Path(downloads['diagram']).read_text(encoding='utf-8'))
    assert diagram['schema'] == 'akaghef.heegaard.punctured/1'
    assert diagram['joining_disk']['boundary_darts'] == [1, 8]
    assert diagram == json.loads((OUT / 'S3.punctured.hgd.json').read_text(encoding='utf-8'))
    # Reopen the files actually offered by this page, using the public dialog.
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files(downloads['diagram'])
    page.locator('#layout-file').set_input_files(downloads['layout'])
    with page.expect_response('**/api/render', timeout=300_000) as response:
        page.locator('#render').click()
    assert response.value.ok, response.value.text()
    expect(page.locator('#import-dialog')).not_to_be_visible()
    expect(page.locator('#title')).to_have_text('Imported Heegaard diagram + D')
    assert page.locator('#drawing').inner_html() == original_svg
    expect(page.locator('#joining-section')).to_be_visible()

    # Invalid D must be reported, never silently dropped from the input.
    diagram['joining_disk']['boundary_darts'] = [1, 7]
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files({
        'name': 'invalid.punctured.hgd.json', 'mimeType': 'application/json',
        'buffer': json.dumps(diagram).encode('utf-8')})
    with page.expect_response('**/api/render', timeout=300_000) as response:
        page.locator('#render').click()
    assert response.value.status == 422
    expect(page.locator('#import-status')).to_contain_text('joining disk')
    expect(page.locator('#import-dialog')).to_be_visible()
    assert page.locator('#drawing').inner_html() == original_svg
    page.locator('#close-import').click()

    for title, count in [('S³', 1), ('Poincaré homology sphere', 12)]:
        page.get_by_role('button', name=title, exact=True).click()
        expect(page.locator('#drawing [data-intersection-id]')).to_have_count(count)
        expect(page.locator('#joining-section')).not_to_be_visible()
    page.get_by_role('button', name='S³ · abalone + D', exact=True).click()
    expect(page.locator('#joining-section')).to_be_visible()
    assert not errors, errors
    browser.close()
print('PASS: common app, alpha/beta/D, tracing, toggle, downloads, file re-import, invalid-D rejection, baseline navigation; owner visual review pending')
