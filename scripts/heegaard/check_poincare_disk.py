"""Thin app check of the requested Poincare example with joining disk D0."""
from pathlib import Path
import json
import sys
import xml.etree.ElementTree as ET
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/heegaard/full_curves'
base = sys.argv[1] if len(sys.argv) > 1 else 'http://localhost:8741'
errors = []

def svg_with_canonical_metadata(svg):
    """Ignore only JSON object-key order, not geometry or provenance values."""
    root = ET.fromstring(svg)
    metadata = root.find('{http://www.w3.org/2000/svg}metadata')
    assert metadata is not None
    metadata.text = json.dumps(json.loads(metadata.text), sort_keys=True, ensure_ascii=False)
    return ET.tostring(root, encoding='unicode')

with sync_playwright() as p:
    browser = p.chromium.launch(channel='msedge', headless=True)
    page = browser.new_page(viewport={'width': 1450, 'height': 1250}, device_scale_factor=1)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(base + '/surface_PoincareD.html')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#title')).to_have_text('Poincaré + D₀')
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(12)
    expect(page.locator('#drawing circle[data-curve-family="alpha"]')).to_have_count(2)
    expect(page.locator('#drawing path[data-curve-family="beta"]')).to_have_count(24)
    expect(page.locator('#drawing path[data-disk-piece="base"]')).to_have_count(1)
    expect(page.locator('#drawing path[data-disk-piece="annulus"]')).to_have_count(2)
    expect(page.locator('#joining-boundary')).to_have_text('x4 → β1 → y1 → α2⁻¹ → y4 → β2 → x6 → α1⁻¹ → x4')
    expect(page.locator('#joining-note')).to_contain_text('逆にたどります')
    expect(page.locator('#cycles')).to_contain_text('x6+')
    expect(page.locator('#cycles')).to_contain_text('y4−')
    original_svg = page.locator('#drawing').inner_html()
    page.screenshot(path=str(OUT / 'PoincareD-page.png'), full_page=True)
    page.locator('#drawing').screenshot(path=str(OUT / 'PoincareD-chart.png'))
    page.locator('#joining-boundary').get_by_role('button', name='x6', exact=True).click()
    expect(page.locator('#status')).to_contain_text('交点 x6+')
    assert page.locator('#drawing [data-link-id="5"].selected').count() >= 7
    page.locator('#clear-selection').click()
    page.locator('#toggle-disk').click()
    assert page.locator('#drawing [data-joining-disk]').evaluate_all(
        '(items)=>items.every(e=>getComputedStyle(e).display==="none")')
    page.locator('#toggle-disk').click()
    page.locator('#legend button[data-family="beta"][data-id="2"]').click()
    expect(page.locator('#status')).to_contain_text('β2')
    page.locator('#clear-selection').click()

    downloads = {}
    for selector, key in [('#download-core', 'diagram'), ('#download-feasible', 'layout')]:
        with page.expect_download() as download:
            page.locator(selector).click()
        downloads[key] = download.value.path()
    diagram = json.loads(Path(downloads['diagram']).read_text(encoding='utf-8'))
    assert diagram['joining_disk']['boundary_darts'] == [31, 24, 43, 8]
    assert diagram['diagram']['intersections'][4]['source_label'] == 'x6'
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files(downloads['diagram'])
    page.locator('#layout-file').set_input_files(downloads['layout'])
    with page.expect_response('**/api/render', timeout=300_000) as response:
        page.locator('#render').click()
    assert response.value.ok, response.value.text()
    expect(page.locator('#import-dialog')).not_to_be_visible()
    expected = json.loads((OUT / 'PoincareD.json').read_text(encoding='utf-8'))['svg']
    actual = response.value.json()['svg']
    assert svg_with_canonical_metadata(actual) == svg_with_canonical_metadata(expected)
    assert svg_with_canonical_metadata(page.locator('#drawing').inner_html()) == svg_with_canonical_metadata(original_svg)
    expect(page.locator('#joining-note')).to_contain_text('逆にたどります')
    page.get_by_role('button', name='S³ · abalone + D', exact=True).click()
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(3)
    expect(page.locator('#joining-description')).to_contain_text('オーナー確認済み')
    page.get_by_role('button', name='Poincaré + D₀', exact=True).click()
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(12)
    assert not errors, errors
    browser.close()
print('PASS: Poincare labels, alpha/beta/D0, boundary trace, toggles, persisted input, same SVG after re-import, Abalone preserved; owner visual review pending')
