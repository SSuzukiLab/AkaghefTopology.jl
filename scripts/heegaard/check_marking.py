"""Actual UI marking/export/import/redraw, plus incompatible-word rejection."""
from pathlib import Path
import copy
import json
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
    page.goto(base + '/surface_AbaloneAuto.html')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#marking-status')).to_contain_text('Unmarked')
    with page.expect_response('**/api/render', timeout=900_000) as response:
        page.locator('#mark-current').click()
    assert response.value.ok, response.value.text()
    result = response.value.json()
    marked = result['diagram']
    assert marked['schema'] == 'akaghef.heegaard.marked/1'
    assert marked['diagram']['joining_disk']['boundary_darts'] == [1, 8]
    expect(page.locator('#marking-status')).to_contain_text('Marked')
    expect(page.locator('#mark-current')).to_be_hidden()
    page.locator('#toggle-seams').click()
    expect(page.locator('#reference-seams')).to_be_visible()
    expect(page.locator('#reference-seams line')).to_have_count(2)
    with page.expect_download() as download:
        page.locator('#download-marked-project').click()
    saved = OUT / 'Abalone-ui.marked-project.json'
    download.value.save_as(str(saved))
    project = json.loads(saved.read_text(encoding='utf-8'))
    assert project['diagram'] == marked
    assert 'layout' not in marked['marking']
    # Import the one-file project through the user's dialog, without CLI steps.
    page.locator('#open-import').click()
    page.locator('#core-file').set_input_files(str(saved))
    with page.expect_response('**/api/render', timeout=900_000) as response:
        page.locator('#render').click()
    assert response.value.ok, response.value.text()
    loaded = response.value.json()
    assert loaded['diagram'] == marked
    assert loaded['marking']['hash'] == result['marking']['hash']
    with page.expect_response('**/api/render', timeout=900_000) as response:
        page.locator('#regenerate').click()
    assert response.value.ok, response.value.text()
    assert response.value.json()['diagram'] == marked
    # Same unmarked core and geometrically valid layout are NOT enough.
    bad = copy.deepcopy(marked)
    word = bad['marking']['routes'][0]['word']
    word.append(word[-1] if word else 1)  # reduced, syntactically valid, different
    rejection = page.request.post(base + '/api/render', data={
        'view': 'surface', 'diagram': bad, 'layout': result['layout']}, timeout=900_000)
    assert rejection.status == 422
    assert 'marking mismatch' in rejection.json()['error']
    missing = page.request.post(base + '/api/render', data={
        'view': 'surface', 'diagram': marked}, timeout=900_000)
    assert missing.status == 422
    assert 'requires its realization layout' in missing.json()['error']
    page.goto(base + '/surface_PoincareDMarked.html')
    page.wait_for_load_state('networkidle')
    expect(page.locator('#marking-status')).to_contain_text('Marked')
    expect(page.locator('#drawing [data-intersection-id]')).to_have_count(12)
    page.locator('#toggle-seams').click()
    expect(page.locator('#reference-seams line')).to_have_count(4)
    with page.expect_response('**/api/render', timeout=900_000) as response:
        page.locator('#regenerate').click()
    assert response.value.ok, response.value.text()
    poincare = response.value.json()
    assert poincare['diagram']['diagram']['joining_disk']['boundary_darts'] == [31, 24, 43, 8]
    assert len(poincare['diagram']['marking']['routes']) == 12
    page.locator('#toggle-seams').click()
    page.screenshot(path=str(OUT / 'Poincare-marked-page.png'), full_page=True)
    page.locator('#drawing').screenshot(path=str(OUT / 'Poincare-marked-chart.png'))
    assert not errors, errors
    browser.close()
print('PASS: explicit mark, project download/import, class-preserving redraw, wrong-word and missing-witness rejection; Poincare D unchanged')
