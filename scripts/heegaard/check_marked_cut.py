"""Inspect the actual marked-cut SVG/page with a headless browser."""
import argparse
from pathlib import Path
from urllib.parse import urljoin
from playwright.sync_api import sync_playwright

p=argparse.ArgumentParser();p.add_argument('url');p.add_argument('output',type=Path)
a=p.parse_args()
with sync_playwright() as pw:
    browser=pw.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1400,'height':1300},device_scale_factor=1)
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(a.url);page.wait_for_load_state('networkidle')
    page.locator('svg').wait_for()
    assert page.locator('[data-crossing]').count()==11
    assert page.locator('[data-puncture]').count()==6
    assert page.locator('[data-curve]').count()==6
    assert page.locator('[data-disk]').count()==1
    assert not errors,errors
    if a.url.startswith('http'):
        page.locator('#legend button').first.click()
        assert page.locator('[data-curve][style*="0.12"]').count()==5
        page.locator('#clear-selection').click()
        page.locator('#toggle-disk').click()
        assert not page.locator('[data-disk]').is_visible()
        page.locator('#toggle-disk').click()
        assert page.locator('[data-disk]').is_visible()
        assert not page.locator('#error').is_visible()
        with page.expect_download() as event:
            page.locator('#download-svg').click()
        download=event.value
        response=page.request.get(urljoin(a.url,'/surface-generated/'+download.suggested_filename))
        assert response.ok
        assert Path(download.path()).read_bytes()==response.body(),'Downloaded SVG differs from checked artifact'
    page.locator('svg').screenshot(path=str(a.output))
    print('Visible: 6 punctures, 6 curves, 11 intersections, D. No page errors. HTTP view: selection, D toggle and exact SVG download checked.' if a.url.startswith('http') else 'Visible: 6 punctures, 6 curves, 11 intersections, D. No page errors.')
    browser.close()
