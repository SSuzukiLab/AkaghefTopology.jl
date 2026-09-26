"""Optional preview test. Requires Playwright and a Chromium executable.
Does not fetch remote pages: it renders this bundle's self-contained HTML.
"""
from pathlib import Path
import json,os,shutil
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
exe=os.environ.get('CHROMIUM_PATH') or shutil.which('chromium') or shutil.which('google-chrome')
with sync_playwright() as p:
    browser=p.chromium.launch(headless=True,**({'executable_path':exe} if exe else {}))
    page=browser.new_page(viewport={'width':1440,'height':1060},device_scale_factor=1)
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.set_content((ROOT/'overview.html').read_text())
    page.screenshot(path=str(ROOT/'artifacts'/'overview-top.png'))
    page.locator('#model').select_option('torus_three_crossings');page.locator('#examples').scroll_into_view_if_needed()
    page.screenshot(path=str(ROOT/'artifacts'/'overview-diagram.png'))
    assert page.locator('#figure .region').count()>0
    count=0
    for option in page.locator('#model option').evaluate_all('(els)=>els.map(e=>e.value)'):
        page.locator('#model').select_option(option)
        assert page.locator('#figure svg').count()==1 and page.locator('#words').inner_text().strip();count+=1
    page.locator('#mesh').check();assert 'hide-mesh' not in page.locator('body').get_attribute('class')
    assert not errors,errors
    browser.close()
result={'executed':True,'examples_switched':count,'mesh_toggle':True,'region_fill':True,'javascript_errors':errors}
(ROOT/'artifacts'/'html-test-status.json').write_text(json.dumps(result,indent=2))
print(result)
