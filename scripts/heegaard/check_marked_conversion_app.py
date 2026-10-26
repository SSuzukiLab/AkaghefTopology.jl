"""Exercise marked source -> existing app button -> Koda converter correspondence."""
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/heegaard/full_curves'
oracle=json.loads((ROOT/'vendor/koda_genus3/data/figure6.odata.json').read_text(encoding='utf-8'))
def cyclic(a,b):
    return len(a)==len(b) and any(a[i:]+a[:i]==b for i in range(len(a)))

with sync_playwright() as pw:
    browser=pw.chromium.launch(channel='msedge',headless=True)
    page=browser.new_page(viewport={'width':1450,'height':1100});errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto('http://localhost:8741/surface_KodaGenus3Auto.html')
    page.wait_for_load_state('networkidle')
    assert page.locator('[data-crossing]').count()==11
    assert page.locator('[data-puncture]').count()==6
    before=page.locator('#drawing').inner_html()
    with page.expect_response(lambda r:r.url.endswith('/api/convert'),timeout=300000) as pending:
        page.locator('#convert-koda').click()
    response=pending.value;assert response.ok,response.text()
    result=response.json();graph=result['ograph']
    page.wait_for_function("!document.querySelector('#koda-section').hidden")
    assert cyclic(graph['gauss'][0],oracle[0][0][::-1])
    assert graph['orientation']==oracle[1]
    assert graph['metadata']['display_source_incidence_checked']
    assert result['closedness']['closed'] and all(result['closedness']['criteria'])
    assert page.locator('#drawing').inner_html()==before
    page.locator('#koda-vertices button').first.click()
    assert page.locator('#drawing [data-crossing].selected').count()==1
    selected=page.locator('#drawing [data-crossing].selected').get_attribute('data-crossing')
    assert selected in page.locator('#status').text_content()
    assert not errors,errors
    (OUT/'KodaMarkedConversion.json').write_text(json.dumps(result,ensure_ascii=False),encoding='utf-8')
    (OUT/'marked-conversion-app-check.json').write_text(json.dumps(dict(
        entry='existing marked surface view and /api/convert',
        source_incidence_checked=True,independent_odata_match=True,
        closedness=result['closedness'],original_svg_unchanged=True,
        source_crossing_selection=True,page_errors=errors),indent=2),encoding='utf-8')
    print('Existing Koda view: full marked source converted; independent o-data matched; original SVG retained; source crossing selected.',flush=True)
    browser.close()
