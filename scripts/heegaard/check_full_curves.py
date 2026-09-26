"""Thin UI acceptance check against the already-running local atlas page."""
import base64
import json
from pathlib import Path
import sys
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "artifacts" / "heegaard" / "full_curves"
base = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8741"
errors = []
with sync_playwright() as p:
    browser = p.chromium.launch(channel="msedge", headless=True)
    page = browser.new_page(viewport={"width": 1500, "height": 1080}, device_scale_factor=1)
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(base + "/surface_Poincare.html")
    page.wait_for_load_state("networkidle")
    page.locator('#drawing [data-intersection-id="12"]').wait_for()
    assert page.locator('#drawing [data-intersection-id]').count() == 12
    assert page.locator('#drawing circle[data-curve-family="alpha"]').count() == 2
    page.screenshot(path=str(OUTPUT / "Poincare-page.png"), full_page=True)
    page.locator("#drawing").screenshot(path=str(OUTPUT / "Poincare-chart.png"))

    page.get_by_role("button", name="β1", exact=True).first.click()
    assert "β1" in page.locator("#status").inner_text()
    assert page.locator('#drawing [data-curve-family="beta"][data-curve-id="2"].dim').count() > 0
    page.locator('#cycles button[data-link-id="6"]').first.click()
    assert "交点 6+" in page.locator("#status").inner_text()
    assert page.locator('#drawing [data-link-id="6"].selected').count() >= 7
    page.locator("#drawing").screenshot(path=str(OUTPUT / "Poincare-intersection-6.png"))
    page.get_by_role("button", name="alpha-cutと対照", exact=True).click()
    assert page.locator("#comparison").is_visible()
    assert page.locator('#base-drawing svg[data-heegaard-view="alpha_cut"]').count() == 1
    page.get_by_role("button", name="両曲線", exact=True).click()
    page.get_by_role("button", name="選択解除", exact=True).click()
    assert page.locator("#drawing .dim").count() == 0

    for title, id, count in [("S³", "S3", 1), ("L(5,2)", "L5_2", 5)]:
        page.get_by_role("button", name=title, exact=True).click()
        page.wait_for_url(base + "/surface_" + id + ".html")
        assert page.locator('#drawing [data-intersection-id]').count() == count
        page.locator("#drawing").screenshot(path=str(OUTPUT / (id + "-chart.png")))
    with page.expect_download() as download:
        page.locator("#download-svg").click()
    downloaded = download.value
    assert downloaded.suggested_filename == "L5_2-surface.svg"
    assert Path(downloaded.path()).read_text(encoding="utf-8") == (OUTPUT / "L5_2.svg").read_text(encoding="utf-8")

    # Follow the visible link to the old application and verify that it still renders.
    page.locator("#original-link").click()
    page.wait_for_load_state("networkidle")
    page.locator('#after-svg svg[data-heegaard-view="alpha_cut"]').wait_for()
    assert page.locator("#title").inner_text() == "L(5,2)"
    assert not errors, errors
    # A narrow viewport keeps the atlas scrollable instead of clipping it.
    page.set_viewport_size({"width": 390, "height": 844})
    page.goto(base + "/surface_S3.html")
    page.wait_for_load_state("networkidle")
    page.locator('#drawing [data-intersection-id="1"]').wait_for()
    assert page.locator("#error").is_hidden()
    assert page.locator(".surface-drawing").evaluate("(e) => e.scrollWidth > e.clientWidth")
    page.screenshot(path=str(OUTPUT / "S3-mobile.png"), full_page=True)
    page.goto(base + "/surface_missing_fixture.html")
    page.locator("#error").wait_for(state="visible")
    assert page.locator("#drawing svg").count() == 0
    # Check rendered notebook images when the separate Quarto execution exists.
    qmd_html = ROOT / "examples" / "HeegaardFullCurves.html"
    qmd_checked = False
    if qmd_html.is_file():
        page.set_viewport_size({"width": 1400, "height": 1000})
        page.goto(qmd_html.as_uri())
        images = page.locator('img[src^="data:image/svg+xml;base64,"]')
        assert images.count() == 4
        for index, id in enumerate(["S3", "L5_2", "Poincare"]):
            source = images.nth(index).get_attribute("src")
            assert base64.b64decode(source.split(",", 1)[1]) == (OUTPUT / (id + ".svg")).read_bytes()
            assert images.nth(index).evaluate("e => e.complete && e.naturalWidth > 0")
        images.nth(2).screenshot(path=str(OUTPUT / "Poincare-qmd.png"))
        qmd_checked = True
    assert not errors, errors
    result = {"url": base + "/surface_Poincare.html", "page_errors": errors,
              "quarto_html_svg_equality": qmd_checked,
              "examples": ["Poincare", "S3", "L5_2"],
              "checked": ["both families", "intersection tracking", "beta highlighting",
                          "alpha-cut comparison", "SVG download equality", "old application link",
                          "narrow viewport", "missing fixture fails explicitly"],
              "screenshots": ["Poincare-page.png", "Poincare-chart.png", "Poincare-intersection-6.png",
                              "S3-chart.png", "L5_2-chart.png", "S3-mobile.png"]}
    (OUTPUT / "browser-check.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    browser.close()
