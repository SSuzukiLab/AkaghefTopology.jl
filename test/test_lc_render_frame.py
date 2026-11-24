"""cf: lib/topology.jl/docs/conventions/LickorishCodes.md — SVG screen orientation."""

from pathlib import Path
import sys
import unittest
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "heegaard"))
from lc_render import svg


class LCRenderFrameTest(unittest.TestCase):
    def test_screen_orientation(self):
        for genus in (1, 2, 3):
            with self.subTest(genus=genus):
                holes = {}
                if genus > 1:
                    for i in range(1, genus + 1):
                        holes[f"top{i}"] = ((120 * i, 80), 10)
                        holes[f"bot{i}"] = ((120 * i, 20), 10)
                data = dict(W=480, H=100, holes=holes, square=genus == 1,
                            curves={"b": dict(points=[(20, 80), (30, 50)], closed=False)})
                root = ElementTree.fromstring(svg(data))
                ns = {"s": "http://www.w3.org/2000/svg"}
                text_y = {label.text: float(label.attrib["y"]) for label in root.findall("s:text", ns)}
                if genus == 1:
                    self.assertGreater(text_y["ℓ₁⁺ →"], text_y["ℓ₁⁻ →"])
                else:
                    for i in range(1, genus + 1):
                        self.assertGreater(text_y[f"ℓ{i}⁺"], text_y[f"ℓ{i}⁻"])
                self.assertIn("M20,80 L30,50", [p.attrib["d"] for p in root.findall("s:path", ns)])


if __name__ == "__main__":
    unittest.main()
