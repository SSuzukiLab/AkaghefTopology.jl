# cf:[[Cas65#^cas65-sec-0014]] section III opens; the explanatory paragraph of Figure 2 begins on that page
# cf:[[Cas65#^cas65-sec-0015]] Figure 2, the collapse of the cube onto X, and the sentence that X is a standard spine of the cube; the vertex count is inferred
# cf:[[Mat73#^mat73-sec-0012]] Example — the Bing house is a special spine of the three-dimensional cube
# cf:[[BP97#^bp97-sec-0547]] the local model prints an upper half-plane and a lower half-plane; calling the upper sheet the wall is inferred
# cf:[[BP97#^bp97-sec-0549]] Definition 2.1.4 and Figure 2.5 — calibration B, draft
# cf:[[BP97#^bp97-sec-0566]] Proposition 3.2.1 — this choice is the screw-consistent class whose colours match

using Topology

"""Coloured o-graph of Bing's house in the three-dimensional cube.

Casler's X is the standard spine. The singular graph has two true vertices
and four triple lines: the front equator, the back equator, the upper loop, and
the lower loop. Slots run counter-clockwise in the oriented disc D. At the
upper vertex the over pair is the equator. At the lower vertex the over pair is
the lower loop. Edge colours, written on the smaller half-edge, are -1 on
the front equator, +1 on the back equator, and 0 on both loops.

The opposite-germ pair and the orientation of D are the choice, among the
screw-consistent choices, whose colours satisfy Proposition 3.2.1. The
convention id stays the draft id.
"""
function bing_house_coloured_ograph()
    ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (3, 7, 1, 5, 4, 8, 2, 6),
        ((2, 4), (2, 4)),
        ((4, -1), (2, 1), (1, 0), (6, 0)))
end
