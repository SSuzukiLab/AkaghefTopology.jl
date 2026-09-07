"""Small, explicit taxonomy for the common move graph.

The numeric values are intentionally user-facing: Area 2, Area 2.5, and Area
3 are different mathematical domains, not storage buckets.  A representation
change is never a topological move.
"""
module DomainTaxonomy

export AREA_2, AREA_2_5, AREA_3, VALID_AREAS, normalize_area,
       normalize_representation, is_topological_move, validate_topological_move,
       move_area

const AREA_2 = 2.0
const AREA_2_5 = 2.5
const AREA_3 = 3.0
const VALID_AREAS = (AREA_2, AREA_2_5, AREA_3)

function normalize_area(area)
    x = Float64(area)
    any(a -> isapprox(x, a; atol=0, rtol=0) , VALID_AREAS) ||
        throw(ArgumentError("area must be 2, 2.5, or 3"))
    x
end

normalize_representation(x) = lowercase(strip(String(x)))

const _AREA2 = Set(["r1", "r2", "r3", "reidemeister1", "reidemeister2", "reidemeister3"])
const _AREA25 = Set(["mixed", "mixed_move", "mixed_link", "lo", "lo_graph", "lo-graph",
                     "mixed/lo-graph", "mixed/lo-graph link", "lo-graph link"])
const _AREA3 = Set(["mp", "ps", "cp", "s", "0-2", "02", "0_2",
                    "mp_inverse", "ps_inverse", "cp_inverse", "s_inverse",
                    "bmp", "b02", "h",
                    "3-manifold mp", "3-manifold ps", "3-manifold cp", "3-manifold s"])

"Return whether a move name denotes a permitted topological move." 
function is_topological_move(move)
    s = lowercase(strip(String(move)))
    isempty(s) && return false
    occursin("transform", s) && return false
    occursin("representation", s) && return false
    occursin("serialize", s) && return false
    s in _AREA2 || s in _AREA25 || s in _AREA3
end

function move_area(move)
    s = lowercase(strip(String(move)))
    s in _AREA2 && return AREA_2
    s in _AREA25 && return AREA_2_5
    s in _AREA3 && return AREA_3
    throw(ArgumentError("unknown or non-topological move type: $move"))
end

function validate_topological_move(move, area)
    is_topological_move(move) ||
        throw(ArgumentError("representation changes and non-topological moves cannot enter transitions"))
    normalize_area(area) == move_area(move) ||
        throw(ArgumentError("move $(move) belongs to Area $(move_area(move)), not Area $(area)"))
    true
end

end
