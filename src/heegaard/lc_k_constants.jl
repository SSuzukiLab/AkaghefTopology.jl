# Source: docs/specs/lickorish-code.md §4.2; scripts/heegaard/lc_input.py::k_constants.
"""Constant K words for genus g (no genus branching in callers)."""

function _lc_canon_cyclic(word::Vector{Int})
    isempty(word) && return Tuple{Int}()
    n = length(word)
    minimum(tuple(word[i:n]..., word[1:i - 1]...) for i in 1:n)
end

function _lc_k_constants(g::Integer)
    g < 1 && return Dict{String, Vector{Int}}()
    lm = Dict(i => i for i in 1:g)
    ln = Dict{Tuple{Int, Int}, Int}()
    label = g
    for i in 1:(g - 1)
        ln[i, i] = label + 1
        ln[i + 1, i] = label + 2
        label += 2
    end
    words = Dict{String, Vector{Int}}()
    for i in 1:g
        w = [-lm[i]]
        i < g && push!(w, -ln[i, i])
        i > 1 && push!(w, ln[i, i - 1])
        words["l$i"] = w
        words["m$i"] = [lm[i]]
    end
    for i in 1:(g - 1)
        words["n$i"] = [ln[i, i], -ln[i + 1, i]]
    end
    words
end

function _lc_validate_k_constants!(lc::LickorishConfiguration)
    g = lc.genus
    g < 1 && return true
    ref = _lc_k_constants(g)
    for k in _lc_k_names(g)
        haskey(lc.words, k) || throw(KodaInputError(:realize_lc, "missing K component", Dict("component" => k)))
        _lc_canon_cyclic(lc.words[k]) == _lc_canon_cyclic(ref[k]) ||
            throw(KodaInputError(:realize_lc, "K component does not match reference constants", Dict("component" => k)))
    end
    true
end
