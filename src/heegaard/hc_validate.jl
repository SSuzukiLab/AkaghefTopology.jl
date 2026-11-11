# cf: lib/topology.jl/docs/specs/lickorish-code.md §3–§4 — 各ラベルは + と − が一度ずつ。オイラー標数 2−2g の比較は _lc_validate_surface! のまま
# cf: lib/topology.jl/docs/specs/heegaard-code.md §2.1 — 非 K 成分の役割は alpha または beta を一つ

"""Validate a resolved Heegaard diagram.

The argument is `HeegaardLayer`, the same struct as `Topology.HeegaardDiagram`.
Returns `true` when the word is consistent. Throws `KodaInputError` when a
label has no opposite passage, a non-K curve has no alpha or beta role, or
the Euler characteristic does not match the stated genus.
"""
function validate_diagram(diagram::HeegaardLayer)
    _validate_diagram_passages(diagram.lc.words)
    _validate_diagram_roles(diagram)
    validate_hc_input(diagram)
    _lc_validate_surface!(diagram.lc)
    true
end

function _validate_diagram_passages(words::AbstractDict)
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in words
        for x in word
            x == 0 && throw(KodaInputError(
                :validate_diagram, "label 0 is not a passage",
                Dict{String,Any}("curve" => String(name))))
            lab = abs(Int(x))
            dest = x > 0 ? plus : minus
            haskey(dest, lab) && throw(KodaInputError(
                :validate_diagram, "passage sign is not unique",
                Dict{String,Any}("label" => lab, "curve" => String(name))))
            dest[lab] = String(name)
        end
    end
    for lab in sort!(collect(keys(plus)))
        haskey(minus, lab) || throw(KodaInputError(
            :validate_diagram, "label has no opposite passage",
            Dict{String,Any}("label" => lab, "missing" => "-")))
    end
    for lab in sort!(collect(keys(minus)))
        haskey(plus, lab) || throw(KodaInputError(
            :validate_diagram, "label has no opposite passage",
            Dict{String,Any}("label" => lab, "missing" => "+")))
    end
    nothing
end

function _validate_diagram_roles(diagram::HeegaardLayer)
    knames = Set(_lc_k_names(diagram.lc.genus))
    for name in sort!(collect(String.(keys(diagram.lc.words))))
        name in knames && continue
        haskey(diagram.roles, name) || throw(KodaInputError(
            :validate_diagram, "alpha/beta role missing",
            Dict{String,Any}("curve" => name)))
        role = diagram.roles[name]
        role in (:alpha, :beta) || throw(KodaInputError(
            :validate_diagram, "invalid role",
            Dict{String,Any}("curve" => name, "role" => String(role))))
    end
    for (name, role) in diagram.roles
        sname = String(name)
        sname in knames && continue
        haskey(diagram.lc.words, sname) || throw(KodaInputError(
            :validate_diagram, "alpha/beta role names a missing component",
            Dict{String,Any}("curve" => sname)))
        role in (:alpha, :beta) || throw(KodaInputError(
            :validate_diagram, "invalid role",
            Dict{String,Any}("curve" => sname, "role" => String(role))))
    end
    nothing
end

export validate_diagram
