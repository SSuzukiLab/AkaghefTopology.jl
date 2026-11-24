# cf: lib/topology.jl/docs/specs/lickorish-code.md §3–§4 — 各ラベルは + と − が一度ずつ。オイラー標数 2−2g の比較は _lc_validate_surface! のまま
# cf: lib/topology.jl/docs/specs/heegaard-code.md §2.1 — 非 K 成分の役割は alpha または beta を一つ

"""Validate a resolved Heegaard diagram.

The argument is `HeegaardCode`, the same struct as `Topology.HeegaardDiagram`.
Returns `true` when the word is consistent. Throws `KodaInputError` when a
label has no opposite passage, a non-K curve has no alpha or beta role, or
the Euler characteristic does not match the stated genus.
"""
function validate_diagram(diagram::HeegaardCode)
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

function _validate_diagram_roles(diagram::HeegaardCode)
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

# cf: threads/PJ260922/data/261114-autonomous-plan/REVIEW.md D8 — which role-assigned curve systems are Heegaard diagrams
# Subsets of 0/1 rows whose sum is zero over Z/2, found by elimination. Each
# returned subset is one dependency; an empty result means the rows are independent.
function _hc_f2_dependencies(rows::Vector{Vector{Bool}})
    basis=Tuple{Vector{Bool},Vector{Int}}[];dependencies=Vector{Int}[]
    for (index,row) in enumerate(rows)
        v=copy(row);members=[index]
        for (b,who) in basis
            pivot=findfirst(b)
            if v[pivot]
                v=xor.(v,b);members=symdiff(members,who)
            end
        end
        if any(v)
            push!(basis,(v,members))
        else
            push!(dependencies,sort!(members))
        end
    end
    dependencies
end
function _hc_family_report(hc::HeegaardCode,family::Symbol)
    g=hc.lc.genus;words=hc.lc.words;names=hc_curve_names(hc,family)
    labels=Dict(name=>abs.(words[name]) for name in names)
    self_crossings=sort!([name for name in names if length(unique(labels[name]))<length(labels[name])])
    family_crossings=Tuple{String,String}[]
    for i in eachindex(names),j in i+1:length(names)
        isempty(intersect(labels[names[i]],labels[names[j]])) || push!(family_crossings,(names[i],names[j]))
    end
    basis=vcat(["l$i" for i in 1:g],["m$i" for i in 1:g])
    parity(name,k)=name==k ? false : isodd(count(in(Set(abs.(words[k]))),labels[name]))
    rows=[Bool[parity(name,k) for k in basis] for name in names]
    dependencies=_hc_f2_dependencies(rows)
    rank=length(names)-length(dependencies)
    simple=isempty(self_crossings) && isempty(family_crossings)
    (curves=names,count_ok=length(names)==g,self_crossings=self_crossings,family_crossings=family_crossings,
     rank=rank,dependent=[names[d] for d in dependencies],
     ok=length(names)==g && simple && rank==g)
end

"""Report whether the alpha and the beta curves of `hc` each form one half of a Heegaard diagram.

`validate_diagram` accepts every role-assigned curve system. This report adds
the conditions that make a family a system of attaching curves on the closed
surface of genus g, read from the words alone:

1. the family has exactly g curves;
2. no curve crosses itself and no two curves of the family cross;
3. the classes of the curves are independent in H₁(Σ; Z/2).

For pairwise disjoint simple curves, condition 3 holds exactly when the
complement of the family is connected, which is then a sphere with 2g holes.
The rank is therefore a statement about the complement only when condition 2
holds. The class of a curve modulo 2 is read from the parities of its
crossings with l₁…l_g and m₁…m_g: these curves form a symplectic basis and
the intersection form is nondegenerate, so the parities determine the class.
No sign or orientation convention enters.

`dependent` lists sets of curves whose classes sum to zero: one curve is
separating or bounds a disk, two curves are parallel, three or more bound a
subsurface together. There is no condition on how alpha meets beta.
"""
function heegaard_family_report(hc::HeegaardCode)
    validate_hc_input(hc)
    alpha=_hc_family_report(hc,:alpha);beta=_hc_family_report(hc,:beta)
    (ok=alpha.ok && beta.ok,genus=hc.lc.genus,alpha=alpha,beta=beta)
end
"""`true` when `heegaard_family_report(hc).ok`: both families are systems of attaching curves."""
is_heegaard_diagram(hc::HeegaardCode)=heegaard_family_report(hc).ok

export heegaard_family_report, is_heegaard_diagram
