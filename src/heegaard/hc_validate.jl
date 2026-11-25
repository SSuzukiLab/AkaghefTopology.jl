# cf: lib/topology.jl/docs/specs/lickorish-code.md §3–§4 — 各ラベルは + と − が一度ずつ。オイラー標数 2−2g の比較は _lc_validate_surface! のまま
# cf: lib/topology.jl/docs/specs/heegaard-code.md §2.1 — 非 K 成分の役割は alpha または beta を一つ

"""
    validate_diagram(diagram::HeegaardCode) -> Bool

Check the implemented passage, role and surface consistency of a resolved HC.

# Arguments
`diagram` is a resolved HC. The checks follow [HC input validation](../specs/heegaard-code.md#4-hc-の入力検査). Missing opposite passages, invalid/missing non-K role assignments and surface inconsistencies raise `KodaInputError`; LC validation and reconstruction errors propagate.

# Returns
`true` if all implemented checks succeed; rejected input throws.

# Guarantees
1. Successful return certifies the implemented passage and non-K role checks and the LC surface check. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not certify each role family as a complete cut system, validate every mark or establish a 3-manifold. Role entries on K curves are skipped by the role-specific check; hand-built records need not meet every requirement enforced by `resolve_hc`.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §4 — distinction between input validation and geometric conversion conditions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = H.resolve_hc(lc, Dict("l1"=>:beta, "m1"=>:alpha));

julia> H.validate_diagram(h)
true
```

Role: condition
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

"""
    heegaard_family_report(hc::HeegaardCode) -> NamedTuple

Report the two families' minimal cut-system checks in a resolved HC.

# Arguments
`hc` must pass `validate_hc_input`; its exceptions propagate. The mathematical interpretation assumes the reference-marked closed surface of [HC input validation](../specs/heegaard-code.md#4-hc-の入力検査). No additional geometric reconstruction is performed here.

# Returns
A named tuple `(ok, genus, alpha, beta)`. Each family record has `curves`, `count_ok`, `self_crossings`, `family_crossings`, `rank`, `dependent` and `ok`; the rank and dependencies are over Z/2 in the reference-surface homology. `ok` is the conjunction of the family flags.

# Guarantees
1. On a valid embedded reference-marked input, a successful flag identifies both families as the minimal complete cut systems of [the conversion domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体), hence gives the closed manifold interpretation there. [proved/agent: homological independence of disjoint simple curves on a closed oriented surface]
2. Reported dependency subsets have zero sum of mod-2 homology classes; they need not enumerate every dependent subset. The input is not modified. [proved/agent: by construction]

# Not guaranteed
A failed flag does not exclude the broader Heegaard-diagram domain of the specification. The report does not certify joining disks, standardness, flow structure or manifold identity. Rank alone does not establish the topology of a family with crossings.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — complete cut systems and the manifold obtained by compression and spherical capping.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1"=>:beta, "m1"=>:alpha));

julia> r = H.heegaard_family_report(hc); (r.ok, r.genus, r.alpha.rank, r.beta.rank)
(true, 1, 1, 1)
```

Role: condition
"""
function heegaard_family_report(hc::HeegaardCode)
    validate_hc_input(hc)
    alpha=_hc_family_report(hc,:alpha);beta=_hc_family_report(hc,:beta)
    (ok=alpha.ok && beta.ok,genus=hc.lc.genus,alpha=alpha,beta=beta)
end
"""
    is_heegaard_diagram(hc::HeegaardCode) -> Bool

Return the success flag of the HC minimal cut-system report.

# Arguments
The input and exceptions are those of `heegaard_family_report`; invalid input can throw rather than return `false`.

# Returns
`heegaard_family_report(hc).ok`.

# Guarantees
1. On a valid reference-marked embedded input, `true` has the minimal complete cut-system and closed-manifold meaning of `heegaard_family_report`. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Despite its name, this predicate does not cover every diagram in [the broader specification domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体). It does not decide equivalence of manifolds or select a joining disk.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — broader diagrams versus minimal complete cut systems.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1"=>:beta, "m1"=>:alpha));

julia> H.is_heegaard_diagram(hc)
true
```

Role: condition
"""
is_heegaard_diagram(hc::HeegaardCode)=heegaard_family_report(hc).ok

export heegaard_family_report, is_heegaard_diagram
