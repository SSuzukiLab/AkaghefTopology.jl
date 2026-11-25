# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.2–§6.4 — omitted selections mean every finite candidate
# cf: lib/topology.jl/docs/specs/module-structure.md §3.1 — movable(x) and movable(x, kind) with omitted choices
# q is the one-step generator ±1. A larger |q| is a composite of these steps, and q = 0 is the identity.
# Role-free K curves stay out of acted_components: the annulus applicator rejects a reference curve that is not the core.
# Role-bearing K curves are included. The existing preparation separates them and then twists the copies.

const _HC_TWIST_STEPS = (1, -1)

function _hc_nonempty_curves(lc)
    [n for n in _component_names(lc) if !isempty(get(lc.words, n, Int[]))]
end

function _hc_twist_acted(input)
    lc = _lc_of(input)
    names = String[]
    for n in _hc_nonempty_curves(lc)
        if _is_k(lc, n)
            input isa HeegaardCode && haskey(input.roles, n) || continue
        end
        push!(names, n)
    end
    names
end

function _enumerate_catch(build)
    out = LCMoveWitness[]
    try
        append!(out, build())
    catch error_value
        error_value isa LCMoveError || rethrow()
    end
    out
end

function _enumerate_affine(input, row, selections)
    lc = _lc_of(input)
    supports = haskey(selections, "support") ? [String(selections["support"])] : _hc_nonempty_curves(lc)
    steps = haskey(selections, "q") ? [Int(selections["q"])] : collect(_HC_TWIST_STEPS)
    acted = haskey(selections, "acted_components") ? selections["acted_components"] : _hc_twist_acted(input)
    isempty(acted) && return LCMoveWitness[]
    out = LCMoveWitness[]
    for support in supports, q in steps
        # q = 0 is absent from the omitted set. An explicit q = 0 still reaches the applicator, which returns the identity.
        next = Dict{String,Any}(
            "q" => q,
            "support" => String(support),
            "acted_components" => String[string(name) for name in acted],
        )
        for (key, value) in selections
            key in ("q", "support", "acted_components") && continue
            next[key] = value
        end
        append!(out, _enumerate_catch(() -> _expand_affine(input, row, next)))
    end
    out
end

function _enumerate_glue(input, row, selections)
    input isa HeegaardCode || return LCMoveWitness[]
    families = Dict{Symbol,Vector{String}}()
    for (name, role) in input.roles
        role in (:alpha, :beta) || continue
        push!(get!(families, role, String[]), String(name))
    end
    out = LCMoveWitness[]
    for members in values(families)
        sort!(members)
        for left in members, right in members
            left == right && continue
            next = Dict{String,Any}("curves" => [left, right], "band" => Any[])
            for (key, value) in selections
                key in ("curves", "band") && continue
                next[key] = value
            end
            append!(out, _enumerate_catch(() -> _expand_glue(input, row, next)))
        end
    end
    out
end

function _destabilize_disk_name(hc)
    names = sort!(collect(keys(hc.disks)))
    length(names) == 1 && return names[1]
    for name in names
        name == "joining" && return name
    end
    nothing
end

# A KodaInputError here means that pair is not a destabilization. Count it
# and keep the other pairs. The first record is the diagram the walk writes out.
const DESTABILIZE_SKIPS = Ref(0)
const DESTABILIZE_SKIP_FIRST = Ref{Any}(nothing)

function destabilize_skip_reset!()
    DESTABILIZE_SKIPS[] = 0
    DESTABILIZE_SKIP_FIRST[] = nothing
    nothing
end

function _jsonish(value)
    value isa AbstractDict && return Dict{String,Any}(string(k) => _jsonish(v) for (k, v) in pairs(value))
    value isa AbstractVector && return Any[_jsonish(v) for v in value]
    value isa Integer && !(value isa Bool) && return Int(value)
    value isa AbstractString && return String(value)
    value isa Symbol && return String(value)
    value isa Bool && return value
    value isa Nothing && return nothing
    string(value)
end

function _note_destabilize_skip!(input, err; alpha=0, beta=0, basis=:none, disk="", site="pair", row_id="hc.stabilization.inverse")
    err isa KodaInputError || return
    DESTABILIZE_SKIPS[] += 1
    DESTABILIZE_SKIP_FIRST[] === nothing || return
    lc = ""
    hc = ""
    try
        lc = lc_text(input.lc)
        hc = hc_text(input; lc_filename="destab-skip-first.lc")
    catch text_err
        lc = "text-failed: " * sprint(showerror, text_err)
    end
    DESTABILIZE_SKIP_FIRST[] = Dict{String,Any}(
        "row_id" => String(row_id),
        "move" => String(row_id),
        "stage" => string(err.stage),
        "message" => String(err.message),
        "site" => String(site),
        "alpha" => Int(alpha),
        "beta" => Int(beta),
        "alpha_basis" => string(basis),
        "disk" => String(disk),
        "detail" => _jsonish(err.detail),
        "lc_text" => lc,
        "hc_text" => hc,
    )
    nothing
end

# The forward row's local pieces are an empty disk and a once-punctured torus
# whose new α and β meet once and meet no older input curve. Those words do not
# occur on a census diagram: the cancelling β still crosses other curves, and
# the surgery rewrites them. §6.4 records that broader excision as a separate
# search. Naming this row runs that search. The one-argument pool does not.
function _enumerate_destabilization(input, row)
    input isa HeegaardCode || return LCMoveWitness[]
    genus = input.lc.genus
    genus > 1 || return LCMoveWitness[]
    disk = _destabilize_disk_name(input)
    disk === nothing && return LCMoveWitness[]
    row_id = string(get(row, "id", "hc.stabilization.inverse"))
    local forgot
    try
        forgot = forget_reference(realize_hc(input))
    catch err
        err isa LCMoveError || err isa HCTriangulationError || err isa ArgumentError || err isa KodaInputError || rethrow()
        _note_destabilize_skip!(input, err; disk=disk, site="forget", row_id=row_id)
        return LCMoveWitness[]
    end
    forgot.diagram === nothing && return LCMoveWitness[]
    code = _undecorated(forgot.diagram)
    out = LCMoveWitness[]
    for alpha in eachindex(code.alpha), beta in eachindex(code.beta)
        length(intersect(code.alpha[alpha].visits, code.beta[beta].visits)) == 1 || continue
        for basis in (:l, :m)
            local reduced
            try
                reduced = hc_destabilize(input; alpha=alpha, beta=beta, alpha_basis=basis, disk_name=disk)
            catch err
                err isa LCMoveError || err isa HCTriangulationError || err isa ArgumentError || err isa KodaInputError || rethrow()
                _note_destabilize_skip!(input, err; alpha=alpha, beta=beta, basis=basis, disk=disk, site="pair", row_id=row_id)
                continue
            end
            reduced.lc.genus == genus - 1 || continue
            push!(out, _destabilization_witness(input, row, reduced.hc, alpha, beta, basis, disk))
        end
    end
    out
end

function _destabilization_witness(input, row, output, alpha, beta, basis, disk)
    names = _component_names(input.lc)
    locations = Tuple{Int,Int}[(i, 1) for i in eachindex(names)]
    before = [copy(get(input.lc.words, name, Int[])) for name in names]
    after = [Int[] for _ in names]
    _witness_from(input, row, names, locations, Int[], before, after, fill(1, length(names));
        extra=Dict{String,Any}(
            "destabilization_apply" => true,
            "destabilized" => output,
            "alpha" => alpha,
            "beta" => beta,
            "alpha_basis" => String(basis),
            "disk_name" => disk,
            "support_direction" => -1,
        ))
end

function _apply_recorded_destabilization(input, witness)
    _input_fingerprint(input) == witness.fingerprint || throw(LCMoveError(
        :witness_mismatch, "witness was not built for these words, places and roles",
        Dict("row" => witness.row_id)))
    output = witness.extra["destabilized"]
    output isa HeegaardCode || throw(LCMoveError(
        :destabilization_support, "recorded destabilization is not an HC",
        Dict("row" => witness.row_id)))
    returning = witness.row_id == "hc.stabilization.inverse"
    inverse = LCMoveWitness(
        returning ? "hc.stabilization.once" : "hc.stabilization.inverse",
        witness.row_id,
        witness.digest,
        witness.relation,
        witness.layer,
        witness.names,
        witness.locations,
        Int[],
        witness.after_frags,
        witness.before_frags,
        witness.capacity,
        _input_fingerprint(output),
        witness.label_policy,
        -witness.genus_delta,
        Dict{String,Any}(
            "destabilization_apply" => true,
            "destabilized" => input,
            "origin_row_id" => returning ? "hc.stabilization.once" : "hc.stabilization.inverse",
            "support_direction" => -get(witness.extra, "support_direction", -1),
        ),
    )
    LCMoveResult(output, inverse, Dict{String,Any}(
        "relation" => witness.relation,
        "row" => witness.row_id,
        "alpha" => get(witness.extra, "alpha", nothing),
        "beta" => get(witness.extra, "beta", nothing),
        "alpha_basis" => get(witness.extra, "alpha_basis", nothing),
        "disk_name" => get(witness.extra, "disk_name", nothing),
        "support_direction" => get(witness.extra, "support_direction", -1),
    ))
end

function _enumerate_k_inclusion(input, row, selections)
    delta = Int(get(get(row, "expansion", Dict()), "genus_delta", get(row, "genus_delta", 0)))
    delta == -1 && return _enumerate_destabilization(input, row)
    delta == 1 || return LCMoveWitness[]
    input isa HeegaardCode || return LCMoveWitness[]
    faces, _ = _lc_faces(_lc_nonempty_words(input.lc.words))
    out = LCMoveWitness[]
    for face in faces
        face.id in values(input.lc.places) && continue
        next = Dict{String,Any}("disk" => face.id, "k_identification" => "ladder_end")
        for (key, value) in selections
            key in ("disk", "k_identification", "support") && continue
            next[key] = value
        end
        append!(out, _enumerate_catch(() -> _expand_k_inclusion(input, row, next)))
    end
    out
end

function _enumerate_role_parallel(input, row, selections)
    input isa HeegaardCode || return LCMoveWitness[]
    lc = input.lc
    curves = String[]
    for name in _component_names(lc)
        _is_k(lc, name) || continue
        get(input.roles, name, nothing) in (:alpha, :beta) || continue
        word = get(lc.words, name, Int[])
        isempty(word) && continue
        length(unique(abs.(word))) == length(word) || continue
        push!(curves, name)
    end
    out = LCMoveWitness[]
    for curve in curves
        next = Dict{String,Any}("curve" => curve)
        for (key, value) in selections
            key == "curve" && continue
            next[key] = value
        end
        append!(out, _enumerate_catch(() -> _move_separate_role(input, row, next)))
    end
    out
end

function _enumerate_mark_edit(input, row, selections)
    input isa HeegaardCode || return LCMoveWitness[]
    field = String(row["expansion"]["field"])
    marks = getproperty(_marks(input), Symbol(field))
    isempty(marks) && return LCMoveWitness[]
    cells = _move_mark_cells(input.lc)
    targets = if field == "points"
        vcat(cells.faces, cells.edges, cells.empties)
    elseif field == "weights"
        vcat(cells.faces, cells.weight_sides)
    else
        [face for face in cells.faces if !(face in values(input.lc.places))]
    end
    out = LCMoveWitness[]
    for name in sort!(collect(keys(marks)))
        haskey(selections, "name") && String(selections["name"]) != String(name) && continue
        for target in targets
            haskey(selections, "target") && String(selections["target"]) != String(target) && continue
            next = Dict{String,Any}("name" => String(name), "target" => field == "disks" ? [String(target)] : String(target))
            append!(out, _enumerate_catch(() -> [_move_mark_edit(input, row, next)]))
        end
    end
    out
end

function _affine_selection_complete(selections)
    haskey(selections, "q") && haskey(selections, "support") && haskey(selections, "acted_components")
end

"""
    movable(input::Union{LickorishCode,HeegaardCode}) -> Vector{LCMoveWitness}

Collect the default finite pool of supported LC or HC move witnesses.

# Arguments
`input` must satisfy the LC/HC checks and
[move applicability requirements](../specs/heegaard-moves.md).
Catalog-loading and non-`LCMoveError` failures propagate. Row-level
`LCMoveError` failures are omitted from this aggregate pool; query a named
row to obtain its diagnostic and supply explicit choices.

# Returns
A vector of `LCMoveWitness` records from the default catalog enumeration.
HC-only rows are omitted for an LC. The pool deliberately excludes
`hc.stabilization.inverse`; query that row explicitly when required.

# Guarantees
1. Each returned record belongs to a successful named-row query in the
   implemented default selection range. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This pool is not complete for geometric moves, destabilizations or arbitrary
twist powers. Empty output can reflect unsupported rows or suppressed
row-level diagnostics. It does not prove rigidity or manifold equivalence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = H.resolve_hc(lc, Dict("l1"=>:beta); standard=true);

julia> count(w -> w.row_id == "lc.twist.annulus", H.movable(h))
4
```

See also [`Topology.HeegaardCodes.move`](@ref), [`heegaard_move_ids`](@ref).

Role: access
"""
function movable(input::Union{LickorishCode,HeegaardCode})
    out = LCMoveWitness[]
    for id in heegaard_move_ids()
        # Once-crossing excision is reached by naming hc.stabilization.inverse.
        # Putting it in this pool would change the recorded walks.
        id == "hc.stabilization.inverse" && continue
        row = heegaard_move_row(id)
        String(row["layer"]) == "HC" && !(input isa HeegaardCode) && continue
        try
            append!(out, movable(input, row))
        catch error_value
            error_value isa LCMoveError || rethrow()
        end
    end
    out
end

# One-argument movable belongs to the parent generic (MoveAPI.movable). The body calls this module's method.
const _HC_MOVABLE_ONE = movable
let child = @__MODULE__
    signature_lc = child.LickorishCode
    signature_hc = child.HeegaardCode
    one_argument = child._HC_MOVABLE_ONE
    Core.eval(parentmodule(child), :(movable(input::Union{$signature_lc,$signature_hc}) = $one_argument(input)))
end
