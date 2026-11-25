"""
    VirtualLink(gauss, orientation, weights, real_only, is_weighted, strand_types, weight_mode) -> VirtualLink

Create mutable storage for a virtual-link or o-graph presentation.

# Arguments
The keyword constructor checks `weight_mode` against `:integer` and `:z2`,
throwing `ArgumentError` for another mode. Positional construction converts
fields to the types listed below but performs no coherence or mode validation;
the caller must ensure the [storage coherence requirements](../specs/virtual-link-data-model.md#22-current-julia-storage-and-coherence).
Conversion errors propagate. The six-argument form selects integer mode.

# Returns
A `VirtualLink` with `gauss::Vector{Vector{Int}}`, `orientation::Vector{Int}`,
`weights::Vector{Vector{Float64}}`, `real_only::Bool`, `is_weighted::Bool`,
`strand_types::Vector{Char}` and `weight_mode::Symbol`. Field meanings, including
the distinction between passage signs, real vertex signs and virtual crossings,
are owned by [Roles of stored data](../specs/virtual-link-data-model.md#2-roles-of-stored-data).

# Guarantees
1. Keyword construction returns independent empty arrays, false `real_only`
   and `is_weighted` flags and the selected mode. [proved/agent: by construction]
2. Positional fields retain supplied values after conversion; arrays already
   of the declared types are shared without copying. [proved/agent: by construction]
3. Coherent data represent the supplied diagram presentation. Construction
   itself neither identifies a link equivalence class nor certifies an
   o-graph's represented 3-manifold. [proved/agent: by construction]

# Conventions
Use the [Gauss/real-Gauss data model](../specs/virtual-link-data-model.md) and
[o-graph terminology](../specs/ograph-terminology.md), including their separate
roles for vertex signs, virtual status and local PD frame.

# Not guaranteed
Raw positional construction and subsequent field mutation can violate every
representation constraint. Assigned weights need not solve any geometric
equation. Empty storage is not a certified MST crossing-free circuit.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-data-model.md §2 — storage roles and the representation boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); weight_mode(v)
:integer
```

See also [`set_data!`](@ref), [`real_gauss_code`](@ref), [`certify_closed_ograph`](@ref).

Role: construct
"""
mutable struct VirtualLink
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    weights::Vector{Vector{Float64}}
    real_only::Bool
    is_weighted::Bool
    strand_types::Vector{Char}
    weight_mode::Symbol
end

"""
    RealEdge(id, crossing, weight, strand, arc) -> RealEdge

Store the traversal data of one directed real arc in a virtual-link presentation.

# Arguments
`id::Int`, `crossing::NTuple{2,Int}`, `weight::Float64`, `strand::Int` and
`arc::Vector{Int}`, or values convertible to those field types. Conversion
errors propagate. The constructor has no link argument and validates no
membership, labels or weight constraint.

# Returns
`RealEdge` with the supplied row identifier, pair of signed endpoint labels,
weight, strand index and raw arc visits. For computed records the interpretation
is that of [`real_edge_table`](@ref), [Edges and weights](../specs/ograph-terminology.md#5-edges-and-weights)
and [Open strands and cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).

# Guarantees
1. Fields equal the supplied values after conversion; an input `Vector{Int}`
   for `arc` is retained without copying. [proved/agent: by construction]
2. Records returned by `real_edge_table` describe arcs of that input
   presentation, including its boundary arcs when present. They do not
   identify arcs across different presentations or classify a space.
   [proved/agent: by construction]

# Not guaranteed
The constructor does not check consistency between `crossing` and `arc`,
verify integer or spin weights, or distinguish a supplied real edge from
malformed data. Its identifiers are not closed PD edge identifiers.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-data-model.md §2.1 — real-arc and boundary-arc interpretation at the VirtualLink boundary.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> rows = real_edge_table(vl_example("s3_abalone"));

julia> copies = [RealEdge(e.id, e.crossing, e.weight, e.strand, e.arc) for e in rows];

julia> Tuple(e.crossing for e in copies)
((-1, 1), (1, -1))
```

See also [`real_edge_table`](@ref).

Role: construct
"""
struct RealEdge
    id::Int
    crossing::NTuple{2,Int}
    weight::Float64
    strand::Int
    arc::Vector{Int}
end

function _check_weight_mode(mode)
    mode in (:integer,:z2) || throw(ArgumentError("weight mode must be :integer or :z2"))
    mode
end
"""
    VirtualLink(; weight_mode=:integer) -> VirtualLink

Create empty mutable storage for a virtual-link presentation.

# Arguments
`weight_mode` must be `:integer` or `:z2`; another value throws `ArgumentError`
naming the accepted modes. See the [data model](../specs/virtual-link-data-model.md).

# Returns
A `VirtualLink` with independent empty Gauss, sign, weight and strand-type
arrays, `real_only == false`, `is_weighted == false`, and the selected mode.

# Guarantees
1. The returned object is empty representation storage ready for `set_data!`. [proved/agent: by construction]
2. Repeated construction allocates independent mutable arrays. [proved/agent: by construction]

# Not guaranteed
Empty storage is not a certified link, closed normal o-graph or MST circuit.

# Sources
cf:docs/specs/virtual-link-data-model.md — representation storage and initialization boundary.

# Examples
```jldoctest
julia> weight_mode(VirtualLink())
:integer
```

Role: construct
"""
VirtualLink(;weight_mode=:integer) = VirtualLink(Vector{Vector{Int}}(), Int[], Vector{Vector{Float64}}(), false, false, Char[], _check_weight_mode(weight_mode))
# Preserve the historical positional constructor as integral mode.
"""
    VirtualLink(g, o, w, r, i, s) -> VirtualLink

Construct integer-mode presentation storage from the six historical fields.

# Arguments
`g`, `o`, `w`, `r`, `i`, `s` supply Gauss data, vertex signs, weights,
real-only and assigned-weight flags, and strand types. They must convert to
the field types of `VirtualLink`; conversion errors propagate. The caller
must satisfy the [data model](../specs/virtual-link-data-model.md), which this
raw constructor does not check.

# Returns
A `VirtualLink` with those fields and integer coefficient mode.

# Guarantees
1. The result has the same field values as the seven-argument constructor with mode `:integer`. [proved/agent: by construction]
2. Already correctly typed arrays are shared without copying; they are not modified by construction. [proved/agent: by construction]

# Not guaranteed
No incidence, closedness, weight equation or framing is validated.

# Sources
cf:docs/specs/virtual-link-data-model.md — raw fields and their representation meaning.

# Examples
```jldoctest
julia> v = VirtualLink(); w = VirtualLink(v.gauss, v.orientation, v.weights, v.real_only, v.is_weighted, v.strand_types); weight_mode(w)
:integer
```

Role: construct
"""
VirtualLink(g,o,w,r,i,s) = VirtualLink(g,o,w,r,i,s,:integer)
"""
    weight_mode(v::VirtualLink) -> Symbol

Return the stored coefficient-mode tag of a virtual-link presentation.

# Arguments
`v` is a `VirtualLink`. No validation is performed; callers who construct or
modify fields directly are responsible for a coherent mode and weight payload.

# Returns
The stored `Symbol`, normally `:integer` or `:z2` after a validated constructor
or setter.

# Guarantees
1. The result equals `v.weight_mode`; the link is unchanged.
   [proved/agent: by construction]
2. The tag concerns storage of the presentation's
   [edge weights](../specs/ograph-terminology.md#5-edges-and-weights), not a
   certification of a framing, spin structure or weight equation.
   [proved/agent: by construction]

# Not guaranteed
A valid tag does not imply assigned or valid weights. Direct field mutation
can create unsupported tags, which this accessor returns unchanged. See the
[coefficient-mode storage policy](../specs/virtual-link-data-model.md#24-weight-assignment-and-endpoint-coordinates).

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §5 and §10 — weights attached to the presentation and their geometric scope.

# Examples
```jldoctest
julia> weight_mode(VirtualLink())
:integer
```

See also [`set_weight_mode!`](@ref), [`set_weight!`](@ref), [`validate_weight`](@ref).

Role: access
"""
weight_mode(v::VirtualLink)=v.weight_mode

"""
    set_weight_mode!(v::VirtualLink, mode::Symbol) -> VirtualLink

Set the coefficient mode of a virtual-link weight assignment.

# Arguments
`mode` must be `:integer` or `:z2`, otherwise `ArgumentError` names the allowed
modes. Changing a weighted Z2 value to integer mode throws `ArgumentError`
explaining that no canonical lift is supplied. Reducing assigned integer-mode
weights invokes the finite-integral checks of [`set_weight!`](@ref), whose
errors propagate. Coherent input storage is the caller's responsibility.

# Returns
The same `VirtualLink` object.

# Guarantees
1. An unchanged mode is a no-op. An unweighted value changes only its mode;
   a weighted value switched to Z2 receives the reductions of its weights.
   [proved/agent: by construction]
2. Gauss data and vertex signs are unchanged. This changes coefficient storage
   for the [edge weights](../specs/ograph-terminology.md#5-edges-and-weights),
   not a certified framing or spin structure. [proved/agent: by construction]

# Not guaranteed
No weight equation or spin admissibility is checked. The mode tag alone does
not establish validity, and unchanged-mode calls do not revalidate weights.
Use an explicit `set_weight!(...; mode=:integer)` assignment for a chosen lift.
The [coefficient-mode storage policy](../specs/virtual-link-data-model.md#24-weight-assignment-and-endpoint-coordinates)
defines the assignment boundary.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §5 and §10 — edge weights and the distinction between equations and geometric structures.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_weight!(v, [3,-2]); set_weight_mode!(v, :z2);

julia> (weight_mode(v), Tuple(only(v.weights)))
(:z2, (1.0, 0.0))
```

See also [`set_weight!`](@ref), [`weight_mode`](@ref), [`validate_weight`](@ref).

Role: convert
"""
function set_weight_mode!(v::VirtualLink, mode::Symbol)
    _check_weight_mode(mode)
    mode==v.weight_mode && return v
    v.is_weighted && mode==:integer && throw(ArgumentError("Z2 weights have no canonical integer lift; assign explicit integer weights with set_weight!(...;mode=:integer)"))
    v.is_weighted ? set_weight!(v,v.weights;mode=mode) : (v.weight_mode=mode; v)
end

"""
    lens_space_rgauss(s::Integer, t::Integer) -> NamedTuple

Return the real-Gauss data of the repository's parameterized lens-space presentation.

# Arguments
- `s`, `t`: lens parameters; `s >= 2` and `gcd(s,t) == 1` are checked with
  `ArgumentError` identifying the failed requirement. The caller must ensure
  all intermediate integer arithmetic and output labels are representable;
  fixed-width overflow is not checked, and conversion to `Int` can fail.
  See [example routing](../specs/virtual-link-examples.md#parameterized-unweighted-or-combed-lens-space).

# Returns
A named tuple with one real-Gauss component `rgauss::Vector{Vector{Int}}`
and its `orientation::Vector{Int}` vertex-sign vector.

# Guarantees
1. The result is the input data for `lens_space_link(s,t)`, the family named by the linked lens-space convention. [proved/agent: by construction]
2. Under the arithmetic precondition, it has `s` true vertices, all with positive vertex sign, and one component of length `2s`. [proved/agent: by construction]
3. Arguments are not modified and fresh vectors are returned. [proved/agent: by construction]

# Conventions
Use [LensSpaceNames](../conventions/LensSpaceNames.md) for oriented names and
translation to other conventions, and the
[real-Gauss data model](../specs/virtual-link-data-model.md) for coordinates.

# Not guaranteed
No integer framing is provided. The external surgery/quotient identification
retains the D1 and D2 dependencies stated in LensSpaceNames; this constructor
does not prove that identification for arbitrary parameters.

# Sources
cf:docs/specs/virtual-link-examples.md — parameterized lens-space family and framing boundary.
cf:docs/conventions/LensSpaceNames.md — naming rule and conditional comparison with surgery and quotient conventions.

# Examples
```jldoctest
julia> input = lens_space_rgauss(2, 1); v = VirtualLink(); set_data!(v; input...);

julia> is_closed(v)
(true, (true, true, true))
```

See also [`lens_space_link`](@ref).

Role: construct
"""
function lens_space_rgauss(s::Integer,t::Integer)
    s>=2 || throw(ArgumentError("s must be at least 2"))
    gcd(s,t)==1 || throw(ArgumentError("s and t must be coprime"))
    code=vcat(-mod.((s:-1:1).*t.-1,s).-1,1:s)
    (rgauss=[Int.(code)],orientation=ones(Int,s))
end

"""
    lens_space_link(s::Integer, t::Integer) -> VirtualLink

Construct the repository's unweighted normal o-graph presentation named L(s,t).

# Arguments
- `s`, `t`: parameters accepted by [`lens_space_rgauss`](@ref); its checked
  lower-bound and coprimality failures throw `ArgumentError`. Its caller
  precondition on representable integer arithmetic also applies here.

# Returns
A fresh unweighted integer-mode `VirtualLink` with the returned family's
real-Gauss data and vertex signs.

# Guarantees
1. `real_gauss_code(result)` agrees with `lens_space_rgauss(s,t)`. [proved/agent: by construction]
2. The presentation is the lens-space family of the linked example specification; its oriented name is fixed by LensSpaceNames. [proved/agent: by construction]
3. Arguments are not modified. [proved/agent: by construction]

# Conventions
Use [LensSpaceNames](../conventions/LensSpaceNames.md) for the orientation
and its dictionary of alternative lens-space names.

# Not guaranteed
No framing is assigned or C1–C3 witness returned. Identification with an
externally oriented surgery or quotient retains the convention page's D1 and
D2 dependencies; the constructor does not establish a general comparison proof.

# Sources
cf:docs/specs/virtual-link-examples.md — unweighted lens-space entrypoint.
cf:docs/conventions/LensSpaceNames.md — repository name and orientation comparison.

# Examples
```jldoctest
julia> is_closed(lens_space_link(2, 1))
(true, (true, true, true))
```

See also [`certify_closed_ograph`](@ref), [`framed_lens_space_link`](@ref).

Role: construct
"""
function lens_space_link(s::Integer,t::Integer)
    input=lens_space_rgauss(s,t)
    link=VirtualLink()
    set_data!(link;rgauss=input.rgauss,orientation=input.orientation)
end

function _normalize_gauss(code)
    code isa AbstractVector{<:Integer} && return [Int.(code)]
    [Int.(c) for c in code]
end

_is_cut_component(component::AbstractVector{<:Integer}) =
    !isempty(component) && component[end] == 0

function _component_visits(component::AbstractVector{<:Integer})
    _is_cut_component(component) ? Int.(component[begin:end-1]) : Int.(component)
end

function _validate_cut_markers(code::Vector{Vector{Int}})
    for (strand,component) in enumerate(code)
        length(component) > 1 && any(==(0),@view component[1:end-1]) &&
            throw(ArgumentError("cut marker 0 is allowed only at the end of strand $strand"))
    end
    true
end

function _validate_cut_flags(flags, strand_count::Integer)
    flags isa AbstractVector || throw(ArgumentError("cut flags must be a vector of Bool"))
    all(flag -> flag isa Bool, flags) ||
        throw(ArgumentError("cut flags must contain actual Bool values"))
    length(flags) == strand_count ||
        throw(DimensionMismatch("one cut flag per Gauss strand is required"))
    Bool[flags...]
end

"""
    strand_cuts(v::VirtualLink) -> Vector{Bool}

Return the open-strand flags of a link presentation.

# Arguments
`v` must have valid [cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
Malformed cut markers throw `ArgumentError`. The remaining link data are not
validated by this query.

# Returns
A fresh Boolean vector in strand order, with `true` for an open strand and
`false` for a closed strand.

# Guarantees
1. The flags record exactly the boundary status in the linked data model,
   with one entry per stored strand. [proved/agent: by construction]
2. The link is unchanged and the result can be modified independently.
   [proved/agent: by construction]

# Not guaranteed
Closed strands alone do not certify a closed normal o-graph or a closed
3-manifold. Crossing coherence and framing admissibility are not tested.

# Sources
cf:docs/specs/virtual-link-data-model.md#21-open-strands-and-cut-markers — open and closed strand status.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); set_cuts!(v, [true]); Tuple(strand_cuts(v))
(true,)
```

See also [`set_cuts!`](@ref), [`is_closed`](@ref).

Role: access
"""
function strand_cuts(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    Bool[_is_cut_component(component) for component in v.gauss]
end

function _strand_arc_count(v::VirtualLink, strand::Integer)
    visits = _component_visits(v.gauss[strand])
    real_count = count(entry -> entry != 0 &&
                              abs(entry) <= length(v.orientation) &&
                              v.orientation[abs(entry)] != 0, visits)
    real_count + (_is_cut_component(v.gauss[strand]) ? 1 : 0)
end

"""
    set_cuts!(v::VirtualLink, flags) -> VirtualLink

Set the open or closed boundary status of each strand in a link presentation.

# Arguments
- `v`: a coherent presentation with valid
  [cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
  Invalid markers throw `ArgumentError`. Changing flags on a weighted link
  throws `ArgumentError`; requesting its existing flags is accepted.
- `flags`: a vector of actual Boolean values, one per strand. Invalid
  container/elements throw `ArgumentError`; wrong length throws `DimensionMismatch`.

# Returns
The same `VirtualLink` with the requested boundary status.

# Guarantees
1. On successful return `strand_cuts(v) == flags`, and the nonboundary Gauss visits and vertex signs are retained. [proved/agent: by construction]
2. When flags change, each weight array has the new real-arc length, retaining the available old prefix and filling additional entries with zeros. [proved/agent: by construction]
3. Changed flags modify `v` in place; identical flags return it unchanged. [proved/agent: by construction]

# Not guaranteed
Changing a boundary is not an equivalence move or a framing transport.
Retained unassigned weight entries are storage values, not a solution of
the new relative boundary equations.

# Sources
cf:docs/specs/virtual-link-data-model.md#21-open-strands-and-cut-markers — strand boundaries and real-arc counts.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); set_cuts!(v, [true]);

julia> last(only(reverse_strands(v).gauss))
0
```

Role: move
"""
function set_cuts!(v::VirtualLink, flags)
    new_flags = _validate_cut_flags(flags, length(v.gauss))
    _validate_cut_markers(v.gauss)
    old_flags = Bool[_is_cut_component(component) for component in v.gauss]
    old_flags == new_flags && return v
    v.is_weighted &&
        throw(ArgumentError("changing cuts on a weighted link requires clearing weights and assigning explicit open-arc weights"))
    old_weights = deepcopy(v.weights)
    for (strand,cut) in enumerate(new_flags)
        component = v.gauss[strand]
        currently_cut = _is_cut_component(component)
        if cut && !currently_cut
            push!(component, 0)
        elseif !cut && currently_cut
            pop!(component)
        end
    end
    v.weights = [begin
        count_ = _strand_arc_count(v, strand)
        values = zeros(Float64, count_)
        old = strand <= length(old_weights) ? old_weights[strand] : Float64[]
        copy_count = min(length(old), count_)
        copy_count > 0 && (values[1:copy_count] = Float64.(old[1:copy_count]))
        values
    end for strand in eachindex(v.gauss)]
    v
end

function _validate_gauss(code::Vector{Vector{Int}}, ori::Vector{Int}; real_only=false)
    _validate_cut_markers(code)
    labels = filter(!=(0), reduce(vcat, code; init=Int[]))
    isempty(labels) && return true
    maximum(abs.(labels)) <= length(ori) || throw(ArgumentError("orientation is missing a crossing"))
    for i in 1:length(ori)
        count(==(i), abs.(labels)) == 2 || throw(ArgumentError("crossing $i must occur exactly twice"))
        real_only && ori[i] == 0 && throw(ArgumentError("real Gauss code cannot contain virtual orientation"))
    end
    true
end

"""
    set_data!(v::VirtualLink; gauss=nothing, rgauss=nothing, orientation=nothing, odata=nothing, pd=nothing, table=nothing, cut=Bool[], strand_types=nothing) -> VirtualLink

Replace a link's stored presentation from one supported input representation.

# Arguments
- The caller supplies one of `gauss`, `rgauss`, `odata`, `pd` or `table` in
  the [input formats](../specs/virtual-link-data-model.md#31-broad-ingress).
  Supplying no representation throws `ArgumentError`; exclusivity is not checked.
  Multiple representations follow the specified
  [ingress precedence](../specs/virtual-link-data-model.md#23-ingress-precedence).
- Gauss inputs require `orientation`; missing signs, invalid crossing counts,
  missing crossing labels, malformed cut markers or virtual status in real-only
  data throw `ArgumentError`. Entries must convert to `Int`. Full mathematical
  incidence and valid sign values remain the caller's responsibility.
- Nonempty PD inputs require one sign per crossing when signs are supplied;
  mismatch throws `DimensionMismatch`, and invalid PD sign values throw
  `ArgumentError`. The caller supplies coherent PD incidence.
- `table` supports the trefoil table entry only; others throw `ArgumentError`.
- `cut`: optional Boolean flags checked by [`set_cuts!`](@ref).
- `strand_types`: optional values convertible to the strand-type field;
  default o-graph types are assigned. The caller supplies a coherent count.
- Three-part `odata` supplies weights accepted by [`set_weight!`](@ref);
  its domain and shape errors propagate. The existing coefficient mode is used.

# Returns
The same `VirtualLink`, containing the selected representation's data.

# Guarantees
1. A successful single-representation call stores the selected input's Gauss presentation with the requested boundary and strand-type data. [proved/agent: by construction]
2. Three-part O-data assigns its supplied weights; other forms clear assigned-weight status and create zero-valued real-arc storage. [proved/agent: by construction]
3. The coefficient mode is retained and `v` is modified in place. [proved/agent: by construction]

# Conventions
Use the [Gauss and PD data model](../specs/virtual-link-data-model.md) and
[PD pipeline](../../../PLOT_PIPELINE.md) for representation-local edge identities.

# Not guaranteed
PD edge IDs need not round-trip; attach edge-indexed auxiliary data to the
result of `pd_code(v)`. This is not a link-isotopy, normal closedness or framing
certificate. Failure can leave `v` partly changed; replacement is not transactional.

# Sources
cf:docs/specs/virtual-link-data-model.md — input formats, stored data roles and representation-local identity.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; gauss=[[1,-2],[-1,2]], orientation=[1,1]);

julia> size(pd_code(v))
(2, 4)
```

Role: convert
"""
function set_data!(v::VirtualLink; gauss=nothing, rgauss=nothing, orientation=nothing,
                   odata=nothing, pd=nothing, table=nothing, cut=Bool[], strand_types=nothing)
    supplied_weights = odata !== nothing && length(odata) == 3
    if odata !== nothing
        gauss, orientation = odata[1], odata[2]
        length(odata) == 3 && (v.weights = [Float64.(w) for w in odata[3]]; v.is_weighted = true)
        rgauss = gauss; gauss = nothing
    elseif table !== nothing
        tuple(table...) == (3,1) || throw(ArgumentError("only the trefoil table example [3,1] is defined without SageMath"))
        pd = [1 2 4 3; 3 4 6 5; 5 6 2 1]
        orientation = [-1,-1,-1]
    end
    if pd !== nothing
        if isempty(pd)
            v.gauss = [Int[]]; v.orientation = Int[]; v.real_only = false
        else
            v.orientation = orientation === nothing ? fill(1, size(pd,1)) : Int.(orientation)
            length(v.orientation) == size(pd,1) ||
                throw(DimensionMismatch("one orientation value per PD crossing is required"))
            all(sign -> sign in (-1,0,1), v.orientation) ||
                throw(ArgumentError("PD orientation values must be -1, 0, or 1"))
            v.gauss = _pd_to_gauss(Int.(pd); orientation=v.orientation)
            v.real_only = false
        end
    elseif gauss !== nothing || rgauss !== nothing
        code = _normalize_gauss(gauss === nothing ? rgauss : gauss)
        ori = orientation === nothing ? throw(ArgumentError("orientation is required")) : Int.(orientation)
        _validate_gauss(code, ori; real_only=rgauss !== nothing)
        v.gauss, v.orientation, v.real_only = code, ori, rgauss !== nothing
    else
        throw(ArgumentError("one input representation is required"))
    end
    v.is_weighted = false
    isempty(cut) || set_cuts!(v, cut)
    v.strand_types = strand_types === nothing ? fill('o', length(v.gauss)) : collect(strand_types)
    if supplied_weights
        v.is_weighted = false
        set_weight!(v,odata[3])
    else
        v.is_weighted = false
        v.weights = [zeros(Float64,_strand_arc_count(v,strand))
                     for strand in eachindex(v.gauss)]
    end
    v
end

function _pd_to_gauss(pd::Matrix{Int}; orientation=nothing)
    # PD edges form oriented components. Crossing sign is encoded by under/over slots.
    signs = orientation === nothing ? fill(1, size(pd,1)) : Int.(orientation)
    length(signs) == size(pd,1) ||
        throw(DimensionMismatch("one orientation value per PD crossing is required"))
    next = Dict{Int,Tuple{Int,Int}}()
    for (i,row) in enumerate(eachrow(pd))
        next[row[1]] = (row[3], -i)
        # The over-strand transition follows the local PD frame.  Positive
        # (and virtual) rows enter at slot 4 and leave at slot 2; negative
        # rows enter at slot 2 and leave at slot 4.  Keeping this directed
        # incidence avoids self-crossing edge IDs overwriting the transition.
        if signs[i] == -1
            next[row[2]] = (row[4], i)
        else
            next[row[4]] = (row[2], i)
        end
    end
    result = Vector{Vector{Int}}(); unseen = Set(keys(next))
    outgoing = Set(first(value) for value in values(next))
    starts = sort!(collect(setdiff(Set(keys(next)), outgoing)))
    function walk!(start, open)
        edge = start; code = Int[]
        while edge in unseen
            delete!(unseen, edge)
            edge, crossing = next[edge]
            push!(code, crossing)
        end
        open && push!(code, 0)
        push!(result, code)
    end
    for start in starts
        walk!(start, true)
    end
    while !isempty(unseen)
        walk!(minimum(unseen), false)
    end
    result
end

"""
    gauss_code(v::VirtualLink) -> Vector{Vector{Int}}

Return a detached copy of the stored Gauss components of a virtual-link presentation.

# Arguments
`v` must have valid [cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
A malformed marker throws `ArgumentError` identifying its placement.
Other fields are not validated by this accessor.

# Returns
A fresh `Vector{Vector{Int}}` containing the stored component words.

# Guarantees
1. The result equals `v.gauss`, including its order, virtual visits and cut
   markers. Neither the outer vector nor its component vectors alias `v.gauss`.
   [proved/agent: by construction]
2. This is the stored [presentation view](../specs/virtual-link-data-model.md#2-roles-of-stored-data),
   not a canonical representative of a virtual-link equivalence class.
   [proved/agent: by construction]
3. The input is unchanged. Mutating the returned vectors does not mutate it.
   [proved/agent: by construction]

# Not guaranteed
The returned words do not include vertex signs, weights or strand types and
are not a complete `VirtualLink` by themselves. Use `real_gauss_code` to
project out virtual crossings. This call checks no realizability or closedness.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-data-model.md §2 — the stored Gauss view and cut-marker retention.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_cuts!(v, [true]);

julia> Tuple(only(gauss_code(v)))
(1, -1, 0)
```

See also [`real_gauss_code`](@ref), [`ograph_state`](@ref).

Role: access
"""
function gauss_code(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    deepcopy(v.gauss)
end
"""
    real_gauss_code(v::VirtualLink) -> Tuple{Vector{Vector{Int}},Vector{Int}}

Return the real-crossing Gauss presentation and its vertex signs.

# Arguments
`v` must have coherent data under the
[VirtualLink data model](../specs/virtual-link-data-model.md).
Malformed cut markers throw `ArgumentError`. Complete consistency of visits
with the orientation array is the caller's responsibility; absent labels can
be omitted instead of diagnosed.

# Returns
`(code, signs)`: newly allocated component words and the corresponding real
vertex signs in the [real-Gauss view](../specs/virtual-link-data-model.md#2-roles-of-stored-data).

# Guarantees
1. Virtual visits are omitted; real crossing labels are consecutive in
   increasing original-label order. Component order and surviving traversal
   order are retained, including a trailing marker on each cut component.
   [proved/agent: by construction]
2. The result represents the real-crossing projection of the input diagram,
   with its passage signs. The input and its weights are unchanged.
   [proved/agent: by construction]

# Not guaranteed
The projection does not encode virtual-crossing locations, weights or strand
types, and is not invertible to the full stored diagram. It neither certifies
closed-normality nor canonicalizes cyclic roots or manifold equivalence.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-data-model.md §2 — raw Gauss storage and the real-Gauss view.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> code, signs = real_gauss_code(vl_example("s3_abalone")); (Tuple(only(code)), Tuple(signs))
((1, -1), (1,))
```

See also [`gauss_code`](@ref), [`ograph_state`](@ref).

Role: convert
"""
function real_gauss_code(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    real = findall(!=(0), v.orientation); mapping = Dict(x=>i for (i,x) in enumerate(real))
    code = Vector{Vector{Int}}()
    for component in v.gauss
        visits = _component_visits(component)
        filtered = [sign(x)*mapping[abs(x)] for x in visits
                    if x != 0 && haskey(mapping, abs(x))]
        _is_cut_component(component) && push!(filtered, 0)
        push!(code, filtered)
    end
    code, v.orientation[real]
end

"""
    real_edge_table(v::VirtualLink) -> Vector{RealEdge}

Return the directed real arcs of a virtual-link presentation in traversal order.

# Arguments
`v` must contain coherent raw Gauss data under the
[data model](../specs/virtual-link-data-model.md). Malformed cut markers are
checked with `ArgumentError`; other coherence is the caller's responsibility.
Weights may be assigned or unassigned and are not validated here.

# Returns
`Vector{RealEdge}` with consecutive `id`, signed endpoint pair `crossing`,
stored `weight::Float64`, strand index `strand`, and raw-visit vector `arc`.
Endpoint labels use the real-Gauss projection; arc labels retain raw visits.
Boundary endpoints use the linked cut convention. An absent weight entry
is represented by `NaN`.

# Guarantees
1. Every cut component contributes one more arc than its number of real visits; a closed component with real visits contributes one per visit. [proved/agent: by construction]
2. The records describe the real arcs of the supplied presentation, including virtual visits between real endpoints. [proved/agent: by construction]
3. The input is not modified and the returned arc vectors are independent. [proved/agent: by construction]

# Conventions
Use [open strands and cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers)
and [Gauss presentation coordinates](../specs/ograph-terminology.md).

# Not guaranteed
A closed component without real visits contributes no record. The table does
not certify a framing, a normal o-graph, or any represented manifold; stored
zero weights are returned even when `is_weighted` is false.

# Sources
cf:docs/specs/virtual-link-data-model.md — raw presentation and real-arc boundary semantics.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); set_cuts!(v, [true]); edges = real_edge_table(v);

julia> (first(edges).crossing[1], last(edges).crossing[2])
(0, 0)
```

Role: access
"""
function real_edge_table(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    real_indices=findall(!=(0),v.orientation)
    mapping=Dict(index=>mapped for (mapped,index) in enumerate(real_indices))
    result=RealEdge[]
    edge_id=0
    for (strand,component0) in enumerate(v.gauss)
        component=_component_visits(component0)
        cut=_is_cut_component(component0)
        real_positions=findall(position -> haskey(mapping,abs(component[position])),eachindex(component))
        weights=strand <= length(v.weights) ? v.weights[strand] : Float64[]
        convert_crossing(x)=x == 0 ? 0 : sign(x)*mapping[abs(x)]
        push_edge!(crossing,arc,weight_index)=begin
            edge_id += 1
            weight=weight_index <= length(weights) ? weights[weight_index] : NaN
            push!(result,RealEdge(edge_id,crossing,weight,strand,Int.(arc)))
        end
        if cut
            if isempty(real_positions)
                push_edge!((0,0),component,1)
            else
                first_position=first(real_positions)
                push_edge!((0,convert_crossing(component[first_position])),
                           component[1:first_position],1)
                for weight_index in 2:length(real_positions)
                    previous_position=real_positions[weight_index-1]
                    target_position=real_positions[weight_index]
                    push_edge!((convert_crossing(component[previous_position]),
                                convert_crossing(component[target_position])),
                               component[previous_position:target_position],weight_index)
                end
                last_position=last(real_positions)
                push_edge!((convert_crossing(component[last_position]),0),
                           component[last_position:end],length(real_positions)+1)
            end
        elseif !isempty(real_positions)
            for (weight_index,target_position) in enumerate(real_positions)
                previous_position=real_positions[mod1(weight_index-1,length(real_positions))]
                arc=Int[component[previous_position]]
                cursor=previous_position
                while cursor != target_position
                    cursor=mod1(cursor+1,length(component))
                    push!(arc,component[cursor])
                end
                crossing=(convert_crossing(first(arc)),convert_crossing(last(arc)))
                push_edge!(crossing,arc,weight_index)
            end
        end
    end
    result
end

"""
    head_map(v::VirtualLink) -> Matrix{Int}

Return the signed successor visits of the stored Gauss presentation.

# Arguments
`v` must contain coherent raw crossing labels and valid
[cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
Malformed markers throw `ArgumentError`; a visit outside the orientation
array can throw `BoundsError`. Repeated or missing signed visits are not
independently validated.

# Returns
A `Matrix{Int}` in the specified
[`head_map` coordinates](../specs/virtual-link-data-model.md#22-current-julia-storage-and-coherence).
The raw view includes virtual crossings.

# Guarantees
1. The result has size `(length(v.orientation), 2)` and does not alias input
   storage. Its nonzero successors come from the input component order.
   [proved/agent: by construction]
2. A closed component includes its cyclic last-to-first successor. A cut
   component's final visit has no successor through its boundary.
   [proved/agent: by construction]
3. This is adjacency data of the supplied
   [Gauss presentation](../specs/ograph-terminology.md#2-standard-gauss-code-displays),
   not a graph or manifold classification. The input is unchanged.
   [proved/agent: by construction]

# Not guaranteed
The call does not check that each signed visit occurs exactly once, identify
which zeros arose from missing visits, or remove virtual crossings. It is
not the real-edge table and does not return a closedness certificate.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-data-model.md §2.1 — raw component order and cut boundaries.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_cuts!(v, [true]);

julia> (head_map(v)[1, 2], head_map(v)[1, 1])
(-1, 0)
```

See also [`real_edge_table`](@ref), [`pd_code`](@ref).

Role: access
"""
function head_map(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    n = length(v.orientation); hm = fill(0, n, 2)
    for component in v.gauss
        clean = filter(!=(0), component); isempty(clean) && continue
        closed = isempty(component) || last(component) != 0
        stop = closed ? length(clean) : length(clean)-1
        for i in 1:max(stop,0)
            current = clean[i]; successor = clean[mod1(i+1, length(clean))]
            hm[abs(current), current > 0 ? 2 : 1] = successor
        end
    end
    hm
end

"""
    pd_code(v::VirtualLink) -> Matrix{Int}

Return the four-slot PD presentation of the stored virtual-link diagram.

# Arguments
`v` must be coherent under the
[data model](../specs/virtual-link-data-model.md) and its
[cut-marker convention](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
Cut-marker errors throw `ArgumentError`; missing signed crossing visits can
raise `KeyError`. Complete crossing-incidence validation is the caller's
responsibility.

# Returns
A four-column integer matrix, one row per stored crossing when visits are
present, or a `0 × 4` matrix for an empty presentation. Entries are edge IDs
in the [PD port convention](../specs/ograph-terminology.md#6-local-ports).

# Guarantees
1. Closed component ends are joined and cut component ends remain distinct.
   Virtual crossings have rows as well as real crossings. The input is unchanged.
   [proved/agent: by construction]
2. The matrix encodes incidence and local cyclic order of the supplied
   diagram; it adds no assertion of planar realizability or a represented
   closed 3-manifold. [proved/agent: by construction]

# Conventions
Use the [PD and Gauss conventions](../specs/virtual-link-data-model.md) and
[local ports](../specs/ograph-terminology.md#6-local-ports).

# Not guaranteed
Virtual status is not recoverable from the PD matrix alone; retain the
orientation vector. Weights, strand labels and empty crossing-free components
are not encoded. The result is not canonical under diagram moves.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §6 — the four-port PD interpretation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = VirtualLink(); set_data!(v; gauss=[[1,-2],[-1,2]], orientation=[1,1]);

julia> size(pd_code(v))
(2, 4)
```

See also [`real_gauss_code`](@ref), [`edge_directions`](@ref).

Role: convert
"""
function pd_code(v::VirtualLink)
    # Sage 10.1+ convention: start at the entering under-strand and walk
    # counter-clockwise. This is the convention used by the MATLAB oracle.
    _validate_cut_markers(v.gauss)
    code = Vector{Vector{Int}}()
    cuts = Bool[]
    for component in v.gauss
        visits = _component_visits(component)
        isempty(visits) && continue
        push!(code,visits); push!(cuts,_is_cut_component(component))
    end
    isempty(code) && return zeros(Int,0,4)

    edge_at_visit = Dict{Int,NTuple{2,Int}}()
    edge = 0
    for (component,cut) in zip(code,cuts)
        first_edge = edge+1
        cut && (edge = first_edge)
        for crossing in component
            edge += 1
            edge_at_visit[crossing] = (edge,edge+1)
        end
        if !cut
            last_crossing = last(component)
            incoming,outgoing = edge_at_visit[last_crossing]
            edge_at_visit[last_crossing] = (incoming,first_edge)
        else
            edge += 1 # reserve the outgoing boundary edge
        end
    end

    pd = zeros(Int,length(v.orientation),4)
    for crossing in eachindex(v.orientation)
        under = edge_at_visit[-crossing]
        over = edge_at_visit[crossing]
        if v.orientation[crossing] == -1
            pd[crossing,:] = [under[1],over[1],under[2],over[2]]
        else
            # `0` is topologically virtual, but a PD row still needs a cyclic
            # slot order. Sage/MATLAB use the +1 local frame for that purely
            # combinatorial choice; virtual status remains in `orientation`.
            pd[crossing,:] = [under[1],over[2],under[2],over[1]]
        end
    end
    pd
end

"""Return edge => crossing-index maps for tails and heads of a PD code."""
function _pd_has_boundary(pd::AbstractMatrix{<:Integer})
    isempty(pd) && return false
    values=vec(Int.(pd))
    any(count(==(edge),values)==1 for edge in unique(values))
end

"""
    edge_directions(pd::AbstractMatrix{<:Integer}) -> Tuple{Dict{Int,Int},Dict{Int,Int}}

Return endpoint directions for the edges of a closed PD presentation.

# Arguments
`pd` must be a closed four-slot PD matrix in the
[repository port convention](../specs/ograph-terminology.md#6-local-ports),
with integer edge labels convertible to `Int`. An edge label occurring only
once is rejected with `ArgumentError` identifying boundary IDs. Other malformed
matrices may raise indexing, conversion or `ArgumentError` failures instead;
this is not a complete PD validator.

# Returns
`(tails, heads)`, two `Dict{Int,Int}` values mapping edge labels to the row
indices of their tail and head crossings. A loop can have equal endpoint rows.

# Guarantees
1. For a valid closed PD input, every edge label has a tail and head, both
   incident to that edge. The input matrix is unchanged. [proved/agent: by construction]
2. The directions concern the strands of the supplied PD presentation, with
   opposite slots continued according to the linked port convention.
   [proved/agent: by construction]
3. An empty matrix returns two empty dictionaries. [proved/agent: by construction]

# Not guaranteed
PD data alone need not fix an orientation for a component with no seeded
undercrossing. No canonical choice across reordered PD presentations is
promised for such components. The result does not establish classical
realizability, a spine, or a 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §6 — the Sage/MATLAB four-slot PD port convention.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); tails, heads = edge_directions(pd_code(v));

julia> Tuple((tails[e.id], heads[e.id]) for e in real_edge_table(v))
((1, 1), (1, 1))
```

See also [`pd_code`](@ref), [`real_edge_table`](@ref).

Role: convert
"""
function edge_directions(pd::AbstractMatrix{<:Integer})
    _pd_has_boundary(pd) &&
        throw(ArgumentError("edge directions require a closed PD code; boundary edge IDs are present"))
    tails = Dict{Int,Int}()
    heads = Dict{Int,Int}()
    rows = [Int.(collect(row)) for row in eachrow(pd)]
    for crossing_index in eachindex(rows)
        crossing = rows[crossing_index]
        tails[crossing[3]] = crossing_index
        edge = crossing[3]
        current = crossing_index
        while !haskey(heads,edge)
            next_crossings = [i for i in eachindex(rows) if i != current && edge in rows[i]]
            if isempty(next_crossings)
                heads[edge] = current
                tails[edge] = current
                current_row = rows[current]
                edge = current_row[1] == edge ? current_row[3] :
                       current_row[4] == edge ? current_row[2] : current_row[4]
            else
                next_crossing = first(next_crossings)
                heads[edge] = next_crossing
                tails[edge] = current
                current = next_crossing
                current_row = rows[current]
                slot = findfirst(==(edge),current_row)
                edge = current_row[mod1(slot+2,4)]
            end
        end
    end

    unassigned = setdiff(Set(vec(Int.(pd))),Set(keys(tails)))
    while !isempty(unassigned)
        edge = pop!(unassigned)
        current = findfirst(row -> edge in row,rows)
        while !haskey(heads,edge)
            tails[edge] = current
            next_crossing = findfirst(i -> i != current && edge in rows[i],eachindex(rows))
            next_crossing === nothing && throw(ArgumentError("invalid PD edge $edge"))
            heads[edge] = next_crossing
            current = next_crossing
            current_row = rows[current]
            slot = findfirst(==(edge),current_row)
            edge = current_row[mod1(slot+2,4)]
            delete!(unassigned,edge)
        end
    end
    tails,heads
end

"""
    regions(pd::AbstractMatrix{<:Integer}) -> Vector{Vector{Int}}

Return the signed face-boundary walks of a closed PD rotation system.

# Arguments
`pd` must be coherent closed PD data in the
[PD convention](../specs/virtual-link-data-model.md).
Boundary edge IDs are checked and rejected with `ArgumentError`. The caller
supplies four entries per crossing and valid paired edge occurrences;
invalid data can raise `ArgumentError`, indexing or key errors, and are not
covered by a total validation guarantee. Integers must convert to `Int`.

# Returns
Vectors of signed edge labels, one per face-boundary walk. Empty PD data
give an empty vector.

# Guarantees
1. For valid input the walks are the left-turning face boundaries of the PD rotation system, with edge signs indicating traversal relative to its directed edges. [proved/agent: by construction]
2. Each directed edge side occurs in the returned boundaries; loops and repeated edge incidence are retained. [proved/agent: by construction]
3. `pd` is not modified. [proved/agent: by construction]

# Conventions
Use [rotation-system faces](../specs/lickorish-code.md)
and the repository PD frame. The historical planar ordering is used by the
[plot pipeline](../../../PLOT_PIPELINE.md).

# Not guaranteed
The function does not establish planarity. On a positive-genus rotation
system these walks are not regions of a sphere diagram. No general equivalence
of ordering with every Sage input is asserted.

# Sources
cf:docs/specs/lickorish-code.md §4 — faces of the PD rotation system.

# Examples
```jldoctest
julia> Tuple(Tuple(r) for r in regions([3 2 4 1; 2 3 1 4]))
((4, 1), (3, -1), (2, -4), (-2, -3))
```

Role: access
"""
function regions(pd::AbstractMatrix{<:Integer})
    _pd_has_boundary(pd) &&
        throw(ArgumentError("regions require a closed PD code; boundary edge IDs are present"))
    rows = [Int.(collect(row)) for row in eachrow(pd)]
    isempty(rows) && return Vector{Vector{Int}}()
    # Opposite-slot loops are not monogons (e.g. the torus [1 2 1 2]).
    # Use occurrences, not crossing indices, to distinguish their endpoints.
    # Source: lickorish-code.md section 4, rotation-system faces.
    # Keep the Sage planar path below unchanged, including its face ordering.
    if any(row -> row[1] == row[3] || row[2] == row[4], rows)
        Dart = Tuple{Int,Int}
        occurrences = Dict{Int,Vector{Dart}}()
        for i in eachindex(rows), slot in 1:4
            push!(get!(occurrences, rows[i][slot], Dart[]), (i,slot))
        end
        twin = Dict{Dart,Dart}()
        for (edge, ends) in occurrences
            length(ends) == 2 || throw(ArgumentError("invalid PD edge $edge"))
            a,b = ends
            twin[a] = b; twin[b] = a
        end
        signed = Dict{Dart,Int}()
        # Slot 3 is outgoing on the under-strand. A component containing
        # only over-strands has no such seed; choose slot 2 consistently.
        for seedslot in (3,2), i in eachindex(rows)
            dart = (i,seedslot)
            while !haskey(signed,dart)
                edge = rows[dart[1]][dart[2]]
                signed[dart] = edge
                arrival = twin[dart]
                signed[arrival] = -edge
                dart = (arrival[1],mod1(arrival[2]+2,4))
            end
        end
        bysigned = Dict(value=>dart for (dart,value) in signed)
        available = Set(keys(bysigned))
        result = Vector{Vector{Int}}()
        while !isempty(available)
            start = bysigned[maximum(available)]
            dart = start
            region = Int[]
            while true
                push!(region,signed[dart])
                delete!(available,signed[dart])
                arrival = twin[dart]
                dart = (arrival[1],mod1(arrival[2]-1,4))
                dart == start && break
            end
            push!(result,region)
        end
        return result
    end
    if length(rows) == 1
        crossing = rows[1]
        return crossing[1] == crossing[4] ?
            [[-crossing[3]],[crossing[1]],[crossing[3],-crossing[1]]] :
            [[crossing[3]],[-crossing[1]],[-crossing[3],crossing[1]]]
    end

    tails,heads = edge_directions(pd)
    unsigned = Set(vec(Int.(pd)))
    loops = sort([edge for edge in unsigned if heads[edge] == tails[edge]])
    available = union(unsigned,Set(-edge for edge in unsigned))
    result = Vector{Vector{Int}}()
    for edge in loops
        crossing = rows[heads[edge]]
        push!(result,crossing[4] == edge ? [edge] : [-edge])
        delete!(available,edge)
        delete!(available,-edge)
    end

    while !isempty(available)
        edge = maximum(available)
        delete!(available,edge)
        region = Int[]
        while !(edge in region)
            push!(region,edge)
            crossing_index = edge > 0 ? heads[edge] : tails[-edge]
            crossing = rows[crossing_index]
            slot = findfirst(==(abs(edge)),crossing)
            next_edge = crossing[mod1(slot-1,4)]
            if [next_edge] in result
                push!(region,-next_edge)
                next_edge = crossing[mod1(slot+1,4)]
            elseif [-next_edge] in result
                push!(region,next_edge)
                next_edge = crossing[mod1(slot+1,4)]
            end
            edge = tails[next_edge] == crossing_index ? next_edge : -next_edge
            delete!(available,edge)
        end
        push!(result,region)
    end
    result
end

"""
    hmove_matrix(v::VirtualLink) -> Matrix{Int}

Return the real-edge weight changes of unit H moves at the real vertices.

# Arguments
`v` must have coherent Gauss and vertex data with defined real-edge and
real-Gauss views. Cut markers must obey the
[data model](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
An incomplete real-edge table throws `ArgumentError`; malformed input may
propagate marker, indexing or conversion errors. Assigned weights are not
required for obtaining this matrix.

# Returns
A `Matrix{Int}` whose rows follow `real_edge_table(v)` and whose columns
follow the real vertices in `real_gauss_code(v)`. Column `j` is the additive
weight change for a unit H move at that real vertex.

# Guarantees
1. The matrix has one row per real edge and one column per real vertex.
   [proved/agent: by construction]
2. Its columns express the local weight change specified by the
   [H-move rule](../specs/virtual-link-moves.md#124-h-move), with repeated
   incidences contributing to the same edge entry. The graph presentation
   and any assigned weights remain unchanged. [proved/agent: by construction]
3. The result depends on the real-edge incidences, not the values of the
   supplied weight vectors. [proved/agent: by construction]

# Not guaranteed
Obtaining the matrix does not establish that assigned weights satisfy a
framing equation, choose a solution, or apply an H move. Application has
the additional weighted-input and parameter requirements on the linked page.
The matrix is presentation-dependent and is not a 3-manifold identifier.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-moves.md §12.4 — the H-move weight update represented by each column.

# Examples
```jldoctest
julia> h = hmove_matrix(lens_space_link(2, 1));

julia> (Tuple(h[:, 1]), Tuple(h[:, 2]))
((0, -1, 0, 1), (0, 1, 0, -1))
```

See also [`calculate_weight`](@ref), [`move!`](@ref).

Role: convert
"""
function hmove_matrix(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    head=zeros(Int,vertex_count,2)
    for (strand,component0) in enumerate(code)
        component=_component_visits(component0)
        isempty(component) && continue
        cut=_is_cut_component(v.gauss[strand])
        for index in eachindex(component)
            vertex=component[index]
            successor = cut && index == lastindex(component) ? 0 :
                component[mod1(index+1,length(component))]
            head[abs(vertex),vertex > 0 ? 2 : 1]=successor
        end
    end
    predecessor(vertex)=begin
        negative=findfirst(==(vertex),view(head,:,1))
        positive=findfirst(==(vertex),view(head,:,2))
        something(positive,0)-something(negative,0)
    end
    vertices=collect(1:vertex_count)
    signed=hcat(-vertices,vertices)
    previous=predecessor.(signed)
    pairs=vcat(hcat(vec(previous),vcat(-vertices,vertices)),
               hcat(vcat(-vertices,vertices),vec(head)))
    crossings=[edge.crossing for edge in real_edge_table(v)]
    indices=[something(findfirst(==(Tuple(pair)),crossings),0) for pair in eachrow(pairs)]
    any(==(0),indices) && throw(ArgumentError("real edge table is incomplete"))
    indices=reshape(indices,vertex_count,4)
    matrix=zeros(Int,length(crossings),vertex_count)
    for vertex in vertices
        matrix[indices[vertex,1],vertex]=-1
        matrix[indices[vertex,2],vertex]=-1
        matrix[indices[vertex,3],vertex]+=1
        matrix[indices[vertex,4],vertex]+=1
    end
    matrix
end

function _disk_component_path(nodes::Vector{Int},
                              edges::Vector{Tuple{Int,Int,Int}},
                              is_cycle::Bool)
    isempty(nodes) && return (Int[],Int[])
    outgoing=Dict{Int,Vector{Tuple{Int,Int}}}()
    incoming=Dict(node=>0 for node in nodes)
    for (source,target,coefficient) in edges
        push!(get!(outgoing,source,Tuple{Int,Int}[]),(target,coefficient))
        incoming[target]=get(incoming,target,0)+1
    end
    all(length(get(outgoing,node,Tuple{Int,Int}[])) <= 1 for node in nodes) ||
        throw(ArgumentError("disk component has branching outgoing edges"))
    all(get(incoming,node,0) <= 1 for node in nodes) ||
        throw(ArgumentError("disk component has branching incoming edges"))

    if is_cycle
        all(length(get(outgoing,node,Tuple{Int,Int}[])) == 1 for node in nodes) ||
            throw(ArgumentError("disk cycle has a missing outgoing edge"))
        all(get(incoming,node,0) == 1 for node in nodes) ||
            throw(ArgumentError("disk cycle has a missing incoming edge"))
        start=minimum(nodes)
    else
        starts=sort([node for node in nodes if get(incoming,node,0) == 0])
        length(starts) == 1 ||
            throw(ArgumentError("open disk component must have one initial node"))
        start=only(starts)
    end

    path=Int[]; edge_signs=Int[]; current=start
    while true
        push!(path,current)
        next_edges=get(outgoing,current,Tuple{Int,Int}[])
        isempty(next_edges) && break
        target,coefficient=only(next_edges)
        push!(edge_signs,coefficient)
        current=target
        is_cycle && current == start && break
        current in path && throw(ArgumentError("disk component repeats before closing"))
    end
    length(path) == length(nodes) && length(edge_signs) == length(edges) ||
        throw(ArgumentError("disk component path does not cover its graph"))
    path,edge_signs
end

"""Return the MATLAB `getDiskTable` data for a virtual link.

The MATLAB columns `ID`, `paths`, `edges_s`, `Ndots`, `isCycle`, `delta`, and
`cp` are retained verbatim as NamedTuple fields.  `cp2 == 2cp` is the exact
integer-scaled form used by Julia's weight solver.  For open paths it uses the
implementation convention `2*isCycle - Ndots`; boundary lane identities are
reported separately and are not folded into a closed framing claim.
"""
function _disk_node_dictionary()
    dp=zeros(Int,12); dm=zeros(Int,12)
    dp[[1,2,4,5,9,12,7,6,8,11,10,3]]=repeat(1:6,2)
    dm[[1,2,4,5,9,12,7,10,3,11,6,8]]=repeat(1:6,2)
    Dict(1=>dp,-1=>dm)
end

function _empty_disk_table()
    (ID=Int[],paths=Vector{Vector{Int}}(),edges_s=Vector{Vector{Int}}(),
     Ndots=Int[],isCycle=Bool[],delta=zeros(Int,0,0),
     cp=Rational{Int}[],cp2=Int[],path_lengths=Int[],
     boundary_start=Union{Nothing,NamedTuple}[],
     boundary_end=Union{Nothing,NamedTuple}[],
     boundary_nodes=NamedTuple[],boundary_lanes=NamedTuple[],
     open_paths=Vector{Vector{Int}}(),
     open_path_rows=Int[])
end

"""
    disk_table(v::VirtualLink) -> NamedTuple

Return the Method-B circuit and open-path data of a virtual-link presentation.

# Arguments
`v` must have coherent Gauss visits, vertex signs and
[cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
The real-Gauss and real-edge views must be defined. Malformed markers or
branching, missing, repeated or uncovered path incidences throw `ArgumentError`
identifying the failed check. Invalid indices and unsupported vertex signs
can instead propagate `BoundsError` or `KeyError`.

# Returns
A `NamedTuple` with one row per circuit or open path:
- `ID`, `paths`, `edges_s`, `Ndots`, `isCycle`: row identifiers, ordered node
  paths, signed edge traversals, dot counts and cycle flags.
- `delta::Matrix{Int}`: signed edge-incidence totals, with columns in
  `real_edge_table` order; `cp` and `cp2` are rational and integer-scaled
  obstruction data; `path_lengths` gives the node-path lengths.
- `boundary_start`, `boundary_end`: optional boundary endpoint records per row.
- `boundary_nodes`, `boundary_lanes`: boundary-end and individual lane records.
- `open_paths`, `open_path_rows`: the open rows and their indices.
The local Method-B model is owned by the
[move-verifier specification](../specs/ograph-move-verifier.md#ps-gate), and the
closedness use by [Enforcement](../specs/ograph-terminology.md#enforcement).

# Guarantees
1. Row-valued fields have the same length as `ID`, and `delta` has that many
   rows and one column per returned real edge. Each row's `delta` entries
   equal the signed occurrence totals in its `edges_s`. [proved/agent: by construction]
2. `path_lengths == length.(paths)`, `cp2 == 2 .* cp`, and `open_paths`
   selects exactly the rows whose `isCycle` flag is false. [proved/agent: by construction]
3. The paths describe the repository's Method-B local model of this
   presentation. Its use as the C3 circuit count has the evidence scope of
   the linked specification. The input is unchanged. [proved/agent: by construction]
4. An input with no real vertices and no real edges returns an empty table.
   Cut boundary arcs can still produce open rows with no real vertices.
   [proved/agent: by construction]

# Not guaranteed
This is not a closed-normal certificate. Open rows are not closed disc
boundaries, and their obstruction data alone do not certify a closed framing.
The general identification with BP97 Figure 1.3 is not re-derived here;
the specification records that dependency. Node labels and row order are
presentation coordinates, not manifold invariants.

# Sources
cf: lib/topology.jl/docs/specs/ograph-move-verifier.md §PS gate — Method-B local-node and lane model, including the source locators and index correspondence.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> d = disk_table(vl_example("s3_abalone"));

julia> Tuple(d.path_lengths)
(5, 1)
```

See also [`closedness_witness`](@ref), [`compare_disc_partition`](@ref).

Role: convert
"""
function disk_table(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    edges=real_edge_table(v)
    edge_count=length(edges)
    vertex_count == 0 && isempty(edges) && return _empty_disk_table()
    dictionary=_disk_node_dictionary()

    endpoint_records=NamedTuple[]
    endpoint_by_edge=Dict{Tuple{Int,Symbol},Int}()
    next_boundary_node=6vertex_count+1
    for edge in edges
        for (side,endpoint) in ((:start,edge.crossing[1]),(:end,edge.crossing[2]))
            endpoint == 0 || continue
            public_side=side == :start ? :in : :out
            ids=(next_boundary_node,next_boundary_node+1,next_boundary_node+2)
            next_boundary_node += 3
            record=(strand=edge.strand,side=public_side,lane=0,node_ids=ids,edge_id=edge.id)
            record_index=length(endpoint_records)+1
            push!(endpoint_records,record)
            endpoint_by_edge[(edge.id,side)]=record_index
        end
    end
    node_count=next_boundary_node-1
    parent=collect(1:node_count)
    find_root(node)=begin
        while parent[node] != node
            parent[node]=parent[parent[node]]
            node=parent[node]
        end
        node
    end
    unite(a,b)=begin
        ra,rb=find_root(a),find_root(b)
        ra != rb && (parent[rb]=ra)
    end
    graph_edges=Tuple{Int,Int,Int}[]
    for edge in edges
        source,target=edge.crossing
        source_nodes = if source == 0
            endpoint_records[endpoint_by_edge[(edge.id,:start)]].node_ids
        else
            source_slots=(1:3) .+ 6 .+ 3*(source > 0)
            Tuple(6*(abs(source)-1) .+ dictionary[orientation[abs(source)]][source_slots])
        end
        target_nodes = if target == 0
            endpoint_records[endpoint_by_edge[(edge.id,:end)]].node_ids
        else
            target_slots=(1:3) .+ 3*(target > 0)
            Tuple(6*(abs(target)-1) .+ dictionary[orientation[abs(target)]][target_slots])
        end
        for index in 1:3
            coefficient=index == 3 ? -edge.id : edge.id
            a,b=index == 3 ?
                (target_nodes[index],source_nodes[index]) :
                (source_nodes[index],target_nodes[index])
            push!(graph_edges,(a,b,coefficient)); unite(a,b)
        end
    end
    roots=unique(find_root.(1:node_count))
    root_index=Dict(root=>index for (index,root) in enumerate(roots))
    delta=zeros(Int,length(roots),edge_count)
    for (a,_,coefficient) in graph_edges
        row=root_index[find_root(a)]
        delta[row,abs(coefficient)] += sign(coefficient)
    end
    component_nodes=[[node for node in 1:node_count if find_root(node) == root]
                     for root in roots]
    component_edges=[[(a,b,coefficient) for (a,b,coefficient) in graph_edges
                      if find_root(a) == root] for root in roots]
    is_cycle=[length(edges_) == length(nodes_)
              for (nodes_,edges_) in zip(component_nodes,component_edges)]
    paths=Vector{Vector{Int}}(undef,length(roots))
    edge_signs=Vector{Vector{Int}}(undef,length(roots))
    for row in eachindex(roots)
        paths[row],edge_signs[row]=_disk_component_path(
            component_nodes[row],component_edges[row],is_cycle[row])
    end

    positive_dots=[0,1,0,0,1,0]
    negative_dots=[0,0,1,0,0,1]
    node_dots=zeros(Int,node_count)
    for vertex in 1:vertex_count, local_node in 1:6
        dots=orientation[vertex] > 0 ? positive_dots : negative_dots
        node_dots[6*(vertex-1)+local_node]=dots[local_node]
    end
    ndots=[sum(node_dots[nodes_]) for nodes_ in component_nodes]
    cp2=2 .* Int.(is_cycle) .- ndots
    cp=cp2 .// 2
    path_lengths=length.(paths)
    ids=collect(eachindex(roots))
    endpoint_at(node) = begin
        for record in endpoint_records
            lane=findfirst(==(node),record.node_ids)
            lane === nothing || return (strand=record.strand,side=record.side,lane=lane)
        end
        nothing
    end
    boundary_start=Union{Nothing,NamedTuple}[
        is_cycle[row] ? nothing : endpoint_at(first(paths[row])) for row in ids]
    boundary_end=Union{Nothing,NamedTuple}[
        is_cycle[row] ? nothing : endpoint_at(last(paths[row])) for row in ids]
    boundary_nodes=[(strand=record.strand,side=record.side,
                     node_ids=record.node_ids,edge_id=record.edge_id)
                    for record in endpoint_records]
    boundary_lanes=[(strand=record.strand,side=record.side,lane=lane,
                     node_id=record.node_ids[lane],node_ids=record.node_ids,
                     edge_id=record.edge_id)
                    for record in endpoint_records for lane in 1:3]
    open_path_rows=findall(!,is_cycle)
    open_paths=paths[open_path_rows]
    (ID=ids,paths=paths,edges_s=edge_signs,Ndots=ndots,isCycle=is_cycle,
     delta=delta,cp=cp,cp2=cp2,path_lengths=path_lengths,
     boundary_start=boundary_start,boundary_end=boundary_end,
     boundary_nodes=boundary_nodes,boundary_lanes=boundary_lanes,
     open_paths=open_paths,
     open_path_rows=open_path_rows)
end

function _nullspace_integer(A::Matrix{Int})
    rows,columns=size(A)
    reduced=Rational{BigInt}.(A)
    pivots=Int[]
    pivot_row=1
    for column in 1:columns
        source=findfirst(row -> reduced[row,column] != 0,pivot_row:rows)
        source === nothing && continue
        source=source+pivot_row-1
        reduced[pivot_row,:],reduced[source,:]=copy(reduced[source,:]),copy(reduced[pivot_row,:])
        reduced[pivot_row,:] ./= reduced[pivot_row,column]
        for row in 1:rows
            row == pivot_row && continue
            reduced[row,:] .-= reduced[row,column].*reduced[pivot_row,:]
        end
        push!(pivots,column); pivot_row += 1
        pivot_row > rows && break
    end
    free=setdiff(collect(1:columns),pivots)
    isempty(free) && return zeros(Int,0,columns)
    basis=Vector{Vector{Int}}()
    for free_column in free
        vector=zeros(Rational{BigInt},columns); vector[free_column]=1
        for (row,pivot) in enumerate(pivots)
            vector[pivot]=-reduced[row,free_column]
        end
        common_denominator=foldl(lcm,(value.den for value in vector);init=BigInt(1))
        integers=BigInt[value.num for value in vector.*common_denominator]
        divisor=foldl(gcd,abs.(integers);init=BigInt(0))
        divisor != 0 && (integers .÷= divisor)
        push!(basis,Int.(integers))
    end
    reduce(vcat,(permutedims(vector) for vector in basis))
end

function _rational_particular(A::Matrix{Int},b::Vector{Int})
    rows,columns=size(A)
    augmented=hcat(Rational{BigInt}.(A),Rational{BigInt}.(b))
    pivot_row=1
    pivots=Int[]
    for column in 1:columns
        source=findfirst(row -> augmented[row,column] != 0,pivot_row:rows)
        source === nothing && continue
        source += pivot_row-1
        augmented[pivot_row,:],augmented[source,:]=copy(augmented[source,:]),copy(augmented[pivot_row,:])
        augmented[pivot_row,:] ./= augmented[pivot_row,column]
        for row in 1:rows
            row == pivot_row && continue
            augmented[row,:] .-= augmented[row,column].*augmented[pivot_row,:]
        end
        push!(pivots,column); pivot_row += 1
        pivot_row > rows && break
    end
    any(row -> all(==(0),augmented[row,1:columns]) && augmented[row,end] != 0,1:rows) &&
        throw(ArgumentError("weight equations are inconsistent"))
    result=zeros(Rational{BigInt},columns)
    for (row,pivot) in enumerate(pivots)
        result[pivot]=augmented[row,end]
    end
    result
end

function _integer_particular(A::Matrix{Int},b::Vector{Int})
    rows,columns=size(A)
    model=Model(GLPK.Optimizer); set_silent(model)
    @variable(model,weight[1:columns],Int)
    for row in 1:rows
        @constraint(model,sum(A[row,column]*weight[column] for column in 1:columns)==b[row])
    end
    @objective(model,Min,0)
    optimize!(model)
    status=termination_status(model)
    status == MOI.INFEASIBLE && return nothing
    status == MOI.OPTIMAL ||
        throw(ArgumentError("integer weight feasibility did not terminate: $status"))
    candidate=round.(Int,value.(weight))
    all(A*candidate .== b) ||
        throw(ArgumentError("integer weight solver returned a non-solution"))
    candidate
end

"""MATLAB `minimizeWeightSolution`: nearest lattice representative in a kernel coset.

The particular solution may be rational.  MATLAB's `Matrix(ZZ).solve_right`
returns that rational solution when the Euler cochain is not in the integral
image; only the homogeneous kernel and its H-move quotient are integral.
"""
function _minimize_weight_solution(weight::AbstractVector{<:Number},kernel::Matrix{Int})
    weight=collect(weight)
    size(kernel,1)==0 && return weight
    # MATLAB's `-W / H`: minimize ||xH+W|| over real row vectors x.
    coefficients=vec(round.(Int,(-permutedims(Float64.(weight))) / Float64.(kernel)))
    current=weight+vec(transpose(coefficients)*kernel)
    value=sum(abs2,current)
    improved=true
    while improved
        improved=false
        for index in eachindex(coefficients)
            for delta in (1,-1)
                candidate_coefficients=copy(coefficients); candidate_coefficients[index]+=delta
                candidate=weight+vec(transpose(candidate_coefficients)*kernel)
                candidate_value=sum(abs2,candidate)
                if candidate_value<value
                    coefficients=candidate_coefficients; current=candidate; value=candidate_value
                    improved=true
                    break
                end
            end
            improved && break
        end
    end
    current
end

"""
    calculate_weight(v::VirtualLink; assign=false) -> Tuple{Vector,Matrix{Int},Matrix{Int}}

Return the legacy particular weight solution, nullspace rows and H-move complement diagnostic.

# Arguments
`v` must be coherent, in integer mode and have closed strands. Z2 mode and
cuts throw `ArgumentError` identifying the unsupported domain. `assign` is a
Boolean; when true, the exact integer solver must report a solution and the
legacy representative must pass exact equation and safe-storage checks, or
`ArgumentError` reports the obstruction or storage failure. Disk-table, linear
algebra and solver errors propagate.

# Returns
`(base, kernel, hperp)`: a particular solution vector, integer-valued rows
spanning the rational nullspace, and a heuristic complement to H moves, with
the scope in [the legacy diagnostic specification](../specs/virtual-link-moves.md).
The base may be rational when `assign=false`.

# Guarantees
1. With `assign=false`, the link is unchanged. With `assign=true`, successful
   return assigns the returned base to real arcs, marks the link weighted,
   and preserves the diagram data. [proved/agent: by construction]
2. An assigned base satisfies the current exact integer weight equation.
   The returned tuple concerns this presentation's weight system, not an
   identification of framing or spin equivalence classes.
   [proved/agent: by construction]

# Not guaranteed
The kernel rows need not generate the full integer lattice; `hperp` is not a
certified integral quotient basis, and representative minimization is not a
global optimum guarantee. Rational-only systems may return diagnostics but
cannot be assigned. Use `solve_integer_weights` for certified lattice data.

# Sources
cf: lib/topology.jl/docs/specs/virtual-link-moves.md §12 — the legacy diagnostic tuple distinguished from exact integer-weight certificates.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); base, kernel, hperp = calculate_weight(v); d = disk_table(v);

julia> all(iszero, 2d.delta * base + d.cp2)
true
```

See also [`solve_integer_weights`](@ref), [`set_weight!`](@ref), [`validate_weight`](@ref).

Role: invariant
"""
function calculate_weight(v::VirtualLink; assign=false)
    v.weight_mode==:integer || throw(ArgumentError("calculate_weight requires integer mode; use solve_spin_weights for Z2"))
    any(strand_cuts(v)) &&
        throw(ArgumentError("calculate_weight currently requires closed strands; use disk_table boundary metadata for open data"))
    if assign
        exact=solve_integer_weights(v)
        exact.status==:integer_solution ||
            throw(ArgumentError("no integral weight assignment ($(exact.status)); use solve_integer_weights for the obstruction"))
    end
    disks=disk_table(v)
    edge_count=size(disks.delta,2)
    base=_rational_particular(2disks.delta,-disks.cp2)
    if any(value -> value.den != 1,base)
        integer_base=_integer_particular(2disks.delta,-disks.cp2)
        integer_base !== nothing && (base=integer_base)
    else
        base=Int[value.num for value in base]
    end
    kernel=_nullspace_integer(disks.delta)
    base=_minimize_weight_solution(base,kernel)
    hmoves=permutedims(hmove_matrix(v))
    selected_rows=Vector{Vector{Int}}()
    span=Float64.(hmoves)
    old_rank=rank(span)
    for row in eachrow(kernel)
        candidate=vcat(span,permutedims(Float64.(row)))
        new_rank=rank(candidate)
        if new_rank > old_rank
            push!(selected_rows,collect(row)); span=candidate; old_rank=new_rank
        end
    end
    selected=isempty(selected_rows) ? zeros(Int,0,edge_count) :
        reduce(vcat,(permutedims(row) for row in selected_rows))
    for selected_index in axes(selected,1)
        current=collect(selected[selected_index,:])
        improved=true
        while improved
            improved=false
            for move in eachrow(hmoves),direction in (-1,1)
                candidate=current+direction.*move
                if sum(abs2,candidate) < sum(abs2,current) ||
                   (sum(abs2,candidate) == sum(abs2,current) && Tuple(candidate) > Tuple(current))
                    current=collect(candidate); improved=true
                end
            end
        end
        selected[selected_index,:]=current
    end
    if assign
        all(x -> isinteger(x) && abs(x)<=big(2)^53-1,base) ||
            throw(ArgumentError("legacy representative is not a safely storable integer; use solve_integer_weights"))
        all(iszero,2 .* BigInt.(disks.delta)*BigInt.(base)+BigInt.(disks.cp2)) ||
            throw(ArgumentError("legacy representative fails the exact integer weight equation"))
        weights = Vector{Vector{Float64}}(); offset = 1
        code,_=real_gauss_code(v)
        for g in code
            n=length(g)
            push!(weights,Float64.(base[offset:offset+n-1])); offset += n
        end
        set_weight!(v, weights)
    end
    base,kernel,selected
end

"""
    set_weight!(v::VirtualLink, weight; mode=v.weight_mode) -> VirtualLink

Assign real-arc weights to a virtual-link presentation in the chosen coefficient mode.

# Arguments
`v` must have coherent [real-arc data](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
`weight` is a numeric vector or matrix in real-arc order, or one vector per
Gauss component. Wrong component or arc counts throw `DimensionMismatch`
identifying the required count. `mode` must be `:integer` or `:z2`, otherwise
`ArgumentError` identifies the mode. Z2 input values must be finite integral
reals, checked with `ArgumentError`; integer mode only requires conversion to
`Float64`. Representation and conversion errors propagate.

# Returns
The same `VirtualLink`, with newly allocated per-component `Float64` weight
vectors, the selected mode and `is_weighted == true`.

# Guarantees
1. Successful assignment preserves Gauss data, signs and strand metadata.
   Shape and value checks finish before replacement of the weight fields.
   [proved/agent: by construction]
2. The weights label the real arcs of this presentation. In Z2 mode their
   stored representatives are 0 or 1; integer mode retains the converted input
   values. The call does not assert a framing or a spin structure.
   [proved/agent: by construction]

# Not guaranteed
Integer mode does not check integrality, finiteness or the weight equation,
and conversion to `Float64` can lose precision. Z2 reduction also does not
check spin admissibility. The
[weight assignment policy](../specs/virtual-link-data-model.md#24-weight-assignment-and-endpoint-coordinates)
specifies storage and flattening, rather than a mathematical equivalence rule.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §5 and §10 — edge weights and the scope of integer-weight equations.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_weight!(v, [3,-2]; mode=:z2);

julia> Tuple(only(v.weights))
(1.0, 0.0)
```

See also [`set_weight_mode!`](@ref), [`real_edge_table`](@ref), [`validate_weight`](@ref).

Role: move
"""
function set_weight!(v::VirtualLink, weight; mode=v.weight_mode)
    _check_weight_mode(mode)
    real_counts=[_strand_arc_count(v,strand) for strand in eachindex(v.gauss)]
    convert_value = mode==:z2 ? (x -> begin
        x isa Real && isfinite(x) && isinteger(x) || throw(ArgumentError("Z2 weights require finite integer-valued inputs"))
        Float64(mod(BigInt(x),2))
    end) : Float64
    weights=Vector{Vector{Float64}}()
    if weight isa AbstractMatrix || weight isa AbstractVector{<:Number}
        flat=vec(weight)
        length(flat)==sum(real_counts) || throw(DimensionMismatch("one weight per real arc is required"))
        idx=1
        for n in real_counts
            push!(weights,convert_value.(flat[idx:idx+n-1])); idx+=n
        end
    else
        length(weight)==length(real_counts) || throw(DimensionMismatch("one weight vector per Gauss component is required"))
        for (values,count_) in zip(weight,real_counts)
            length(values)==count_ || throw(DimensionMismatch("one weight per real arc is required"))
            push!(weights,convert_value.(values))
        end
    end
    v.weights=weights; v.weight_mode=mode; v.is_weighted=true; v
end

"""
    set_weight_by_crossing!(v::VirtualLink, crossings, weights) -> VirtualLink

Assign real-arc weights selected by their signed crossing endpoints.

# Arguments
- `v`: a coherent link accepted by [`real_edge_table`](@ref) and
  [`set_weight!`](@ref); their domain and shape errors propagate.
- `crossings`: a matrix whose rows are signed endpoint pairs in the real-edge
  table's coordinates. Entries must convert to `Int`; the caller supplies
  two columns and intended pair identities, which are not separately checked.
- `weights`: one value per row, otherwise `DimensionMismatch`. Integer mode
  converts values to `Float64`; Z2 mode requires finite integer-valued real
  values, otherwise `ArgumentError`, and uses the mode's normalized values.

# Returns
The same `VirtualLink` with assigned weights in its existing coefficient mode.

# Guarantees
1. Every real edge with a selected endpoint pair receives the associated value; the last repeated assignment to that pair wins. [proved/agent: by construction]
2. Unselected edges retain their real-edge-table weights, subject to the selected mode's conversion; `is_weighted` is true on successful return. [proved/agent: by construction]
3. The link's geometry is retained and its weights are modified in place. [proved/agent: by construction]

# Conventions
Use [real arcs and their weights](../specs/ograph-terminology.md#5-edges-and-weights)
and the [weight-assignment coordinates](../specs/virtual-link-data-model.md#24-weight-assignment-and-endpoint-coordinates).

# Not guaranteed
Unmatched endpoint pairs are ignored. Several arcs with the same endpoint
pair receive the same assignment; the pair need not select a unique occurrence.
No integer framing equation, spin condition or finiteness check in integer
mode is supplied by assignment alone.

# Sources
cf:docs/specs/virtual-link-data-model.md — real-arc weight storage and the distinction between data assignment and mathematical validity.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); knot_complement!(v);

julia> set_weight_by_crossing!(v, [2 -2; -3 1; 3 4; 4 -4], [1,-1,-1,1]);

julia> Tuple(e.weight for e in real_edge_table(v))
(-1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, -1.0)
```

Role: move
"""
function set_weight_by_crossing!(v::VirtualLink,crossings,weights)
    pairs=[Tuple(Int.(pair)) for pair in eachrow(Int.(crossings))]
    length(pairs) == length(weights) || throw(DimensionMismatch("one weight per crossing pair is required"))
    convert_value = v.weight_mode==:z2 ? (x -> begin
        x isa Real && isfinite(x) && isinteger(x) || throw(ArgumentError("Z2 weights require finite integer-valued inputs"))
        Float64(mod(BigInt(x),2))
    end) : Float64
    assignments=Dict(pair=>convert_value(weight) for (pair,weight) in zip(pairs,weights))
    table=real_edge_table(v)
    by_strand=[Float64[] for _ in v.gauss]
    for edge in table
        push!(by_strand[edge.strand],get(assignments,edge.crossing,edge.weight))
    end
    set_weight!(v,by_strand)
end

function _graph_a_connected(v::VirtualLink)
    all(!=(0), v.orientation) || return false
    any(word -> any(==(0), word), v.gauss) && return false
    code, orientation = real_gauss_code(v)
    isempty(orientation) && return length(code) <= 1
    any(isempty, code) && return false
    length(connected_components(ograph_graph_a(v))) == 1
end

"""
    is_closed(v::VirtualLink) -> Tuple{Bool,NTuple{3,Bool}}

Report the C1–C3 closed-normal checks at the mutable VirtualLink boundary.

# Arguments
`v` must have coherent stored Gauss and vertex data. Virtual crossings, cut
markers and empty components with real vertices are outside the supported
closed branch and return failed checks. Malformed cut markers, graph data
or disk paths can propagate the errors of `real_gauss_code`, `ograph_graph_a`
or `disk_table` instead of a Boolean result.

# Returns
`(closed, criteria)` with `closed::Bool` and `criteria::NTuple{3,Bool}` in
C1, C2, C3 order, using [the enforcement specification](../specs/ograph-terminology.md#enforcement).

# Guarantees
1. On return, `closed == all(criteria)`. The input is unchanged.
   [proved/agent: by construction]
2. For coherent supported nonempty presentations these are the repository's
   closed-normal criteria; their geometric scope and the C3 correspondence
   dependency are stated on the linked specification page. [proved/agent: by construction]
3. Unsupported virtual or cut inputs report false C2 and C3 flags. The
   empty single circuit reports `(false, (true, true, false))`, consistently
   with the [MST boundary distinction](../specs/ograph-terminology.md#mst-crossing-free-circuit).
   [proved/agent: by construction]

# Not guaranteed
False C2 or C3 on an unsupported input is a domain rejection, not a computation
of that condition for a different mathematical class. This predicate does not
validate arbitrary manually corrupted `VirtualLink` data, certify weights,
or return a replayable `ClosedNormalOGraph`. Its empty-input behavior differs
from the MST closedness-witness extension.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4, C1–C3 — closed-normal conditions; the repository enforcement and MST boundary are specified separately.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> is_closed(vl_example("s3_abalone"))
(true, (true, true, true))
```

See also [`certify_closed_ograph`](@ref), [`closedness_witness`](@ref).

Role: condition
"""
function is_closed(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    # C1, C2, C3 at the VirtualLink boundary ([BP97, §1.4]).
    # C1 is the Gauss-component count.  The zero-vertex MST circuit is not
    # this predicate: `disk_table` has no rows there, so the third flag fails.
    c1 = length(code) == 1
    unsupported = any(==(0), v.orientation) ||
                  any(word -> any(==(0), word), v.gauss) ||
                  (!isempty(orientation) && any(isempty, code))
    unsupported && return false, (c1, false, false)
    c2 = _graph_a_connected(v)
    c3 = size(disk_table(v).delta,1) == length(orientation)+1
    all((c1,c2,c3)), (c1,c2,c3)
end

include("virtual_link_symmetries.jl")
"""
    cyclic_shift!(v::VirtualLink, strand::Integer, amount::Integer) -> VirtualLink

Change the cyclic origin of one stored closed Gauss component in place.

# Arguments
`strand` must index an existing component of `v.gauss`; otherwise `BoundsError`
propagates. Cut markers must obey
[the data model](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers),
otherwise `ArgumentError` identifies the malformed marker. For a cut component
only `amount == 0` is accepted; any nonzero amount throws `ArgumentError`
identifying the undefined cyclic shift. `amount` is an integer shift.

# Returns
The same `VirtualLink` object `v`, with the selected Gauss component shifted.

# Guarantees
1. The new selected word equals `circshift` of its old value by `amount`.
   Other components and every field other than `gauss` are unchanged.
   [proved/agent: by construction]
2. On a closed unweighted presentation this changes a word's origin, not
   its oriented cyclic word, under the
   [Gauss-code interpretation](../specs/ograph-terminology.md#2-standard-gauss-code-displays).
   [proved/agent: by construction]
3. A zero shift of a valid cut component is accepted and returns `v`.
   [proved/agent: by construction]

# Not guaranteed
Weights are not rotated with the word. This operation therefore does not
promise to preserve the weighted presentation. It does not canonicalize,
relabel crossings, certify closedness or change the cyclic direction.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §2 — cyclic component words and their origins; cut components use the separate data-model boundary rule.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_cuts!(v, [true]);

julia> cyclic_shift!(v, 1, 0) === v
true
```

See also [`gauss_code`](@ref).

Role: move
"""
function cyclic_shift!(v::VirtualLink, strand::Integer, amount::Integer)
    _validate_cut_markers(v.gauss)
    _is_cut_component(v.gauss[strand]) && amount != 0 &&
        throw(ArgumentError("cyclic_shift! on a cut strand is undefined; only amount=0 is allowed"))
    v.gauss[strand] = circshift(v.gauss[strand], amount); v
end

"""
    disjoint_sum!(a::VirtualLink, b::VirtualLink; preserve_weights=true) -> VirtualLink

Append a disjoint virtual-link presentation to the first operand in place.

# Arguments
`a` and `b` must contain coherent presentations with valid
[cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
Malformed markers throw `ArgumentError`. `preserve_weights` must be an actual
`Bool`, otherwise `ArgumentError` names that requirement. When preservation
is requested and both operands are weighted, their weight modes must agree;
otherwise `ArgumentError` identifies the mismatch. Use independent operands
if the second must remain unchanged through aliasing.

# Returns
The same object `a`, containing the first operand's components followed by
the second operand's relabelled components.

# Guarantees
1. The appended presentation uses fresh crossing labels and preserves its
   passage signs, vertex signs, strand types and cut markers. The result
   represents the disjoint union of the two presented virtual-link diagrams.
   [proved/agent: by construction]
2. Weight preservation follows the
   [disjoint-sum rule](../specs/ograph-connected-sum.md#preconditions-and-postconditions).
   If preservation is disabled or either operand is unweighted, the result
   is marked unweighted with zero-filled arc weight vectors. [proved/agent: by construction]
3. The first operand's weight mode is retained. `real_only` is true exactly
   when it was true on both inputs. An independent `b` is unchanged.
   [proved/agent: by construction]

# Not guaranteed
This is not a connected sum and does not produce a closed-normal o-graph
certificate. It does not validate supplied weights or claim that a disjoint
diagram union is the connected sum of represented 3-manifolds. Aliased inputs
are not isolated by a copy of the whole operand.

# Sources
cf: lib/topology.jl/docs/specs/ograph-connected-sum.md §Preconditions and postconditions — disjoint-sum weight retention before the separate connected-sum operation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> a = vl_example("s3_abalone"); set_cuts!(a, [true]); b = deepcopy(a);

julia> disjoint_sum!(a, b);

julia> Tuple(strand_cuts(a))
(true, true)
```

See also [`connected_sum!`](@ref).

Role: convert
"""
function disjoint_sum!(a::VirtualLink, b::VirtualLink; preserve_weights=true)
    _validate_cut_markers(a.gauss); _validate_cut_markers(b.gauss)
    preserve_weights isa Bool || throw(ArgumentError("preserve_weights must be an explicit Boolean"))
    weighted=preserve_weights && a.is_weighted && b.is_weighted
    weighted && a.weight_mode != b.weight_mode &&
        throw(ArgumentError("weighted disjoint sum requires matching weight modes"))
    a_weights=deepcopy(a.weights)
    b_weights=deepcopy(b.weights)
    offset=length(a.orientation)
    append!(a.gauss, [[sign(x)*offset+x for x in g] for g in b.gauss])
    append!(a.orientation,b.orientation)
    append!(a.strand_types,b.strand_types)
    a.real_only = a.real_only && b.real_only
    if weighted
        a.weights=vcat(a_weights,b_weights)
        a.is_weighted=true
    else
        a.is_weighted=false
        a.weights=[zeros(Float64,_strand_arc_count(a,strand))
                   for strand in eachindex(a.gauss)]
    end
    a
end

function _connected_sum_selection(selection)
    if selection isa AbstractMatrix
        size(selection,2) == 2 && size(selection,1) in (1,2) ||
            throw(ArgumentError("connected-sum selection must have shape 1×2 or 2×2"))
        all(value -> value isa Integer,selection) ||
            throw(ArgumentError("connected-sum selection entries must be integers"))
        strands=(Int(selection[1,1]),Int(selection[1,2]))
        anchors=size(selection,1) == 2 ?
            (Int(selection[2,1]),Int(selection[2,2])) : nothing
        return strands,anchors
    end

    applicable(iterate,selection) ||
        throw(ArgumentError("connected-sum selection must contain two strand indices"))
    values=collect(selection)
    length(values) == 2 ||
        throw(ArgumentError("connected-sum selection must contain two strand indices"))
    all(value -> value isa Integer,values) ||
        throw(ArgumentError("connected-sum strand indices must be integers"))
    (Int(values[1]),Int(values[2])),nothing
end

function _rotate_connected_sum_component(component::Vector{Int},anchor::Int,strand::Int)
    positions=findall(==(anchor),component)
    length(positions) == 1 ||
        throw(ArgumentError("connected-sum anchor $anchor must occur exactly once on strand $strand"))
    index=only(positions)
    vcat(component[index:end],component[begin:index-1])
end

function _rotate_connected_sum_weights(weights::Vector{Float64},component::Vector{Int},
                                       anchor::Int,strand::Int)
    positions=findall(==(anchor),component)
    length(positions) == 1 ||
        throw(ArgumentError("connected-sum anchor $anchor must occur exactly once on strand $strand"))
    index=only(positions)
    length(weights)==length(component) ||
        throw(DimensionMismatch("connected-sum weights must match real-Gauss component length"))
    vcat(weights[index:end],weights[begin:index-1])
end

const _CONNECTED_SUM_MAX_INTEGER=big(2)^53-1
function _connected_sum_integer(value)
    value isa Real && isfinite(value) && isinteger(value) ||
        throw(ArgumentError("connected-sum weight lift requires finite integer values"))
    integer=BigInt(value)
    abs(integer)<=_CONNECTED_SUM_MAX_INTEGER ||
        throw(ArgumentError("connected-sum weight exceeds the exact Float64 integer range"))
    integer
end

function _connected_sum_lift(first_weights::Vector{Float64},second_weights::Vector{Float64})
    isempty(first_weights) && throw(ArgumentError("first connected-sum component has no edge"))
    isempty(second_weights) && throw(ArgumentError("second connected-sum component has no edge"))
    foreach(_connected_sum_integer,first_weights)
    foreach(_connected_sum_integer,second_weights)
    u=Float64(_connected_sum_integer(first(first_weights)))
    v=Float64(_connected_sum_integer(first(second_weights)))
    # Local edges x1,...,x10 satisfy
    # x1+x10=u, x5+x6=v, x3=x8=1, x2+x4+x7+x9=-1.
    # The stored cyclic edge order starts with x10.
    vcat(u,first_weights[2:end],0.0,0.0,1.0,0.0,0.0,
         second_weights[2:end],v,0.0,1.0,-1.0)
end

"""
    connected_sum!(v::VirtualLink, selection=(1,2); transport_weights=true) -> VirtualLink

Replace two selected normal o-graph summands by their connected-sum presentation.

# Arguments
- `v`: the caller supplies disjoint closed normal summands satisfying the
  [connected-sum preconditions](../specs/ograph-connected-sum.md#preconditions-and-postconditions).
  These topological preconditions are not independently certified here; without
  them the manifold interpretation is not guaranteed. Open strands are rejected
  with `ArgumentError` explaining that open boundaries would be discarded.
- `selection`: two distinct strand indices, optionally with anchors, in the
  [selection convention](../specs/ograph-connected-sum.md#real-gauss-splice).
  Accepted containers are a two-element iterable, a `1×2` matrix, or a `2×2`
  matrix with anchors in its second row. `ArgumentError` reports invalid shape,
  noninteger entries, equal indices or a nonunique/missing anchor;
  out-of-range strand indices throw `BoundsError`.
- `transport_weights`: must be `Bool`, otherwise `ArgumentError`.
  Transport of assigned weights requires integer mode and framed-valid weights;
  `ArgumentError` identifies unsupported Z2 transport or invalid weights.
  Input and output weight validation is delegated to [`validate_weight`](@ref).

# Returns
The same `VirtualLink` object after replacement of its data.

# Guarantees
1. On successful return the strand count decreases by one and the true vertex count increases by four. [proved/agent: by construction]
2. For two summands meeting the linked preconditions, the specified local operation represents their oriented connected sum. [proved/literature: [[K10#^k10-sec-0033]] Figure 6 and Theorem 1.1]
3. When assigned weights are transported, the result passes framed weight validation and follows the specification's fixed-exterior integer lift. [proved/agent: by construction]
4. `transport_weights=false` returns an unweighted presentation; `v` is modified in place. [proved/agent: by construction]

# Conventions
The [real-Gauss splice](../specs/ograph-connected-sum.md#real-gauss-splice)
and [integer-weight lift](../specs/ograph-connected-sum.md#fixed-exterior-integer-weight-lift)
fix the labels, seams and weight representative.

# Not guaranteed
The weight lift is not established as a canonical connected sum of geometric
framing classes. Arbitrary pairs of strands need not be separate manifold
summands. No independent C1–C3 certification is performed on the result.

# Sources
cf:[[K10#^k10-sec-0033]] Figure 6 and Theorem 1.1 — oriented connected sum, under the linked local-operation specification.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> solve_integer_weights(v; assign=true); other = deepcopy(v);

julia> disjoint_sum!(v, other); connected_sum!(v);

julia> (v.is_weighted, is_framed_weight(v), is_closed(v))
(true, true, (true, (true, true, true)))
```

See also [`disjoint_sum!`](@ref), [`certify_closed_ograph`](@ref).

Role: move
"""
function connected_sum!(v::VirtualLink, selection=(1,2);transport_weights=true)
    transport_weights isa Bool ||
        throw(ArgumentError("transport_weights must be an explicit Boolean"))
    any(strand_cuts(v)) &&
        throw(ArgumentError("connected_sum! currently requires closed strands; open boundaries would be discarded"))
    transport=transport_weights && v.is_weighted
    if transport
        v.weight_mode==:integer ||
            throw(ArgumentError("connected-sum integer lift does not support Z2 mode"))
        validate_weight(v;framed=true)
    end
    code,orientation=real_gauss_code(v)
    (i,j),anchors=_connected_sum_selection(selection)
    i==j && throw(ArgumentError("distinct strands required"))
    1 <= i <= length(code) && 1 <= j <= length(code) || throw(BoundsError(code,(i,j)))
    first_component=code[i]
    second_component=code[j]
    first_weights=transport ? copy(v.weights[i]) : Float64[]
    second_weights=transport ? copy(v.weights[j]) : Float64[]
    if anchors !== nothing
        transport && begin
            first_weights=_rotate_connected_sum_weights(first_weights,first_component,anchors[1],i)
            second_weights=_rotate_connected_sum_weights(second_weights,second_component,anchors[2],j)
        end
        first_component=_rotate_connected_sum_component(first_component,anchors[1],i)
        second_component=_rotate_connected_sum_component(second_component,anchors[2],j)
    end
    vertex_count=length(orientation)
    a=vertex_count+1
    b=vertex_count+2
    c=vertex_count+3
    d=vertex_count+4
    bridge1=[-b,-a,a,d]
    bridge2=[-d,-c,c,b]
    code[i]=vcat(first_component,bridge1,second_component,bridge2)
    deleteat!(code,j)
    types=copy(v.strand_types)
    types[i]='o'
    deleteat!(types,j)
    output_weights=transport ? deepcopy(v.weights) : Vector{Vector{Float64}}()
    if transport
        output_weights[i]=_connected_sum_lift(first_weights,second_weights)
        deleteat!(output_weights,j)
    end
    target=VirtualLink(weight_mode=v.weight_mode)
    set_data!(target;rgauss=code,orientation=vcat(orientation,[-1,1,-1,1]),
              strand_types=types)
    if transport
        set_weight!(target,output_weights)
        validate_weight(target;framed=true)
    end
    for field in fieldnames(VirtualLink)
        setfield!(v,field,getfield(target,field))
    end
    v
end

"""
    knot_complement!(v::VirtualLink) -> VirtualLink

Return the version-2 combinatorial expansion associated with the knot-complement construction.

# Arguments
- `v`: a coherent closed-strand presentation in the
  [VirtualLink data model](../specs/virtual-link-data-model.md).
  Open strands are checked and rejected with `ArgumentError`. The caller
  supplies valid real-Gauss incidence and vertex signs; no knot recognition,
  planarity or single-component certification is performed.

# Returns
The same `VirtualLink` with expanded real-Gauss data, unassigned zero weights,
and o-graph strand types. Its coefficient mode is retained.

# Guarantees
1. The output has four times the input's true vertex count and twice its number of real-Gauss components. [proved/agent: by construction]
2. The result is unweighted and `real_only` is `true`. [proved/agent: by construction]
3. `v` is modified in place. [proved/agent: by construction]

# Not guaranteed
The [knot-complement word specification](../specs/ur-diagrams.md#knot-complement-word)
records the unverified identification with a spine of the input knot exterior,
its peripheral marking and its applicability beyond knots. This function does
not certify a knot-complement homeomorphism or C1–C3. Assigned weights and
their framing interpretation are discarded.

# Sources
cf:docs/specs/virtual-link-data-model.md — meaning of the input and output presentation data; no knot-exterior identification is supplied there.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); knot_complement!(v);

julia> Tuple(v.orientation)
(-1, 1, -1, 1)
```

Role: convert
"""
function knot_complement!(v::VirtualLink)
    any(strand_cuts(v)) &&
        throw(ArgumentError("knot_complement! currently requires closed strands"))
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    expand(vertex,inverted)=begin
        base=4(abs(vertex)-1)+1
        offset=(vertex < 0 ? 1+2*(orientation[abs(vertex)] < 0) : 0)+2inverted
        sign(vertex).*(base .+ mod.([0,1].+offset,4))
    end
    forward=[reduce(vcat,(expand(vertex,false) for vertex in component);init=Int[]) for component in code]
    inverse=[reduce(vcat,(expand(vertex,true) for vertex in reverse(component));init=Int[]) for component in code]
    v.gauss=vcat(forward,inverse)
    v.orientation=repeat([-1,1,-1,1],vertex_count)
    v.real_only=true
    v.is_weighted=false
    v.weights=[zeros(Float64,length(g)) for g in v.gauss]
    v.strand_types=fill('o',length(v.gauss))
    v
end

"""
    vl_to_urdiagram(v::VirtualLink) -> URDiagram

Return the legacy U(R) diagram associated with a link's real-Gauss projection.

# Arguments
`v` must have integer weight mode and closed strands; Z2 mode and open
strands throw `ArgumentError` identifying the unsupported domain. The caller
supplies coherent [link data](../specs/virtual-link-data-model.md).
Malformed cut markers throw `ArgumentError`; inconsistent crossing data can
raise indexing or conversion errors. Assigned edge weights are not required.

# Returns
A new `URDiagram` in the legacy `VL2URD` ordering, with coefficient one and
the formal perturbations supplied by [`set_data_cvs!`](@ref).

# Guarantees
1. The result represents the real-crossing data in the legacy UR-diagram
   interface; explicit virtual crossings and assigned edge weights do not
   enter that result. [proved/agent: by construction]
2. The source link is unchanged, and the returned diagram has independent
   mutable containers. [proved/agent: by construction]

# Not guaranteed
This conversion does not certify a closed normal o-graph, identify a
3-manifold invariant, take a regulator limit or return a trace value.
The [U(R) specification](../specs/ur-diagrams.md#weighted-endpoint-data) owns
the endpoint data and links their ordering convention. This contract does not
assert a correspondence with the distinct exact U(R) MST backend. Multi-component input is permitted but no
general manifold interpretation is established here.

# Sources
cf:lib/topology.jl/src/virtual_link.jl vl_to_urdiagram; lib/topology.jl/src/ur_diagram.jl set_data_cvs! — the legacy conversion and formal-perturbation interfaces.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_weight!(v, [3,-2]); set_weight_mode!(v, :z2);

julia> try vl_to_urdiagram(v); false catch e; e isa ArgumentError end
true
```

See also [`ur_invariant`](@ref), [`ur_prime_invariant`](@ref).

Role: convert
"""
function vl_to_urdiagram(v::VirtualLink)
    v.weight_mode==:integer || throw(ArgumentError("UR invariants require integer mode; Z2 weights are not integer lifts"))
    any(strand_cuts(v)) &&
        throw(ArgumentError("vl_to_urdiagram currently requires closed strands"))
    code,orientation=real_gauss_code(v)
    vertices=Int[]
    for (index,component) in enumerate(code)
        index > 1 && push!(vertices,0)
        append!(vertices,reverse(component))
    end
    diagram=URDiagram()
    set_data_cvs!(diagram,1,vertices,orientation)
end

"""
    ur_invariant(v::VirtualLink) -> Any

Return the legacy symbolic trace of the link's UR-diagram projection.

# Arguments
`v` must have integer mode and closed strands; Z2 mode or cuts throw
`ArgumentError` naming the unsupported domain. Its real-Gauss data must be
coherent under the [data model](../specs/virtual-link-data-model.md).
Conversion, symbolic simplification and trace errors propagate. The function
does not require or check assigned framing weights or closed-normality.

# Returns
The scalar or symbolic expression produced by `trace_value` for
`vl_to_urdiagram(v)`. Its concrete type depends on the algebraic result;
`Any` describes that polymorphic return, not a certified numeric value.

# Guarantees
1. The result is the trace expression of the converted UR diagram, with the
   same meaning as [`vl_to_urdiagram`](@ref) and `trace_value`. The input link
   is unchanged. [proved/agent: by construction]
2. The input's assigned edge weights do not enter this expression.
   [proved/agent: by construction]

# Not guaranteed
This legacy entry is not the exact MST backend exposed by `ur_prime_invariant`.
The [legacy trace specification](../specs/ur-diagrams.md#legacy-trace-expressions)
does not establish general move invariance or identify this expression with
an invariant of a closed 3-manifold. It does
not promise a finite numeric value, removal of regulators, or a framed invariant.

# Sources
cf: lib/topology.jl/src/ur_diagram.jl trace_value — the legacy UR trace operation; this is an implementation reference, not a manifold-invariance theorem.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_weight!(v, [3,-2]); set_weight_mode!(v, :z2);

julia> try ur_invariant(v); false catch e; e isa ArgumentError end
true
```

See also [`vl_to_urdiagram`](@ref), [`ur_prime_invariant`](@ref).

Role: invariant
"""
function ur_invariant(v::VirtualLink)
    trace_value(vl_to_urdiagram(v))
end

"""
    ur_prime_invariant(v::VirtualLink) -> Rational{BigInt}

Return the exact U(R) MST trace value of a link's real-Gauss presentation.

# Arguments
`v` must have closed strands and one nonempty real-Gauss component, with
coherent signed visits under the [link data model](../specs/virtual-link-data-model.md).
Cuts, an empty or multiple-component real projection, and invalid visit
multiplicities throw `ArgumentError`. Incoherent underlying link data can
raise indexing errors. A singular defining trace matrix throws `DomainError`.
Assigned edge weights and their mode do not enter this evaluator.

# Returns
An exact rational scalar from `Z_MST(v, URAlgebra())`.

# Guarantees
1. The result is the U(R) specialization of the library's MST trace on the
   supplied real-Gauss presentation, with precisely the semantics of that
   `Z_MST` overload. [proved/agent: by construction]
2. The input's Gauss data and weights are unchanged, and assigned weights
   do not affect the value. [proved/agent: by construction]

# Not guaranteed
This entry does not certify closed-normality. The
[exact MST specification](../specs/ur-diagrams.md#exact-mst-expression) owns
the trace convention and records the absent general identification of this
scalar with an invariant of the represented 3-manifold.
The tested finite examples do not establish that theorem. No finite value
is supplied for singular trace data, and no framing dependence is asserted.

# Sources
cf:lib/topology.jl/src/ur_mst_link.jl invariant_expression and Z_MST; lib/topology.jl/src/ur_mst.jl evaluate_invariant — the exact U(R) trace interface, not a general manifold-invariance theorem.

# Examples
```jldoctest
julia> positive = VirtualLink(); set_data!(positive; rgauss=[[1,-1]], orientation=[1]); negative = VirtualLink(); set_data!(negative; rgauss=[[1,-1]], orientation=[-1]);

julia> (ur_prime_invariant(positive) == -1, ur_prime_invariant(negative) == 1)
(true, true)
```

See also [`Z_MST`](@ref), [`ur_invariant`](@ref).

Role: invariant
"""
ur_prime_invariant(v::VirtualLink) = Z_MST(v, URAlgebra())
