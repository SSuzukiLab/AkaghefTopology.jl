"""
    EdgeLayout(id, points, tail, head, tail_slot, head_slot) -> EdgeLayout

Store the routed polyline and incidence metadata of a PD edge.

# Arguments
Arguments must have or convert to the field types: integer identifiers and
slots, and `Vector{ComplexF64}` points. Conversion failures propagate.
The caller supplies the routing conditions of
[plot specification §6.5](../specs/virtual-link-plot.md#65-edge-routing);
this constructor checks none of them. Invalid fields do not certify a route.

# Returns
`EdgeLayout` with edge `id`, ordered `points`, crossing indices `tail` and
`head`, and corresponding PD slots `tail_slot` and `head_slot`.

# Guarantees
1. Fields retain the supplied data after conversion, representing the
   declared edge route when the caller's routing conditions hold.
   [proved/agent: by construction]
2. Inputs are not modified; an already typed point vector is shared.
   [proved/agent: by construction]

# Conventions
Layout coordinates: [plot specification §6.6](../specs/virtual-link-plot.md#66-linklayout).

# Not guaranteed
Construction does not check endpoint incidence, orthogonality, self-intersections
or global layout validity. It is not a knot-type or manifold certificate.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §6.5 — routed PD-edge metadata and required geometric properties.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; gauss=[[1,-2],[-1,2]], orientation=[1,1]); layout = orthogonal_layout(v); e = first(layout.edges);

julia> stored = EdgeLayout(e.id, e.points, e.tail, e.head, e.tail_slot, e.head_slot);

julia> first(stored.points) == layout.crossings[stored.tail] && last(stored.points) == layout.crossings[stored.head]
true
```

Role: construct
"""
struct EdgeLayout
    id::Int
    points::Vector{ComplexF64}
    tail::Int
    head::Int
    tail_slot::Int
    head_slot::Int
end

"""
    LinkLayout(pd, crossings, edges, bending_numbers, orientation) -> LinkLayout

Store a routed geometric representation of a PD link diagram.

# Arguments
Supply data convertible to the field types in
[plot specification §6.6](../specs/virtual-link-plot.md#66-linklayout).
Conversion failures propagate. Incidence, coordinate, bending and orientation
consistency are the caller's responsibility; no geometric check or layout
optimization is performed by this constructor.

# Returns
`LinkLayout` containing the PD matrix, crossing coordinates, routed
`EdgeLayout` records, signed bending numbers and crossing orientations.

# Guarantees
1. Fields retain the supplied data after conversion. When the caller's
   data satisfy the plot contract, the result stores that routed diagram.
   [proved/agent: by construction]
2. Inputs are not modified and already typed containers are shared.
   [proved/agent: by construction]

# Conventions
Layout coordinates and orientation metadata:
[plot specification §6.6](../specs/virtual-link-plot.md#66-linklayout) and
[§4.2](../specs/virtual-link-plot.md#42-crossing-orientation).

# Not guaranteed
The type alone does not establish the geometric properties of
[§11.2](../specs/virtual-link-plot.md#112-mandatory-geometric-properties),
bend optimality, or topological equivalence to a source diagram. Shared
mutable containers can invalidate previously checked geometry.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §6.6 and §11.2 — layout data and separate geometric acceptance.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; gauss=[[1,-2],[-1,2]], orientation=[1,1]); l = orthogonal_layout(v);

julia> stored = LinkLayout(l.pd, l.crossings, l.edges, l.bending_numbers, l.orientation);

julia> (length(stored.crossings), sort([e.id for e in stored.edges])) == (2, collect(1:4))
true
```

See also [`orthogonal_layout`](@ref), [`check_layout`](@ref).

Role: convert
"""
struct LinkLayout
    pd::Matrix{Int}
    crossings::Vector{ComplexF64}
    edges::Vector{EdgeLayout}
    bending_numbers::Vector{Int}
    orientation::Vector{Int}
end

"""
    minimal_bending_numbers(pd::AbstractMatrix{<:Integer}) -> Vector{Int}

Return signed edge bend counts minimizing the specified region-flow model.

# Arguments
Checked: conversion to machine integers and region construction must succeed;
a nonoptimal solver termination raises `ArgumentError`.
Caller: `pd` must satisfy L2 of the
[layout input conditions](../specs/virtual-link-plot.md#6-layout-stage-contracts) for the
[minimum-bend model](../specs/virtual-link-plot.md#63-minimum-bend-milp).

# Returns
Signed bend counts in increasing edge-label order; an empty matrix gives `Int[]`.

# Guarantees
1. On successful nonempty input, the counts come from an optimal solution of
   the [specified model](../specs/virtual-link-plot.md#63-minimum-bend-milp)
   for this PD embedding. [proved/agent: by construction and optimal termination check]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The optimum need not be unique, and a particular bend vector or drawing is
not canonical. Optimality is for this region model, not all embeddings of
the represented link; numerical optimizer tolerances are not an exact certificate.

# Sources
cf: [plot specification §6.3](../specs/virtual-link-plot.md#63-minimum-bend-milp) — model and outer-face convention.

# Examples
```jldoctest
julia> isempty(minimal_bending_numbers(zeros(Int, 0, 4)))
true
```

Role: construct
"""
function minimal_bending_numbers(pd::AbstractMatrix{<:Integer})
    isempty(pd) && return Int[]
    edge_ids = sort(unique(vec(Int.(pd))))
    edge_index = Dict(edge => index for (index,edge) in enumerate(edge_ids))
    link_regions = sort(regions(pd); by=length, alg=MergeSort)
    model = Model(GLPK.Optimizer)
    set_silent(model)
    @variable(model,bend[1:2length(edge_ids)] >= 0,Int)

    source_flow(edge) = edge > 0 ? bend[2edge_index[edge]-1] : bend[2edge_index[-edge]]
    sink_flow(edge) = source_flow(-edge)
    for (index,region) in enumerate(link_regions)
        capacity = index < length(link_regions) ? length(region)-4 : length(region)+4
        @constraint(model,sum(sink_flow(edge)-source_flow(edge) for edge in region) == capacity)
    end
    @objective(model,Min,sum(bend))
    optimize!(model)
    termination_status(model) == MOI.OPTIMAL ||
        throw(ArgumentError("orthogonal bend model is $(termination_status(model))"))
    values = round.(Int,value.(bend))
    [values[2index-1]-values[2index] for index in eachindex(edge_ids)]
end

function _piece_lengths(pd::Matrix{Int},bending::Vector{Int};
                        DEBUG_MODE::Bool=false,debug_io::IO=stderr)
    edge_ids = sort(unique(vec(pd)))
    edge_index = Dict(edge => index for (index,edge) in enumerate(edge_ids))
    segments_by_edge = Dict(edge => [(edge,k) for k in 0:abs(bending[edge_index[edge]])]
                            for edge in edge_ids)
    pieces = Dict{Tuple{Int,Int},Vector{Any}}(
        segment => Any[segment] for segments in values(segments_by_edge) for segment in segments)
    link_regions = sort(regions(pd); by=length, alg=MergeSort)
    new_regions = Vector{Vector{Tuple{Any,Int}}}()
    for region in link_regions[1:end-1]
        new_region = Tuple{Any,Int}[]
        for signed_edge in region
            if signed_edge > 0
                segments = segments_by_edge[signed_edge]
                sign_of_bend = sign(bending[edge_index[signed_edge]])
                append!(new_region,[(segment,sign_of_bend) for segment in segments[1:end-1]])
                push!(new_region,(last(segments),1))
            else
                segments = reverse(segments_by_edge[-signed_edge][2:end])
                sign_of_bend = sign(bending[edge_index[-signed_edge]])
                append!(new_region,[(segment,-sign_of_bend) for segment in segments])
                push!(new_region,(first(segments_by_edge[-signed_edge]),1))
            end
        end
        push!(new_regions,new_region)
    end

    next_segment = maximum(edge_ids)+1
    while true
        bad_index = findfirst(region -> any(item -> item[2] == -1,region),new_regions)
        bad_index === nothing && break
        bad_region = new_regions[bad_index]
        a = findfirst(item -> item[2] == -1,bad_region)
        bend_sum = -1
        b = a
        while bend_sum != 2
            b = b == length(bad_region) ? 1 : b+1
            bend_sum += bad_region[b][2]
        end
        shared = bad_region[b][1]
        other_index = findfirst(i -> i != bad_index && any(item -> item[1] == shared,new_regions[i]),
                                eachindex(new_regions))
        first_new = next_segment
        second_new = next_segment+1
        next_segment += 2
        if bad_region[b][1] isa Int
            piece_keys=collect(keys(pieces))
            piece_index=findfirst(key -> bad_region[b][1] in pieces[key],piece_keys)
            piece_index !== nothing && push!(pieces[piece_keys[piece_index]],second_new)
        else
            push!(pieces[bad_region[b][1]],second_new)
        end

        if a < b
            first_region = vcat(bad_region[1:a-1],[(bad_region[a][1],0),(first_new,1)],bad_region[b:end])
            second_region = vcat(bad_region[a+1:b-1],[(second_new,1),(first_new,1)])
        else
            first_region = vcat(bad_region[b:a-1],[(bad_region[a][1],0),(first_new,1)])
            second_region = vcat(bad_region[1:b-1],[(second_new,1),(first_new,1)],bad_region[a+1:end])
        end
        if other_index !== nothing
            other = new_regions[other_index]
            c = findfirst(item -> item[1] == shared,other)
            original_turn = other[c][2]
            other[c] = (other[c][1],0)
            insert!(other,c+1,(second_new,original_turn))
        end
        deleteat!(new_regions,bad_index)
        push!(new_regions,first_region,second_region)
    end

    all_segments = Any[segment for segments in values(segments_by_edge) for segment in segments]
    for region in new_regions, item in region
        item[1] in all_segments || push!(all_segments,item[1])
    end
    model = Model(GLPK.Optimizer)
    set_silent(model)
    variables = Dict{Any,VariableRef}()
    for segment in all_segments
        variables[segment] = @variable(model,lower_bound=1,integer=true)
    end
    max_direction=0
    legacy_dropped=Tuple{Int,Any,Int}[]
    for (region_index,region) in enumerate(new_regions)
        positive_horizontal = VariableRef[]
        negative_horizontal = VariableRef[]
        positive_vertical = VariableRef[]
        negative_vertical = VariableRef[]
        direction = 0
        for (segment,turn) in region
            max_direction=max(max_direction,direction)
            direction>3 && push!(legacy_dropped,(region_index,segment,direction))
            # Sage applies the modulo only to the first branch and compares the
            # rest bare, so a segment reached with direction in 5..7 is dropped
            # from every bucket. This normalisation is a deliberate divergence,
            # latent on all diagrams measured so far. See `PLOT_PIPELINE.md` IS5.
            bucket = direction % 4 == 0 ? positive_horizontal :
                     direction % 4 == 1 ? positive_vertical :
                     direction % 4 == 2 ? negative_horizontal : negative_vertical
            push!(bucket,variables[segment])
            turn == 1 && (direction += 1)
        end
        @constraint(model,sum(positive_horizontal)-sum(negative_horizontal) == 0)
        @constraint(model,sum(positive_vertical)-sum(negative_vertical) == 0)
    end
    DEBUG_MODE && begin
        println(debug_io,"[piece-lengths] max_direction=",max_direction)
        println(debug_io,"[piece-lengths] MATLAB legacy-dropped-directions=",legacy_dropped)
    end
    @objective(model,Min,sum(values(variables)))
    optimize!(model)
    termination_status(model) == MOI.OPTIMAL ||
        throw(ArgumentError("orthogonal length model is $(termination_status(model))"))
    segment_lengths = Dict(segment => round(Int,value(variable)) for (segment,variable) in variables)
    Dict(piece => sum(segment_lengths[segment] for segment in segments)
         for (piece,segments) in pieces)
end

_step(point::ComplexF64,length::Real,direction::Int) =
    point + length*(direction == 0 ? 1 : direction == 1 ? im : direction == 2 ? -1 : -im)

function _route_edges(pd::Matrix{Int},orientation::Vector{Int},bending::Vector{Int},lengths)
    edge_ids = sort(unique(vec(pd)))
    edge_index = Dict(edge => index for (index,edge) in enumerate(edge_ids))
    # A stored zero means "virtual crossing": no real over/under crossing and
    # no rendered vertex. Routing nevertheless needs a local frame to interpret
    # PD slots. Sage `Link.orientation()` returns +1 for these planarisation
    # crossings, and MATLAB explicitly applies `ori_(ori_==0)=1` before its
    # slot-dependent calculations. Normalize only this local routing view;
    # `LinkLayout.orientation` keeps the original zero as topological/rendering
    # metadata. This is a frame convention, not a change of crossing type.
    routing_orientation=replace(orientation,0=>1)
    crossing_positions = fill(complex(NaN,NaN),size(pd,1))
    crossing_rotations = fill(0,size(pd,1))
    crossing_positions[1] = 0.0+0.0im
    available_crossings = collect(2:size(pd,1))
    discovered = [1]
    used_edges = Set{Int}()
    layouts = EdgeLayout[]

    while length(used_edges) < length(edge_ids)
        crossing_index = 0
        slot = 0
        for candidate in discovered
            found = findfirst(edge -> !(edge in used_edges),pd[candidate,:])
            if found !== nothing
                crossing_index = candidate
                slot = found
                break
            end
        end
        crossing_index == 0 && throw(ArgumentError("disconnected PD component is not yet positioned"))
        edge = pd[crossing_index,slot]
        push!(used_edges,edge)
        target = findfirst(i -> i != crossing_index && edge in pd[i,:],axes(pd,1))
        if target === nothing
            target = crossing_index
            target_slot = findfirst(candidate -> candidate != slot &&
                                     pd[target,candidate] == edge,axes(pd,2))
        else
            target_slot = findfirst(==(edge),pd[target,:])
        end

        route_crossing = crossing_index
        route_slot = slot
        route_target = target
        route_target_slot = target_slot
        direction = mod(crossing_rotations[route_crossing]+route_slot-1,4)
        turn = bending[edge_index[edge]] < 0 ? -1 : 1
        edge_lengths = [lengths[(edge,k)] for k in 0:abs(bending[edge_index[edge]])]
        if route_slot == 1 || (route_slot == 4 && routing_orientation[route_crossing] == 1) ||
                        (route_slot == 2 && routing_orientation[route_crossing] == -1)
            turn = -turn
            reverse!(edge_lengths)
        end

        points = ComplexF64[crossing_positions[route_crossing]]
        point = first(points)
        for length in edge_lengths
            point = _step(point,length,direction)
            push!(points,point)
            direction = mod(direction+turn,4)
        end
        last_direction = mod(direction-turn,4)

        # The minimum-bend lengths close the regions, but a cycle edge can be
        # encountered from a crossing whose rotation was fixed by a different
        # spanning-tree edge.  If that forward traversal misses an already
        # positioned target, traverse the same edge from the target slot.  The
        # reverse traversal is the direction constrained by the stored target
        # position and preserves the chosen bend count.
        if target != crossing_index && !(target in available_crossings) &&
                        last(points) != crossing_positions[target]
            route_crossing = target
            route_slot = target_slot
            route_target = crossing_index
            route_target_slot = slot
            direction = mod(crossing_rotations[route_crossing]+route_slot-1,4)
            turn = bending[edge_index[edge]] < 0 ? -1 : 1
            edge_lengths = [lengths[(edge,k)] for k in 0:abs(bending[edge_index[edge]])]
            if route_slot == 1 ||
                    (route_slot == 4 && routing_orientation[route_crossing] == 1) ||
                    (route_slot == 2 && routing_orientation[route_crossing] == -1)
                turn = -turn
                reverse!(edge_lengths)
            end
            points = ComplexF64[crossing_positions[route_crossing]]
            point = first(points)
            for length in edge_lengths
                point = _step(point,length,direction)
                push!(points,point)
                direction = mod(direction+turn,4)
            end
            last_direction = mod(direction-turn,4)
        end

        if target in available_crossings
            deleteat!(available_crossings,findfirst(==(target),available_crossings))
            push!(discovered,target)
            crossing_positions[target] = last(points)
            crossing_rotations[target] = mod(last_direction-target_slot+3,4)
        end
        push!(layouts,EdgeLayout(edge,points,route_crossing,route_target,
                                 route_slot,route_target_slot))
    end
    crossing_positions,layouts
end

"""
    orthogonal_layout(v::VirtualLink; bending_numbers=nothing, planarization=:matlab) -> LinkLayout

Construct an orthogonal layout of the closed virtual-link presentation.

# Arguments
Checked: cut components, an unsupported planarization route when `v.real_only`,
and a supplied bend-count length mismatch raise `ArgumentError`; conversion,
planarization, optimization and routing errors propagate.
Caller: satisfy L1–L3 of the
[layout stages](../specs/virtual-link-plot.md#6-layout-stage-contracts),
as applicable to real-Gauss conversion, the PD embedding and supplied bends.
`planarization` selects `:matlab` or `:layer8` only for a real-only input.

# Returns
A [LinkLayout](../specs/virtual-link-plot.md#66-linklayout) for the working
presentation, including any inserted virtual crossings. An empty PD produces
an empty layout.

# Guarantees
1. The result records the PD embedding, routed edge polylines and chosen bend
   counts in the [specified coordinates](../specs/virtual-link-plot.md#66-linklayout).
   Default bend counts solve the specified minimum-bend model.
   [proved/agent: by construction]
2. The input is not modified; the returned orientation can share its storage
   with a non-real-only input. [proved/agent: by construction]

# Not guaranteed
The drawing is not canonical, and the renderer does not run the full geometric
acceptance checks. An optimizer-selected layout is not a proof of all P1–P7
properties, nor a minimum over every embedding of the represented link.

# Sources
cf: [plot specification §6](../specs/virtual-link-plot.md#6-layout-stage-contracts) and [§11](../specs/virtual-link-plot.md#11-conformance-and-acceptance) — layout versus acceptance.

# Examples
```jldoctest
julia> layout = orthogonal_layout(VirtualLink()); (isempty(layout.edges), isempty(layout.crossings))
(true, true)
```

Role: convert
"""
function orthogonal_layout(v::VirtualLink;bending_numbers=nothing,planarization=:matlab)
    any(!isempty(component) && last(component)==0 for component in v.gauss) &&
        throw(ArgumentError("orthogonal_layout does not support open cut components; use plot_svg"))
    source=v
    if v.real_only
        gauss,orientation=planarization === :layer8 ?
            planarize_virtual_gauss_layer8(v.gauss,v.orientation) :
            planarization === :matlab ? planarize_virtual_gauss(v.gauss,v.orientation) :
            throw(ArgumentError("unsupported planarization route: $planarization"))
        source=VirtualLink()
        set_data!(source;gauss=gauss,orientation=orientation)
    end
    pd = pd_code(source)
    isempty(pd) && return LinkLayout(pd,ComplexF64[],EdgeLayout[],Int[],Int[])
    bending = bending_numbers === nothing ? minimal_bending_numbers(pd) : Int.(bending_numbers)
    length(bending) == length(unique(vec(pd))) ||
        throw(ArgumentError("one bending number is required for each PD edge"))
    lengths = _piece_lengths(pd,bending)
    crossings,edges = _route_edges(pd,source.orientation,bending,lengths)
    LinkLayout(pd,crossings,edges,bending,source.orientation)
end
