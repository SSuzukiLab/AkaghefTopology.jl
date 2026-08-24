struct EdgeLayout
    id::Int
    points::Vector{ComplexF64}
    tail::Int
    head::Int
    tail_slot::Int
    head_slot::Int
end

struct LinkLayout
    pd::Matrix{Int}
    crossings::Vector{ComplexF64}
    edges::Vector{EdgeLayout}
    bending_numbers::Vector{Int}
    orientation::Vector{Int}
end

"""Solve the same minimum-bend region-flow model used by Sage's Link.plot.

This model is degenerate: several distinct bending vectors attain the same
optimal objective, so the MILP backend decides which embedding is drawn.
Matching the backend family does not remove it -- Sage/GLPK and this port
disagree on the Borromean rings.

As of 2026-08-20 that difference is **accepted**: plots are no longer accepted
by comparison with MATLAB, so any optimum is a valid answer. Tests must compare
the objective value and never the bending vector. See `PLOT_PIPELINE.md` IS1
and section 4 part A.
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

function _piece_lengths(pd::Matrix{Int},bending::Vector{Int})
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
    for region in new_regions
        positive_horizontal = VariableRef[]
        negative_horizontal = VariableRef[]
        positive_vertical = VariableRef[]
        negative_vertical = VariableRef[]
        direction = 0
        for (segment,turn) in region
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
        if route_slot == 1 || (route_slot == 4 && orientation[route_crossing] == 1) ||
                        (route_slot == 2 && orientation[route_crossing] == -1)
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
                    (route_slot == 4 && orientation[route_crossing] == 1) ||
                    (route_slot == 2 && orientation[route_crossing] == -1)
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

"""Compute a Julia-native orthogonal layout for a real Gauss diagram.

Given the same `bending_numbers`, the stages below reproduce the corner
sequence of `allocate_pos.py` exactly on the hopf link and on both C251020
figures. That agreement does **not** imply that the degenerate minimum-bend
MILP selects the same optimum for every diagram.

`bending_numbers` index the labelling returned by `pd_code`, not any PD matrix
passed to `set_data!` (IS4). Pinning them is a debugging aid, not a fix: the
hand-tuned overrides in `examples/` were chosen against the MATLAB path and
trigger IS7 here. See `PLOT_PIPELINE.md` section 4 part B.
"""
function orthogonal_layout(v::VirtualLink;bending_numbers=nothing)
    source=v
    if v.real_only
        gauss,orientation=planarize_virtual_gauss(v.gauss,v.orientation)
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
