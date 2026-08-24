mutable struct VirtualLink
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    weights::Vector{Vector{Float64}}
    real_only::Bool
    is_weighted::Bool
    strand_types::Vector{Char}
end

struct RealEdge
    id::Int
    crossing::NTuple{2,Int}
    weight::Float64
    strand::Int
    arc::Vector{Int}
end

VirtualLink() = VirtualLink(Vector{Vector{Int}}(), Int[], Vector{Vector{Float64}}(), false, false, Char[])

function _normalize_gauss(code)
    code isa AbstractVector{<:Integer} && return [Int.(code)]
    [Int.(c) for c in code]
end

function _validate_gauss(code::Vector{Vector{Int}}, ori::Vector{Int}; real_only=false)
    labels = filter(!=(0), reduce(vcat, code; init=Int[]))
    isempty(labels) && return true
    maximum(abs.(labels)) <= length(ori) || throw(ArgumentError("orientation is missing a crossing"))
    for i in 1:length(ori)
        count(==(i), abs.(labels)) == 2 || throw(ArgumentError("crossing $i must occur exactly twice"))
        real_only && ori[i] == 0 && throw(ArgumentError("real Gauss code cannot contain virtual orientation"))
    end
    true
end

"""Set Gauss, real-Gauss, O-data, PD-code, or tabular-knot input.

A `pd` matrix does not round-trip: it is rebuilt through `_pd_to_gauss` and
then re-derived by `pd_code`, which may relabel every edge. Per-edge data
carried alongside a PD matrix -- bending numbers, replay positions -- must be
indexed to what `pd_code` returns. See `PLOT_PIPELINE.md` IS4.
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
            v.gauss = _pd_to_gauss(Int.(pd))
            v.orientation = orientation === nothing ? fill(1, size(pd,1)) : Int.(orientation)
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
    for i in eachindex(cut)
        cut[i] && !isempty(v.gauss[i]) && v.gauss[i][end] != 0 && push!(v.gauss[i], 0)
    end
    v.strand_types = strand_types === nothing ? fill('o', length(v.gauss)) : collect(strand_types)
    if supplied_weights
        v.is_weighted = true
    else
        v.is_weighted = false
        v.weights = [zeros(Float64, length(filter(!=(0), g))) for g in v.gauss]
    end
    v
end

function _pd_to_gauss(pd::Matrix{Int})
    # PD edges form oriented components. Crossing sign is encoded by under/over slots.
    next = Dict{Int,Tuple{Int,Int}}()
    for (i,row) in enumerate(eachrow(pd))
        next[row[1]] = (row[3], -i)
        next[row[2]] = (row[4], i)
    end
    result = Vector{Vector{Int}}(); unseen = Set(keys(next))
    while !isempty(unseen)
        start = first(unseen); edge = start; code = Int[]
        while edge in unseen
            delete!(unseen, edge); edge, crossing = next[edge]; push!(code, crossing)
        end
        push!(result, code)
    end
    result
end

gauss_code(v::VirtualLink) = deepcopy(v.gauss)
function real_gauss_code(v::VirtualLink)
    real = findall(!=(0), v.orientation); mapping = Dict(x=>i for (i,x) in enumerate(real))
    code = [[sign(x)*mapping[abs(x)] for x in g if x != 0 && haskey(mapping, abs(x))] for g in v.gauss]
    code, v.orientation[real]
end

"""Return the stable fields of MATLAB `REdgeTable` in traversal order."""
function real_edge_table(v::VirtualLink)
    real_indices=findall(!=(0),v.orientation)
    mapping=Dict(index=>mapped for (mapped,index) in enumerate(real_indices))
    result=RealEdge[]
    edge_id=0
    for (strand,component0) in enumerate(v.gauss)
        component=filter(!=(0),component0)
        isempty(component) && continue
        real_positions=findall(position -> haskey(mapping,abs(component[position])),eachindex(component))
        for (weight_index,target_position) in enumerate(real_positions)
            previous_position=real_positions[mod1(weight_index-1,length(real_positions))]
            arc=Int[component[previous_position]]
            cursor=previous_position
            while cursor != target_position
                cursor=mod1(cursor+1,length(component))
                push!(arc,component[cursor])
            end
            convert_crossing(x)=sign(x)*mapping[abs(x)]
            crossing=(convert_crossing(first(arc)),convert_crossing(last(arc)))
            weights=strand <= length(v.weights) ? v.weights[strand] : Float64[]
            weight=weight_index <= length(weights) ? weights[weight_index] : NaN
            edge_id += 1
            push!(result,RealEdge(edge_id,crossing,weight,strand,arc))
        end
    end
    result
end

function head_map(v::VirtualLink)
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

function pd_code(v::VirtualLink)
    # Sage 10.1+ convention: start at the entering under-strand and walk
    # counter-clockwise. This is the convention used by the MATLAB oracle.
    code = [filter(!=(0), component) for component in v.gauss]
    code = filter(!isempty, code)
    isempty(code) && return zeros(Int,0,4)

    edge_at_visit = Dict{Int,NTuple{2,Int}}()
    edge = 0
    for component in code
        first_edge = edge+1
        for crossing in component
            edge += 1
            edge_at_visit[crossing] = (edge,edge+1)
        end
        last_crossing = last(component)
        incoming,outgoing = edge_at_visit[last_crossing]
        edge_at_visit[last_crossing] = (incoming,first_edge)
    end

    pd = zeros(Int,length(v.orientation),4)
    for crossing in eachindex(v.orientation)
        under = edge_at_visit[-crossing]
        over = edge_at_visit[crossing]
        if v.orientation[crossing] == -1
            pd[crossing,:] = [under[1],over[1],under[2],over[2]]
        else
            pd[crossing,:] = [under[1],over[2],under[2],over[1]]
        end
    end
    pd
end

"""Return edge => crossing-index maps for tails and heads of a PD code."""
function edge_directions(pd::AbstractMatrix{<:Integer})
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

"""Return signed, left-turning region boundaries of a PD code.

Output and ordering match Sage `Link.regions()`, which the minimum-bend MILP
consumes after a stable sort by length. Verified on the hopf link; see
`PLOT_PIPELINE.md` section 4.
"""
function regions(pd::AbstractMatrix{<:Integer})
    rows = [Int.(collect(row)) for row in eachrow(pd)]
    isempty(rows) && return Vector{Vector{Int}}()
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

function hmove_matrix(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    head=zeros(Int,vertex_count,2)
    for component in code
        isempty(component) && continue
        for index in eachindex(component)
            vertex=component[index]
            head[abs(vertex),vertex > 0 ? 2 : 1]=component[mod1(index+1,length(component))]
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

function disk_table(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    vertex_count == 0 && return (delta=zeros(Int,0,0),cp2=Int[],path_lengths=Int[])
    edges=real_edge_table(v)
    edge_count=length(edges)
    dp=zeros(Int,12); dm=zeros(Int,12)
    dp[[1,2,4,5,9,12,7,6,8,11,10,3]]=repeat(1:6,2)
    dm[[1,2,4,5,9,12,7,10,3,11,6,8]]=repeat(1:6,2)
    dictionary=Dict(1=>dp,-1=>dm)

    node_count=6vertex_count
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
        source_slots=(1:3) .+ 6 .+ 3*(source > 0)
        target_slots=(1:3) .+ 3*(target > 0)
        source_nodes=6*(abs(source)-1) .+ dictionary[orientation[abs(source)]][source_slots]
        target_nodes=6*(abs(target)-1) .+ dictionary[orientation[abs(target)]][target_slots]
        for index in 1:3
            coefficient=index == 3 ? -edge.id : edge.id
            a,b=source_nodes[index],target_nodes[index]
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
    positive_dots=[0,1,0,0,1,0]
    negative_dots=[0,0,1,0,0,1]
    cp2=fill(2,length(roots))
    for vertex in 1:vertex_count, local_node in 1:6
        row=root_index[find_root(6*(vertex-1)+local_node)]
        dots=orientation[vertex] > 0 ? positive_dots : negative_dots
        cp2[row] -= dots[local_node]
    end
    path_lengths=[count(==(root),find_root.(1:node_count)) for root in roots]
    (delta=delta,cp2=cp2,path_lengths=path_lengths)
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
    all(value -> value.den == 1,result) ? Int[value.num for value in result] : nothing
end

"""MATLAB `minimizeWeightSolution`: nearest lattice representative in a kernel coset."""
function _minimize_weight_solution(weight::Vector{Int},kernel::Matrix{Int})
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

function calculate_weight(v::VirtualLink; assign=false)
    disks=disk_table(v)
    edge_count=size(disks.delta,2)
    # All audited Live-Script cases have an integral exact particular solution.
    # Obtain it before constructing an integer-programming model: GLPK's model
    # compilation dominates this small problem and made the full regression
    # suite allocate gigabytes per repeated calculation.  The lattice
    # minimizer below is unchanged, so this preserves the selected coset
    # representative.  Keep the solver only for genuinely non-integral cases.
    base=_rational_particular(2disks.delta,-disks.cp2)
    if base === nothing
        model=Model(GLPK.Optimizer); set_silent(model)
        @variable(model,weight[1:edge_count],Int)
        @variable(model,magnitude[1:edge_count] >= 0)
        @constraint(model,magnitude .>= weight)
        @constraint(model,magnitude .>= -weight)
        for row in axes(disks.delta,1)
            @constraint(model,2sum(disks.delta[row,column]*weight[column] for column in 1:edge_count) == -disks.cp2[row])
        end
        @objective(model,Min,sum(magnitude))
        optimize!(model)
        termination_status(model) == MOI.OPTIMAL ||
            throw(ArgumentError("no integral o-graph weight exists"))
        base=round.(Int,value.(weight))
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
        weights = Vector{Vector{Float64}}(); offset = 1
        code,_=real_gauss_code(v)
        for g in code
            n=length(g); push!(weights,Float64.(base[offset:offset+n-1])); offset += n
        end
        set_weight!(v, weights)
    end
    base,kernel,selected
end

function set_weight!(v::VirtualLink, weight)
    if weight isa AbstractMatrix
        flat = vec(weight); idx = 1; weights = Vector{Vector{Float64}}()
        for g in v.gauss
            n = length(filter(!=(0),g)); push!(weights, Float64.(flat[idx:min(idx+n-1,end)])); idx += n
        end
        v.weights = weights
    elseif weight isa AbstractVector{<:Number}
        set_weight!(v, reshape(weight,1,:))
    else
        v.weights = [Float64.(w) for w in weight]
    end
    v.is_weighted = true; v
end

"""Assign selected MATLAB `REdgeTable.Crossing` weights without relying on row order."""
function set_weight_by_crossing!(v::VirtualLink,crossings,weights)
    pairs=[Tuple(Int.(pair)) for pair in eachrow(Int.(crossings))]
    length(pairs) == length(weights) || throw(DimensionMismatch("one weight per crossing pair is required"))
    assignments=Dict(pair=>Float64(weight) for (pair,weight) in zip(pairs,weights))
    table=real_edge_table(v)
    by_strand=[Float64[] for _ in v.gauss]
    for edge in table
        push!(by_strand[edge.strand],get(assignments,edge.crossing,edge.weight))
    end
    set_weight!(v,by_strand)
end

function _graph_a_connected(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    vertex_count == 0 && return length(code) <= 1
    source=VirtualLink(); set_data!(source;rgauss=code,orientation=orientation)
    pd=pd_code(source)
    parent=collect(1:3vertex_count+1)
    used=Set{Int}()
    find_root(node)=begin
        while parent[node] != node
            parent[node]=parent[parent[node]]; node=parent[node]
        end
        node
    end
    unite(a,b)=begin
        ra,rb=find_root(a),find_root(b); ra != rb && (parent[rb]=ra)
    end
    order=[4,1,2,3,4,1]
    for edge in unique(vec(pd))
        occurrences=findall(==(edge),pd)
        length(occurrences) == 2 || throw(ArgumentError("PD edge $edge must occur twice"))
        endpoints=Vector{Vector{Int}}()
        for occurrence in occurrences
            crossing=occurrence[1]; slot=occurrence[2]
            push!(endpoints,3*(crossing-1) .+ order[slot:slot+2])
            union!(used,last(endpoints))
        end
        for index in 1:3
            unite(endpoints[1][index],endpoints[2][index])
        end
    end
    length(unique(find_root.(collect(used)))) == 1
end

function is_closed(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    c1 = length(code) == 1
    c2 = _graph_a_connected(v)
    c3 = size(disk_table(v).delta,1) == length(orientation)+1
    all((c1,c2,c3)), (c1,c2,c3)
end

function mirror_manifold!(v::VirtualLink)
    v.gauss = [-g for g in v.gauss]; v.orientation .*= -1; v
end
function cyclic_shift!(v::VirtualLink, strand::Integer, amount::Integer)
    v.gauss[strand] = circshift(v.gauss[strand], amount); v
end

function _cyclic_match(component::Vector{Int},pattern::Vector{Int})
    length(pattern)<=length(component) || return nothing
    for start in eachindex(component)
        all(component[mod1(start+offset,length(component))]==pattern[offset+1]
            for offset in 0:length(pattern)-1) && return start
    end
    nothing
end

"""Apply a non-weighted 2-to-3 Matveev--Piergallini move.

This is the direct half of MATLAB `VirtualLink.move("MP",...)` used by
`C250811MSTinvUR`: the two selected real vertices are replaced by the
three-vertex template of type `A1` through `D4`.  The inverse 3-to-2 path
and weighted edge transport remain intentionally separate because MATLAB's
own implementation uses a different relabelling branch for those cases.
"""
function mp_move!(v::VirtualLink,parameter::AbstractString,vertices)
    v.is_weighted && throw(ArgumentError("weighted MP edge transport is not implemented"))
    selected=Int.(collect(vertices)); length(selected)==2 ||
        throw(ArgumentError("only direct 2-to-3 MP moves are implemented"))
    all(vertex -> 1<=vertex<=length(v.orientation),selected) || throw(BoundsError(v.orientation,selected))
    templates=make_move_data()
    before=only(filter(pattern -> pattern.parameter==parameter,templates["MP-L"]))
    after=only(filter(pattern -> pattern.parameter==parameter,templates["MP-R"]))
    v.orientation[selected]==before.orientation ||
        throw(ArgumentError("MP parameter $parameter does not match selected crossing orientations"))
    new_vertex=length(v.orientation)+1
    mapped_before=[[sign(vertex)*selected[abs(vertex)] for vertex in fragment] for fragment in before.gauss]
    mapped_after=[[sign(vertex)*(abs(vertex)==3 ? new_vertex : selected[abs(vertex)]) for vertex in fragment]
                  for fragment in after.gauss]
    matches=Tuple{Int,Int}[]
    for fragment in mapped_before
        candidates=Tuple{Int,Int}[]
        for (component_index,component) in enumerate(v.gauss)
            start=_cyclic_match(component,fragment)
            start===nothing || push!(candidates,(component_index,start))
        end
        length(candidates)==1 || throw(ArgumentError("MP fragment has $(length(candidates)) matches; move is ambiguous or invalid"))
        push!(matches,only(candidates))
    end
    # This direct splice is valid when the three template fragments occupy
    # disjoint, non-wrapping source intervals.  It covers the C250811 calls
    # and deliberately fails closed for the more general reconnection case.
    occupied=Set{Tuple{Int,Int}}()
    for (fragment,(component_index,start)) in zip(mapped_before,matches)
        component=v.gauss[component_index]
        start+length(fragment)-1<=length(component) ||
            throw(ArgumentError("wrapping MP fragments require the general reconnection path"))
        for position in start:start+length(fragment)-1
            key=(component_index,position)
            key in occupied && throw(ArgumentError("overlapping MP fragments require the general reconnection path"))
            push!(occupied,key)
        end
    end
    replacements=Dict{Tuple{Int,Int},Tuple{Int,Vector{Int}}}()
    for (index,(component_index,start)) in enumerate(matches)
        replacements[(component_index,start)]=(length(mapped_before[index]),mapped_after[index])
    end
    result=Vector{Vector{Int}}()
    for (component_index,component) in enumerate(v.gauss)
        rewritten=Int[]; position=1
        while position<=length(component)
            key=(component_index,position)
            if haskey(replacements,key)
                consumed,replacement=replacements[key]
                append!(rewritten,replacement); position+=consumed
            else
                push!(rewritten,component[position]); position+=1
            end
        end
        push!(result,rewritten)
    end
    orientation=vcat(v.orientation,0)
    for local_vertex in 1:3
        actual=local_vertex==3 ? new_vertex : selected[local_vertex]
        orientation[actual]=after.orientation[local_vertex]
    end
    v.real_only ? set_data!(v;rgauss=result,orientation=orientation) :
                  set_data!(v;gauss=result,orientation=orientation)
end
function disjoint_sum!(a::VirtualLink, b::VirtualLink)
    offset=length(a.orientation)
    append!(a.gauss, [[sign(x)*offset+x for x in g] for g in b.gauss])
    append!(a.orientation,b.orientation)
    a.real_only = true
    a.is_weighted = false
    a.weights = [zeros(Float64,length(g)) for g in a.gauss]
    a
end
function connected_sum!(v::VirtualLink, strands=(1,2))
    code,orientation=real_gauss_code(v)
    i,j=strands; i==j && throw(ArgumentError("distinct strands required"))
    1 <= i <= length(code) && 1 <= j <= length(code) || throw(BoundsError(code,(i,j)))
    vertex_count=length(orientation)
    bridge1=[vertex_count+2,vertex_count+1,-vertex_count-1,-vertex_count-4]
    bridge2=[vertex_count+4,vertex_count+3,-vertex_count-3,-vertex_count-2]
    code[i]=vcat(code[i],bridge1,code[j],bridge2)
    deleteat!(code,j)
    v.gauss=code
    v.orientation=vcat(orientation,[-1,1,-1,1])
    v.real_only=true
    v.is_weighted=false
    v.weights=[zeros(Float64,length(g)) for g in code]
    v
end

"""Convert a knot diagram to the version-2 knot-complement o-graph."""
function knot_complement!(v::VirtualLink)
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

"""Convert a VirtualLink to the URDiagram ordering used by MATLAB `VL2URD`."""
function vl_to_urdiagram(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertices=Int[]
    for (index,component) in enumerate(code)
        index > 1 && push!(vertices,0)
        append!(vertices,reverse(component))
    end
    diagram=URDiagram()
    set_data_cvs!(diagram,1,vertices,orientation)
end

function ur_invariant(v::VirtualLink)
    trace_value(vl_to_urdiagram(v))
end

"""Direct determinant form of MATLAB `calcInvariant("UR_")`."""
function ur_prime_invariant(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    length(code) == 1 || throw(ArgumentError("UR_ invariant currently requires one strand"))
    component=code[1]
    start=findfirst(>(0),component)
    start === nothing && throw(ArgumentError("the Gauss code has no positive crossing"))
    component=circshift(component,1-start)
    cumulative=cumsum(component .> 0)
    index=Dict(vertex=>cumulative[position] for (position,vertex) in enumerate(component))
    vertex_count=length(orientation)
    matrix=zeros(Int,vertex_count,vertex_count)
    for vertex in 1:vertex_count
        matrix[index[vertex],index[-vertex]]=orientation[vertex]
    end
    adjusted=matrix+Matrix{Int}(I,vertex_count,vertex_count)-circshift(Matrix{Int}(I,vertex_count,vertex_count),(1,0))
    -1/det(adjusted)
end
