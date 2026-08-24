mutable struct URDiagram
    C::Any
    V::Vector{Int}
    E::Matrix{Int}
    W::Vector{Any}
end

URDiagram() = URDiagram(1, Int[], zeros(Int, 5, 0), Any[])
LinearAlgebra.rank(d::URDiagram) = count(==(0), d.V) + 1
next_index(d::URDiagram) = isempty(d.V) ? 1 : maximum(abs.(d.V)) + 1

function set_data_cvw!(d::URDiagram, C, V, W)
    vertices = Int.(collect(V))
    edges = sort(unique(filter(>(0), vertices)))
    weights = W isa AbstractArray ? Any[W...] : Any[W]
    length(weights) == length(edges) || throw(ArgumentError("one weight is required for each positive edge"))
    d.C, d.V, d.E, d.W = C, vertices, zeros(Int, 5, length(edges)), weights
    d.E[1, :] = edges
    d
end

function set_data_cvs!(d::URDiagram, C, V, S; epsilon=nothing)
    set_data_cvw!(d, C, V, S)
    eps = epsilon === nothing ? Any[sym("eps$(e)") for e in d.E[1, :]] :
          epsilon isa AbstractArray ? Any[epsilon...] : fill(epsilon, length(d.W))
    length(eps) == length(d.W) || throw(ArgumentError("epsilon length does not match edge count"))
    d.W = Any[d.W[i] + eps[i] for i in eachindex(d.W)]
    d
end

set_data_v!(d::URDiagram, V) = set_data_cvw!(d, 1, V, fill(NaN, length(unique(filter(>(0), V)))))

"""Construct the standard matrix-form diagram used by the trace formula."""
function set_data_cm!(d::URDiagram,C,matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    positive=[Int[] for _ in 1:rows]
    negative=[Int[] for _ in 1:rows]
    edge=0
    for column in 1:columns, row in 1:rows
        edge += 1
        push!(positive[row],edge)
        push!(negative[column],-edge)
    end
    vertices=Int[]
    for section in 1:rows
        append!(vertices,positive[section]); append!(vertices,negative[section])
        section < rows && push!(vertices,0)
    end
    set_data_cvw!(d,C,vertices,vec(matrix))
end

"""The cyclic correction matrix B_N used in the Weyl trace formula."""
function standard_matrix_adjustment(matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    rows==0 && return zeros(eltype(matrix),0,0)
    prototype=matrix[1]
    adjustment=fill(zero(prototype),rows,rows)
    rows==1 && return adjustment
    for index in 1:rows
        adjustment[index,index]=one(prototype)
        adjustment[index==1 ? rows : index-1,index]-=one(prototype)
    end
    adjustment
end

"""Closed trace formula for a standard matrix-form Weyl operator.

The cofactor expansion is deliberate: it preserves exact Laurent-polynomial
coefficients without requiring a symbolic linear-algebra package.
"""
function standard_matrix_trace(C,matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    -C/_det(matrix+standard_matrix_adjustment(matrix))
end

function edge_index(d::URDiagram, e::Integer)
    i = findfirst(==(abs(e)), vec(d.E[1, :]))
    i === nothing ? 0 : i
end
vertex_index(d::URDiagram, v::Integer) = something(findfirst(==(v), d.V), 0)
get_weight(d::URDiagram, e::Integer) = d.W[edge_index(d, e)]

function validate_data(d::URDiagram)
    positives = sort(filter(>(0), d.V))
    negatives = sort(-filter(<(0), d.V))
    positives == unique(positives) || throw(ArgumentError("positive edge labels must be unique"))
    negatives == unique(negatives) || throw(ArgumentError("negative edge labels must be unique"))
    mismatch = union(setdiff(Set(positives), Set(negatives)), setdiff(Set(negatives), Set(positives)))
    positives == negatives || throw(ArgumentError("vertices are not paired: $mismatch"))
    positives == vec(d.E[1, :]) || throw(ArgumentError("edge table does not match vertices"))
    length(d.W) == length(positives) || throw(ArgumentError("weight count does not match edges"))
    true
end

function check_data(d::URDiagram, C, V, W)
    d.C == C || throw(AssertionError("C differs"))
    d.V == collect(V) || throw(AssertionError("V differs"))
    d.W == Any[W...] || throw(AssertionError("W differs"))
    true
end

function swap!(d::URDiagram, pair; skip=false)
    a, b = Int.(pair)
    p1, p2 = vertex_index(d, a), vertex_index(d, b)
    !skip && (p1 == 0 || p2 == 0) && throw(ArgumentError("vertices must exist"))
    !skip && abs(p1 - p2) != 1 && throw(ArgumentError("vertices must be adjacent"))
    s1, s2 = sign(d.V[p1]), sign(d.V[p2])
    d.V[p1], d.V[p2] = d.V[p2], d.V[p1]
    if a == -b
        wi = edge_index(d, a); w = d.W[wi]
        if !xor(p1 < p2, s1 < 0)
            d.W[wi] = w / (1 - w); d.C = d.C / (1 - w)
        else
            d.W[wi] = w / (1 + w); d.C = d.C / (1 + w)
        end
    elseif s1 != s2
        w1, w2 = get_weight(d, a), get_weight(d, b)
        q1, q2 = vertex_index(d, -a), vertex_index(d, -b)
        nw = xor(s1 > 0, p2 < p1) ? -w1 * w2 : w1 * w2
        n = next_index(d)
        lo, hi = min(q1, q2), max(q1, q2)
        vlo = q1 < q2 ? -s1*n : -s2*n
        vhi = q1 < q2 ? -s2*n : -s1*n
        insert!(d.V, lo, vlo); insert!(d.V, hi + 1, vhi)
        d.E = hcat(d.E, zeros(Int, 5)); d.E[1, end] = n; push!(d.W, nw)
    end
    d
end

function delete_edge!(d::URDiagram, edge::Integer)
    i = edge_index(d, edge); i == 0 && throw(ArgumentError("edge does not exist"))
    d.V = filter(v -> abs(v) != abs(edge), d.V)
    d.E = d.E[:, setdiff(axes(d.E, 2), i)]; deleteat!(d.W, i); d
end

function add_edges!(d::URDiagram, pair; skip=false)
    a, b = abs.(Int.(pair))
    pp = vertex_index.(Ref(d), (a, b)); pn = vertex_index.(Ref(d), (-a, -b))
    !skip && (any(==(0), pp) || any(==(0), pn) || abs(pp[1]-pp[2]) != 1 || abs(pn[1]-pn[2]) != 1) &&
        throw(ArgumentError("both signed pairs must be adjacent"))
    for p in sort([pp[2], pn[2]], rev=true); deleteat!(d.V, p); end
    ia, ib = edge_index(d, a), edge_index(d, b)
    d.W[ia] = d.W[ia] + d.W[ib]
    d.E = d.E[:, setdiff(axes(d.E, 2), ib)]; deleteat!(d.W, ib); d
end

function dilation!(d::URDiagram, dil::Integer, arg::Integer; skip=false)
    dil = abs(dil); pos = [vertex_index(d, dil), vertex_index(d, -dil), vertex_index(d, arg)]
    !skip && (any(==(0), pos) || sort(pos) != collect(minimum(pos):maximum(pos))) &&
        throw(ArgumentError("dilation vertices and target must be adjacent"))
    e1, e2, e3 = pos[2]-pos[1], sign(arg), sign(pos[3]-pos[1])
    w = get_weight(d, dil); mul = e1 < 0 ? w/(1-w) : w+1
    exponent=e2*e3
    powered=exponent<0 && mul isa Integer ? 1/(mul^(-exponent)) : mul^exponent
    d.W[edge_index(d, arg)] = get_weight(d, arg) * powered
    old = copy(d.V[pos]); shift = e3 .* [1, 1, -2]
    d.V[pos .+ shift] = old; d
end

function compose_dilation!(d::URDiagram, pair; skip=false)
    a, b = abs.(Int.(pair)); pos = collect(vertex_index.(Ref(d), (a, b, -a, -b)))
    !skip && (any(==(0), pos) || sort(pos) != collect(minimum(pos):maximum(pos))) &&
        throw(ArgumentError("dilations must form one adjacent block"))
    # MATLAB `composeDilation`: first normalize the two possible crossed
    # orderings, then choose the weight law from the endpoint directions.
    direction=(pos[1]>pos[3],pos[2]>pos[4])
    if direction[1] && !direction[2]
        swap!(d,(d.V[pos[1]],d.V[pos[3]]);skip=skip)
    elseif direction[2] && !direction[1]
        swap!(d,(d.V[pos[2]],d.V[pos[4]]);skip=skip)
    end
    ia, ib = edge_index(d, a), edge_index(d, b); w1, w2 = d.W[ia], d.W[ib]
    d.W[ia] = all(direction) ? 1-(1-w1)*(1-w2) : (1+w1)*(1+w2)-1
    for p in sort([pos[2], pos[4]], rev=true); deleteat!(d.V, p); end
    d.E = d.E[:, setdiff(axes(d.E, 2), ib)]; deleteat!(d.W, ib); d
end

"""Insert a zero separator and its dilation edge after two vertex positions.

`positive_position` and `negative_position` are one-based positions in the
current vertex word, matching MATLAB `URDiagram.put0`.
"""
function put_zero!(d::URDiagram, positive_position::Integer, negative_position::Integer)
    0 <= positive_position <= length(d.V) || throw(BoundsError(d.V,positive_position))
    0 <= negative_position <= length(d.V) || throw(BoundsError(d.V,negative_position))
    new_edge=next_index(d)
    if positive_position <= negative_position
        insert!(d.V,positive_position+1,new_edge)
        insert!(d.V,negative_position+2,-new_edge)
    else
        insert!(d.V,negative_position+1,-new_edge)
        insert!(d.V,positive_position+2,new_edge)
    end
    d.E=hcat(d.E,zeros(Int,5)); d.E[1,end]=new_edge; push!(d.W,0)
    d
end

function _reduct_add!(d::URDiagram)
    changed = false
    for i in 1:max(length(d.V)-1, 0)
        a, b = d.V[i], d.V[i+1]
        if a > 0 && b > 0
            na, nb = vertex_index(d, -a), vertex_index(d, -b)
            if abs(na-nb) == 1
                add_edges!(d, (a, b)); changed = true; break
            end
        end
    end
    changed
end

function _groups(d::URDiagram)
    boundaries=vcat(0,findall(!=(0),diff(sign.(d.V))))
    group(position)=count(boundary -> position > boundary,boundaries)
    ([group(vertex_index(d,edge)) for edge in vec(d.E[1,:])],
     [group(vertex_index(d,-edge)) for edge in vec(d.E[1,:])])
end

"""MATLAB reduct3: merge edges whose positive and negative ends share groups."""
function reduct3!(d::URDiagram)
    positive_groups,negative_groups=_groups(d)
    deleted=Set{Int}()
    pairs=Tuple{Int,Int}[]
    edges=vec(d.E[1,:])
    for left in eachindex(edges), right in left+1:length(edges)
        right in deleted && continue
        if positive_groups[left] == positive_groups[right] &&
           negative_groups[left] == negative_groups[right]
            push!(pairs,(edges[left],edges[right])); push!(deleted,right)
        end
    end
    for pair in pairs
        add_edges!(d,pair;skip=true)
    end
    d
end

"""Push the supplied separator dilations through adjacent chords (MATLAB reduct4)."""
function reduct4!(d::URDiagram, dilation_edges)
    dilations=Int.(collect(dilation_edges))
    positive_groups,negative_groups=_groups(d)
    edge_labels=vec(d.E[1,:])
    dilation_indices=Set(edge_index.(Ref(d),abs.(dilations)))
    # MATLAB `reduct4` computes both the candidates and their directions once,
    # before any swap.  Reclassifying later edges after an earlier swap skips
    # source-selected reductions and changes the resulting trace.
    candidates=[index for index in eachindex(edge_labels)
                if index ∉ dilation_indices && abs(positive_groups[index]-negative_groups[index])==1]
    composed=Tuple{Int,Int}[]
    for index in candidates
        edge=edge_labels[index]
        direction=positive_groups[index]-negative_groups[index]
        positive,negative=vertex_index(d,edge),vertex_index(d,-edge)
        between=direction==1 ? d.V[negative+1:positive-1] : d.V[positive+1:negative-1]
        positive>negative && (between=reverse(between))
        for other in filter(>(0),between)
            swap!(d,(edge,other))
        end
        for other in reverse(filter(<(0),between))
            swap!(d,(-edge,other))
        end
        direction==1 && swap!(d,(edge,-edge))
        section=count(==(0),d.V[1:vertex_index(d,edge)])+1
        section<=length(dilations) || throw(ArgumentError("missing separator dilation for rank block $section"))
        dilation=dilations[section]
        negative_position,separator_position=vertex_index(d,-edge),vertex_index(d,dilation)
        for other in d.V[negative_position+1:separator_position-1]
            dilation!(d,edge,other)
        end
        push!(composed,(dilation,edge))
    end
    for pair in reverse(composed)
        compose_dilation!(d,pair)
    end
    d
end

reduct4!(d::URDiagram,dilation_edge::Integer)=reduct4!(d,[dilation_edge])

"""One-rank separator-dilation reduction (MATLAB reduct5)."""
function reduct5!(d::URDiagram)
    put_zero!(d,length(d.V),length(d.V))
    reduct4!(d,d.V[end-1])
end

"""Resolve one sign inversion, matching MATLAB reduct6."""
function reduct6!(d::URDiagram)
    seen_negative=false
    for index in eachindex(d.V)
        if d.V[index] < 0
            seen_negative=true
        elseif d.V[index] > 0 && seen_negative
            swap!(d,d.V[[index,index-1]])
            break
        end
    end
    d
end

"""Resolve at most one inversion in every zero-separated rank block."""
function reduct7!(d::URDiagram)
    seen_negative=false; changed=false; complete=Bool[]
    # MATLAB uses `for idx=1:length(obj.V)`: the upper bound is captured
    # before a swap inserts its auxiliary edge pair.  Traversing a growing
    # Julia vector in the same pass performs extra inversions and makes the
    # multi-rank reduction diverge from the source implementation.
    initial_length=length(d.V)
    for index in 1:initial_length
        vertex=d.V[index]
        if vertex < 0
            seen_negative=true
        elseif vertex > 0 && seen_negative && !changed
            swap!(d,d.V[[index,index-1]])
            changed=true
        elseif vertex == 0
            push!(complete,!changed); seen_negative=false; changed=false
        end
    end
    push!(complete,!changed)
    d,complete
end

"""Multi-rank separator-dilation reduction (MATLAB reduct8)."""
function reduct8!(d::URDiagram)
    positions=vcat(findall(==(0),d.V).-1,length(d.V))
    for index in eachindex(positions)
        positions[index]+=2*(index-1)
        put_zero!(d,positions[index],positions[index])
    end
    dilations=Int[d.V[position+1] for position in positions]
    reduct4!(d,dilations)
    previous=vcat(-2,vertex_index.(Ref(d),dilations))
    snapshot=copy(d.V)
    for index in eachindex(dilations)
        start=previous[index]+3
        stop=previous[index+1]-1
        start>stop && continue
        for vertex in reverse(snapshot[start:stop])
            vertex>0 && break
            dilation!(d,dilations[index],vertex)
        end
    end
    d
end

function simplify1!(d::URDiagram;maxiter=50)
    count=0
    while length(d.V)>2 && count<maxiter
        count+=1
        reduct3!(d); reduct5!(d); reduct6!(d)
    end
    d
end

function _simplify2!(d::URDiagram;maxiter=50,history=nothing)
    count=0
    while count<maxiter
        count+=1
        reduct3!(d); reduct8!(d)
        _,complete=reduct7!(d)
        history === nothing || push!(history,(iteration=count,vertices=copy(d.V),complete=copy(complete)))
        all(complete) && break
    end
    reduct3!(d)
    d
end

function simplify2!(d::URDiagram;maxiter=50)
    _simplify2!(d;maxiter=maxiter)
end

"""Run `simplify2!` while retaining structural state after each source iteration."""
function simplify2_history!(d::URDiagram;maxiter=50)
    history=NamedTuple[]
    _simplify2!(d;maxiter=maxiter,history=history)
    d,history
end

function simplify!(d::URDiagram; maxiter=50)
    simplify2!(d;maxiter=maxiter)
end

function _det(A::AbstractMatrix)
    n = size(A, 1); n == 0 && return 1; n == 1 && return A[1,1]
    acc = 0
    for j in 1:n
        minor = A[2:end, setdiff(1:n, j)]
        acc = acc + (isodd(j) ? 1 : -1) * A[1,j] * _det(minor)
    end
    acc
end

function _matrix_display(d::URDiagram)
    n = rank(d); M = Matrix{Any}(undef, n, n); fill!(M, 0)
    section(v) = count(==(0), d.V[1:vertex_index(d, v)]) + 1
    for (i, e) in enumerate(vec(d.E[1, :]))
        M[section(e), section(-e)] = d.W[i]
    end
    M
end

"""Rank-one block matrix used by MATLAB `URDiagram.trace3` and `simplify3`."""
function _block_matrix(d::URDiagram)
    rank(d)==1 || throw(ArgumentError("block matrix is defined only for rank-one diagrams"))
    block=1; active=true; blocks=Int[]
    for vertex in d.V
        !active && vertex>0 && (active=true; block+=1)
        vertex<0 && (active=false)
        push!(blocks,block)
    end
    matrix=Matrix{Any}(undef,block,block); fill!(matrix,0)
    for (index,edge) in enumerate(vec(d.E[1,:]))
        matrix[blocks[vertex_index(d,edge)],blocks[vertex_index(d,-edge)]]=d.W[index]
    end
    matrix
end

trace2(d::URDiagram) = d.C / _det(-_matrix_display(d))
trace1(d::URDiagram) = length(d.W) == 1 ? -d.C/d.W[1] : throw(ArgumentError("trace1 requires one edge"))
trace3(d::URDiagram) = rank(d) == 1 ? -d.C/_det(_block_matrix(d)+standard_matrix_adjustment(_block_matrix(d))) :
    throw(ArgumentError("trace3 requires rank one"))
trace_value(d::URDiagram) = (simplify!(d); trace2(d))

function odata(d::URDiagram)
    all(w -> w isa Number && abs(w) == 1, d.W) || throw(ArgumentError("O-data requires weights ±1"))
    parts = String[]
    for v in reverse(d.V)
        v == 0 ? push!(parts, "⊗") : push!(parts, "v_{$(abs(v))}^{" * (v > 0 ? "+" : "-") * (get_weight(d,v)>0 ? "r}" : "l}"))
    end
    join(parts)
end

function weight_expression(d::URDiagram)
    parts = String[]
    for v in d.V
        v == 0 ? push!(parts, "⊗") : push!(parts, (get_weight(d,v)>0 ? "" : "bar(") * "WO_$(abs(v))$(v>0 ? 2 : 1)" * (get_weight(d,v)>0 ? "" : ")"))
    end
    join(parts)
end

function unicode_art(d::URDiagram)
    io = IOBuffer()
    println(io, "C=$(d.C), rank=$(rank(d))")
    println(io, join((v == 0 ? "⊗" : v > 0 ? "+$(v)" : "$(v)" for v in d.V), "─"))
    for (e,w) in zip(vec(d.E[1,:]), d.W); println(io, "e$(e): $(w)"); end
    String(take!(io))
end

function trans_knot_complement(code)
    a = Int.(collect(code)); pos = Int[]; neg = Int[]
    for v in a
        push!(pos, 4v - (v>0 ? 2 : 0) + (v<0 ? 3 : 0))
        push!(pos, 4v - (v>0 ? 3 : 0) + (v<0 ? 1 : 0))
        push!(neg, 4v - (v>0 ? 0 : 0) + (v<0 ? 2 : 0))
        push!(neg, 4v - (v>0 ? 1 : 0) + (v<0 ? 0 : 0))
    end
    vcat(pos, 0, reverse(neg))
end
