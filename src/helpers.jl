mutable struct SmallObject{T}
    value::T
end

function bernstein(n::Integer, t)
    n >= 1 || throw(ArgumentError("degree must be positive"))
    ts = t isa Number ? [t] : collect(t)
    all(x -> 0 <= x <= 1, ts) || throw(ArgumentError("parameters must lie in [0,1]"))
    [binomial(n,i) * (1-x)^(n-i) * x^i for i in 0:n, x in ts]
end

function bezier(t, control_points::AbstractMatrix)
    size(control_points,2) >= 2 || throw(ArgumentError("at least two control points are required"))
    control_points * bernstein(size(control_points,2)-1,t)
end

"""One row of MATLAB's generated move-data tables.

`gauss_first`, `gauss_last`, `vertex_count`, `weights`, and `capacity` are the
fields added by `makeMovesData.m/modify`.  The four-argument constructor is
retained for source compatibility with the first Julia migration snapshot.
"""
struct MovePattern
    parameter::String
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    capacity::Vector{Int}
    gauss_first::Vector{Int}
    gauss_last::Vector{Int}
    vertex_count::Int
    weights::Vector{Vector{Float64}}
end

function MovePattern(parameter::AbstractString, gauss, orientation, capacity=Int[]; weights=nothing)
    fragments=[Int.(collect(fragment)) for fragment in gauss]
    caps=isempty(capacity) ? collect(1:length(fragments)) : Int.(collect(capacity))
    length(caps)==length(fragments) || throw(DimensionMismatch("one capacity label per Gauss fragment is required"))
    firsts=[isempty(fragment) ? 0 : first(fragment) for fragment in fragments]
    lasts=[isempty(fragment) ? 0 : last(fragment) for fragment in fragments]
    labels=abs.(reduce(vcat,fragments;init=Int[]))
    vertex_count=isempty(labels) ? 0 : maximum(labels)
    fragment_weights=weights===nothing ?
        [zeros(Float64,length(fragment)+1) for fragment in fragments] :
        [Float64.(collect(weight)) for weight in weights]
    length(fragment_weights)==length(fragments) ||
        throw(DimensionMismatch("one weight rule per Gauss fragment is required"))
    MovePattern(String(parameter),fragments,Int.(collect(orientation)),caps,
                firsts,lasts,vertex_count,fragment_weights)
end

"""Apply MATLAB `VirtualLink.calcReverseStrand` to a local move template."""
function _reverse_template(gauss,orientation,flags)
    reversed=[Int.(copy(fragment)) for fragment in gauss]
    changed=Int[]
    for index in eachindex(reversed)
        flags[index] || continue
        reversed[index]=reverse(reversed[index])
        append!(changed,reversed[index])
    end
    result=Int.(copy(orientation))
    positive=filter(>(0),changed)
    paired=filter(vertex -> -vertex in changed,positive)
    affected=setdiff(unique(abs.(changed)),unique(paired))
    virtual=filter(vertex -> 1<=vertex<=length(result) && result[vertex]==0,affected)
    for vertex in affected
        result[vertex]*=-1
    end
    # MATLAB keeps a virtual orientation at zero and swaps the two signed visits
    # when only one side of that virtual crossing belongs to a reversed strand.
    for vertex in virtual
        for fragment in reversed, index in eachindex(fragment)
            if fragment[index]==vertex
                fragment[index]=-vertex
            elseif fragment[index]==-vertex
                fragment[index]=vertex
            end
        end
    end
    reversed,result
end

function _mp_patterns(side::Symbol)
    labels=["A","B","C","D"]
    gauss=side===:left ? [
        [[-1],[1,2],[-2]],
        [[1],[-1,-2],[2]],
        [[-1,2],[1],[-2]],
        [[2,-1],[1],[-2]],
    ] : [
        [[1,-3,2],[3],[-1,-2]],
        [[-1,3,-2],[-3],[1,2]],
        [[2,-3],[1,3],[-1,-2]],
        [[-3,1],[3,2],[-1,-2]],
    ]
    orientations=side===:left ? [[1,1],[-1,-1],[-1,1],[1,-1]] :
                                 [[1,-1,1],[-1,1,-1],[1,1,-1],[-1,-1,1]]
    flags=[[false,false,false],[false,false,true],[false,true,false],[false,true,true]]
    patterns=MovePattern[]
    for (label,fragments,orientation) in zip(labels,gauss,orientations)
        for (number,reverse_flags) in enumerate(flags)
            transformed,transformed_orientation=_reverse_template(fragments,orientation,reverse_flags)
            parameter=label*string(number)
            weights=[zeros(Float64,length(fragment)+1) for fragment in transformed]
            if parameter in ("A4","B4")
                weights=[Float64[NaN] for _ in transformed]
            elseif side===:right && parameter in ("A2","B2")
                weights[3]=Float64[1,-1,0]
            elseif side===:right && parameter=="C3"
                weights[2]=Float64[0,1,-1]
            elseif side===:right && parameter=="D3"
                weights[2]=Float64[-1,1,0]
            end
            push!(patterns,MovePattern(parameter,transformed,transformed_orientation;weights=weights))
        end
    end
    patterns
end

"""Julia-native form of `makeMovesData.m`; returns data instead of writing MAT files."""
function make_move_data()
    Dict(
        "MP-L" => _mp_patterns(:left),
        "MP-R" => _mp_patterns(:right),
        "02-L" => [MovePattern("dir",[Int[],Int[]],Int[];weights=[Float64[0],Float64[0]])],
        "02-R" => [MovePattern("inv",[[-1,-2,2,1]],[1,-1];weights=[Float64[0,0,1,-1,0]])],
        "CP-L" => [MovePattern("dir",[[1,-2,2,-3,3,-1]],[-1,-1,-1])],
        "CP-R" => [MovePattern("inv",[[1,2,-5,3,-2,-3,-4,4,-1,5]],[1,1,1,-1,-1])],
        # INS24 Fig.16, independently transcribed in the source visual record.
        # S is not an alias for BP's 3-to-5 CP: its local patch is 0-to-4.
        "S-L" => [MovePattern("INS",[Int[]],Int[])],
        "S-R" => [MovePattern("INS",[[1,2,-3,-4,4,3,-2,-1]],[-1,-1,-1,-1])],
        "PS-L" => [MovePattern("",[Int[],Int[]],Int[])],
        "PS-R" => [
            MovePattern("1",[[1,2],[-1,-2]],[1,-1]),
            MovePattern("2",[[1,2],[-1,-2]],[-1,1]),
            MovePattern("3",[[1,2],[-2,-1]],[1,-1]),
            MovePattern("4",[[1,2],[-2,-1]],[-1,1]),
        ],
        "BMP-L" => [
            MovePattern("E1",[[1,2],[-1],[-2]],[-1,1],[1,2,3]),
            MovePattern("F1",[[-1,-2],[1],[2]],[1,-1],[1,2,3]),
            MovePattern("E2",[[1,2],[-1],[-2]],[-1,1],[1,2,3]),
            MovePattern("F2",[[-1,-2],[1],[2]],[1,-1],[1,2,3]),
        ],
        "BMP-R" => [
            MovePattern("E1",[[3],[2],[1],[-1,-2,-3]],[-1,1,1],[1,3,2,0]),
            MovePattern("F1",[[-3],[-2],[-1],[1,2,3]],[1,-1,-1],[1,3,2,0]),
            MovePattern("E2",[[3],[2],[1],[-3,-2,-1]],[1,-1,-1],[1,3,2,0]),
            MovePattern("F2",[[-3],[-2],[-1],[3,2,1]],[-1,1,1],[1,3,2,0]),
        ],
        "B02-L" => [MovePattern("",[Int[],Int[]],[1,-1],[1,2])],
        "B02-R" => [
            MovePattern("1",[[1],[2],[-1,-2]],[1,-1],[2,1,3]),
            MovePattern("2",[[1],[2],[-2,-1]],[-1,1],[2,1,3]),
        ],
    )
end
