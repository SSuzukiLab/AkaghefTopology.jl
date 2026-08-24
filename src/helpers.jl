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

struct MovePattern
    parameter::String
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    capacity::Vector{Int}
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
    for vertex in setdiff(unique(abs.(changed)),unique(paired))
        result[vertex]*=-1
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
            push!(patterns,MovePattern(label*string(number),transformed,
                transformed_orientation,collect(1:length(transformed))))
        end
    end
    patterns
end

"""Julia-native replacement for `makeMovesData.m`; no MAT-file side effects."""
function make_move_data()
    Dict(
        "MP-L" => _mp_patterns(:left),
        "MP-R" => _mp_patterns(:right),
        "02-L" => [MovePattern("dir",[Int[],Int[]],Int[],Int[])],
        "02-R" => [MovePattern("inv",[[-1,-2,2,1]],[1,-1],Int[])],
        "CP-L" => [MovePattern("dir",[[1,-2,2,-3,3,-1]],[-1,-1,-1],Int[])],
        "CP-R" => [MovePattern("inv",[[1,2,-5,3,-2,-3,-4,4,-1,5]],[1,1,1,-1,-1],Int[])],
        "PS-L" => [MovePattern("",[Int[],Int[]],Int[],Int[])],
        "PS-R" => [
            MovePattern("1",[[1,2],[-1,-2]],[1,-1],Int[]),
            MovePattern("2",[[1,2],[-2,-1]],[1,-1],Int[]),
            MovePattern("3",[[1,2],[-1,-2]],[-1,1],Int[]),
            MovePattern("4",[[1,2],[-2,-1]],[-1,1],Int[]),
        ],
        "BMP-L" => [
            MovePattern("E1",[[1,2],[-1],[-2]],[-1,1],[1,2,3]),
            MovePattern("F1",[[-1,-2],[1],[2]],[1,-1],[1,2,3]),
        ],
        "BMP-R" => [
            MovePattern("E1",[[3],[2],[1],[-1,-2,-3]],[-1,1,1],[1,3,2,0]),
            MovePattern("F1",[[-3],[-2],[-1],[1,2,3]],[1,-1,-1],[1,3,2,0]),
        ],
        "B02-L" => [MovePattern("",[Int[],Int[]],[1,-1],[1,2])],
        "B02-R" => [
            MovePattern("1",[[1],[2],[-1,-2]],[1,-1],[2,1,3]),
            MovePattern("2",[[1],[2],[-2,-1]],[-1,1],[2,1,3]),
        ],
    )
end
