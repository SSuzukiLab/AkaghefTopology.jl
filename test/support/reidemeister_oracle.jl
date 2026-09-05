"""Independent signed-Gauss Reidemeister oracle.

The production move implementation is intentionally not used in this file.
In particular, no production move dispatcher, witness checker, or private
move helper is called while deciding whether a local witness is applicable.
The oracle works directly from the `KnotGaussDiagram` representation:

* positive visits are heads and negative visits are tails;
* a position edge is `(i,next(i))`, with `(0,0)` the unique empty-code gap;
* R1 is one adjacent head/tail pair;
* R2 is an adjacent HH pair and an adjacent TT pair, with opposite crossing
  orientations;
* R3 is a matched HH/TT/HT triple with equal sign*parity*direction.

The witness retains rooted source positions.  Consequently two witnesses that
produce the same cyclic target still have different keys when their source
edges differ.
"""

struct OracleReidemeisterWitness
    type::String
    direction::Symbol
    vertices::Tuple{Vararg{Int}}
    edge_positions::Tuple{Vararg{NTuple{2,Int}}}
    before::Tuple{Vararg{Int}}
    after::Tuple{Vararg{Int}}
    variant::Symbol
    source::KnotGaussDiagram

    function OracleReidemeisterWitness(type, direction, vertices, edges,
                                       before, after, variant,
                                       source::KnotGaussDiagram)
        normalized_edges = Tuple((Int(edge[1]), Int(edge[2])) for edge in edges)
        new(String(type), Symbol(direction), Tuple(Int.(vertices)),
            normalized_edges, Tuple(Int.(before)), Tuple(Int.(after)),
            Symbol(variant), source)
    end
end

oracle_next_position(i::Int, length_::Int) =
    length_ == 0 ? 0 : mod1(i + 1, length_)

oracle_gap_edges(length_::Int) =
    length_ == 0 ? [(0, 0)] :
    [(i, oracle_next_position(i, length_)) for i in 1:length_]

function oracle_edge_start(edge::NTuple{2,Int}, length_::Int)
    length_ == 0 && return edge == (0, 0) ? 0 : nothing
    i, j = edge
    1 <= i <= length_ && j == oracle_next_position(i, length_) ? i : nothing
end

function oracle_occurrences(d::KnotGaussDiagram, label::Int)
    findall(value -> abs(value) == label, d.gauss)
end

function oracle_type_direction(param, type::String)
    param === nothing && return :both
    direction = Symbol(lowercase(String(param)))
    direction in (:insert, :delete, :both) ||
        throw(ArgumentError("invalid $type direction parameter $param"))
    direction
end

function oracle_delete_diagram(d::KnotGaussDiagram,
                              removed_positions,
                              removed_labels)
    removed_position_set = Set(Int.(collect(removed_positions)))
    removed_label_set = Set(Int.(collect(removed_labels)))
    old_to_new = Dict{Int,Int}()
    next_label = 1
    for old_label in 1:length(d.orientation)
        old_label in removed_label_set || begin
            old_to_new[old_label] = next_label
            next_label += 1
        end
    end
    code = Int[sign(d.gauss[position]) * old_to_new[abs(d.gauss[position])]
               for position in eachindex(d.gauss)
               if position ∉ removed_position_set]
    orientation = Int[d.orientation[label]
                       for label in 1:length(d.orientation)
                       if label ∉ removed_label_set]
    KnotGaussDiagram(code, orientation)
end

function oracle_insert_code(d::KnotGaussDiagram,
                            inserts::Dict{Int,Vector{Int}})
    isempty(d.gauss) && return copy(get(inserts, 0, Int[]))
    code = Int[]
    for position in eachindex(d.gauss)
        push!(code, d.gauss[position])
        append!(code, get(inserts, position, Int[]))
    end
    code
end

function oracle_insert_diagram(d::KnotGaussDiagram,
                               inserts::Dict{Int,Vector{Int}},
                               new_orientation)
    KnotGaussDiagram(oracle_insert_code(d, inserts),
                     vcat(collect(d.orientation), Int.(new_orientation)))
end

function oracle_r1_candidates(d::KnotGaussDiagram; direction::Symbol=:delete)
    result = OracleReidemeisterWitness[]
    length_ = length(d.gauss)
    n = crossing_count(d)
    if direction == :delete
        for i in 1:length_
            j = oracle_next_position(i, length_)
            x, y = d.gauss[i], d.gauss[j]
            abs(x) == abs(y) && sign(x) != sign(y) || continue
            label = abs(x)
            push!(result, OracleReidemeisterWitness(
                "R1", :delete, [label], [(i, j)], [x, y], Int[],
                x > 0 ? :head_tail : :tail_head, d))
        end
    elseif direction == :insert
        for gap in oracle_gap_edges(length_)
            label = n + 1
            for (order, order_name) in (([label, -label], :head_tail),
                                        ([-label, label], :tail_head)),
                crossing_sign in (-1, 1)
                variant = Symbol(string(order_name, crossing_sign == 1 ?
                                        "_positive" : "_negative"))
                push!(result, OracleReidemeisterWitness(
                    "R1", :insert, [label], [gap], Int[], order, variant, d))
            end
        end
    else
        throw(ArgumentError("oracle R1 direction must be :insert or :delete"))
    end
    result
end

function oracle_r2_candidates(d::KnotGaussDiagram; direction::Symbol=:delete)
    result = OracleReidemeisterWitness[]
    length_ = length(d.gauss)
    n = crossing_count(d)
    if direction == :delete
        for i in 1:length_
            j = oracle_next_position(i, length_)
            head_values = (d.gauss[i], d.gauss[j])
            head_values[1] > 0 && head_values[2] > 0 || continue
            a, b = abs.(head_values)
            a != b || continue
            d.orientation[a] == -d.orientation[b] || continue
            a_occurrence = only(filter(position -> position != i,
                                       oracle_occurrences(d, a)))
            b_occurrence = only(filter(position -> position != j,
                                       oracle_occurrences(d, b)))
            d.gauss[a_occurrence] < 0 && d.gauss[b_occurrence] < 0 || continue
            (oracle_next_position(a_occurrence, length_) == b_occurrence ||
             oracle_next_position(b_occurrence, length_) == a_occurrence) ||
                continue
            tail_start = oracle_next_position(a_occurrence, length_) == b_occurrence ?
                         a_occurrence : b_occurrence
            tail_edge = (tail_start, oracle_next_position(tail_start, length_))
            variant = d.gauss[tail_start] == -a ? :parallel : :antiparallel
            before = [head_values[1], head_values[2],
                      d.gauss[tail_edge[1]], d.gauss[tail_edge[2]]]
            push!(result, OracleReidemeisterWitness(
                "R2", :delete, [a, b], [(i, j), tail_edge], before, Int[],
                variant, d))
        end
    elseif direction == :insert
        gaps = oracle_gap_edges(length_)
        for head_gap in gaps, tail_gap in gaps,
            head_order in ((1, 2), (2, 1)),
            tail_order in ((-1, -2), (-2, -1)),
            crossing_sign in (-1, 1)
            a, b = n + 1, n + 2
            head = [head_order[1] == 1 ? a : b,
                    head_order[2] == 1 ? a : b]
            tail = [tail_order[1] == -1 ? -a : -b,
                    tail_order[2] == -1 ? -a : -b]
            head_name = head_order == (1, 2) ? :head_ab : :head_ba
            tail_name = tail_order == (-1, -2) ? :tail_ab : :tail_ba
            variant = Symbol(string(head_name, "_", tail_name, "_",
                                    crossing_sign == 1 ? "positive" : "negative"))
            push!(result, OracleReidemeisterWitness(
                "R2", :insert, [a, b], [head_gap, tail_gap], Int[],
                vcat(head, tail), variant, d))
            # If both arcs occupy one source gap, the other concatenation is
            # a second position-specific witness.  Its local `after` field is
            # tail then head; signs determine the split during application.
            head_gap == tail_gap || continue
            push!(result, OracleReidemeisterWitness(
                "R2", :insert, [a, b], [head_gap, tail_gap], Int[],
                vcat(tail, head), variant, d))
        end
    else
        throw(ArgumentError("oracle R2 direction must be :insert or :delete"))
    end
    result
end

oracle_between(a::Int, b::Int, x::Int) =
    a < b ? (a < x < b) : (x > a || x < b)

function oracle_chords_cross(first_endpoints::NTuple{2,Int},
                             second_endpoints::NTuple{2,Int})
    a, b = first_endpoints
    c, d = second_endpoints
    oracle_between(a, b, c) != oracle_between(a, b, d)
end

"""Return independent matched-triple data, or `nothing` for a near miss."""
function oracle_r3_data(d::KnotGaussDiagram, arcs)
    length(arcs) == 3 || return nothing
    length_ = length(d.gauss)
    normalized_arcs = Tuple((Int(arc[1]), Int(arc[2])) for arc in arcs)
    positions = Int[]
    values = Int[]
    for arc in normalized_arcs
        start = oracle_edge_start(arc, length_)
        start === nothing && return nothing
        push!(positions, start)
        push!(positions, oracle_next_position(start, length_))
        push!(values, d.gauss[start])
        push!(values, d.gauss[oracle_next_position(start, length_)])
    end
    length(unique(positions)) == 6 || return nothing
    labels = abs.(values)
    unique_labels = unique(labels)
    length(unique_labels) == 3 || return nothing
    all(count(==(label), labels) == 2 for label in unique_labels) || return nothing
    all(abs(values[2k - 1]) != abs(values[2k]) for k in 1:3) || return nothing
    kinds = map(1:3) do k
        first_value, second_value = values[2k - 1], values[2k]
        first_value > 0 && second_value > 0 ? :HH :
            first_value < 0 && second_value < 0 ? :TT : :HT
    end
    sort(kinds) == [:HH, :HT, :TT] || return nothing

    chord_arcs = Dict{Int,Vector{Int}}()
    for arc_index in 1:3
        for label in (abs(values[2arc_index - 1]), abs(values[2arc_index]))
            push!(get!(chord_arcs, label, Int[]), arc_index)
        end
    end
    all(length(arc_indices) == 2 && arc_indices[1] != arc_indices[2]
        for arc_indices in Base.values(chord_arcs)) || return nothing

    ordered_labels = sort(collect(keys(chord_arcs)))
    three_signs = Int[]
    for label in ordered_labels
        endpoint_positions = Tuple(positions[index]
                                   for (index, value) in enumerate(values)
                                   if abs(value) == label)
        parity = 1
        for other in ordered_labels
            other == label && continue
            other_positions = Tuple(positions[index]
                                    for (index, value) in enumerate(values)
                                    if abs(value) == other)
            oracle_chords_cross(endpoint_positions, other_positions) &&
                (parity = -parity)
        end
        occurrences = findall(value -> abs(value) == label, values)
        tail_arc = cld(only(filter(index -> values[index] < 0, occurrences)), 2)
        head_arc = cld(only(filter(index -> values[index] > 0, occurrences)), 2)
        direction = mod(head_arc - tail_arc, 3) == 1 ? 1 : -1
        push!(three_signs, d.orientation[label] * parity * direction)
    end
    all(three_sign == first(three_signs) for three_sign in three_signs) ||
        return nothing
    (values=values, labels=ordered_labels, kinds=kinds,
     positions=positions, three_signs=three_signs)
end

function oracle_r3_candidates(d::KnotGaussDiagram)
    result = OracleReidemeisterWitness[]
    length_ = length(d.gauss)
    length_ < 6 && return result
    for i in 1:length_, j in i + 1:length_, k in j + 1:length_
        arcs = [(i, oracle_next_position(i, length_)),
                (j, oracle_next_position(j, length_)),
                (k, oracle_next_position(k, length_))]
        data = oracle_r3_data(d, arcs)
        data === nothing && continue
        before = Tuple(d.gauss[position] for position in data.positions)
        after = collect(before)
        for arc_index in 1:3
            left, right = 2arc_index - 1, 2arc_index
            after[left], after[right] = after[right], after[left]
        end
        push!(result, OracleReidemeisterWitness(
            "R3", :self_inverse, data.labels, arcs, before, after,
            Symbol(join(string.(data.kinds), "_")), d))
    end
    result
end

function oracle_candidates(d::KnotGaussDiagram, type; direction=nothing)
    move_type = uppercase(String(type))
    requested = oracle_type_direction(direction, move_type)
    if move_type == "R1"
        requested == :both ?
            vcat(oracle_r1_candidates(d; direction=:delete),
                 oracle_r1_candidates(d; direction=:insert)) :
            oracle_r1_candidates(d; direction=requested)
    elseif move_type == "R2"
        requested == :both ?
            vcat(oracle_r2_candidates(d; direction=:delete),
                 oracle_r2_candidates(d; direction=:insert)) :
            oracle_r2_candidates(d; direction=requested)
    elseif move_type == "R3"
        requested in (:both, :self_inverse) ||
            throw(ArgumentError("invalid oracle R3 direction parameter $direction"))
        oracle_r3_candidates(d)
    else
        throw(ArgumentError("unknown oracle Reidemeister move type $type"))
    end
end

function oracle_selector_match(c::OracleReidemeisterWitness, vertices, edges)
    selected_vertices = vertices === nothing ? Int[] : Int.(collect(vertices))
    isempty(selected_vertices) ||
        (length(selected_vertices) == length(c.vertices) &&
         Tuple(abs.(selected_vertices)) == c.vertices) || return false
    selected_edges = if edges === nothing
        NTuple{2,Int}[]
    elseif edges isa AbstractMatrix
        [Tuple(Int.(row)) for row in eachrow(edges)]
    else
        [begin
             values = Int.(collect(edge))
             length(values) == 2 ||
                 throw(DimensionMismatch("each edge must have two positions"))
             (values[1], values[2])
         end for edge in edges]
    end
    isempty(selected_edges) || Tuple(selected_edges) == c.edge_positions
end

function oracle_selected_candidates(d::KnotGaussDiagram, type;
                                    direction=nothing, vertices=Int[], edges=nothing)
    filter(candidate -> oracle_selector_match(candidate, vertices, edges),
           oracle_candidates(d, type; direction=direction))
end

function oracle_r2_fragments(c::OracleReidemeisterWitness)
    length(c.after) == 4 || throw(ArgumentError("invalid oracle R2 fragment"))
    if c.edge_positions[1] == c.edge_positions[2] && c.after[1] < 0
        (collect(c.after[3:4]), collect(c.after[1:2]))
    else
        (collect(c.after[1:2]), collect(c.after[3:4]))
    end
end

function oracle_sign_suffix(c::OracleReidemeisterWitness)
    text = String(c.variant)
    endswith(text, "_positive") ? 1 :
        endswith(text, "_negative") ? -1 :
        throw(ArgumentError("invalid oracle crossing-sign variant"))
end

"""Apply an oracle witness without consulting the production implementation."""
function oracle_apply(d::KnotGaussDiagram,
                       c::OracleReidemeisterWitness)
    c.source == d || throw(ArgumentError("stale oracle witness"))
    length_ = length(d.gauss)
    move_type = uppercase(c.type)
    if move_type == "R1" && c.direction == :delete
        length(c.vertices) == 1 && length(c.edge_positions) == 1 ||
            throw(ArgumentError("invalid oracle R1 deletion witness"))
        edge = c.edge_positions[1]
        start = oracle_edge_start(edge, length_)
        start === nothing && throw(ArgumentError("invalid oracle R1 edge"))
        positions = (start, oracle_next_position(start, length_))
        Tuple(d.gauss[collect(positions)]) == c.before ||
            throw(ArgumentError("oracle R1 before mismatch"))
        abs(c.before[1]) == abs(c.before[2]) &&
            sign(c.before[1]) != sign(c.before[2]) ||
            throw(ArgumentError("oracle R1 local rule mismatch"))
        abs(c.before[1]) == c.vertices[1] ||
            throw(ArgumentError("oracle R1 label mismatch"))
        return oracle_delete_diagram(d, positions, c.vertices)
    elseif move_type == "R1" && c.direction == :insert
        length(c.vertices) == 1 && length(c.edge_positions) == 1 ||
            throw(ArgumentError("invalid oracle R1 insertion witness"))
        gap = oracle_edge_start(c.edge_positions[1], length_)
        gap === nothing && throw(ArgumentError("invalid oracle R1 gap"))
        c.vertices == (crossing_count(d) + 1,) ||
            throw(ArgumentError("oracle R1 label is not fresh"))
        length(c.after) == 2 &&
            Tuple(abs.(c.after)) == (c.vertices[1], c.vertices[1]) &&
            sign(c.after[1]) != sign(c.after[2]) ||
            throw(ArgumentError("oracle R1 fragment mismatch"))
        return oracle_insert_diagram(d, Dict(gap => collect(c.after)),
                                     [oracle_sign_suffix(c)])
    elseif move_type == "R2" && c.direction == :delete
        length(c.vertices) == 2 && length(c.edge_positions) == 2 ||
            throw(ArgumentError("invalid oracle R2 deletion witness"))
        first_start = oracle_edge_start(c.edge_positions[1], length_)
        second_start = oracle_edge_start(c.edge_positions[2], length_)
        (first_start === nothing || second_start === nothing) &&
            throw(ArgumentError("invalid oracle R2 edge"))
        positions = (first_start, oracle_next_position(first_start, length_),
                     second_start, oracle_next_position(second_start, length_))
        Tuple(d.gauss[collect(positions)]) == c.before ||
            throw(ArgumentError("oracle R2 before mismatch"))
        c.before[1] > 0 && c.before[2] > 0 &&
            c.before[3] < 0 && c.before[4] < 0 ||
            throw(ArgumentError("oracle R2 endpoint kinds mismatch"))
        Set(abs.(c.before)) == Set(c.vertices) ||
            throw(ArgumentError("oracle R2 labels mismatch"))
        c.vertices[1] != c.vertices[2] ||
            throw(ArgumentError("oracle R2 labels must differ"))
        d.orientation[c.vertices[1]] == -d.orientation[c.vertices[2]] ||
            throw(ArgumentError("oracle R2 orientations mismatch"))
        return oracle_delete_diagram(d, positions, c.vertices)
    elseif move_type == "R2" && c.direction == :insert
        length(c.vertices) == 2 && length(c.edge_positions) == 2 ||
            throw(ArgumentError("invalid oracle R2 insertion witness"))
        c.vertices == (crossing_count(d) + 1, crossing_count(d) + 2) ||
            throw(ArgumentError("oracle R2 labels are not fresh"))
        head_start = oracle_edge_start(c.edge_positions[1], length_)
        tail_start = oracle_edge_start(c.edge_positions[2], length_)
        (head_start === nothing || tail_start === nothing) &&
            throw(ArgumentError("invalid oracle R2 gap"))
        head, tail = oracle_r2_fragments(c)
        all(value -> value > 0, head) && all(value -> value < 0, tail) ||
            throw(ArgumentError("oracle R2 fragment endpoint kinds mismatch"))
        Set(abs.(head)) == Set(c.vertices) && Set(abs.(tail)) == Set(c.vertices) ||
            throw(ArgumentError("oracle R2 fragment labels mismatch"))
        inserts = if head_start == tail_start
            Dict(head_start => collect(c.after))
        else
            Dict(head_start => head, tail_start => tail)
        end
        crossing_sign = oracle_sign_suffix(c)
        return oracle_insert_diagram(d, inserts, [crossing_sign, -crossing_sign])
    elseif move_type == "R3" && c.direction == :self_inverse
        length(c.edge_positions) == 3 && length(c.before) == 6 &&
            length(c.after) == 6 ||
            throw(ArgumentError("invalid oracle R3 witness"))
        positions = Int[]
        for edge in c.edge_positions
            start = oracle_edge_start(edge, length_)
            start === nothing && throw(ArgumentError("invalid oracle R3 edge"))
            append!(positions, (start, oracle_next_position(start, length_)))
        end
        length(unique(positions)) == 6 ||
            throw(ArgumentError("oracle R3 edges overlap"))
        Tuple(d.gauss[positions]) == c.before ||
            throw(ArgumentError("oracle R3 before mismatch"))
        data = oracle_r3_data(d, c.edge_positions)
        data === nothing && throw(ArgumentError("oracle R3 local rule mismatch"))
        Tuple(data.labels) == c.vertices ||
            throw(ArgumentError("oracle R3 labels mismatch"))
        code = collect(d.gauss)
        for arc_index in 1:3
            left, right = positions[2arc_index - 1], positions[2arc_index]
            code[left], code[right] = code[right], code[left]
        end
        Tuple(code[positions]) == c.after ||
            throw(ArgumentError("oracle R3 after mismatch"))
        return KnotGaussDiagram(code, collect(d.orientation))
    end
    throw(ArgumentError("invalid oracle witness type/direction"))
end

function oracle_witness_text(c)
    edges = join((string(edge[1], ":", edge[2]) for edge in c.edge_positions), ",")
    "$((c.type))|$((c.direction))|v=$(join(c.vertices, ','))|e=$edges|" *
    "b=$(join(c.before, ','))|a=$(join(c.after, ','))|$(c.variant)"
end

function oracle_multiset(candidates)
    counts = Dict{String,Int}()
    for candidate in candidates
        key = oracle_witness_text(candidate)
        counts[key] = get(counts, key, 0) + 1
    end
    counts
end
