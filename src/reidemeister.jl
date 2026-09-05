"""Immutable Reidemeister moves for `KnotGaussDiagram`.

The local rules used here are the signed Gauss-diagram rules of Ganzell et
al., *Reidemeister moves for virtual knots*, arXiv:2005.12957, section 2:
R1 has adjacent head/tail ends, R2 has adjacent heads and adjacent tails with
opposite crossing signs, and R3 has the three matched HH, TT and HT arcs with
equal 3-signs.  `KnotGaussDiagram` stores positive visits as over endpoints
and negative visits as under endpoints.  Thus the internal names `head` and
`tail` below reverse the paper's arrow terminology (Ganzell's arrow points
over-to-under, i.e. positive-to-negative).  This is a global arrow reversal:
every direction factor in the 3-sign changes sign, so the condition that all
three 3-signs are equal is unchanged.

Candidates are source-bound witnesses.  `apply_reidemeister` does not search
again, and refuses a candidate whose source diagram or local witness is stale.
"""

struct ReidemeisterCandidate
    type::String
    direction::Symbol
    vertices::Tuple{Vararg{Int}}
    edge_positions::Tuple{Vararg{NTuple{2,Int}}}
    before::Tuple{Vararg{Int}}
    after::Tuple{Vararg{Int}}
    variant::Symbol
    source::KnotGaussDiagram

    function ReidemeisterCandidate(type, direction, vertices, edge_positions,
                                   before, after, variant, source::KnotGaussDiagram)
        edges = Tuple((Int(e[1]), Int(e[2])) for e in edge_positions)
        new(String(type), Symbol(direction), Tuple(Int.(vertices)), edges,
            Tuple(Int.(before)), Tuple(Int.(after)), Symbol(variant), source)
    end
end

"""Compatibility aliases for the witness fields used by older callers."""
function Base.getproperty(c::ReidemeisterCandidate, name::Symbol)
    name === :edges && return getfield(c, :edge_positions)
    name === :local_before && return getfield(c, :before)
    name === :local_after && return getfield(c, :after)
    name === :diagram && return getfield(c, :source)
    getfield(c, name)
end

_rm_type(type) = uppercase(String(type))

function _rm_next(i::Int, length_::Int)
    length_ == 0 ? 0 : mod1(i + 1, length_)
end

function _rm_gap_edges(length_::Int)
    length_ == 0 ? [(0, 0)] : [ (i, _rm_next(i, length_)) for i in 1:length_ ]
end

function _rm_edge_start(edge::NTuple{2,Int}, length_::Int)
    length_ == 0 && return edge == (0, 0) ? 0 : nothing
    i, j = edge
    1 <= i <= length_ && j == _rm_next(i, length_) ? i : nothing
end

function _rm_clean_vertices(vertices)
    vertices === nothing && return Int[]
    Int.(collect(vertices))
end

function _rm_coerce_edges(edges)
    edges === nothing && return NTuple{2,Int}[]
    edges isa AbstractMatrix && return [Tuple(Int.(row)) for row in eachrow(edges)]
    result = NTuple{2,Int}[]
    for edge in edges
        x = Int.(collect(edge))
        length(x) == 2 || throw(DimensionMismatch("each edge must have two positions"))
        push!(result, (x[1], x[2]))
    end
    result
end

function _rm_selector_match(c::ReidemeisterCandidate, vertices, edges)
    selected = _rm_clean_vertices(vertices)
    isempty(selected) || (length(selected) == length(c.vertices) &&
        Tuple(abs.(selected)) == c.vertices) || return false
    selected_edges = _rm_coerce_edges(edges)
    isempty(selected_edges) || Tuple(selected_edges) == c.edge_positions || return false
    true
end

function _rm_direction_param(param, type::String)
    param === nothing && return :both
    direction = Symbol(lowercase(String(param)))
    direction in (:insert, :delete, :both) ||
        throw(ArgumentError("invalid $type direction parameter $param"))
    direction
end

function _rm_local(c::KnotGaussDiagram, starts::Vector{Int})
    l = length(c.gauss)
    [c.gauss[mod1(start + offset, l)] for start in starts for offset in 0:0]
end

function _rm_occurrences(d::KnotGaussDiagram, label::Int)
    findall(x -> abs(x) == label, d.gauss)
end

function _rm_relabel_after_delete(code::Vector{Int}, removed::Vector{Int}, orientation)
    removed_set = Set(removed)
    old_to_new = Dict{Int,Int}()
    next = 1
    for old in 1:length(orientation)
        old in removed_set || (old_to_new[old] = next; next += 1)
    end
    mapped = Int[sign(x) * old_to_new[abs(x)] for x in code if abs(x) ∉ removed_set]
    new_orientation = [orientation[old] for old in 1:length(orientation) if old ∉ removed_set]
    KnotGaussDiagram(mapped, new_orientation)
end

function _rm_insert_at_gaps(d::KnotGaussDiagram,
                            inserts::Dict{Int,Vector{Int}})
    l = length(d.gauss)
    l == 0 && return copy(get(inserts, 0, Int[]))
    code = Int[]
    for i in 1:l
        push!(code, d.gauss[i])
        append!(code, get(inserts, i, Int[]))
    end
    # This helper is used only with crossing labels already allocated by the
    # caller; the caller supplies the corresponding orientation separately.
    code
end

function _rm_make_insert_r1(d::KnotGaussDiagram, gap::Int, fragment::Vector{Int}, sign_::Int)
    code = _rm_insert_at_gaps(d, Dict(gap => fragment))
    KnotGaussDiagram(code, vcat(collect(d.orientation), sign_))
end

function _rm_make_insert_r2(d::KnotGaussDiagram, gap1::Int, frag1::Vector{Int},
                            gap2::Int, frag2::Vector{Int}, signs::Vector{Int})
    inserts = gap1 == gap2 ? Dict(gap1 => vcat(frag1, frag2)) :
        Dict(gap1 => frag1, gap2 => frag2)
    code = _rm_insert_at_gaps(d, inserts)
    KnotGaussDiagram(code, vcat(collect(d.orientation), signs))
end

function _rm_r1_candidates(d::KnotGaussDiagram; direction=:delete)
    result = ReidemeisterCandidate[]
    l, n = length(d.gauss), crossing_count(d)
    if direction == :delete
        l == 0 && return result
        for i in 1:l
            j = _rm_next(i, l); x, y = d.gauss[i], d.gauss[j]
            abs(x) == abs(y) && sign(x) != sign(y) || continue
            a = abs(x)
            push!(result, ReidemeisterCandidate("R1", :delete, [a], [(i, j)],
                [x, y], Int[], x > 0 ? :head_tail : :tail_head, d))
        end
    else
        for gap in _rm_gap_edges(l)
            i = first(gap)
            for order in ((1, -1), (-1, 1)), crossing_sign in (-1, 1)
                a = n + 1
                fragment = [order[1] * a, order[2] * a]
                variant = order == (1, -1) ?
                    (crossing_sign == 1 ? :head_tail_positive : :head_tail_negative) :
                    (crossing_sign == 1 ? :tail_head_positive : :tail_head_negative)
                push!(result, ReidemeisterCandidate("R1", :insert, [a], [gap],
                    Int[], fragment, variant, d))
            end
        end
    end
    result
end

function _rm_r2_candidates(d::KnotGaussDiagram; direction=:delete)
    result = ReidemeisterCandidate[]
    l, n = length(d.gauss), crossing_count(d)
    if direction == :delete
        l < 4 && return result
        for i in 1:l
            j = _rm_next(i, l); x, y = d.gauss[i], d.gauss[j]
            x > 0 && y > 0 && abs(x) != abs(y) || continue
            a, b = abs(x), abs(y)
            d.orientation[a] == -d.orientation[b] || continue
            # The other arc must be a directed adjacent tail pair.  Its order
            # may agree with, or oppose, the head arc order.
            ta = only(filter(!=(i), _rm_occurrences(d, a)))
            tb = only(filter(!=(j), _rm_occurrences(d, b)))
            (d.gauss[ta] < 0 && d.gauss[tb] < 0) || continue
            (tb == _rm_next(ta, l) || ta == _rm_next(tb, l)) || continue
            tail_start = tb == _rm_next(ta, l) ? ta : tb
            tail_edge = (tail_start, _rm_next(tail_start, l))
            variant = d.gauss[tail_start] == -a ? :parallel : :antiparallel
            push!(result, ReidemeisterCandidate("R2", :delete, [a, b],
                [(i, j), tail_edge], [x, y, d.gauss[tail_start], d.gauss[_rm_next(tail_start,l)]],
                Int[], variant, d))
        end
    else
        gaps = _rm_gap_edges(l)
        for g1i in eachindex(gaps), g2i in eachindex(gaps)
            g1, g2 = gaps[g1i], gaps[g2i]
            # The head and tail bigon arcs are ordered, so retaining both
            # ordered gap choices records both local placements.
            for head_order in ((1, 2), (2, 1)), tail_order in ((-1, -2), (-2, -1)), eps in (-1, 1)
                a, b = n + 1, n + 2
                head = [head_order[1] == 1 ? a : b, head_order[2] == 1 ? a : b]
                tail = [tail_order[1] == -1 ? -a : -b, tail_order[2] == -1 ? -a : -b]
                target = _rm_make_insert_r2(d, first(g1), head, first(g2), tail, [eps, -eps])
                variant = (head_order == (1,2) ? :head_ab : :head_ba,
                           tail_order == (-1,-2) ? :tail_ab : :tail_ba)
                push!(result, ReidemeisterCandidate("R2", :insert, [a,b], [g1,g2],
                    Int[], vcat(head, tail), Symbol(string(variant[1], "_", variant[2], eps > 0 ? "_positive" : "_negative")), d))
                # Keep target construction live for validation of the local
                # fragments without storing a second mutable object.
                target isa KnotGaussDiagram || error("internal R2 insertion failure")
                # When both arcs are placed in one source gap, the opposite
                # concatenation is a distinct rooted local witness.  The
                # after-field carries this order; its signs let application
                # and verification recover the HH/TT split without changing
                # the long-standing variant names.
                g1 == g2 || continue
                reverse_target = _rm_make_insert_r2(d, first(g1), tail,
                                                     first(g2), head, [eps, -eps])
                reverse_target isa KnotGaussDiagram ||
                    error("internal same-gap R2 insertion failure")
                push!(result, ReidemeisterCandidate(
                    "R2", :insert, [a,b], [g1,g2], Int[], vcat(tail, head),
                    Symbol(string(variant[1], "_", variant[2],
                                  eps > 0 ? "_positive" : "_negative")), d))
            end
        end
    end
    result
end

_rm_between(a::Int, b::Int, x::Int, l::Int) =
    a < b ? (a < x < b) : (x > a || x < b)

_rm_ccw_direction(tail_arc::Int, head_arc::Int) =
    mod(head_arc - tail_arc, 3) == 1 ? 1 : -1

function _rm_chords_cross(p::NTuple{2,Int}, q::NTuple{2,Int}, l::Int)
    a,b = p; c,e = q
    (_rm_between(a,b,c,l) != _rm_between(a,b,e,l))
end

"""Return (kind, chord, parity, direction, 3-sign) for a matched triple."""
function _rm_r3_data(d::KnotGaussDiagram, arcs::Vector{NTuple{2,Int}})
    l = length(d.gauss)
    values = Int[]
    for e in arcs
        append!(values, (d.gauss[e[1]], d.gauss[e[2]]))
    end
    endpoint_positions = vcat([[e[1], e[2]] for e in arcs]...)
    labels = abs.(values)
    length(unique(labels)) == 3 || return nothing
    all(count(==(x), labels) == 2 for x in unique(labels)) || return nothing
    all(abs(values[2k-1]) != abs(values[2k]) for k in 1:3) || return nothing
    kinds = map(1:3) do k
        a,b = values[2k-1], values[2k]
        a > 0 && b > 0 ? :HH : a < 0 && b < 0 ? :TT : :HT
    end
    sort(kinds) == [:HH,:HT,:TT] || return nothing
    # For each selected chord, identify the two selected arcs containing its
    # ends.  Increasing rooted position is counterclockwise.
    chord_arcs = Dict{Int,Vector{Int}}()
    for k in 1:3
        for label in (abs(values[2k-1]), abs(values[2k]))
            push!(get!(chord_arcs,label,Int[]), k)
        end
    end
    all(length(v) == 2 && v[1] != v[2] for v in Base.values(chord_arcs)) || return nothing
    verts = sort(collect(keys(chord_arcs)))
    data = NamedTuple[]
    for label in verts
        pair = Tuple(chord_arcs[label])
        endpoints = Tuple(endpoint_positions[k] for (k,v) in enumerate(values) if abs(v) == label)
        parity = 1
        for other in verts
            other == label && continue
            op = Tuple(endpoint_positions[k] for (k,v) in enumerate(values) if abs(v) == other)
            _rm_chords_cross((endpoints[1],endpoints[2]), (op[1],op[2]), l) && (parity = -parity)
        end
        occurrence = findall(v -> abs(v) == label, values)
        tail_index = only(filter(k -> values[k] < 0, occurrence))
        head_index = only(filter(k -> values[k] > 0, occurrence))
        tail_arc = cld(tail_index,2); head_arc = cld(head_index,2)
        direction = _rm_ccw_direction(tail_arc, head_arc)
        threesign = d.orientation[label] * parity * direction
        push!(data, (label=label, parity=parity, direction=direction, threesign=threesign, arcs=pair))
    end
    all(x -> x.threesign == first(data).threesign, data) || return nothing
    (values=values, labels=verts, kinds=kinds, data=data)
end

function _rm_r3_candidates(d::KnotGaussDiagram)
    result = ReidemeisterCandidate[]
    l = length(d.gauss)
    l < 6 && return result
    # Select any three disjoint adjacent endpoint pairs; unrelated endpoints
    # may lie between the selected arcs.
    for i in 1:l, j in i+1:l, k in j+1:l
        starts = [i,j,k]
        arcs = [(s,_rm_next(s,l)) for s in starts]
        positions = vcat([[a[1],a[2]] for a in arcs]...)
        length(unique(positions)) == 6 || continue
        data = _rm_r3_data(d, arcs)
        data === nothing && continue
        # Canonicalize the same unordered set of arcs once.
        before = Tuple(d.gauss[p] for p in positions)
        after = collect(before)
        for q in 1:3; after[2q-1],after[2q] = after[2q],after[2q-1]; end
        push!(result, ReidemeisterCandidate("R3", :self_inverse, data.labels, arcs,
            before, after, Symbol(join(string.(data.kinds), "_")), d))
    end
    result
end

"""Check that a public witness carries the variant implied by local data."""
function _rm_variant_matches(d::KnotGaussDiagram, c::ReidemeisterCandidate)
    t = _rm_type(c.type)
    if t == "R1" && c.direction == :delete
        length(c.before) == 2 || return false
        return c.variant == (c.before[1] > 0 ? :head_tail : :tail_head)
    elseif t == "R1" && c.direction == :insert
        length(c.after) == 2 || return false
        order = c.after[1] > 0 ? :head_tail : :tail_head
        suffix = occursin("positive", String(c.variant)) ? :positive :
            occursin("negative", String(c.variant)) ? :negative : nothing
        suffix === nothing && return false
        return sign(c.after[1]) != sign(c.after[2]) &&
            c.variant == Symbol(string(order, "_", suffix))
    elseif t == "R2" && c.direction == :delete
        length(c.before) == 4 || return false
        return c.variant == (c.before[3] == -c.vertices[1] ? :parallel : :antiparallel)
    elseif t == "R2" && c.direction == :insert
        length(c.after) == 4 && length(c.vertices) == 2 || return false
        a,b = c.vertices
        # Distinct gaps have a fixed role order (HH at edge 1, TT at edge 2).
        # For one gap, both concatenations are valid; a negative first pair
        # identifies the TT-then-HH witness.
        head, tail = if length(c.edge_positions) == 2 &&
                         c.edge_positions[1] == c.edge_positions[2] &&
                         c.after[1] < 0
            c.after[3:4], c.after[1:2]
        else
            c.after[1:2], c.after[3:4]
        end
        ho = Tuple(abs.(head)) == (a,b) ? :head_ab : Tuple(abs.(head)) == (b,a) ? :head_ba : nothing
        to = Tuple(abs.(tail)) == (a,b) ? :tail_ab : Tuple(abs.(tail)) == (b,a) ? :tail_ba : nothing
        (ho === nothing || to === nothing) && return false
        suffix = occursin("positive", String(c.variant)) ? :positive :
            occursin("negative", String(c.variant)) ? :negative : nothing
        suffix === nothing && return false
        c.variant == Symbol(string(ho, "_", to, "_", suffix))
    elseif t == "R3" && c.direction == :self_inverse
        length(c.edge_positions) == 3 || return false
        values = Int[]
        for e in c.edge_positions
            i = _rm_edge_start(e, length(d.gauss)); i === nothing && return false
            append!(values, (d.gauss[i], d.gauss[_rm_next(i,length(d.gauss))]))
        end
        kinds = map(1:3) do k
            values[2k-1] > 0 && values[2k] > 0 ? :HH :
                values[2k-1] < 0 && values[2k] < 0 ? :TT : :HT
        end
        c.variant == Symbol(join(string.(kinds), "_"))
    else
        false
    end
end

"""Return all source-bound R1/R2/R3 witnesses matching the selectors."""
function movable(d::KnotGaussDiagram, type; param=nothing, vertices=Int[], edges=nothing, strict=false)
    t = _rm_type(type)
    candidates = if t == "R1"
        direction = _rm_direction_param(param,t)
        direction == :insert ? _rm_r1_candidates(d; direction=:insert) :
            direction == :delete ? _rm_r1_candidates(d; direction=:delete) :
            vcat(_rm_r1_candidates(d; direction=:delete), _rm_r1_candidates(d; direction=:insert))
    elseif t == "R2"
        direction = _rm_direction_param(param,t)
        direction == :insert ? _rm_r2_candidates(d; direction=:insert) :
            direction == :delete ? _rm_r2_candidates(d; direction=:delete) :
            vcat(_rm_r2_candidates(d; direction=:delete), _rm_r2_candidates(d; direction=:insert))
    elseif t == "R3"
        param === nothing || Symbol(lowercase(String(param))) in (:both,:self_inverse) ||
            throw(ArgumentError("invalid R3 direction parameter $param"))
        _rm_r3_candidates(d)
    else
        throw(ArgumentError("unknown Reidemeister move type $type"))
    end
    filter(c -> _rm_selector_match(c, vertices, edges), candidates)
end

function movable(d::KnotGaussDiagram, type, arg::NamedTuple)
    param = haskey(arg, :param) ? arg.param : haskey(arg, :p) ? arg.p : nothing
    vertices = haskey(arg, :vertices) ? arg.vertices : haskey(arg, :v) ? arg.v : Int[]
    edges = haskey(arg, :edges) ? arg.edges : haskey(arg, :e) ? arg.e : nothing
    strict = haskey(arg, :strict) ? arg.strict : false
    movable(d, type; param=param, vertices=vertices, edges=edges, strict=strict)
end

function _rm_check_source(c::ReidemeisterCandidate, d::KnotGaussDiagram)
    c.source == d || throw(ArgumentError("stale Reidemeister candidate: source diagram differs"))
    _rm_type(c.type) in ("R1", "R2", "R3") || throw(ArgumentError("invalid candidate type"))
    c.direction in (:insert, :delete, :self_inverse) || throw(ArgumentError("invalid candidate direction"))
    all(length(e) == 2 for e in c.edge_positions) || throw(ArgumentError("invalid candidate edge witness"))
    true
end

function _rm_apply_raw(d::KnotGaussDiagram, c::ReidemeisterCandidate)
    t = _rm_type(c.type)
    if t == "R1" && c.direction == :delete
        length(c.vertices) == 1 && length(c.edge_positions) == 1 || throw(ArgumentError("invalid R1 witness"))
        _rm_variant_matches(d,c) || throw(ArgumentError("invalid R1 variant witness"))
        i, j = c.edge_positions[1]; l = length(d.gauss)
        j == _rm_next(i,l) || throw(ArgumentError("R1 edge is not adjacent"))
        x, y = d.gauss[i], d.gauss[j]
        abs(x) == abs(y) && sign(x) != sign(y) || throw(ArgumentError("R1 ends are not a head/tail pair"))
        abs(x) == c.vertices[1] || throw(ArgumentError("R1 vertex mismatch"))
        _rm_relabel_after_delete(Int[d.gauss[k] for k in 1:l if k != i && k != j],
            [c.vertices[1]], collect(d.orientation))
    elseif t == "R1" && c.direction == :insert
        length(c.vertices) == 1 && length(c.edge_positions) == 1 || throw(ArgumentError("invalid R1 insertion witness"))
        _rm_variant_matches(d,c) || throw(ArgumentError("invalid R1 variant witness"))
        gap = _rm_edge_start(c.edge_positions[1], length(d.gauss)); gap === nothing && throw(ArgumentError("invalid R1 insertion gap"))
        a = c.vertices[1]; a == crossing_count(d)+1 || throw(ArgumentError("R1 insertion label is not fresh"))
        length(c.after) == 2 && Tuple(abs.(c.after)) == (a,a) && sign(c.after[1]) != sign(c.after[2]) || throw(ArgumentError("invalid R1 insertion fragment"))
        crossing_sign = occursin("positive", String(c.variant)) ? 1 : -1
        _rm_make_insert_r1(d, gap, collect(c.after), crossing_sign)
    elseif t == "R2" && c.direction == :delete
        length(c.vertices) == 2 && length(c.edge_positions) == 2 || throw(ArgumentError("invalid R2 witness"))
        _rm_variant_matches(d,c) || throw(ArgumentError("invalid R2 variant witness"))
        a,b = c.vertices; a != b || throw(ArgumentError("R2 vertices must differ"))
        i,j = c.edge_positions[1]; k,m = c.edge_positions[2]; l=length(d.gauss)
        j == _rm_next(i,l) && m == _rm_next(k,l) || throw(ArgumentError("R2 edges are not adjacent"))
        d.gauss[i] > 0 && d.gauss[j] > 0 && d.gauss[k] < 0 && d.gauss[m] < 0 || throw(ArgumentError("R2 arcs have wrong head/tail signs"))
        [abs(d.gauss[i]),abs(d.gauss[j])] == [a,b] || throw(ArgumentError("R2 head witness mismatch"))
        Set(abs.([d.gauss[k],d.gauss[m]])) == Set([a,b]) || throw(ArgumentError("R2 tail witness mismatch"))
        d.orientation[a] == -d.orientation[b] || throw(ArgumentError("R2 crossing signs are not opposite"))
        remove = Set([i,j,k,m]); code=Int[d.gauss[q] for q in 1:l if q ∉ remove]
        _rm_relabel_after_delete(code, [a,b], collect(d.orientation))
    elseif t == "R2" && c.direction == :insert
        length(c.vertices) == 2 && length(c.edge_positions) == 2 || throw(ArgumentError("invalid R2 insertion witness"))
        _rm_variant_matches(d,c) || throw(ArgumentError("invalid R2 variant witness"))
        a,b = c.vertices; [a,b] == [crossing_count(d)+1,crossing_count(d)+2] || throw(ArgumentError("R2 labels are not fresh"))
        g1,g2=c.edge_positions; i=_rm_edge_start(g1,length(d.gauss)); k=_rm_edge_start(g2,length(d.gauss))
        (i===nothing || k===nothing) && throw(ArgumentError("invalid R2 insertion gaps"))
        length(c.after)==4 || throw(ArgumentError("invalid R2 insertion fragment"))
        head, tail = if i == k && c.after[1] < 0
            collect(c.after[3:4]), collect(c.after[1:2])
        else
            collect(c.after[1:2]), collect(c.after[3:4])
        end
        all(>(0),head) && all(<(0),tail) && Set(abs.(head))==Set([a,b]) && Set(abs.(tail))==Set([a,b]) || throw(ArgumentError("invalid R2 insertion fragment"))
        eps = occursin("positive", String(c.variant)) ? 1 : -1
        if i == k
            KnotGaussDiagram(_rm_insert_at_gaps(d, Dict(i => collect(c.after))),
                             vcat(collect(d.orientation), [eps, -eps]))
        else
            _rm_make_insert_r2(d,i,head,k,tail,[eps, -eps])
        end
    elseif t == "R3" && c.direction == :self_inverse
        length(c.edge_positions)==3 && length(c.before)==6 && length(c.after)==6 || throw(ArgumentError("invalid R3 witness"))
        _rm_variant_matches(d,c) || throw(ArgumentError("invalid R3 variant witness"))
        l=length(d.gauss); positions=Int[]
        for edge in c.edge_positions
            i=_rm_edge_start(edge,l); i===nothing && throw(ArgumentError("R3 edge is not adjacent")); push!(positions,i); push!(positions,_rm_next(i,l))
        end
        length(unique(positions)) == 6 || throw(ArgumentError("R3 arcs overlap"))
        actual=collect(d.gauss)[positions]; Tuple(actual) == c.before || throw(ArgumentError("R3 local-before witness mismatch"))
        r3data=_rm_r3_data(d, collect(c.edge_positions)); r3data === nothing && throw(ArgumentError("R3 triple fails matched/equal-3-sign criterion"))
        Tuple(r3data.labels) == c.vertices || throw(ArgumentError("R3 vertex witness mismatch"))
        labels=abs.(actual); unique(labels) |> length == 3 && all(count(==(x),labels)==2 for x in unique(labels)) || throw(ArgumentError("R3 labels do not match"))
        kinds=map(1:3) do q
            actual[2q-1]>0 && actual[2q]>0 ? :HH : actual[2q-1]<0 && actual[2q]<0 ? :TT : :HT
        end
        sort(kinds)==[:HH,:HT,:TT] || throw(ArgumentError("R3 arcs are not HH/TT/HT"))
        verts=sort(unique(labels)); Tuple(verts) == c.vertices || throw(ArgumentError("R3 vertex witness mismatch"))
        transformed=collect(d.gauss); for q in 1:3; transformed[positions[2q-1]],transformed[positions[2q]]=transformed[positions[2q]],transformed[positions[2q-1]]; end
        Tuple(transformed[positions]) == c.after || throw(ArgumentError("R3 local-after witness mismatch"))
        KnotGaussDiagram(transformed,collect(d.orientation))
    else
        throw(ArgumentError("candidate direction does not match type"))
    end
end

# Verification-side constructions.  These intentionally duplicate the small
# sequence edits instead of calling the application dispatcher or its helpers;
# a verifier must remain useful when the application path is faulty.
function _rm_verify_delete(d::KnotGaussDiagram, removed_positions::Set{Int}, removed_labels)
    old_to_new = Dict{Int,Int}(); next = 1
    removed_set = Set(Int.(removed_labels))
    for old in 1:crossing_count(d)
        old in removed_set || (old_to_new[old] = next; next += 1)
    end
    code = Int[sign(d.gauss[i]) * old_to_new[abs(d.gauss[i])] for i in eachindex(d.gauss)
               if i ∉ removed_positions]
    orientation = Int[d.orientation[i] for i in 1:crossing_count(d) if i ∉ removed_set]
    KnotGaussDiagram(code, orientation)
end

function _rm_verify_r3_condition(d::KnotGaussDiagram, arcs)
    l=length(d.gauss); values=Int[]; positions=Int[]
    for e in arcs
        1 <= e[1] <= l && e[2] == _rm_next(e[1],l) || return nothing
        append!(values,(d.gauss[e[1]],d.gauss[e[2]])); append!(positions,(e[1],e[2]))
    end
    length(unique(positions)) == 6 || return nothing
    labels=abs.(values); unique_labels=unique(labels)
    length(unique_labels)==3 && all(count(==(x),labels)==2 for x in unique_labels) || return nothing
    all(abs(values[2q-1]) != abs(values[2q]) for q in 1:3) || return nothing
    kinds=map(1:3) do q
        a,b=values[2q-1],values[2q]
        a>0 && b>0 ? :HH : a<0 && b<0 ? :TT : :HT
    end
    sort(kinds)==[:HH,:HT,:TT] || return nothing
    chord_arcs=Dict{Int,Vector{Int}}()
    for q in 1:3
        for label in (abs(values[2q-1]),abs(values[2q]))
            push!(get!(chord_arcs,label,Int[]),q)
        end
    end
    all(length(x)==2 && x[1]!=x[2] for x in Base.values(chord_arcs)) || return nothing
    labels=sort(collect(keys(chord_arcs))); three_signs=Int[]
    between(a,b,x)=a<b ? a<x<b : x>a || x<b
    for label in labels
        ep=Tuple(positions[q] for (q,v) in enumerate(values) if abs(v)==label)
        parity=1
        for other in labels
            other==label && continue
            op=Tuple(positions[q] for (q,v) in enumerate(values) if abs(v)==other)
            (between(ep[1],ep[2],op[1]) != between(ep[1],ep[2],op[2])) && (parity=-parity)
        end
        occ=findall(v->abs(v)==label,values)
        tail_arc=cld(only(filter(q->values[q]<0,occ)),2)
        head_arc=cld(only(filter(q->values[q]>0,occ)),2)
        direction = mod(head_arc-tail_arc,3)==1 ? 1 : -1
        push!(three_signs,d.orientation[label]*parity*direction)
    end
    all(==(first(three_signs)),three_signs) || return nothing
    (values=values,labels=labels,kinds=kinds)
end

function _rm_verify_insert(d::KnotGaussDiagram, inserts::Dict{Int,Vector{Int}}, new_signs)
    code = Int[]
    if isempty(d.gauss)
        append!(code, get(inserts, 0, Int[]))
    else
        for i in eachindex(d.gauss)
            push!(code, d.gauss[i])
            append!(code, get(inserts, i, Int[]))
        end
    end
    KnotGaussDiagram(code, vcat(collect(d.orientation), Int.(new_signs)))
end

"""Apply a previously returned witness to exactly its source diagram."""
function apply_reidemeister(d::KnotGaussDiagram, c::ReidemeisterCandidate)
    _rm_check_source(c,d)
    result = _rm_apply_raw(d,c)
    verify_reidemeister_transition(d,result,c) || throw(ArgumentError("candidate failed structural verification"))
    result
end

"""Apply the unique selected move (or reject ambiguity in strict mode)."""
function move(d::KnotGaussDiagram, type; param=nothing, vertices=Int[], edges=nothing, strict=true)
    candidates=movable(d,type;param=param,vertices=vertices,edges=edges,strict=strict)
    isempty(candidates) && throw(ArgumentError("No valid Reidemeister moves found"))
    strict && length(candidates) != 1 && throw(ArgumentError("Reidemeister move is ambiguous"))
    apply_reidemeister(d,first(candidates))
end

move(d::KnotGaussDiagram, type, arg::NamedTuple) = begin
    param=haskey(arg,:param) ? arg.param : haskey(arg,:p) ? arg.p : nothing
    vertices=haskey(arg,:vertices) ? arg.vertices : haskey(arg,:v) ? arg.v : Int[]
    edges=haskey(arg,:edges) ? arg.edges : haskey(arg,:e) ? arg.e : nothing
    strict=haskey(arg,:strict) ? arg.strict : true
    move(d,type;param=param,vertices=vertices,edges=edges,strict=strict)
end

function _rm_verify_variant(d::KnotGaussDiagram, c::ReidemeisterCandidate)
    t = _rm_type(c.type)
    if t == "R1" && c.direction == :delete
        length(c.before) == 2 && c.variant == (c.before[1] > 0 ? :head_tail : :tail_head)
    elseif t == "R1" && c.direction == :insert
        length(c.after) == 2 && sign(c.after[1]) != sign(c.after[2]) || return false
        order = c.after[1] > 0 ? :head_tail : :tail_head
        suffix = endswith(String(c.variant), "_positive") ? :positive :
                 endswith(String(c.variant), "_negative") ? :negative : nothing
        suffix !== nothing && c.variant == Symbol(string(order, "_", suffix))
    elseif t == "R2" && c.direction == :delete
        length(c.before) == 4 && c.variant == (c.before[3] == -c.vertices[1] ? :parallel : :antiparallel)
    elseif t == "R2" && c.direction == :insert
        length(c.after) == 4 && length(c.vertices) == 2 || return false
        a,b = c.vertices
        head, tail = if c.edge_positions[1] == c.edge_positions[2] && c.after[1] < 0
            c.after[3:4], c.after[1:2]
        else
            c.after[1:2], c.after[3:4]
        end
        ho = Tuple(abs.(head)) == (a,b) ? :head_ab : Tuple(abs.(head)) == (b,a) ? :head_ba : nothing
        to = Tuple(abs.(tail)) == (a,b) ? :tail_ab : Tuple(abs.(tail)) == (b,a) ? :tail_ba : nothing
        suffix = endswith(String(c.variant), "_positive") ? :positive :
                 endswith(String(c.variant), "_negative") ? :negative : nothing
        ho !== nothing && to !== nothing && suffix !== nothing &&
            c.variant == Symbol(string(ho, "_", to, "_", suffix))
    elseif t == "R3" && c.direction == :self_inverse
        length(c.edge_positions) == 3 || return false
        values = Int[]
        for e in c.edge_positions
            i = _rm_edge_start(e, length(d.gauss)); i === nothing && return false
            append!(values, (d.gauss[i], d.gauss[_rm_next(i,length(d.gauss))]))
        end
        kinds = map(1:3) do k
            values[2k-1] > 0 && values[2k] > 0 ? :HH :
                values[2k-1] < 0 && values[2k] < 0 ? :TT : :HT
        end
        c.variant == Symbol(join(string.(kinds), "_"))
    else
        false
    end
end

function _rm_verify_expected(d::KnotGaussDiagram,c::ReidemeisterCandidate)
    # This is deliberately a separate structural construction from the
    # application dispatcher above: it checks the recorded local rule first.
    t=_rm_type(c.type); l=length(d.gauss)
    if t=="R1" && c.direction==:delete
        i,j=c.edge_positions[1]; j==_rm_next(i,l) || return nothing
        d.gauss[i]==c.before[1] && d.gauss[j]==c.before[2] || return nothing
        Tuple(abs.(c.before))==(c.vertices[1],c.vertices[1]) || return nothing
        sign(c.before[1]) != sign(c.before[2]) || return nothing
        return _rm_verify_delete(d, Set([i,j]), [c.vertices[1]])
    elseif t=="R1" && c.direction==:insert
        i=_rm_edge_start(c.edge_positions[1],l); i===nothing && return nothing
        c.vertices==(crossing_count(d)+1,) || return nothing
        code = l == 0 ? collect(c.after) : begin
            out=Int[]
            for q in 1:l; push!(out,d.gauss[q]); q==i && append!(out,c.after); end
            out
        end
        crossing_sign = occursin("positive", String(c.variant)) ? 1 : -1
        return KnotGaussDiagram(code,vcat(collect(d.orientation),crossing_sign))
    elseif t=="R2" && c.direction==:delete
        length(c.before)==4 || return nothing
        i,j=c.edge_positions[1]; k,m=c.edge_positions[2]
        j==_rm_next(i,l) && m==_rm_next(k,l) || return nothing
        Tuple(collect(d.gauss)[[i,j,k,m]])==c.before || return nothing
        Set(abs.(c.before))==Set(c.vertices) || return nothing
        c.before[1] > 0 && c.before[2] > 0 && c.before[3] < 0 && c.before[4] < 0 || return nothing
        d.orientation[c.vertices[1]] == -d.orientation[c.vertices[2]] || return nothing
        remove=Set([i,j,k,m]); return _rm_verify_delete(d,remove,c.vertices)
    elseif t=="R2" && c.direction==:insert
        i=_rm_edge_start(c.edge_positions[1],l); k=_rm_edge_start(c.edge_positions[2],l)
        (i===nothing || k===nothing) && return nothing
        c.vertices==(crossing_count(d)+1,crossing_count(d)+2) || return nothing
        code = l == 0 ? collect(c.after) : begin
            out=Int[]
            for q in 1:l
                push!(out,d.gauss[q])
                q==i && append!(out,c.after[1:2])
                q==k && append!(out,c.after[3:4])
            end
            out
        end
        eps = occursin("positive", String(c.variant)) ? 1 : -1
        return KnotGaussDiagram(code,vcat(collect(d.orientation),eps,-eps))
    elseif t=="R3" && c.direction==:self_inverse
        length(c.edge_positions)==3 || return nothing
        positions=Int[]
        for e in c.edge_positions; i=_rm_edge_start(e,l); i===nothing && return nothing; append!(positions,[i,_rm_next(i,l)]); end
        length(unique(positions))==6 || return nothing
        Tuple(collect(d.gauss)[positions])==c.before || return nothing
        r3data=_rm_verify_r3_condition(d, collect(c.edge_positions)); r3data === nothing && return nothing
        Tuple(r3data.labels) == c.vertices || return nothing
        code=collect(d.gauss); for q in 1:3; code[positions[2q-1]],code[positions[2q]]=code[positions[2q]],code[positions[2q-1]]; end
        Tuple(code[positions]) == c.after || return nothing
        return KnotGaussDiagram(code,collect(d.orientation))
    end
    nothing
end

"""Independently certify a transition modulo cyclic shift and relabelling."""
function verify_reidemeister_transition(before::KnotGaussDiagram,
                                        after::KnotGaussDiagram,
                                        c::ReidemeisterCandidate)
    c.source == before || return false
    t = _rm_type(c.type)
    t in ("R1", "R2", "R3") || return false
    c.direction == :self_inverse && t != "R3" && return false
    c.direction == :delete && !isempty(c.after) && return false
    c.direction == :insert && !isempty(c.before) && return false
    _rm_verify_variant(before,c) || return false
    expected = _rm_verify_expected(before,c)
    expected === nothing && return false
    canonicalize(expected) == canonicalize(after) || return false
    if c.direction == :self_inverse
        c.after == Tuple(c.before[i] for i in (2,1,4,3,6,5))
    else
        true
    end
end

verify_reidemeister_transition(c::ReidemeisterCandidate,
                               before::KnotGaussDiagram, after::KnotGaussDiagram) =
    verify_reidemeister_transition(before,after,c)
