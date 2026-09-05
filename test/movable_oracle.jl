using Test

if !isdefined(Main, :KnotGaussDiagram)
    include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
end
if !isdefined(Main, :ReidemeisterCandidate)
    include(joinpath(@__DIR__, "..", "src", "reidemeister.jl"))
end
if !isdefined(Main, :OracleReidemeisterWitness)
    include(joinpath(@__DIR__, "support", "reidemeister_oracle.jl"))
end

"""Convert a production witness to the independent witness representation."""
function movable_oracle_as_oracle(candidate::ReidemeisterCandidate)
    OracleReidemeisterWitness(candidate.type, candidate.direction,
                              candidate.vertices, candidate.edge_positions,
                              candidate.before, candidate.after,
                              candidate.variant, candidate.source)
end

function movable_oracle_production_candidates(d::KnotGaussDiagram,
                                              move_type::AbstractString,
                                              direction)
    direction === nothing ? movable(d, move_type) :
        movable(d, move_type; param=direction)
end

function movable_oracle_compare_candidates(d::KnotGaussDiagram,
                                           move_type::AbstractString,
                                           direction)
    production = movable_oracle_production_candidates(d, move_type, direction)
    oracle = oracle_candidates(d, move_type; direction=direction)
    production_keys = oracle_multiset(movable_oracle_as_oracle.(production))
    oracle_keys = oracle_multiset(oracle)
    production_keys == oracle_keys
end

function movable_oracle_apply_and_compare(d::KnotGaussDiagram,
                                          move_type::AbstractString,
                                          direction)
    production = movable_oracle_production_candidates(d, move_type, direction)
    oracle = oracle_candidates(d, move_type; direction=direction)
    by_key = Dict{String,Vector{ReidemeisterCandidate}}()
    for candidate in production
        key = oracle_witness_text(movable_oracle_as_oracle(candidate))
        push!(get!(by_key, key, ReidemeisterCandidate[]), candidate)
    end
    for oracle_candidate in oracle
        key = oracle_witness_text(oracle_candidate)
        choices = get(by_key, key, ReidemeisterCandidate[])
        isempty(choices) && return false
        production_candidate = popfirst!(choices)
        production_target = apply_reidemeister(d, production_candidate)
        oracle_target = oracle_apply(d, oracle_candidate)
        production_target == oracle_target || return false
        verify_reidemeister_transition(d, production_target,
                                       production_candidate) || return false
    end
    all(isempty, values(by_key))
end

function movable_oracle_inverse_exists(d::KnotGaussDiagram,
                                       move_type::AbstractString,
                                       candidate::ReidemeisterCandidate)
    target = apply_reidemeister(d, candidate)
    inverse_direction = candidate.direction == :insert ? :delete :
                        candidate.direction == :delete ? :insert : nothing
    inverse_candidates = movable_oracle_production_candidates(
        target, move_type, inverse_direction)
    any(inverse_candidate -> begin
            try
                # Deletion densely relabels surviving crossings, and insertion
                # appends fresh labels.  The inverse therefore returns the
                # same rooted signed Gauss diagram only modulo the repository's
                # cyclic/root and crossing-label equivalence.
                canonicalize(apply_reidemeister(target, inverse_candidate)) ==
                    canonicalize(d)
            catch
                false
            end
        end, inverse_candidates)
end

function movable_oracle_mapped_value(value::Int, permutation)
    sign(value) * permutation[abs(value)]
end

"""Map a witness under a cyclic root shift and/or label permutation."""
function movable_oracle_map_witness(candidate::ReidemeisterCandidate,
                                    transformed_source::KnotGaussDiagram;
                                    shift::Int=0, permutation=nothing)
    length_ = length(candidate.source.gauss)
    map_position = position -> length_ == 0 ? 0 : mod1(position + shift, length_)
    map_edge = edge -> (map_position(edge[1]), map_position(edge[2]))
    edges = [map_edge(edge) for edge in candidate.edge_positions]
    labels = permutation === nothing ? collect(candidate.vertices) :
             [permutation[label] for label in candidate.vertices]
    labels = uppercase(candidate.type) == "R3" ? sort(labels) : labels
    map_value = value -> permutation === nothing ? value :
                         movable_oracle_mapped_value(value, permutation)
    before = map(map_value, collect(candidate.before))
    after = map(map_value, collect(candidate.after))
    variant = candidate.variant

    # R3's three arcs are semantically an unordered selected set.  Production
    # stores them in rooted-position order, so restore that order after a root
    # shift and carry each local pair and kind along with its edge.
    if uppercase(candidate.type) == "R3"
        order = sortperm(first.(edges))
        edge_pairs = [(edges[q], before[2q - 1], before[2q],
                       after[2q - 1], after[2q]) for q in 1:3]
        edge_pairs = edge_pairs[order]
        edges = [pair[1] for pair in edge_pairs]
        before = reduce(vcat, [[pair[2], pair[3]] for pair in edge_pairs];
                        init=Int[])
        after = reduce(vcat, [[pair[4], pair[5]] for pair in edge_pairs];
                       init=Int[])
        kinds = split(String(candidate.variant), "_")
        variant = Symbol(join(kinds[order], "_"))
    end
    OracleReidemeisterWitness(candidate.type, candidate.direction, labels, edges,
                              before, after, variant, transformed_source)
end

function movable_oracle_transformed_multiset(candidates, transformed_source;
                                             shift::Int=0, permutation=nothing)
    oracle_multiset(movable_oracle_map_witness.(candidates,
                                                 Ref(transformed_source);
                                                 shift=shift,
                                                 permutation=permutation))
end

@testset "independence and exact witness multisets" begin
    oracle_source = read(joinpath(@__DIR__, "support", "reidemeister_oracle.jl"),
                         String)
    @test !occursin("_rm_", oracle_source)
    @test !occursin("movable(", oracle_source)
    @test !occursin("apply_reidemeister", oracle_source)
    @test !occursin("verify_reidemeister_candidate", oracle_source)

    # Every rooted normalized signed code and every orientation decoration up
    # to n=3 is small enough to compare the full candidate multiset directly.
    for n in 0:3, diagram in enumerate_knot_gauss(n)
        for move_type in ("R1", "R2")
            for direction in (:delete, :insert, :both)
                @test movable_oracle_compare_candidates(diagram, move_type,
                                                        direction)
            end
        end
        @test movable_oracle_compare_candidates(diagram, "R3", nothing)
    end
end

@testset "empty, cyclic-boundary, and same-gap witnesses" begin
    empty_diagram = KnotGaussDiagram(Int[], Int[])
    @test length(movable(empty_diagram, "R1"; param=:insert)) == 4
    @test movable_oracle_compare_candidates(empty_diagram, "R1", :insert)
    @test movable_oracle_compare_candidates(empty_diagram, "R2", :insert)

    # The empty code has one gap.  HH-TT and TT-HH are different local
    # concatenations, even though some of their cyclic targets coincide.
    same_gap = movable(empty_diagram, "R2"; param=:insert,
                       edges=[(0, 0), (0, 0)])
    @test length(same_gap) == 16
    @test length(Set(c.after for c in same_gap)) == 8
    @test any(c -> c.after[1] > 0, same_gap)
    @test any(c -> c.after[1] < 0, same_gap)
    @test movable_oracle_apply_and_compare(empty_diagram, "R2", :insert)
    @test length(Set(canonical_hash(apply_reidemeister(empty_diagram, c))
                     for c in same_gap)) < length(same_gap)

    one_crossing = KnotGaussDiagram([1, -1], [1])
    nonempty_same_gap = movable(one_crossing, "R2"; param=:insert,
                                edges=[(1, 2), (1, 2)])
    @test length(nonempty_same_gap) == 16
    @test movable_oracle_compare_candidates(one_crossing, "R2", :insert)
    @test movable_oracle_apply_and_compare(one_crossing, "R2", :insert)

    r1_boundary = KnotGaussDiagram([1, 2, -2, -1], [1, 1])
    r1_boundary_candidates = movable(r1_boundary, "R1"; param=:delete)
    @test any(c -> c.edge_positions == ((4, 1),), r1_boundary_candidates)
    @test movable_oracle_apply_and_compare(r1_boundary, "R1", :delete)

    r2_boundary = KnotGaussDiagram([1, -1, -2, 2], [1, -1])
    r2_boundary_candidates = movable(r2_boundary, "R2"; param=:delete)
    @test any(c -> c.edge_positions == ((4, 1), (2, 3)),
              r2_boundary_candidates)
    @test movable_oracle_apply_and_compare(r2_boundary, "R2", :delete)

    r3_boundary = cyclic_shift(
        KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1]), 1)
    @test any(c -> any(==( (6, 1)), c.edge_positions),
              movable(r3_boundary, "R3"))
    @test movable_oracle_apply_and_compare(r3_boundary, "R3", nothing)
end

@testset "matched R3 signs and near misses" begin
    valid_code = [-1, -2, -3, 1, 2, 3]
    valid = KnotGaussDiagram(valid_code, [1, -1, 1])
    @test length(movable(valid, "R3")) == 2
    @test movable_oracle_apply_and_compare(valid, "R3", nothing)
    data = oracle_r3_data(valid, [(1, 2), (3, 4), (5, 6)])
    @test data !== nothing
    @test sort(data.kinds) == [:HH, :HT, :TT]
    @test length(unique(data.three_signs)) == 1
    @test length(unique(valid.orientation)) == 2 # equal 3-sign is not raw sign equality

    near_misses = [
        KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, 1, 1]), # no HH/TT arc
        KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, 1, 1]), # unequal 3-sign
        KnotGaussDiagram([1, -1, 2, -2, 3, -3], [1, -1, 1]), # all HT
        KnotGaussDiagram([1, 2, -1, 3, -2, -3], [1, -1, 1]), # matched shape, wrong 3-signs
    ]
    for diagram in near_misses
        @test isempty(oracle_candidates(diagram, "R3"))
        @test isempty(movable(diagram, "R3"))
    end

    # Explicit selected triples cover overlap, non-adjacency, and wrong local
    # endpoint data without asking either matcher to discover a replacement.
    overlap = ReidemeisterCandidate("R3", :self_inverse, [1, 2, 3],
                                    [(1, 2), (2, 3), (5, 6)],
                                    [-1, -2, -2, -3, 2, 3],
                                    [-2, -1, -3, -2, 3, 2],
                                    :TT_HT_HH, valid)
    @test_throws ArgumentError apply_reidemeister(valid, overlap)
    @test !verify_reidemeister_transition(valid,
                                           KnotGaussDiagram(valid_code,
                                                            [1, -1, 1]),
                                           overlap)
    @test_throws ArgumentError oracle_apply(valid,
        OracleReidemeisterWitness(overlap.type, overlap.direction,
                                  overlap.vertices, overlap.edge_positions,
                                  overlap.before, overlap.after,
                                  overlap.variant, valid))
end

@testset "source binding, stale candidates, and immutable witnesses" begin
    diagram = KnotGaussDiagram([1, -1], [1])
    candidate = only(movable(diagram, "R1"; param=:delete,
                             edges=[(1, 2)]))
    stale = KnotGaussDiagram([-1, 1], [1])
    @test_throws ArgumentError apply_reidemeister(stale, candidate)
    @test !verify_reidemeister_transition(stale,
                                           KnotGaussDiagram(Int[], Int[]),
                                           candidate)
    @test_throws ArgumentError oracle_apply(stale,
        movable_oracle_as_oracle(candidate))

    forged_before = ReidemeisterCandidate(candidate.type, candidate.direction,
        candidate.vertices, candidate.edge_positions, [1, 1], Int[],
        candidate.variant, diagram)
    forged_variant = ReidemeisterCandidate(candidate.type, candidate.direction,
        candidate.vertices, candidate.edge_positions, candidate.before,
        candidate.after, :not_a_generated_variant, diagram)
    @test_throws ArgumentError apply_reidemeister(diagram, forged_before)
    @test_throws ArgumentError apply_reidemeister(diagram, forged_variant)
    @test !verify_reidemeister_transition(diagram,
                                           KnotGaussDiagram(Int[], Int[]),
                                           forged_before)
    @test !verify_reidemeister_transition(diagram,
                                           KnotGaussDiagram(Int[], Int[]),
                                           forged_variant)
    @test_throws MethodError setindex!(candidate.vertices, 2, 1)
    @test_throws MethodError setindex!(candidate.edge_positions, (2, 1), 1)
end

@testset "cyclic-shift and relabel equivariance" begin
    fixtures = [
        KnotGaussDiagram([1, 2, -1, -2], [1, -1]),
        KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1]),
        KnotGaussDiagram([-1, -2, 4, -3, 1, -4, 2, 3], [1, -1, 1, 1]),
    ]
    for diagram in fixtures
        move_types = diagram == fixtures[1] ? ("R2",) :
                     diagram == fixtures[2] ? ("R3",) : ("R1", "R2")
        for move_type in move_types
            direction = move_type == "R3" ? nothing : :delete
            base = movable_oracle_production_candidates(diagram, move_type,
                                                        direction)
            for shift in 0:length(diagram.gauss)-1
                shifted = cyclic_shift(diagram, shift)
                expected = movable_oracle_transformed_multiset(
                    base, shifted; shift=shift)
                actual = oracle_multiset(movable_oracle_as_oracle.(
                    movable_oracle_production_candidates(shifted, move_type,
                                                          direction)))
                @test actual == expected
            end
        end
    end

    diagram = fixtures[2]
    permutation = [3, 1, 2]
    relabelled = reindex(diagram, permutation)
    for move_type in ("R3",)
        base = movable_oracle_production_candidates(diagram, move_type, nothing)
        expected = movable_oracle_transformed_multiset(
            base, relabelled; permutation=permutation)
        actual = oracle_multiset(movable_oracle_as_oracle.(
            movable_oracle_production_candidates(relabelled, move_type, nothing)))
        @test actual == expected
    end
end

@testset "inverse closure on representative placements" begin
    fixtures = [
        KnotGaussDiagram(Int[], Int[]),
        KnotGaussDiagram([1, -1], [1]),
        KnotGaussDiagram([1, 2, -1, -2], [1, -1]),
        KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1]),
    ]
    for diagram in fixtures
        for move_type in ("R1", "R2")
            for direction in (:delete, :insert)
                candidates = movable_oracle_production_candidates(
                    diagram, move_type, direction)
                @test all(candidate -> movable_oracle_inverse_exists(
                    diagram, move_type, candidate), candidates)
            end
        end
        candidates = movable_oracle_production_candidates(diagram, "R3", nothing)
        @test all(candidate -> movable_oracle_inverse_exists(
            diagram, "R3", candidate), candidates)
    end
end
