using Test
using Random

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
include(joinpath(@__DIR__, "..", "src", "reidemeister.jl"))

@testset "Reidemeister candidate witnesses" begin
    empty = KnotGaussDiagram(Int[], Int[])
    r1i = movable(empty, "R1"; param=:insert)
    @test length(r1i) == 4
    for c in r1i
        d = apply_reidemeister(empty, c)
        @test crossing_count(d) == 1
        @test verify_reidemeister_transition(empty, d, c)
    @test all(c -> c.vertices == (1,), movable(d, "R1"; param=:delete))
    end

    r1d = KnotGaussDiagram([1, -1], [1])
    dels = movable(r1d, "R1"; param=:delete)
    @test length(dels) == 2                 # both cyclic rooted arcs
    for c in dels
        @test apply_reidemeister(r1d, c) == empty
        @test verify_reidemeister_transition(r1d, empty, c)
    end

    # R2: opposite crossing signs are required; the same-sign near miss is not.
    r2 = KnotGaussDiagram([1, 2, -1, -2], [1, -1])
    r2c = only(movable(r2, "R2"; param=:delete))
    @test r2c.vertices == (1, 2)
    @test apply_reidemeister(r2, r2c) == empty
    @test verify_reidemeister_transition(r2, empty, r2c)
    @test isempty(movable(KnotGaussDiagram([1, 2, -1, -2], [1, 1]), "R2"; param=:delete))
    r2i = first(movable(empty, "R2"; param=:insert))
    r2inserted = apply_reidemeister(empty, r2i)
    @test any(c -> canonicalize(apply_reidemeister(r2inserted, c)) == canonicalize(empty),
              movable(r2inserted, "R2"; param=:delete))
    # Same-gap R2 has both concatenation orders as distinct witnesses.
    same_gap = movable(empty, "R2"; param=:insert, edges=[(0, 0), (0, 0)])
    @test length(same_gap) == 16
    @test any(c -> c.after[1] > 0, same_gap)
    @test any(c -> c.after[1] < 0, same_gap)
    for c in same_gap
        @test verify_reidemeister_transition(empty, apply_reidemeister(empty, c), c)
    end

    # R3 uses matched arcs and equal 3-sign, not equal raw crossing signs.
    r3 = KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1])
    r3c = first(movable(r3, "R3"))
    r3after = apply_reidemeister(r3, r3c)
    @test r3c.direction == :self_inverse
    @test r3c.before == (-1, -2, -3, 1, 2, 3)
    @test verify_reidemeister_transition(r3, r3after, r3c)
    @test canonicalize(apply_reidemeister(r3after, first(movable(r3after, "R3")))) == canonicalize(r3)
    @test isempty(movable(KnotGaussDiagram([-1, 2, -3, 1, -2, 3], [1, 1, 1]), "R3"))
    @test isempty(movable(KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1]), "R3"; vertices=[1, 2]))
    separated = KnotGaussDiagram([-1, -2, 4, -3, 1, -4, 2, 3], [1, -1, 1, 1])
    separated_candidates = movable(separated, "R3")
    @test any(c -> c.vertices == (1, 2, 3) && c.edge_positions == ((1,2),(4,5),(7,8)), separated_candidates)
    separated_candidate = first(separated_candidates)
    @test verify_reidemeister_transition(separated,
        apply_reidemeister(separated, separated_candidate), separated_candidate)
    rotated = cyclic_shift(r3, 1)
    @test !isempty(movable(rotated, "R3"))
    @test any(c -> verify_reidemeister_transition(rotated, apply_reidemeister(rotated,c), c), movable(rotated,"R3"))
    paper_arrow_reversal = KnotGaussDiagram([-x for x in r3.gauss], collect(r3.orientation))
    @test length(movable(paper_arrow_reversal, "R3")) == length(movable(r3, "R3"))
    @test Set(c.edge_positions for c in movable(paper_arrow_reversal, "R3")) ==
          Set(c.edge_positions for c in movable(r3, "R3"))

    # NamedTuple aliases retain the MATLAB-shaped p/v/e interface.
    @test length(movable(r1d, "R1", (p=:delete, v=[1], e=[(1, 2)]))) == 1
    @test_throws ArgumentError move(r1d, "R1"; param=:delete) # two rooted arcs are ambiguous
    @test move(r1d, "R1"; param=:delete, vertices=[1], edges=[(1,2)]) == empty
    @test move(r1d, "R1", (p=:delete, v=[1], e=[(2,1)], strict=true)) == empty
    @test_throws ArgumentError movable(r1d, "R1"; param=:bogus)
end

@testset "Reidemeister stale and forged witness rejection" begin
    empty = KnotGaussDiagram(Int[], Int[])
    d = KnotGaussDiagram([1, -1], [1])
    c = first(filter(x -> x.direction == :delete, movable(d, "R1")))
    @test_throws ArgumentError apply_reidemeister(KnotGaussDiagram([-1, 1], [1]), c)
    @test !verify_reidemeister_transition(KnotGaussDiagram([-1, 1], [1]),
                                           KnotGaussDiagram(Int[], Int[]), c)
    forged = ReidemeisterCandidate(c.type, c.direction, [1], [(1,2)], [1,1], Int[], c.variant, d)
    @test_throws ArgumentError apply_reidemeister(d, forged)
    forged_variant = ReidemeisterCandidate(c.type, c.direction, c.vertices,
        c.edge_positions, c.before, c.after, :not_a_generated_variant, d)
    @test_throws ArgumentError apply_reidemeister(d, forged_variant)
    @test !verify_reidemeister_transition(d, KnotGaussDiagram(Int[], Int[]), forged_variant)
    @test_throws MethodError setindex!(c.vertices, 2, 1)
    @test_throws MethodError setindex!(c.edge_positions, (2,1), 1)
    source = read(joinpath(@__DIR__, "..", "src", "reidemeister.jl"), String)
    verifier = source[findfirst("function _rm_verify_expected", source)[1]:end]
    for forbidden in ("_rm_apply_raw", "_rm_relabel_after_delete", "_rm_insert_at_gaps",
                      "_rm_make_", "_rm_r3_data", "_rm_ccw_direction")
        @test !occursin(forbidden, verifier)
    end
end
