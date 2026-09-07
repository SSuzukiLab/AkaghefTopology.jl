using Test
using Topology

# Minimal closed PS-I source: one cyclic component, with the two PS-R
# fragments [1,2] and [-1,-2].  The validator is deliberately tested through
# the public VirtualLink/set_data! surface, not by mutating internals.
function ps_fixture()
    v = VirtualLink()
    # A lens-space presentation (L(4,1)) with the PS-III local pair.
    set_data!(v; rgauss=[[-4,-3,-2,-1,1,2,3,4]],
              orientation=[1,-1,-1,-1])
    v
end

@testset "closed-normal PS validation" begin
    source = ps_fixture()
    candidates = movable(source, "PS"; vertices=[1,2], strict=true)
    @test !isempty(candidates)
    candidate = first(candidates)
    @test candidate.parameter == "3"
    @test isdefined(Topology, :validate_ps_move)
    report = Topology.validate_ps_move(source, candidate)
    @test report.valid
    @test report.source_closed
    @test report.target_closed
    @test report.source_conditions isa NTuple{3,Bool}
    @test report.target_conditions isa NTuple{3,Bool}
    @test report.target !== nothing
    @test move(source, "PS"; vertices=[1,2], strict=true).gauss == report.target.gauss

    before = deepcopy(source)
    Topology.validate_ps_move(source, candidate)
    @test source.gauss == before.gauss
    @test source.orientation == before.orientation
    @test source.weights == before.weights

    open_source = VirtualLink()
    set_data!(open_source; rgauss=[[1,2],[-1,-2]], orientation=[1,-1])
    local_candidate = MoveCandidate("PS", "1", [1,2], NTuple{2,Int}[], :inverse,
                                    Topology.make_move_data()["PS-R"][1],
                                    Topology.make_move_data()["PS-L"][1],
                                    [(1,1),(2,1)])
    open_report = Topology.validate_ps_move(open_source, local_candidate)
    @test !open_report.valid
    @test open_report.reason == :source_not_closed
    @test isempty(movable(open_source, "PS"; vertices=[1,2], strict=true))

    # Under the source-confirmed C2 graph, this former tiny fixture is not a
    # closed source.  Keep the rejection explicit; do not infer target-side
    # behavior from a source that fails the frozen C2 closedness check.
    tiny = VirtualLink()
    set_data!(tiny; rgauss=[[1,2,-1,-2]], orientation=[1,-1])
    tiny_data = Topology.make_move_data()
    tiny_candidate = MoveCandidate("PS", "1", [1,2], NTuple{2,Int}[], :inverse,
                                   tiny_data["PS-R"][1], tiny_data["PS-L"][1],
                                   [(1,1),(1,3)])
    tiny_report = Topology.validate_ps_move(tiny, tiny_candidate)
    @test !tiny_report.source_closed
    @test !tiny_report.target_closed
    @test tiny_report.reason == :source_not_closed

    nonps = MoveCandidate("CP", "", [1,2], NTuple{2,Int}[], :forward,
                          candidate.before, candidate.after, candidate.locations)
    nonps_report = Topology.validate_ps_move(source, nonps)
    @test !nonps_report.valid
    @test nonps_report.reason == :not_ps

    weighted = VirtualLink()
    set_data!(weighted; rgauss=[[-4,-3,-2,-1,1,2,3,4]],
              orientation=[1,-1,-1,-1])
    set_weight!(weighted, [zeros(8)])
    weighted_data = Topology.make_move_data()
    weighted_candidate = MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                                       weighted_data["PS-R"][3], weighted_data["PS-L"][1],
                                       [(1,5),(1,3)])
    weighted_report = Topology.validate_ps_move(weighted, weighted_candidate)
    @test !weighted_report.valid
    @test weighted_report.reason == :weighted_unsupported

    @test all(Topology.validate_ps_move(source, c).valid for c in movable(source, "PS"; vertices=[1,2]))

    data = Topology.make_move_data()
    malformed_cases = [
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                       nothing, candidate.after, candidate.locations), :malformed_candidate),
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                       candidate.before, nothing, candidate.locations), :malformed_candidate),
        (MoveCandidate("PS", "2", [1,2], NTuple{2,Int}[], :inverse,
                       candidate.before, candidate.after, candidate.locations), :source_mismatch),
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :forward,
                       candidate.before, candidate.after, candidate.locations), :source_mismatch),
        (MoveCandidate("PS", "3", [2,1], NTuple{2,Int}[], :inverse,
                       candidate.before, candidate.after, candidate.locations), :source_mismatch),
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                       candidate.before, candidate.after, [(1,4),(1,3)]), :source_mismatch),
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                       candidate.before, candidate.after, [(1,5),(1,5)]), :source_mismatch),
        (MoveCandidate("PS", "3", [1,2], NTuple{2,Int}[], :inverse,
                       data["PS-R"][2], data["PS-L"][1], [(1,5),(1,3)]), :source_mismatch),
    ]
    for (forged, reason) in malformed_cases
        forged_report = Topology.validate_ps_move(source, forged)
        @test !forged_report.valid
        @test forged_report.reason == reason
        @test forged_report.error === nothing
    end

    edge_candidate = MoveCandidate("PS", "3", [1,2], [(99,100)], :inverse,
                                   candidate.before, candidate.after, candidate.locations)
    edge_report = Topology.validate_ps_move(source, edge_candidate)
    @test !edge_report.valid
    @test edge_report.reason == :edge_selected_unsupported
    @test isempty(movable(source, "PS"; vertices=[1,2], edges=[(99,100)], strict=true))
    @test isempty(movable(weighted, "PS"; vertices=[1,2], strict=true))

    # C260110 LHS prefix: weighted initial state, two C3 MP moves, then
    # exactly the QMD clear-weight transition before the closed PS check.
    c260 = VirtualLink()
    set_data!(c260; rgauss=[[-1,2,-2,3,-3,1]], orientation=[1,1,1])
    set_weight!(c260, [[1,-2,1,-1,1,-1]])
    move!(c260, "MP"; param="C3", vertices=[1,2])
    move!(c260, "MP"; param="C3", vertices=[2,3])
    set_data!(c260; rgauss=c260.gauss, orientation=c260.orientation)
    c260_closed, c260_conditions = is_closed(c260)
    @test c260_closed
    @test c260_conditions == (true, true, true)
    c260_candidates = movable(c260, "PS"; vertices=[5,2], strict=true)
    @test !isempty(c260_candidates)
    c260_report = Topology.validate_ps_move(c260, first(c260_candidates))
    @test c260_report.valid
    c260_target = move(c260, "PS"; vertices=[5,2], strict=true)
    @test is_closed(c260_target)[1]
end
