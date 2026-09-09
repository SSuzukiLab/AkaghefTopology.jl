using Test
using Topology

function _cut_spin_cp_fixture(t;extra_empty=false)
    source=spin_cp_fixture(t)
    geometry=Topology._weighted_unweight(source)
    if extra_empty
        push!(geometry.gauss,[0])
    end
    set_cuts!(geometry,fill(true,length(geometry.gauss)))
    values=reduce(vcat,source.weights;init=Float64[])
    weights=extra_empty ? vcat(values,0,0) : vcat(values,0)
    set_weight!(geometry,weights;mode=:z2)
    geometry
end

function _cut_spin_cp_candidate(source)
    geometry=Topology._weighted_unweight(source)
    only(movable(geometry,"CP";vertices=[1,2,3],strict=true))
end

@testset "cut spin CP captures and transports relative boundary data" begin
    for t in 0:1
        source=_cut_spin_cp_fixture(t)
        candidate=_cut_spin_cp_candidate(source)
        before=(deepcopy(source.gauss),deepcopy(source.weights),copy(source.orientation))
        report=spin_cp_transport(source,candidate;normalize=false)
        target=report.target
        @test report.case==t
        @test report.modulus==2 && report.boundary_modulus==4
        @test last(only(target.gauss))==0
        @test length(only(target.weights))==length(only(target.gauss))
        @test length(report.boundary_transport)==2
        @test validate_spin_weight(target;boundary=report.target_boundary)
        @test source.gauss==before[1] && source.weights==before[2] && source.orientation==before[3]

        # The inverse candidate uses the transported target vertex roles.
        inverse=only(movable(Topology._weighted_unweight(target),"CP";
                             vertices=report.target_vertices,strict=true))
        back=spin_cp_transport(target,inverse;normalize=false,
                               boundary=report.target_boundary).target
        @test back.gauss==source.gauss
        @test back.weights==source.weights
    end
end

@testset "cut spin CP accepts explicit beta and rejects tampering atomically" begin
    source=_cut_spin_cp_fixture(0)
    candidate=_cut_spin_cp_candidate(source)
    captured=spin_cp_transport(source,candidate;normalize=false)
    explicit=(source=captured.source_boundary,target=captured.target_boundary)
    checked=spin_cp_transport(source,candidate;normalize=false,boundary=explicit)
    @test checked.target.gauss==captured.target.gauss
    @test checked.target.weights==captured.target.weights

    bad=deepcopy(source)
    snapshot=(deepcopy(bad.gauss),deepcopy(bad.weights),copy(bad.orientation))
    bad_boundary=copy(captured.source_boundary)
    key=first(keys(bad_boundary)); bad_boundary[key]=mod(bad_boundary[key]+1,4)
    @test_throws ArgumentError spin_cp_move(bad,candidate;normalize=false,boundary=bad_boundary)
    @test bad.gauss==snapshot[1] && bad.weights==snapshot[2] && bad.orientation==snapshot[3]

    crossing=MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,
        candidate.edges,candidate.direction,candidate.before,candidate.after,
        [(1,length(only(source.gauss))-1)])
    @test_throws ArgumentError spin_cp_transport(source,crossing;normalize=false)

    # Row-aligned input is normalized to the same directed endpoint map.
    data=spin_boundary_data(source)
    row_targets=Any[data.isCycle[row] ? 2 : data.measured4[row]
                    for row in eachindex(data.isCycle)]
    row_report=spin_cp_transport(source,candidate;normalize=false,boundary=row_targets)
    @test row_report.source_boundary==captured.source_boundary

    # A caller-provided target beta is checked in addition to the transported
    # source beta, so a single altered residue cannot be smuggled through.
    bad_target=copy(captured.target_boundary)
    target_key=first(keys(bad_target)); bad_target[target_key]=mod(bad_target[target_key]+2,4)
    @test_throws ArgumentError spin_cp_transport(source,candidate;normalize=false,
        boundary=(source=captured.source_boundary,target=bad_target))
end

@testset "cut spin CP accepts local H gauges while preserving beta" begin
    source=_cut_spin_cp_fixture(1)
    geometry=Topology._weighted_unweight(source)
    candidate=only(movable(geometry,"CP";vertices=[1,2,3],strict=true))
    for mask in 0:7
        gauged,_=Topology._spin_cp_hgauge(source,[1,2,3],mask)
        measured=spin_boundary_data(gauged)
        open_rows=findall(!,Bool.(measured.isCycle))
        beta=Dict((measured.boundary_start[row],measured.boundary_end[row])=>
                  measured.measured4[row] for row in open_rows)
        report=spin_cp_transport(gauged,candidate;normalize=true,boundary=beta)
        @test validate_spin_weight(report.target;boundary=report.target_boundary)
        @test report.boundary_modulus==4
    end
end

@testset "cut spin CP retains a separate exterior cut component" begin
    source=VirtualLink(;weight_mode=:z2)
    set_data!(source;rgauss=[[1,-2,2,-3,3,-1,4,-5,5,-4,0],[6,-7,7,-6,0]],
              orientation=[-1,-1,-1,-1,1,1,1])
    set_weight!(source,[[0,0,1,1,1,0,1,0,0,1,0],[1,0,0,1,0]];mode=:z2)
    geometry=Topology._weighted_unweight(source)
    candidate=only(movable(geometry,"CP";vertices=[1,2,3],strict=true))
    report=spin_cp_transport(source,candidate;normalize=false)
    @test any(word -> word==[6,-7,7,-6,0],report.target.gauss)
    @test length(report.boundary_transport)==4
    @test validate_spin_weight(report.target;boundary=report.target_boundary)
end

@testset "cut spin CP preserves an exterior empty interval" begin
    source=_cut_spin_cp_fixture(1;extra_empty=true)
    candidate=_cut_spin_cp_candidate(source)
    report=spin_cp_transport(source,candidate;normalize=false)
    @test report.target.gauss[end]==[0]
    @test report.target.weights[end]==[0.0]
    @test length(report.boundary_transport)==4
    @test validate_spin_weight(report.target;boundary=report.target_boundary)
    @test length(unique(report.target_boundary))==length(report.target_boundary)
end
