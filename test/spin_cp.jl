using Test, Topology

function spin_cp_fixture(t)
    v=VirtualLink(;weight_mode=:z2)
    if t==0
        set_data!(v;rgauss=[1,-2,2,-3,3,-1],orientation=[-1,-1,-1])
        set_weight!(v,[1,0,1,0,1,0])
    else
        # The earlier integer-weight-infeasible target example, with its
        # source integer solution explicitly reduced mod 2.
        set_data!(v;rgauss=[1,-2,2,-3,3,-1,4,-5,5,-4],
            orientation=[-1,-1,-1,-1,1])
        set_weight!(v,[0,0,1,-1,1,-2,1,0,0,-1])
    end
    v
end

function spin_cp_same_weighted_cycle(a,b)
    ac,ao=real_gauss_code(a); bc,bo=real_gauss_code(b)
    ao==bo && length(ac)==length(bc) || return false
    all(any(s -> circshift(bc[i],1-s)==ac[i] &&
        circshift(b.weights[i],1-s)==a.weights[i],
        Topology._cyclic_occurrences(bc[i],ac[i])) for i in eachindex(ac))
end

@testset "Spin CP inverse normalization across all RHS H gauges" begin
    source=spin_cp_fixture(1)
    forward=only(movable(source,"CP";vertices=[1,2,3],strict=true))
    initial=spin_cp_transport(source,forward;normalize=false)
    inverse=only(movable(initial.target,"CP";vertices=initial.target_vertices,strict=true))
    for mask in 0:31
        gauged,_=Topology._spin_cp_hgauge(initial.target,initial.target_vertices,mask)
        report=spin_cp_transport(gauged,inverse)
        @test report.case==1 && is_spin_weight(report.target)
        redo=only(movable(report.target,"CP";vertices=report.target_vertices,strict=true))
        back=spin_cp_move(report.target,redo;normalize=false)
        @test spin_cp_same_weighted_cycle(back,report.normalized_source)
    end
end

@testset "Spin CP published edge occurrences and boundary equations" begin
    for t in 0:1
        l,r=Topology._spin_cp_patterns(t)
        w=only(l.weights); q=only(r.weights)
        @test w==[0,0,1,t,1,0,0]
        @test q==[0,0,1-t,0,1,0,0,1,1,0,0]
        # Independent local closed-disk constraints and the three B passages.
        @test all(iszero,mod.([q[4]-q[5]-q[6]+q[7]+q[9]+q[10],
            q[2]-q[4]+q[6]-q[10],q[5]+q[6]+q[8],1-q[8]],2))
        @test mod(q[2]+q[3]-sum(w[2:6])+1,2)==0
        @test mod(q[10]-w[3]-w[5]+2,2)==0
        @test mod(-q[2]+q[4]+q[5]+q[8]-2,2)==0
    end
end

@testset "Spin CP public entry, both cases, exact normalized inverse" begin
    for t in 0:1
        source=spin_cp_fixture(t)
        original=deepcopy(source)
        @test is_spin_weight(source)
        candidate=only(movable(source,"CP";vertices=[1,2,3],strict=true))
        report=spin_cp_transport(source,candidate;normalize=false)
        target=report.target
        @test report.case==t && isempty(report.source_h_vertices)
        @test report.source_already_normalized
        @test weight_mode(target)==:z2 && is_spin_weight(target)
        component,start=only(report.target_locations)
        @test [target.weights[component][mod1(start+k,length(target.weights[component]))]
            for k in 1:9]==[0,1-t,0,1,0,0,1,1,0]
        @test source.gauss==original.gauss && source.weights==original.weights
        @test move(source,candidate).weights==target.weights
        @test move(source,"CP";vertices=[1,2,3],strict=true).weights==target.weights
        inverse=only(movable(target,"CP";vertices=report.target_vertices,strict=true))
        back=spin_cp_transport(target,inverse;normalize=false).target
        @test spin_cp_same_weighted_cycle(back,source)
        @test weight_mode(back)==:z2 && is_spin_weight(back)
    end
end

@testset "Spin CP explicit H normalization and rejection" begin
    source=spin_cp_fixture(1)
    geometry=Topology._weighted_unweight(source)
    candidate=only(movable(geometry,"CP";vertices=[1,2,3],strict=true))
    # Exhaust every source-local H gauge; equality is only claimed after
    # applying the H record, not as a raw-weight inverse of H*CP.
    for mask in 0:7
        gauged,_=Topology._spin_cp_hgauge(source,[1,2,3],mask)
        @test is_spin_weight(gauged)
        report=spin_cp_transport(gauged,candidate)
        @test is_spin_weight(report.target)
        @test report.case==1
        inverse=only(movable(report.target,"CP";vertices=report.target_vertices,strict=true))
        back=spin_cp_move(report.target,inverse;normalize=false)
        @test spin_cp_same_weighted_cycle(back,report.normalized_source)
    end
    gauged,_=Topology._spin_cp_hgauge(source,[1,2,3],1)
    @test_throws ArgumentError spin_cp_move(gauged,candidate;normalize=false)
    bad=deepcopy(source); bad.weights[1][3]=0
    @test_throws ArgumentError spin_cp_move(bad,candidate)
    @test isempty(movable(bad,"CP";vertices=[1,2,3],strict=true))
    @test_throws ArgumentError spin_cp_move(geometry,candidate)
    cut=Topology._weighted_unweight(source); set_cuts!(cut,[true])
    set_weight!(cut,vcat(only(source.weights),0);mode=:z2)
    @test_throws ArgumentError spin_cp_move(cut,candidate)
    @test source.weights==spin_cp_fixture(1).weights
end
