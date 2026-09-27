using Test
using Topology

_symmetry_snapshot(v)=Tuple(deepcopy(getfield(v,f)) for f in fieldnames(VirtualLink))

@testset "diagram primitives and MST semantic operations" begin
    v=lens_space_link(2,1)
    saved=_symmetry_snapshot(v)
    for f in (switch_crossings,reflect_diagram,reverse_strands,mst_sign_change,mst_reflect_reverse)
        @test _symmetry_snapshot(f(f(v)))==saved
        @test _symmetry_snapshot(v)==saved
    end
    @test switch_crossings(v).gauss==[-w for w in v.gauss]
    @test reflect_diagram(v).gauss==v.gauss
    @test switch_crossings(v).orientation==reflect_diagram(v).orientation==-v.orientation
    @test reverse_strands(v).gauss==reverse.(v.gauss)
    @test reverse_strands(v).orientation==v.orientation
    @test mst_reflect_reverse(v).gauss==reverse.(v.gauss)
    @test first(is_closed(mst_sign_change(v)))
    @test first(is_closed(mst_reflect_reverse(v)))
    for (f,mutating) in ((switch_crossings,switch_crossings!),
                          (reflect_diagram,reflect_diagram!),
                          (reverse_strands,reverse_strands!),
                          (mst_sign_change,mst_sign_change!),
                          (mst_reflect_reverse,mst_reflect_reverse!))
        target=deepcopy(v)
        @test mutating(target)===target
        @test _symmetry_snapshot(target)==_symmetry_snapshot(f(v))
    end
    cut=deepcopy(v); set_cuts!(cut,[true])
    @test last(only(reverse_strands(cut).gauss))==0
    @test _symmetry_snapshot(reverse_strands(reverse_strands(cut)))==_symmetry_snapshot(cut)
    @test_throws ArgumentError mst_sign_change(cut)
    @test_throws ArgumentError mst_reflect_reverse(cut)
    virtual=VirtualLink(); set_data!(virtual;gauss=[[1,2,-1,-2]],orientation=[1,0])
    @test switch_crossings(virtual).gauss==[[-1,2,1,-2]]
    @test_throws ArgumentError reflect_diagram(virtual)
    for mode in (:integer,:z2)
        weighted=deepcopy(v); weighted.weight_mode=mode
        weighted.is_weighted=true; weighted.weights=[zeros(4)]
        old=_symmetry_snapshot(weighted)
        for f in (switch_crossings!,reflect_diagram!,reverse_strands!,mst_sign_change!,mst_reflect_reverse!,mirror_manifold!)
            @test_throws ArgumentError f(weighted)
            @test _symmetry_snapshot(weighted)==old
        end
    end
end

@testset "0-2 over/under parameter, not direction reversal" begin
    # Exact finite corpus: standard lens presentations, s=2:7, gcd(s,t)=1;
    # every directed edge. Source is MST25 Fig. 2.9 + Prop. A.1 conjugation.
    for s in 2:7, t in 1:s-1
        gcd(s,t)==1 || continue
        base=lens_space_link(s,t)
        word=only(first(real_gauss_code(base)))
        for i in eachindex(word)
            site=[word[i],word[mod1(i+1,length(word))]]
            standard=move(base,"02";param="dir",vertices=site)
            switched=move(base,"02";param="dir-switched",vertices=site)
            expected=switch_crossings(move(switch_crossings(base),"02";param="dir",vertices=-site))
            @test switched.gauss==expected.gauss
            @test switched.orientation==expected.orientation
            @test first(is_closed(switched))==first(is_closed(base))
            n=length(base.orientation)
            for (target,param) in ((standard,"inv"),(switched,"inv-switched"))
                candidate=only(movable(target,"02";param=param,vertices=[n+1,n+2],strict=true))
                restored=move(target,candidate)
                @test restored.gauss==base.gauss
                @test restored.orientation==base.orientation
            end
            candidate=only(movable(base,"02";param="dir-switched",vertices=site))
            @test move(base,candidate).gauss==switched.gauss
        end
    end
    base=lens_space_link(2,1)
    cut=deepcopy(base); set_cuts!(cut,[true])
    added=move(cut,"02";param="dir-switched",vertices=[-2,-1])
    @test strand_cuts(added)==[true]
    @test move(added,"02";param="inv-switched",vertices=[3,4]).gauss==cut.gauss
    weighted=deepcopy(base); set_weight!(weighted,zeros(4))
    old=_symmetry_snapshot(weighted)
    @test isempty(movable(weighted,"02";param="dir-switched",vertices=[-2,-1]))
    @test_throws ArgumentError move!(weighted,"02";param="dir-switched",vertices=[-2,-1])
    candidate=only(movable(base,"02";param="dir-switched",vertices=[-2,-1]))
    @test_throws ArgumentError move!(weighted,candidate)
    @test _symmetry_snapshot(weighted)==old
    # Preserve standard integer transport; this extension does not remove it.
    added=move(weighted,"02";param="dir",vertices=[-2,-1])
    restored=move(added,"02";param="inv",vertices=[3,4])
    @test restored.gauss==weighted.gauss
    @test restored.weights==weighted.weights
    # Reversing the local arrows gives the same visits as switched, but the
    # WRONG vertex signs. This is not the implemented parameter.
    bad=VirtualLink()
    set_data!(bad;rgauss=[[-2,3,4,-4,-3,-1,1,2]],orientation=[1,1,1,-1])
    @test is_closed(bad)==(false,(true,true,false))
end
