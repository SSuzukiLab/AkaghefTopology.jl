# cf: lib/topology.jl/src/heegaard/hc_validate.jl heegaard_family_report — which role-assigned curve systems are Heegaard diagrams
using Test, Topology
@testset "Heegaard family report" begin
    H=Topology.HeegaardCodes
    # Standard genus 2 diagram: alpha = m1, m2 and two beta curves.
    standard=H._lc_parse_text("genus: 2\nl1 : -1 -3\nl2 : -2 +4\nm1 : +1 +5\nm2 : +2 +7\nn1 : +3 +6 -4 -8\nb1 : -5 -6\nb2 : -7 +8\n")
    ok=H.heegaard_family_report(H.resolve_hc(standard,Dict("b1"=>:beta,"b2"=>:beta);standard=true))
    @test ok.ok && ok.alpha.rank==2 && ok.beta.rank==2 && isempty(ok.beta.dependent)
    @test H.is_heegaard_diagram(H.resolve_hc(standard,Dict("b1"=>:beta,"b2"=>:beta);standard=true))
    # x runs round the opening l1+ and y round l1-: different words, both parallel to l1.
    parallel=H._lc_parse_text("genus: 2\nl1 : -1 -3\nl2 : -2 +4\nm1 : +1 +5 +7\nm2 : +2\nn1 : +3 +6 -4 +8\nx : -5 -6\ny : -7 -8\n")
    pair=H.heegaard_family_report(H.resolve_hc(parallel,Dict("x"=>:alpha,"y"=>:alpha)))
    @test !pair.ok && pair.alpha.count_ok && pair.alpha.rank==1 && pair.alpha.dependent==[["x","y"]]
    @test !pair.beta.count_ok
    # One of them with a genuine second curve is not enough curves on the other side.
    @test !H.is_heegaard_diagram(H.resolve_hc(parallel,Dict("x"=>:alpha,"y"=>:beta)))
    # Crossing-free circles bound disks.
    circles=H._lc_parse_text("genus: 0\na :\nb :\nplace a : outer\nplace b : a+\n")
    empty=H.heegaard_family_report(H.resolve_hc(circles,Dict("a"=>:alpha,"b"=>:beta)))
    @test !empty.ok && empty.alpha.dependent==[["a"]] && !empty.alpha.count_ok
    # A curve of one family that crosses another curve of the same family.
    crossing=H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\nalpha1 : -3 +4 -5 +6\nbeta1 : -4 +5 -6 -2\n")
    same=H.heegaard_family_report(H.resolve_hc(crossing,Dict("alpha1"=>:alpha,"beta1"=>:alpha)))
    @test !same.ok && same.alpha.family_crossings==[("alpha1","beta1")]
    @test H.is_heegaard_diagram(H.resolve_hc(crossing,Dict("alpha1"=>:alpha,"beta1"=>:beta)))
    # Three classes with one relation: the third is the sum of the first two.
    @test H._hc_f2_dependencies([Bool[1,0,0,0],Bool[0,1,0,0],Bool[1,1,0,0]])==[[1,2,3]]
    @test isempty(H._hc_f2_dependencies([Bool[1,0],Bool[0,1]]))
    @test H._hc_f2_dependencies([Bool[0,0]])==[[1]]
end
