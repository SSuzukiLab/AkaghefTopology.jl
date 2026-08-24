using Test, Topology
const ROOT=normpath(joinpath(@__DIR__,".."))
for file in ["VLExample.jl","C240610URD.jl","C240623URD.jl","URDforKnotCmpl.jl",
             "dev_URDiagram.jl","OgraphPlotSuite.jl","C251030WeylAlgTraceFormula.jl",
             "C260319MSTinv_unknot.jl","InvariantWorkflows.jl","C251020ReplayFixture.jl"]
    include(joinpath(ROOT,"examples",file))
end

include(joinpath(ROOT,"examples","C250406TensorNetwork.jl"))
include(joinpath(ROOT,"examples","C250825Dsl2.jl"))
include(joinpath(ROOT,"examples","C260225M_Oinv.jl"))
include(joinpath(ROOT,"examples","C250405Sweedler.jl"))
include(joinpath(ROOT,"examples","C250823MSTinv_extalg.jl"))
include(joinpath(ROOT,"examples","C250818MSTinv_uqsl2.jl"))
include(joinpath(ROOT,"examples","C260125MSTinv_uqsl2.jl"))
include(joinpath(ROOT,"examples","C260220LensSpaces.jl"))
include(joinpath(ROOT,"examples","C251020KnotDsl2bfromUR.jl"))
include(joinpath(ROOT,"examples","C250811MSTinvUR.jl"))

@testset "URDiagram core and Live Scripts" begin
    d=URDiagram(); set_data_cvw!(d,1,[1,-2,-3,3,-1,2],[sym(:x),sym(:y),sym(:z)])
    @test validate_data(d); @test rank(d)==1
    @test trace3(d)==-1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
    block=Topology._block_matrix(d)
    @test block==Any[0 sym(:x) 0; sym(:z) 0 0; sym(:y) 0 0]
    swap!(d,[3,-3]); @test d.V==[1,-2,3,-3,-1,2]
    @test d.C==1/(1-sym(:z)); @test d.W[3]==sym(:z)/(1-sym(:z))
    @test_throws ArgumentError validate_data(set_data_cvw!(URDiagram(),1,[1,-2,-3,4,1,2],[1,-1,1]))
    @test length(run_C240610().W) >= 1
    trace610=run_C240610_trace(); @test trace610.diagram.V==[1,-1]; @test trace610.limit==-1
    diagrams623=run_C240623(); @test length(diagrams623)==4; @test diagrams623[end].V==[2,-2]
    traces623=run_C240623_traces(); @test [item.limit for item in traces623]==[-1//9,-1,-1]
    @test Set(keys(run_URDforKnotCmpl())) == Set(["unknot","trefoil","figure-eight"])
    @test run_URDforKnotCmpl_manual_prefix().V==[-5,1,-2,2,5,-1]
    d2_exploratory=run_URDforKnotCmpl_d2_exploratory_cell()
    @test d2_exploratory.prefix.V==[-5,1,-2,2,5,-1]
    @test all(snapshot -> all(isfinite,snapshot.W),
              (d2_exploratory.after_simplify,d2_exploratory.before_last_dilation,d2_exploratory.final))
    @test isapprox(run_URDforKnotCmpl_figure_eight_trace().trace,-4.526591341323682e-5;rtol=1e-12)
    dev, err=run_dev_urdiagram(); @test err isa ArgumentError; @test validate_data(dev)
    dev_reductions=run_dev_reduction_cells()
    @test dev_reductions.reduct3==[1,-3,3,4,-1,0,-4]
    @test dev_reductions.reduct5_first==[1,2,-1,-2,4,-4]
    @test dev_reductions.reduct5_second==[5,-5]
    @test dev_reductions.simplify1==[4,-4]
    @test run_dev_denominator_loop().denominators==[1,3,4,3,1,Inf,1,3,4,3,1,Inf,1,3,4]
end

@testset "VirtualLink core and VLExample" begin
    links=run_vl_examples(); @test length(links)==19
    hopf=links["hopf"]; @test gauss_code(hopf)==[[1,-2],[-1,2]]
    @test size(head_map(hopf))==(2,2); @test size(pd_code(hopf))==(2,4)
    @test pd_code(hopf)==[3 2 4 1; 2 3 1 4]
    @test regions(pd_code(hopf))==[[4,1],[3,-1],[2,-4],[-2,-3]]
    @test minimal_bending_numbers(pd_code(hopf))==[2,2,4,0]
    layout=orthogonal_layout(hopf)
    @test length(layout.crossings)==2
    @test sort([edge.id for edge in layout.edges])==collect(1:4)
    @test all(edge -> first(edge.points)==layout.crossings[edge.tail],layout.edges)
    planar_gauss,planar_orientation=planarize_virtual_gauss(
        [[2,-1,-2,1,4,-3,-4,3]],[-1,1,-1,1])
    @test planar_gauss==[[2,-5,-1,-2,5,1,4,-6,-3,-4,6,3]]
    @test planar_orientation==[-1,1,-1,1,0,0]
    virtual=VirtualLink()
    set_data!(virtual;rgauss=[[2,-1,-2,1,4,-3,-4,3]],orientation=[-1,1,-1,1])
    virtual_layout=orthogonal_layout(virtual)
    @test length(virtual_layout.crossings)==6
    @test count(==(0),virtual_layout.orientation)==2
    base,kernel,hperp=calculate_weight(links["koda"])
    @test length(base)==18; @test size(kernel,2)==18; @test size(hperp,2)==18
    closed,criterion=is_closed(links["s2xs1"]); @test criterion==(true,true,true); @test closed
    mirror_manifold!(hopf); @test hopf.orientation==[-1,-1]
    cyclic_shift!(hopf,1,1); @test length(hopf.gauss[1])==2
    summand=VirtualLink(); set_data!(summand;rgauss=[[1,-1]],orientation=[1])
    disjoint_sum!(summand,deepcopy(summand)); connected_sum!(summand)
    @test summand.gauss==[[1,-1,4,3,-3,-6,2,-2,6,5,-5,-4]]
    @test summand.orientation==[1,1,-1,1,-1,1]
    complement=VirtualLink(); set_data!(complement;rgauss=[[1,-1]],orientation=[1])
    knot_complement!(complement)
    @test complement.gauss==[[1,2,-2,-3],[-4,-1,3,4]]
    @test complement.orientation==[-1,1,-1,1]
    set_weight!(complement,[[-1,0,1,0],[1,0,0,-1]])
    edges=real_edge_table(complement)
    @test [edge.crossing for edge in edges]==[(-3,1),(1,2),(2,-2),(-2,-3),
                                               (4,-4),(-4,-1),(-1,3),(3,4)]
    @test [edge.weight for edge in edges]==[-1,0,1,0,1,0,0,-1]
    table_weighted=VirtualLink(); set_data!(table_weighted;rgauss=[[1,-1]],orientation=[1])
    knot_complement!(table_weighted)
    set_weight_by_crossing!(table_weighted,[2 -2; -3 1; 3 4; 4 -4],[1,-1,-1,1])
    @test [edge.weight for edge in real_edge_table(table_weighted)]==[-1,0,1,0,1,0,0,-1]
    positive=VirtualLink(); set_data!(positive;rgauss=[[1,-1]],orientation=[1])
    negative=VirtualLink(); set_data!(negative;rgauss=[[1,-1]],orientation=[-1])
    @test ur_prime_invariant(positive)==-1
    @test ur_prime_invariant(negative)==1
    lens=VirtualLink(); set_data!(lens;gauss=[[1,-3,2,-2,3,-1]],orientation=[1,1,0])
    lens_base,_,_=calculate_weight(lens;assign=true)
    @test lens_base==[-1,-1,1,0]
    @test lens.weights==[Float64.([-1,-1,1,0])]
    s2=VirtualLink(); set_data!(s2;rgauss=[[-1,-2,3,4,-4,-3,2,1]],orientation=[1,1,1,1])
    s2_base,_,s2_hperp=calculate_weight(s2)
    @test s2_base==[1,-1,-1,0,1,0,-1,0]
    @test s2_hperp==[0 0 0 1 0 0 0 -1]
    z2=uqsl2_borel_small(1//2;style=:L)
    z2_isolated=uqsl2_borel_small(1//2;style=:L)
    @test z2 !== z2_isolated
    z2_isolated.antipode.terms[(1,1)]=999
    @test z2.antipode.terms[(1,1)]!=999
    abalone_hopf=VirtualLink(); set_data!(abalone_hopf;rgauss=[[1,-1]],orientation=[1])
    set_weight!(abalone_hopf,[[-2,1]])
    @test hopf_invariant(abalone_hopf,z2)≈-1
    set_weight!(s2,[s2_base])
    @test hopf_invariant(s2,z2)≈0 atol=1e-12
    z5=uqsl2_borel_small(2//5;style=:L)
    @test hopf_invariant(lens,z5)≈(-0.381966011250059-1.17557050458494im) atol=1e-10
    sweedler=sweedler_algebra()
    @test hopf_invariant(abalone_hopf,sweedler)≈-1
    @test hopf_invariant(links["koda_fig13"],sweedler)≈0 atol=1e-12
    exterior=exterior_algebra(2)
    @test exterior.dimension==4
    @test exterior.antipode.terms==Dict{Tuple,ComplexF64}(
        (1,1)=>1,(2,2)=>-1,(3,3)=>-1,(4,4)=>1)
    @test all(values(verify_hopf(sweedler)))
    @test all(values(verify_hopf(z2)))
    @test all(values(verify_hopf(exterior;degrees=count_ones.(0:3))))
    lens_w0,lens_kernel,_=calculate_weight(lens)
    @test lens_w0==[-1,-1,1,0]
    @test vec(lens_kernel[1,:])==[0,-1,0,1]
    set_weight!(lens,[lens_w0]); @test hopf_invariant(lens,exterior)≈4
    set_weight!(lens,[lens_w0+3vec(lens_kernel[1,:])]); @test hopf_invariant(lens,exterior)≈4
    @test hopf_invariant(abalone_hopf,exterior)≈5
    set_weight!(s2,[s2_base]); @test hopf_invariant(s2,exterior)≈28
    for multiple in 1:3
        set_weight!(s2,[s2_base+multiple*vec(s2_hperp)])
        @test hopf_invariant(s2,exterior)≈28
    end
    koda_ext=deepcopy(links["koda_fig13"]); calculate_weight(koda_ext;assign=true)
    @test hopf_invariant(koda_ext,exterior)≈40
    invariant_results=run_invariant_workflows()
    @test invariant_results["lens_uqsl2_5"]≈(-0.381966011250059-1.17557050458494im) atol=1e-10
    @test invariant_results["s2xs1_uqsl2_2"]≈0 atol=1e-12
    @test invariant_results["lens_exterior_2"]≈4
    @test invariant_results["s2xs1_exterior_2"]≈28
    @test invariant_results["abalone_exterior_2"]≈5
    @test invariant_results["koda_exterior_2"]≈40
end

@testset "SVG rendering" begin
    out=mktempdir(); d=run_C240610(); v=run_vl_examples()["koda"]
    p1=plot_svg(d,joinpath(out,"ur.svg")); p2=plot_svg(v,joinpath(out,"vl.svg"))
    ur_svg=read(p1,String); vl_svg=read(p2,String)
    @test occursin("<svg",ur_svg); @test occursin("<svg",vl_svg)
    # MATLAB URDiagram.plot is exactly one leftward horizontal quiver.
    @test count("<line",ur_svg)==1
    @test count("<path",ur_svg)==1
    @test occursin("x1=\"780\" y1=\"300\" x2=\"120\" y2=\"300\"",ur_svg)
    @test filesize(p2)>500
    @test !occursin("<polygon",vl_svg)
    @test occursin(" Q ",vl_svg)
    replay_link=VirtualLink(); set_data!(replay_link;rgauss=[[1,-2],[-1,2]],orientation=[1,1])
    replay_positions=[ComplexF64[complex(abs(arc.source),0),complex(abs(arc.target),0)]
                      for arc in Topology._replay_arc_metadata(replay_link)]
    replay_svg=read(plot_svg(replay_link,joinpath(out,"vl-replay.svg");replay_positions=replay_positions),String)
    @test occursin("<svg",replay_svg)
    @test_throws ArgumentError plot_svg(replay_link,joinpath(out,"vl-invalid-replay.svg");replay_positions=replay_positions[1:end-1])
    c251020=run_C251020_knot_dsl2b()
    fixture_records=Any[]
    for (call_id,link) in [(1,c251020.borromean),(2,c251020.whitehead)]
        crossings=[collect(edge.crossing) for edge in real_edge_table(link)]
        positions=[[[abs(crossing[1]),0],[abs(crossing[2]),0]] for crossing in crossings]
        push!(fixture_records,Dict("relativePath"=>"C251020KnotDsl2bfromUR/call_$(lpad(call_id,3,'0'))_after.mat",
            "callId"=>call_id,"objects"=>Dict("vl"=>Dict("properties"=>Dict("REdgeTable"=>Dict(
                "Crossing"=>crossings,"Position"=>positions))))))
    end
    fixture_json=joinpath(out,"c251020-fixture.json")
    open(fixture_json,"w") do io; JSON.print(io,fixture_records); end
    c251020_output=joinpath(out,"c251020-fixture")
    rendered=render_c251020_replay_fixtures(fixture_json,c251020_output)
    @test Set(keys(rendered))==Set([1,2])
    @test all(path -> isfile(path) && filesize(path)>500,values(rendered))
    delete!(fixture_records[1]["objects"]["vl"]["properties"]["REdgeTable"],"Position")
    open(fixture_json,"w") do io; JSON.print(io,fixture_records); end
    @test_throws ArgumentError load_c251020_replay_fixtures(fixture_json)
    unknot_svg=read(plot_svg(run_vl_examples()["unknot"],joinpath(out,"unknot.svg")),String)
    @test occursin("<circle",unknot_svg)
    @test count("<path",unknot_svg)==1
    groups=run_ograph_plot_suite()
    @test sum(length,values(groups))==26
    # 19 VLExample + 26 other VirtualLink + 1 dev_URDiagram quiver.
    @test length(run_vl_plot_calls())+sum(length,values(groups))+1==46
    vl_calls=run_vl_plot_calls()
    @test vl_calls[13].link.weights==[[0,0,0,1,1,1,-1,0,1,0]]
    @test vl_calls[14].link.weights==[[1,-1,-1,0,1,0,-1,0]]
    @test vl_calls[15].link.weights==[[0,0,0,1,0,0,0,-1]]
    c250818=groups["projects/invariant/C250818MSTinv_uqsl2"]
    @test [call.link.weights for call in c250818]==[
        [[-1,-1,1,0]],
        [[1,-1,-1,0,1,0,-1,0]],
        [[-2,1]],
        [[-1,1,0,1,0,0,-1,1,0,0,0,0]],
        [[-1,0,1,0],[1,0,0,-1]],
        [[-1,0,1,2],[1,-2,0,-1]],
    ]
    c260125=groups["projects/invariant/C260125MSTinv_uqsl2"]
    @test c260125[2].link.gauss==[[1,-1,4,3,-3,-6,-2,2,6,5,-5,-4]]
    @test c260125[3].link.gauss==[[1,-1,4,3,-3,-6,-2,2,6,5,-5,-4,
                                  9,8,-8,-11,-7,7,11,10,-10,-9]]
    @test c260125[8].link.weights==[[0,0,0,1,0,0,0,-1]]
end

@testset "Topology helpers" begin
    B=bernstein(3,[0.0,0.5,1.0]); @test size(B)==(4,3); @test all(isapprox.(sum(B,dims=1),1))
    P=[0.0 1 2; 0.0 1 0]; @test bezier(0.5,P)[:,1]≈[1.0,0.5]
    moves=make_move_data(); @test Set(keys(moves))==Set(["MP-L","MP-R","02-L","02-R","CP-L","CP-R","PS-L","PS-R","BMP-L","BMP-R","B02-L","B02-R"])
    @test length(moves["MP-L"])==16
    @test length(moves["MP-R"])==16
    a2=only(filter(pattern -> pattern.parameter=="A2",moves["MP-L"]))
    @test a2.gauss==[[-1],[1,2],[-2]] && a2.orientation==[1,-1]
    s=SmallObject(3); s.value=4; @test s.value==4
    matrix=[2.0 1.0; 3.0 4.0]
    trace_result=only(run_C251030_trace_formula([matrix]))
    @test trace_result.trace≈-1/15
    @test trace_result.diagram.V==[1,3,-1,-2,0,2,4,-3,-4]
    standard_forms=run_C251030_standard_forms()
    @test length(standard_forms)==4
    @test evaluate.(standard_forms[2].adjustment,Ref(fill(1//1,4)))==[1//1 -1//1;-1//1 1//1]
    for form in standard_forms
        values=fill(2//1,form.order^2)
        numeric_matrix=fill(2//1,form.order,form.order)
        correction=Matrix{Rational{Int}}(I,form.order,form.order)-circshift(Matrix{Rational{Int}}(I,form.order,form.order),(1,0))
        expected=det(numeric_matrix+correction)
        @test evaluate(form.determinant,values)==expected
        @test evaluate(form.trace,values)==-1/expected
    end
    symbolic_reductions=run_C251030_symbolic_reductions()
    @test [diagram.V for diagram in symbolic_reductions]==[[4,-4],[3,-3],[2,-2]]
    @test all(validate_data, symbolic_reductions)
    numeric_probe=run_C251030_numeric_reduction_probe(Rational{Int}[2 1 0 3; 0 3 1 1; 1 0 4 2; 2 1 0 5])
    @test numeric_probe.C==-1//51
    @test numeric_probe.V==[13,-13]
    @test numeric_probe.W==Any[-113//17]
    @test numeric_probe.inverse_trace==-339//1
    final_matrix=run_C251030_final_matrix_workflow()
    @test final_matrix.order==4
    @test final_matrix.residual==zero(final_matrix.residual)
    @test final_matrix.observed==final_matrix.expected
    variants=run_C251030_knot_trace_variants()
    @test variants.first_direct.numerator.terms==Dict((-1,)=>1//1)
    @test variants.first_reduced.denominator.terms==Dict((0,)=>1//1,(1,)=>-1//1,(2,)=>1//1)
    @test variants.second_diagram.W[1]==zero(laurent_variable(:t))
    @test variants.second_trace_failure==DivideError
    twisted=run_C260319_unknots()
    @test length(twisted)==3
    @test length(twisted[3].orientation)==12
    @test twisted[2].gauss==[[1,2,-2,-3,5,6,-6,-7],
                              [-8,-5,7,8,-4,-1,3,4]]
    @test all(link -> validate_data(vl_to_urdiagram(link)),twisted)
end

@testset "C251020 rank-one Laurent workflow" begin
    t=laurent_variable(:t)
    @test evaluate(t^-2+3t-1,[2//1])==21//4
    @test evaluate(t^-2+3t-1,[big(2)//big(1)])==21//4
    @test limit_zero((t^2+3t)/(one(t)-t))==0
    @test limit_zero((t+t^2)/(t+t^2))==1
    x,y=laurent_variables(2;prefix="x")
    @test ((x+y)*(x-y))/(x-y)==LaurentRational(x+y)
    @test ((x^-1+y)*(x-y))/(x-y)==x^-1+y
    @test x/(x+y)+y/(x+y)==LaurentRational(one(x))
    @test ((x-1)*(x+y))/((x-1)*(y-1))==(x+y)/(y-1)
    @test ((x*y-1)*(x+y))/((x*y-1)*(y-1))==(x+y)/(y-1)
    @test zero(x)/(x-y)==LaurentRational(zero(x))
    result=run_C251020_knot_dsl2b()
    @test result.figure_eight.reduced.V==[2,-2]
    # `simplify2` cancels the common Laurent factors exactly:
    # C=-t/(t^2-3t+1), hence C^-1=-t+3-t^-1.
    inverse=result.figure_eight.inverse_coefficient
    @test inverse.numerator.terms==Dict((1,)=>-1//1,(0,)=>3//1,(-1,)=>-1//1)
    @test inverse.denominator.terms==Dict((0,)=>1//1)
    @test evaluate(inverse,[2//1])==1//2
    @test evaluate(result.figure_eight.reduced.C,[3//1])==-3//1
    @test evaluate(result.figure_eight.reduced.C,[5//1])==-5//11
    @test [form.numerator.terms for form in result.helper_cells]==[
        Dict((1,)=>-1//1),Dict((1,)=>-1//1,(2,)=>1//1,(3,)=>-1//1),
        Dict((-1,)=>1//1,(0,)=>-3//1,(1,)=>1//1),
    ]
    @test all(form -> form.denominator.terms==Dict((0,)=>1//1),result.helper_cells)
    @test result.hopf.coefficient==1
    @test result.hopf.alexander==LaurentRational(one(laurent_variables(2)[1]))
    hopf_numeric=URDiagram(); set_data_cvw!(hopf_numeric,1,[-2,1,0,2,-1],[1,1])
    @test alexander_calc_inv2_at(hopf_numeric,[2//1,3//1])==1
    multi=result.multi_rank_alexander
    borromean_input=_vl_to_urdiagram_signs(result.borromean)
    @test alexander_calc_inv2_at(borromean_input,[2//1,3//1,5//1])==evaluate(multi.borromean,[2//1,3//1,5//1])==4
    @test alexander_calc_inv2_at(borromean_input,[3//2,5//3,7//4])==evaluate(multi.borromean,[3//2,5//3,7//4])==1//6
    @test alexander_calc_inv2_at(borromean_input,[5//3,7//5,11//7])==evaluate(multi.borromean,[5//3,7//5,11//7])==16//175
    whitehead_input=result.whitehead_urdiagram
    @test alexander_calc_inv2_at(whitehead_input,[2//3,3//5])==evaluate(multi.whitehead,[2//3,3//5])==5//9
    @test alexander_calc_inv2_at(whitehead_input,[3//2,5//3])==evaluate(multi.whitehead,[3//2,5//3])==2//25
    @test alexander_calc_inv2_at(whitehead_input,[5//3,7//4])==evaluate(multi.whitehead,[5//3,7//4])==24//245
end

@testset "C250811 non-weighted MP moves" begin
    result=run_C250811_mp_moves()
    @test result.trials["A2"].gauss==[[-1,-2,2,-5,1,3,4,-4,-3,5]]
    @test result.trials["B2"].orientation==[1,-1,-1,1,-1]
    @test result.trials["B4"].gauss==[[-5,3,4,-4,-3,-2,5,-1,1,2]]
    @test result.trials["A4"].orientation==[-1,1,-1,1,-1]
    @test result.chained.gauss==[[1,-5,-3,2,-2,-1,3,6,4,-4,5,-6]]
    @test length(real_edge_table(result.chained))==12
    weighted=deepcopy(result.trials["A2"]); set_weight!(weighted,[zeros(10)])
    @test_throws ArgumentError mp_move!(weighted,"A2",[2,1])
    reductions=run_C250811_ur_reductions()
    @test all(validate_data,values(reductions))
    @test reductions.cell13.V==[-2,2] && reductions.cell13.C==-1 && reductions.cell13.W==Any[1.0]
    @test reductions.cell15.V==[-2,2] && trace3(reductions.cell15)==-1
    @test reductions.cell16.V==[-1,1] && trace3(reductions.cell16)==1
    sweep=run_C250811_epsilon_sweep()
    # MATLAB C250811 output.xml line 62: limit(value,ep,0) is eight copies of -1/4.
    @test sweep.values==fill(BigInt(-1)//BigInt(4),8)
    @test all(form -> form isa LaurentRational,sweep.rational_forms)
    conversions=run_C250811_conversion_cells()
    @test Set(keys(conversions))==Set([:hopf,:knot,:three_component,:lens,:s2xs1,:two_zero,:virtual])
    @test all(entry -> validate_data(entry.diagram),values(conversions))
    @test conversions[:hopf].diagram.V==[1,5,-5,-2,0,2,6,-6,-1]
    @test conversions[:three_component].trace==sym("(1/(1-(1+eps2)))/0")
    weighted_result=run_C250811_weighted_invariant()
    @test weighted_result.base==[-1,1,-1,0,0,0,0,0,0,0]
    # MATLAB C250811 output.xml line 88: SparseEx scalar -0.38197 + 1.7142e-13i.
    @test weighted_result.hopf≈-0.381966011250105 atol=1e-10
end


@testset "C250406 legacy tensor network" begin
    expression1,expression2,expression3=run_C250406_tensor_network()
    @test expression1==expression2==expression3
    @test size(expression1)==(4,4,4,4)
    stored=reshape(expression1,16,16)
    @test sum(abs,stored)==24
    @test stored[1,16]==-1
    @test stored[4,4]==1
    @test stored[11,[1,11,16]]==[1,1,1]
    @test stored[16,[6,11,16]]==[1,1,1]
    exploratory=run_C250406_exploratory_cells()
    @test sum(abs,exploratory.integral_cointegral)==0
    @test size(exploratory.coproduct_product)==(4,4)
    @test sum(abs,exploratory.coproduct_product)==2
    @test sum(abs,exploratory.integral_pair)==2
    @test sum(abs,exploratory.tp_twist)==sum(abs,exploratory.tm_twist)==2
    @test sum(abs,exploratory.counit_unit)==2
    @test sum(abs,exploratory.mixed)==6
    framing=run_C250406_sweedler_framing()
    @test framing.full==framing.s4_identity
    @test framing.full!=framing.changed
    @test all(sum(abs,tensor)==4 for tensor in values(framing))
    twists=run_C250406_sweedler_twists()
    @test all(size(tensor)==(4,4,4,4) for tensor in values(twists))
    @test [sum(abs,tensor) for tensor in values(twists)]==[24,24,40,40]
    # The source script compares all four; they are pairwise distinct cells.
    @test length(Set(collect(values(twists))))==4
    rings=run_C250406_ring_cells()
    @test sum(abs,rings.sweedler_ring)==0 && sum(abs,rings.sweedler_with_ring)==0
    @test maximum(abs,rings.small_borel_3_ring)<1e-10
    @test maximum(abs,rings.small_borel_5_ring)<1e-10
    @test sum(abs,rings.cyclic_ring)==25 && count(!iszero,rings.cyclic_ring)==25
    @test rings.cyclic_with_ring==rings.cyclic_without_ring
    @test all(values(rings.small_borel_3_checks))
    @test all(values(rings.cyclic_checks))
end


@testset "C250825 finite doubles" begin
    result=run_C250825_dsl2(order=3)
    @test all(values(result.hopf_checks))
    @test all(values(result.dual_checks))
    @test result.algebra.dimension==9
    @test result.dual.dimension==9
    @test result.drinfeld.dimension==81
    @test result.heisenberg.dimension==81
    @test length(result.drinfeld.multiplication.terms)>0
    @test length(result.heisenberg.multiplication.terms)>0
    @test result.pairing≈0
    cyclic=cyclic_group_algebra(5)
    @test all(values(verify_hopf(cyclic)))
    generator=basis_element(cyclic,2)
    @test (generator^5).coefficients≈basis_element(cyclic,1).coefficients
end

@testset "Remaining finite-Hopf Live Scripts" begin
    m_oinv=run_C260225_m_oinv()
    @test all(values(m_oinv.uqsl2_checks))
    @test all(values(m_oinv.cyclic_checks))
    @test m_oinv.weight.base==[1,-1,-1,0,1,0,-1,0]
    @test m_oinv.cyclic_invariant≈3

    sweedler=run_C250405_sweedler()
    @test all(values(sweedler.hopf_checks))
    @test sweedler.abalone_invariant≈-1
    @test sweedler.koda_invariant≈0 atol=1e-12

    exterior=run_C250823_extalg()
    @test all(values(exterior.hopf_checks))
    @test exterior.lens_values≈[4,4]
    @test exterior.abalone_value≈5
    @test exterior.s2xs1_values≈fill(28,4)
    @test exterior.koda_value≈40

    twisted=run_C260319_invariants()
    @test all(values(twisted.cyclic_checks))
    @test all(values(twisted.kac_paljutkin_checks))
    # MATLAB R2026a replay: repeated indices inside Tp/Tm must be diagonalized.
    @test twisted.cyclic≈[3,3,9]
    @test twisted.kac_paljutkin≈[8,48,8]

    kp=kac_paljutkin_algebra()
    x,y,z=basis_element(kp,2),basis_element(kp,3),basis_element(kp,5)
    @test z*x≈y*z
    @test z*y≈x*z
    @test 2*(z*z)≈basis_element(kp,1)+x+y-x*y

    uqsl2=run_C250818_uqsl2()
    @test all(values(uqsl2.hopf_checks))
    @test isapprox(uqsl2.lens_values,fill(-0.3819660112501-1.17557050458495im,2);atol=1e-10)
    @test isapprox(uqsl2.abalone_value,-0.809016994374947-0.587785252292473im;atol=1e-10)
    @test isapprox(uqsl2.s2xs1_values[1],0;atol=1e-10)
    @test isapprox(uqsl2.s2xs1_values[2],98.1762745781211-5.61284970724479im;atol=1e-9)
    @test all(value -> isapprox(value,0;atol=1e-10),uqsl2.s2xs1_values[3:end])
    @test isapprox(uqsl2.koda_value,0;atol=1e-10)
    @test isapprox(uqsl2.connected_values,[-0.809016994374947-0.587785252292473im,
                                            -0.809016994374947+0.587785252292473im];atol=1e-10)
    @test all(value -> isapprox(value,0;atol=1e-10),uqsl2.knot_complement_values)

    uqsl2_2026=run_C260125_uqsl2()
    @test all(values(uqsl2_2026.hopf_checks))
    @test uqsl2_2026.surgery_standard≈[-1,1,1,1]
    @test uqsl2_2026.surgery_alternate≈[-1,1,1,1]
    @test uqsl2_2026.abalone_sums≈[-1,-3,1]
    @test [row.hopf for row in uqsl2_2026.alternate_sums]≈[-1,-3,-9]
    @test [row.ur_prime for row in uqsl2_2026.alternate_sums]≈[-1,1,-1]

    lens_spaces=run_C260220_lens_spaces()
    @test length(lens_spaces.sweep)==31
    @test lens_spaces.weight.base==[1,-1,-1,0]
    @test lens_spaces.disk_edge_sums≈[0,0,-1]
    # MATLAB C260220 L(2,1): cellfun(@length,T.paths)' = [8 3 1].
    @test only(filter(row -> row.s==2 && row.t==1,lens_spaces.sweep)).path_lengths==[8,3,1]
    # Extracted MATLAB Live Script output: these are every sweep item that
    # completed before its local try/catch.  MATLAB printed `error` for the
    # seven omitted pairs; the Julia disk routine handles them instead.
    matlab_paths=Dict(
        (2,1)=>[8,3,1],
        (3,1)=>[10,3,4,1], (3,2)=>[7,3,7,1],
        (4,1)=>[12,3,4,4,1], (4,3)=>[8,3,8,4,1],
        (5,1)=>[14,3,4,4,4,1], (5,3)=>[9,3,9,4,4,1], (5,4)=>[9,3,9,4,4,1],
        (6,1)=>[16,3,4,4,4,4,1], (6,5)=>[10,3,10,4,4,4,1],
        (7,1)=>[18,3,4,4,4,4,4,1], (7,4)=>[11,3,11,4,4,4,4,1],
        (7,5)=>[11,3,11,4,4,4,4,1], (7,6)=>[11,3,11,4,4,4,4,1],
        (8,1)=>[20,3,4,4,4,4,4,4,1], (8,5)=>[12,3,12,4,4,4,4,4,1],
        (8,7)=>[12,3,12,4,4,4,4,4,1],
        (9,1)=>[22,3,4,4,4,4,4,4,4,1], (9,5)=>[13,3,13,4,4,4,4,4,4,1],
        (9,7)=>[13,3,13,4,4,4,4,4,4,1], (9,8)=>[13,3,13,4,4,4,4,4,4,1],
        (10,1)=>[24,3,4,4,4,4,4,4,4,4,1],
        (10,7)=>[14,3,14,4,4,4,4,4,4,4,1], (10,9)=>[14,3,14,4,4,4,4,4,4,4,1],
    )
    computed_paths=Dict((row.s,row.t)=>row.path_lengths for row in lens_spaces.sweep)
    @test all(computed_paths[pair]==paths for (pair,paths) in matlab_paths)
    @test setdiff(Set(keys(computed_paths)),Set(keys(matlab_paths)))==Set([(5,2),(7,2),(7,3),(8,3),(9,2),(9,4),(10,3)])
end
