using Topology

_suite_link(code,orientation;real=false)=begin
    link=VirtualLink()
    real ? set_data!(link;rgauss=code,orientation=orientation) :
           set_data!(link;gauss=code,orientation=orientation)
    link
end

function _knot_complement_abalone(weight)
    link=_suite_link([[1,2,-2,-3],[-4,-1,3,4]],[-1,1,-1,1];real=true)
    set_weight!(link,[weight[1:4],weight[5:8]])
    link
end

"""All non-VLExample o-graph plot instances reached by the scoped Live Scripts."""
function run_ograph_plot_suite(;output_root=nothing)
    groups=Dict{String,Vector{NamedTuple}}()
    add(group,name,link;bending=nothing,labels=:weight)=push!(get!(groups,group,NamedTuple[]),
        (name=name,link=link,bending=bending,labels=labels))

    s2xs1=_suite_link([[-1,-2,3,4,-4,-3,2,1]],[1,1,1,1];real=true)
    set_weight!(s2xs1,[[1,-1,-1,0,1,0,-1,0]])
    add("Docs/experiment/C260225M_Oinv","s2xs1",s2xs1)

    add("projects/invariant/C250811MSTinvUR","hopf",
        _suite_link([[1,-2],[-1,2]],[1,1]))
    poincare=_suite_link([[4,-5,3,-4,2,-3,1,-2,5,-1]],[1,1,1,1,1];real=true)
    add("projects/invariant/C250811MSTinvUR","poincare",poincare)

    lens21=_suite_link([[1,-3,2,-2,3,-1]],[1,1,0])
    set_weight!(lens21,[[-1,-1,1,0]])
    calculation=_suite_link([[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],[1,1,1,1,0,0])
    set_weight!(calculation,[[1,-1,-1,0,1,0,-1,0]])
    abalone=_suite_link([[1,-1]],[1];real=true)
    set_weight!(abalone,[[-2,1]])
    cs2=_suite_link([[1,-1,4,3,-3,-6,2,-2,6,5,-5,-4]],[1,1,-1,1,-1,1];real=true)
    set_weight!(cs2,[[-1,1,0,1,0,0,-1,1,0,0,0,0]])
    kc1=_knot_complement_abalone([-1,0,1,0,1,0,0,-1])
    kc2=_knot_complement_abalone([-1,0,1,2,1,-2,0,-1])
    for group in ["projects/invariant/C250818MSTinv_uqsl2"]
        add(group,"lens_2_1",deepcopy(lens21))
        add(group,"s2xs1_calculation",deepcopy(calculation))
        add(group,"abalone",deepcopy(abalone))
        add(group,"connected_sum_2",deepcopy(cs2))
        add(group,"knot_complement_weight_1",deepcopy(kc1))
        add(group,"knot_complement_weight_2",deepcopy(kc2))
    end
    for group in ["projects/invariant/C250823MSTinv_extalg"]
        add(group,"lens_2_1",deepcopy(lens21))
        add(group,"s2xs1_calculation",deepcopy(calculation))
    end

    add("projects/invariant/C251020KnotDsl2bfromUR","borromean",
        _suite_link([[-1,6,-4,3],[-2,4,-5,1],[-3,5,-6,2]],[1,1,1,-1,-1,-1]))
    add("projects/invariant/C251020KnotDsl2bfromUR","whitehead",
        _suite_link([[1,-4,5,-3],[3,-1,2,-5,4,-2]],[-1,-1,-1,1,1]))

    group="projects/invariant/C260125MSTinv_uqsl2"
    cs1=_suite_link([[1,-1]],[1];real=true)
    cs2_negative=_suite_link([[1,-1,4,3,-3,-6,-2,2,6,5,-5,-4]],
                             [1,1,-1,1,-1,1];real=true)
    set_weight!(cs2_negative,[[0,1,0,1,-1,0,0,0,0,0,0,0]])
    cs3=_suite_link([[1,-1,4,3,-3,-6,-2,2,6,5,-5,-4,9,8,-8,-11,-7,7,11,10,-10,-9]],
                    [1,1,-1,1,-1,1,1,-1,1,-1,1];real=true)
    set_weight!(cs3,[[1,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]])
    add(group,"connected_sum_1",cs1)
    add(group,"connected_sum_2",cs2_negative)
    add(group,"connected_sum_3",cs3)
    add(group,"lens_2_1",deepcopy(lens21))
    add(group,"s2xs1_calculation",deepcopy(calculation))
    add(group,"abalone",deepcopy(abalone))
    add(group,"connected_sum_2_again",deepcopy(cs2))
    weighted_s2xs1=deepcopy(s2xs1); set_weight!(weighted_s2xs1,[[0,0,0,1,0,0,0,-1]])
    add(group,"s2xs1_weighted",weighted_s2xs1)
    add(group,"knot_complement_weight_1",deepcopy(kc1))
    add(group,"knot_complement_weight_2",deepcopy(kc2))
    abalone_positive=_suite_link([[1,-1]],[1];real=true); set_weight!(abalone_positive,[[-2,1]])
    abalone_negative=_suite_link([[1,-1]],[-1];real=true); set_weight!(abalone_negative,[[1,-2]])
    add(group,"abalone_positive",abalone_positive)
    add(group,"abalone_negative",abalone_negative)

    lens_space=_suite_link([[-2,-1,1,2]],[1,1];real=true)
    set_weight!(lens_space,[[1,-1,-1,0]])
    add("projects/spine/example/C260220LensSpaces","lens_2_1",lens_space;labels=:all)

    if output_root !== nothing
        for (group,calls) in groups
            directory=joinpath(output_root,split(group,'/')...)
            mkpath(directory)
            for (index,call) in enumerate(calls)
                plot_svg(call.link,joinpath(directory,"figure_$(lpad(index,3,'0'))_$(call.name).svg");
                         bending_numbers=call.bending,edge_labels=call.labels)
            end
        end
    end
    groups
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_ograph_plot_suite(output_root=joinpath(@__DIR__,"..","artifacts","julia_plots"))
end
