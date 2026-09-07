using Topology

function _link(code,orientation;real=false)
    link=VirtualLink()
    real ? set_data!(link;rgauss=code,orientation=orientation) :
           set_data!(link;gauss=code,orientation=orientation)
    link
end

"""Named Virtual-Gauss inputs used by `VLExample.mlx` and its cited MST examples.

The `lens_3_2` entry records the exact parameters used by `VLExample.mlx`;
`lens_space_rgauss(s,t)` constructs the full coprime lens-space family.
`lens_3_1_mst2_example23` is the two-vertex presentation in MST2 Example 23,
distinct from the three-vertex lens-space-family presentation.
`matlab_input` preserves whether the source cell used `RGauss` or `odata`,
although both map to Julia's `set_data!(...; rgauss=...)` boundary.
"""
function vl_virtual_gauss_inputs()
    lens32=lens_space_rgauss(3,2)
    [
        (name="s3_abalone",rgauss=[[1,-1]],orientation=[1],
         matlab_input=:odata),
        (name="lens_2_1",rgauss=[[-2,-1,1,2]],orientation=[1,1],
         matlab_input=:RGauss),
        (name="lens_3_1_mst2_example23",rgauss=[[-1,1,-2,2]],orientation=[1,1],
         matlab_input=:RGauss,source="MST2 Example 23"),
        (name="lens_3_2",rgauss=lens32.rgauss,orientation=lens32.orientation,
         matlab_input=:RGauss),
        (name="poincare",rgauss=[[4,-5,3,-4,2,-3,1,-2,5,-1]],
         orientation=ones(Int,5),matlab_input=:RGauss),
        (name="s2xs1",rgauss=[[-1,-2,3,4,-4,-3,2,1]],
         orientation=[1,1,1,1],matlab_input=:odata),
        (name="odata_variant",rgauss=[[-1,2,3,4,-4,-3,-2,1]],
         orientation=[1,-1,1,1],matlab_input=:odata),
        (name="koda_fig13",rgauss=[[2,-1,-2,1,4,-3,-4,3]],
         orientation=[-1,1,-1,1],matlab_input=:RGauss),
    ]
end

"""List the names accepted by `vl_example_input`, `vl_example`, and `vl_example_disk_table`."""
vl_example_names()=sort!([input.name for input in vl_virtual_gauss_inputs()])

"""Return the stored raw input for a named `VLExample` link."""
function vl_example_input(name)
    key=String(name)
    for input in vl_virtual_gauss_inputs()
        input.name==key && return input
    end
    throw(KeyError(key))
end

"""Construct fresh `VirtualLink` objects for `vl_virtual_gauss_inputs()`."""
function vl_virtual_gauss_examples()
    Dict(input.name=>_link(input.rgauss,input.orientation;real=true)
         for input in vl_virtual_gauss_inputs())
end

"""Return a fresh named `VLExample` link from either a string or symbol."""
function vl_example(name)
    input=vl_example_input(name)
    _link(input.rgauss,input.orientation;real=true)
end

"""Compute the MATLAB-faithful disk table for a named `VLExample` link."""
vl_example_disk_table(name)=disk_table(vl_example(name))

"""Reproduce the 19 plotting calls in `topology/Manifold/VLExample.mlx`."""
function run_vl_plot_calls(;output_dir=nothing)
    calls=NamedTuple[]
    add(name,link;bending=nothing)=push!(calls,(name=name,link=link,bending=bending))
    virtual_gauss=vl_virtual_gauss_examples()

    add("unknot",_link([Int[]],Int[]))
    trefoil=_link([[-1,2,-3,1,-2,3]],[-1,-1,-1])
    add("trefoil",trefoil); add("trefoil_sage",deepcopy(trefoil))
    abalone=_link([[1,-1]],[1]); add("abalone",abalone); add("abalone_sage",deepcopy(abalone))
    add("hopf",_link([[1,-2],[-1,2]],[1,1]))
    add("borromean",_link([[-1,6,-4,3],[-2,4,-5,1],[-3,5,-6,2]],[1,1,1,-1,-1,-1]))
    add("whitehead",_link([[1,-4,5,-3],[3,-1,2,-5,4,-2]],[-1,-1,-1,1,1]))
    koda=[[-1,-2,-8,3,9,1,-4,6,5,4,-6,-5,-7,8,-3,7,2,-9]]
    add("koda",_link(koda,[-1,1,1,1,-1,1,1,1,1]))
    add("lens_2_1",_link([[1,-2,3,-3,2,-1]],[1,0,1]))
    s=3; arr=collect(1:2s-1); orientation=Int.((1 .- (-1).^arr)./2)
    signed_arr=arr.*(-1).^arr
    arr=vcat(-signed_arr,reverse(signed_arr))
    add("lens_3_1",_link([arr],orientation))
    add("lens_3_2",deepcopy(virtual_gauss["lens_3_2"]))
    koda_virtual=_link(koda,[-1,1,1,1,-1,0,0,0,0])
    set_weight!(koda_virtual,[[0,0,0,1,1,1,-1,0,1,0]])
    add("koda_virtual",koda_virtual)

    s2xs1=deepcopy(virtual_gauss["s2xs1"])
    set_weight!(s2xs1,[[1,-1,-1,0,1,0,-1,0]])
    add("s2xs1",s2xs1)
    weighted=deepcopy(s2xs1); set_weight!(weighted,[[0,0,0,1,0,0,0,-1]])
    add("s2xs1_weighted",weighted)
    calculation=_link([[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],[1,1,1,1,0,0])
    set_weight!(calculation,[[1,-1,-1,0,1,0,-1,0]])
    add("s2xs1_calculation",calculation)
    odata_variant=deepcopy(virtual_gauss["odata_variant"])
    set_weight!(odata_variant,[[1,0,-1,-1,1,0,0,-1]])
    add("odata_variant",odata_variant)
    trefoil_negative=_link([[1,-2,3,-1,2,-3]],[-1,-1,-1])
    set_weight!(trefoil_negative,[[1,-1,0,0,0,-1]])
    add("trefoil_negative",trefoil_negative)
    koda_fig13=deepcopy(virtual_gauss["koda_fig13"])
    set_weight!(koda_fig13,[[-1,0,0,0,0,0,0,0]])
    add("koda_fig13",koda_fig13)

    if output_dir !== nothing
        mkpath(output_dir)
        for (index,call) in enumerate(calls)
            plot_svg(call.link,joinpath(output_dir,"figure_$(lpad(index,3,'0'))_$(call.name).svg");
                     bending_numbers=call.bending)
        end
    end
    calls
end

function run_vl_examples(;output_dir=nothing)
    calls=run_vl_plot_calls(output_dir=output_dir)
    Dict(call.name=>call.link for call in calls)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_vl_plot_calls(output_dir=joinpath(@__DIR__,"..","artifacts","julia_plots","topology","Manifold","VLExample"))
end
