# cf: docs/specs/heegaard-snappy.md — fixed bounded inverse search and evidence boundaries.

"""Return safe native 2->3 neighbors with inverse replay certificates (unoriented)."""
hc_pachner_neighbors(input;kwargs...) = _hc_snappy_call(
    input isa HCClosedTriangulation ? hc_triangulation_data(input) : input;
    extra_args=["--pachner-neighbors"],kwargs...)

"""Search safe positive 2->3 moves and full branchings, retaining the best intermediate.

No resolved HC is claimed: a successful Koda stage returns its checked surface
cut certificate and the explicit remaining PL-identification/overlay stop.
The current native move set preserves vertex count; multiple vertices stop.
"""
function hc_prepare_inverse(input;max_depth::Int=2,max_states::Int=64,max_extra_tetrahedra::Int=2,
                            max_nodes::Int=100_000,max_solutions::Int=1024,
                            runner=nothing,timeout::Real=120)
    max_depth>=0 && max_states>0 && max_extra_tetrahedra>=0 || throw(ArgumentError("invalid move limits"))
    original=input isa HCClosedTriangulation ? hc_triangulation_data(input) : input
    t=_hc_explicit_triangulation(original)
    limits=Dict("max_depth"=>max_depth,"max_states"=>max_states,"max_extra_tetrahedra"=>max_extra_tetrahedra,
        "branching_nodes"=>max_nodes,"branching_solutions"=>max_solutions,"backend_timeout_seconds"=>timeout)
    result=Dict{String,Any}("status"=>"stopped","stop"=>nothing,"hc"=>nothing,"limits"=>limits,
        "attempts"=>Any[],"moves"=>Any[],"triangulation"=>original)
    if t.topology.vertices!=1
        result["stop"]="multiple_finite_vertices_unreduced"
        result["detail"]="safe native 2->3 preserves finite vertices; a certified 4->1 reducer is not available"
        return result
    end
    frontier=Any[(data=original,depth=0,moves=Any[])]
    seen=Set{String}(); processed=0; truncated=false; any_branch_budget=false
    while !isempty(frontier) && processed<max_states
        state=popfirst!(frontier); processed+=1
        accepted=Ref{Any}(nothing); failures=Any[]; enumeration_stop=nothing
        try
            hc_find_branching(state.data;max_nodes=max_nodes,max_solutions=max_solutions,accept=orders->begin
                try
                    accepted[]=hc_flow_certificate(state.data;orders=orders,max_nodes=max_nodes)
                    true
                catch e
                    if e isa HCTriangulationError && e.reason in
                        (:not_normal_ograph,:not_one_suture,:bw_not_two_disks,:suture_not_simple)
                        push!(failures,Dict("reason"=>string(e.reason),"detail"=>e.detail));false
                    else
                        rethrow()
                    end
                end
            end)
        catch e
            if e isa HCTriangulationError && e.reason in
                (:not_branchable_on_input,:branching_filter_exhausted,:branching_search_budget,:branching_solution_budget)
                enumeration_stop=string(e.reason)
                any_branch_budget |= e.reason in (:branching_search_budget,:branching_solution_budget)
            else
                result["stop"]="inverse_certificate_failure"
                result["detail"]=sprint(showerror,e)
                result["triangulation"]=state.data;result["moves"]=state.moves
                return result
            end
        end
        push!(result["attempts"],Dict("state"=>processed,"depth"=>state.depth,
            "tetrahedra"=>state.data["tetrahedra"],"branching_stop"=>enumeration_stop,"rejected_branchings"=>failures))
        if accepted[]!==nothing
            f=accepted[]
            result["status"]="koda_only";result["stop"]="polygon_normalization_and_reference_overlay_required"
            result["triangulation"]=state.data;result["moves"]=state.moves
            result["flow_certificate"]=f.certificate;result["koda_boundary_code"]=f.raw
            result["odata"]=Dict("gauss"=>f.encoded.value.gauss,"orientation"=>f.encoded.value.orientation)
            result["branching"]=f.encoded.branching
            result["correspondence"]=f.encoded.correspondence
            result["joining_region"]=f.region
            result["surface_cut"]=hc_surface_cut_system(f.diagram)
            result["processed_states"]=processed
            return result
        end
        state.depth>=max_depth && continue
        Int(state.data["tetrahedra"])>=t.n+max_extra_tetrahedra && continue
        neighbors=hc_pachner_neighbors(state.data;runner=runner,timeout=timeout)
        push!(seen,neighbors["isosig"])
        for child in neighbors["states"]
            child["isosig"] in seen && continue
            if length(frontier)+processed>=max_states
                truncated=true;break
            end
            push!(seen,child["isosig"])
            push!(frontier,(data=child["triangulation"],depth=state.depth+1,moves=vcat(state.moves,Any[child["move"]])))
        end
    end
    result["processed_states"]=processed
    result["stop"]=truncated ? "pachner_state_budget" : any_branch_budget ? "branching_budget" : "no_flow_branching_within_bounds"
    result
end

export hc_prepare_inverse, hc_pachner_neighbors
