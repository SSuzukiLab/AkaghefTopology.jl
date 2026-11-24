# cf: docs/specs/heegaard-snappy.md — finite branching CSP and exact inverse boundary.
# cf: Benedetti–Petronio, AGT 14 (2014), §1.1 — full branching means compatible edge orders.

_hc_perm_sign(p) = isodd(sum(p[i] > p[j] for i in 1:4 for j in i+1:4)) ? -1 : 1

function _hc_explicit_triangulation(input)
    data = input isa HCClosedTriangulation ? hc_triangulation_data(input) : input
    n = Int(data["tetrahedra"])
    gluings = [(tet=Int(g["tet"]), omitted=Int(g["omitted"]),
                to_tet=Int(g["to_tet"]), to_omitted=Int(g["to_omitted"]),
                permutation=Tuple(Int.(g["permutation"]))) for g in data["gluings"]]
    topology = _hc_triangulation_topology(n, gluings)
    signs = haskey(data,"signs") ? Int.(data["signs"]) : zeros(Int,n)
    length(signs)==n || throw(HCTriangulationError(:invalid_orientation,"one sign per tetrahedron required"))
    supplied = haskey(data,"signs")
    supplied && !all(s->s in (-1,1),signs) && throw(HCTriangulationError(:invalid_orientation,"sign must be +/-1"))
    !supplied && (signs[1]=1)
    reached=Set([1])
    changed=true
    while changed
        changed=false
        for g in gluings
            s,t=g.tet+1,g.to_tet+1
            factor=-_hc_perm_sign(g.permutation)
            if signs[s]!=0 && signs[t]!=0
                signs[t]==factor*signs[s] || throw(HCTriangulationError(:nonorientable,"face signs disagree"))
            elseif signs[s]!=0
                signs[t]=factor*signs[s]; changed=true
            elseif signs[t]!=0
                signs[s]=factor*signs[t]; changed=true
            end
            if s in reached && !(t in reached); push!(reached,t); changed=true; end
            if t in reached && !(s in reached); push!(reached,s); changed=true; end
        end
    end
    length(reached)==n || throw(HCTriangulationError(:disconnected,"face pairing is disconnected"))
    (;n,gluings,topology,signs,orientation_source=supplied ? "input" : "chosen coherent orientation")
end

"""Use the existing Koda reverse after checking its sphere/one-suture flow domain.

The output is a punctured cellular Koda diagram. Marking it on the fixed
Lickorish reference surface remains the separate existing `odata_to_hc` boundary.
"""
function hc_triangulation_to_heegaard_koda(input; kwargs...)
    f=hc_flow_certificate(input;kwargs...)
    correspondence=Dict{String,Any}("genus"=>f.raw["genus"],"chi"=>f.raw["chi"],
        "intersections"=>f.raw["intersections"],"joining_region"=>f.region,
        "table_decode"=>Dict("sector_count"=>f.data["sector_count"],"spine_chi"=>f.data["spine_chi"]))
    conversion=KodaConversionResult(f.punctured,
        Dict{String,Any}("flow_witness"=>f.certificate,"source_identity"=>"labelled"),correspondence,:constructed)
    (;value=f.punctured,conversion,encoding=f.encoded,flow_witness=f.certificate)
end

"""Find full branching on the unchanged face pairing, or check supplied vertex orders.

An exhausted search proves only nonbranchability of this triangulation. Budget
exhaustion is a separate error. This function does not certify manifold links.
"""
function hc_find_branching(input; orders=nothing, max_nodes::Int=100_000,
                           accept=nothing, max_solutions::Int=1024)
    max_nodes>0 || throw(ArgumentError("max_nodes must be positive"))
    max_solutions>0 || throw(ArgumentError("max_solutions must be positive"))
    t=_hc_explicit_triangulation(input)
    permutations=parentmodule(@__MODULE__)._cp_recipe_permutations4()
    compatible(g,a,b) = [g.permutation[x+1] for x in a if x!=g.omitted] ==
                         [x for x in b if x!=g.to_omitted]
    if orders!==nothing
        length(orders)==t.n && all(o->sort(collect(o))==collect(0:3),orders) ||
            throw(HCTriangulationError(:invalid_branching,"one vertex permutation per tetrahedron required"))
        all(g->compatible(g,orders[g.tet+1],orders[g.to_tet+1]),t.gluings) ||
            throw(HCTriangulationError(:incompatible_branching,"face orders disagree"))
        return (;orders=[Int.(collect(o)) for o in orders],signs=t.signs,
                 orientation_source=t.orientation_source,search_nodes=0)
    end
    constraints=[(s=g.tet+1,t=g.to_tet+1,
        allowed=[compatible(g,a,b) for a in permutations,b in permutations]) for g in t.gluings]
    function propagate!(domains)
        changed=true
        while changed
            changed=false
            for c in constraints
                old_s,old_t=domains[c.s],domains[c.t]
                if c.s==c.t
                    new_s=[a for a in old_s if c.allowed[a,a]]
                    isempty(new_s) && return false
                    changed |= length(new_s)!=length(old_s)
                    domains[c.s]=new_s
                else
                    new_s=[a for a in old_s if any(b->c.allowed[a,b],old_t)]
                    new_t=[b for b in old_t if any(a->c.allowed[a,b],new_s)]
                    (isempty(new_s)||isempty(new_t)) && return false
                    changed |= length(new_s)!=length(old_s) || length(new_t)!=length(old_t)
                    domains[c.s],domains[c.t]=new_s,new_t
                end
            end
        end
        true
    end
    nodes=Ref(0); solutions=Ref(0)
    function search(domains)
        nodes[]+=1
        nodes[]>max_nodes && throw(HCTriangulationError(:branching_search_budget,string(max_nodes)))
        propagate!(domains) || return nothing
        unresolved=findall(d->length(d)>1,domains)
        if isempty(unresolved)
            solution=[permutations[only(d)] for d in domains]
            hc_find_branching(input;orders=solution)
            solutions[]+=1
            solutions[]>max_solutions && throw(HCTriangulationError(:branching_solution_budget,string(max_solutions)))
            return (accept===nothing || accept(solution)) ? solution : nothing
        end
        k=unresolved[argmin(length.(domains[unresolved]))]
        for candidate in domains[k]
            child=copy(domains); child[k]=[candidate]
            found=search(child)
            found!==nothing && return found
        end
        nothing
    end
    result=search([collect(1:24) for _ in 1:t.n])
    result===nothing && throw(HCTriangulationError(solutions[]==0 ? :not_branchable_on_input : :branching_filter_exhausted,
        "exhausted full branching CSP; $(solutions[]) solutions checked"))
    hc_find_branching(input;orders=result) # independent verification of the chosen assignment
    (;orders=result,signs=t.signs,orientation_source=t.orientation_source,search_nodes=nodes[],solutions_checked=solutions[])
end

include("hc_triangulation_surface.jl")
include("hc_triangulation_search.jl")

"""Enumerate branchings with explicit completeness; bounded stops retain checked solutions."""
function hc_each_branching(input;max_nodes::Int=100_000,max_solutions::Int=1024)
    solutions=Vector{Vector{Vector{Int}}}()
    stop=nothing; complete=false
    try
        hc_find_branching(input;max_nodes=max_nodes,max_solutions=max_solutions,
            accept=o->(push!(solutions,[collect(x) for x in o]);false))
    catch e
        e isa HCTriangulationError || rethrow()
        if e.reason in (:not_branchable_on_input,:branching_filter_exhausted)
            complete=true
        elseif e.reason in (:branching_search_budget,:branching_solution_budget)
            stop=e.reason
        else
            rethrow()
        end
    end
    (;solutions,complete,stop,max_nodes,max_solutions)
end

"""Recorded GKS branching, if this triangulation carries one. Empty is absent."""
function _hc_recorded_branch_orders(input)
    recorded = if input isa HCClosedTriangulation
        input.branch_orders
    elseif input isa AbstractDict
        get(input, "branch_orders", nothing)
    else
        nothing
    end
    (recorded === nothing || isempty(recorded)) && return nothing
    recorded
end

"""Encode a one-vertex triangulation using existing o-graph port conventions.

Returns the OData, the VirtualLink and explicit face/vertex/edge correspondence.
Acceptance requires the existing closed-o-graph checker and exact re-decoding
of every face map, including orientation. Flow and marked HC are separate.
`orders` is the branching. When omitted, a recorded `branch_orders` field is
used. That field is not an orientation.
"""
function hc_triangulation_to_odata(input; orders=nothing, max_nodes::Int=100_000)
    orders === nothing && (orders = _hc_recorded_branch_orders(input))
    t=_hc_explicit_triangulation(input)
    t.topology.vertices==1 || throw(HCTriangulationError(:multiple_finite_vertices,string(t.topology.vertices)))
    t.topology.euler_characteristic==0 || throw(HCTriangulationError(:nonzero_euler,string(t.topology)))
    branching=hc_find_branching(input;orders=orders,max_nodes=max_nodes)
    rank=[[findfirst(==(x),o)-1 for x in 0:3] for o in branching.orders]
    signs=[t.signs[i]*_hc_perm_sign(branching.orders[i]) for i in 1:t.n]
    # R1, triangulation to o-data. S is the branched tetrahedron sign above.
    # GKS stores the dual tetrahedron as minus the repository vertex sign.
    # K08 §3's left-handed frame is the one further minus against the
    # right-handed tetrahedron (∂Hα: outward normal of Hα first): vertex sign
    # = −S. R2, the standard-form entry, does not repeat it. There l is written
    # with +ε and β with −ε. F2 writes this same right-handed orientation into
    # the vertex numbers and leaves the minus here.
    vertex_signs=-signs
    T=parentmodule(@__MODULE__)
    # Read the existing table: do not copy a second set of local face conventions.
    rho=T._CP_GKS_FACE_OMITTED_VERTEX
    ports=[[only(c for c in 0:3 if rho[vertex_signs[i]][c+1]==rank[i][a+1]) for a in 0:3] for i in 1:t.n]
    successor=Dict{Int,Int}()
    expected=Dict{Tuple{Int,Int},Tuple{Int,Int,NTuple{4,Int}}}()
    pair_ports=NamedTuple[]
    for (index,g) in enumerate(t.gluings)
        s,u=g.tet+1,g.to_tet+1
        a,b=ports[s][g.omitted+1],ports[u][g.to_omitted+1]
        (a>=2)!=(b>=2) || throw(HCTriangulationError(:port_orientation,"both faces point the same way"))
        x,y=(isodd(a) ? s : -s),(isodd(b) ? u : -u)
        tail,head=a>=2 ? (x,y) : (y,x)
        haskey(successor,tail) && throw(HCTriangulationError(:reused_port,string(tail)))
        successor[tail]=head
        p=ntuple(i->rank[u][g.permutation[branching.orders[s][i]+1]+1],4)
        ar,br=rank[s][g.omitted+1],rank[u][g.to_omitted+1]
        invp=ntuple(i->only(j-1 for j in 1:4 if p[j]==i-1),4)
        expected[(s-1,ar)]=(u-1,br,p)
        expected[(u-1,br)]=(s-1,ar,invp)
        push!(pair_ports,(input_pair=index,source_port=a,target_port=b,tail=tail,head=head))
    end
    sort(collect(values(successor)))==sort(collect(keys(successor))) ||
        throw(HCTriangulationError(:invalid_passage_permutation,"incoming and outgoing visits differ"))
    remaining=Set(keys(successor)); words=Vector{Vector{Int}}()
    while !isempty(remaining)
        start=minimum(remaining); x=start; word=Int[]
        while x in remaining
            push!(word,x); delete!(remaining,x); x=successor[x]
        end
        x==start || throw(HCTriangulationError(:invalid_passage_permutation,"orbit does not close"))
        push!(words,word)
    end
    v=T.VirtualLink(); T.set_data!(v;rgauss=words,orientation=vertex_signs)
    decoded=try
        T.disk_tetrahedra(v)
    catch e
        throw(HCTriangulationError(:not_normal_ograph,sprint(showerror,e)))
    end
    decoded.disk_count==t.topology.edges || throw(HCTriangulationError(:correspondence_failure,"edge/sector counts differ"))
    actual=Dict{Tuple{Int,Int},Tuple{Int,Int,NTuple{4,Int}}}()
    for g in decoded.face_gluings
        p=g.permutation
        actual[(g.tet,g.omitted)]=(g.to_tet,g.to_omitted,p)
        actual[(g.to_tet,g.to_omitted)]=(g.tet,g.omitted,ntuple(i->only(j-1 for j in 1:4 if p[j]==i-1),4))
    end
    # disk_tetrahedra keeps sign = -vertex_sign. Face maps must still be the input's.
    expected==actual && [x.sign for x in decoded.tetrahedra]==-vertex_signs ||
        throw(HCTriangulationError(:correspondence_failure,"exact oriented face correspondence failed"))
    real_edge_ids=Dict(Tuple(e.crossing)=>e.id for e in T.real_edge_table(v))
    local_edges=decoded.ordered_edges
    input_edge_to_disk=[[decoded.tetrahedra[i].edges[only(j for j in 1:6
        if local_edges[j]==minmax(rank[i][a+1],rank[i][b+1]))]
        for (a,b) in local_edges] for i in 1:t.n]
    corr=Dict{String,Any}("tetrahedron_to_vertex"=>collect(1:t.n),"branch_orders"=>branching.orders,
        "input_vertex_to_rank"=>rank,"input_local_edge_to_disk"=>input_edge_to_disk,
        "face_pair_to_real_edge"=>[real_edge_ids[(p.tail,p.head)] for p in pair_ports],
        "face_occurrence_to_port"=>ports,"face_pair_ports"=>pair_ports,
        "local_rank_edge_to_disk"=>[collect(x.edges) for x in decoded.tetrahedra],
        "local_rank_edge_to_sheet_germ"=>[collect(x.local_nodes) for x in decoded.tetrahedra],
        "exact_oriented_face_correspondence"=>true,"orientation_source"=>branching.orientation_source)
    (;value=OData(words,vertex_signs,Dict{String,Any}("source"=>"explicit triangulation")),
      virtual_link=v,branching=branching,correspondence=corr,
      koda_status=(status=:not_computed,reason=:flow_witness_required))
end
