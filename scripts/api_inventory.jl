# Read-only inventory and pilot probe; not an acceptance judge.
# cf: lib/topology.jl/docs/specs/api-contracts.md §§2–3 — exported types and pilot examples.
if isdefined(Main, :Topology)
    const Topology = Main.Topology
else
    using Topology
end
include(joinpath(@__DIR__, "api_index.jl"))

let
    problems = String[]
    _, registry = APIIndex.registry(Topology, problems)
    seen = Base.IdSet{Any}()
    for mod in APIIndex.modules(Topology), name in sort!(names(mod); by=string)
        isdefined(mod, name) || continue
        value = getfield(mod, name)
        value isa Type || continue
        value in seen && continue
        push!(seen, value)
        groups = APIIndex.objects_for(value, registry; ancestors=false)
        println("TYPE ", APIIndex.qualified(value), " => ", join(groups, ", "))
        length(groups) == 1 || push!(problems, string(value, " has ", length(groups), " objects"))
    end
    println("REGISTRY PROBLEMS: ", length(problems))
    foreach(println, problems)
end

let
    graph = Topology.certify_closed_ograph(Topology.OGraphState([1, -1], [1]))
    p = Topology.ograph_pi1_presentation(graph)
    println("S3 PRESENTATION")
    show(stdout, (p.generator_count, p.relators)); println()
    p = Topology.ograph_pi1_presentation(Topology.certify_closed_ograph(Topology.lens_space_link(2, 1)))
    h = Topology.ograph_h1(p)
    println("L21 H1")
    show(stdout, (h.free_rank, h.torsion_invariant_factors, h.h1_finite_order)); println()
end
