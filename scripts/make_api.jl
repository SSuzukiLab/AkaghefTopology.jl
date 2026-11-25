# cf: lib/topology.jl/docs/specs/api-contracts.md §§2–5 — contract inventory,
# object membership, source inheritance, and generated Markdown layout.
# cf: docs/COMPUTATIONS.md §Source reference comments — scoped cf: lines.

if !isdefined(@__MODULE__, :Topology)
    if isdefined(Main, :Topology)
        const Topology = Main.Topology
    else
        using Topology
    end
end
using Markdown
using TOML

include(joinpath(@__DIR__, "api_index.jl"))
APIIndex.generate(Topology)
