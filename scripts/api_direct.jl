# Direct driver allowed by docs/ENVIRONMENT.md, "Invoking Julia from a
# delegated or sandboxed shell". Include the real source, without stubs.
empty!(LOAD_PATH)
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "Project.toml")), "@stdlib")
include(joinpath(@__DIR__, "..", "src", "Topology.jl"))
include(joinpath(@__DIR__, isempty(ARGS) ? "make_api.jl" : ARGS[1]))
