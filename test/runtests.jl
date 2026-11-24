# Pkg.test entry. Runs the tier named by TOPOLOGY_TEST_TIER (default "unit") for all
# modules through the tiered runner; policy in test/TESTING.md. The former contents
# of this file are split into test/prelude_top.jl, test/topology_core.jl and the
# entries of test/manifest.toml.
include(joinpath(@__DIR__, "run.jl"))
main(String[]) == 0 || error("tests failed; see FAILED line above")
