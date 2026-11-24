# Heegaard suite entry: runs the heegaard modules at the tier named by
# TOPOLOGY_TEST_TIER (default "unit"). Policy in test/TESTING.md. The former inline
# tests are test/heegaard/core.jl; fixtures are test/heegaard/prelude.jl.
include(joinpath(@__DIR__, "..", "run.jl"))
main(["heegaard"]) == 0 || error("heegaard tests failed; see FAILED line above")
