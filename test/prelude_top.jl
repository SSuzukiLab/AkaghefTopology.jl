# Shared prelude of the top-level suite (formerly the head of test/runtests.jl).
using Test, Topology
const ROOT = normpath(joinpath(@__DIR__, ".."))
