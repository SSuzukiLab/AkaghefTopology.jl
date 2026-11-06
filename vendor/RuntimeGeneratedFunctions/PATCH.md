# Julia 1.12 precompile compatibility patch

This directory is RuntimeGeneratedFunctions.jl 0.5.25, vendored under its MIT
license because Catlab.jl 0.17.6 cannot precompile it on Julia 1.12.6 in this
Windows package environment.

The local delta is restricted to the `PrecompileTools.@compile_workload` smoke
call in `src/RuntimeGeneratedFunctions.jl`: `increment(41)` is invoked through
`Base.invokelatest`. The runtime-generated-function API and cache logic are
unchanged. The unpatched package failed with a `generated_callfunc` world-age
`MethodError`; the patched package was verified by ordinary
`julia --startup-file=no --pkgimages=no --project=. -e "using Catlab"` and by
the Heegaard Catlab layout tests.

Upstream `master` still contained the direct call when checked on 2026-08-31.
Remove this source override after a released RuntimeGeneratedFunctions version
is verified to fix the Julia 1.12 workload.
