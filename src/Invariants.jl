"""
Generic names for presentation-independent 3-manifold invariants.

`first_homology`, `linking_form`, `Z_Kuperberg`, `mirror`, and `connected_sum`
have no methods in this module. `casson_walker.jl` gives Casson–Walker values
from surgery coefficients and Seifert data. Those functions do not read a Heegaard code.
"""
module Invariants

include("casson_walker.jl")

"""
    first_homology(x::HeegaardCode) -> Topology.HeegaardCodes.HCFirstHomology

Return the abelianization of the Heegaard presentation associated with a code.

# Arguments
`x` is a coherent `HeegaardCode` in the
[Heegaard data model](../specs/heegaard-code.md). The curve-family count
condition of `Heegaard.hc_first_homology` is checked; failure raises
`HCClosedDiagramExcluded` with `:alpha_beta_count_ne_genus`. Passage and
presentation errors propagate from that function. Complete cutting systems
are caller obligations for the closed-manifold interpretation and are not
certified by this count check.

# Returns
`HCFirstHomology` with `free_rank`, cyclic `torsion` orders, nonzero Smith
`invariant_factors`, `generator_family`, `relator_family` and `condition`
metadata, as in `Heegaard.HCFirstHomology`.

# Guarantees
1. The result is the abelianization in
   [Heegaard terminology §9.1](../specs/heegaard-terminology.md#91-群と不変量で使う用語),
   with the same value and metadata as `Heegaard.hc_first_homology(x)`.
   [proved/agent: by construction]
2. The input code is not modified. [proved/agent: by construction]

# Not guaranteed
The count gate does not prove that the supplied curves are complete cutting
systems. Without the geometric hypotheses, the result is a presented
abelian group, not a certified manifold invariant. There is no `VirtualLink`
method at this entry point; its closed normal o-graph pipeline is
[`ograph_pi1_presentation`](@ref) followed by [`ograph_h1`](@ref).

# Sources
cf:lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — Heegaard presentation and its abelianization.
cf:lib/topology.jl/src/Topology.jl first_homology; lib/topology.jl/src/heegaard/hc_compare.jl hc_first_homology — forwarding and checked domain.

# Examples
```jldoctest
julia> H = Topology.Heegaard; lc = H.LickorishCode(1, Dict("a"=>[1,2],"b"=>[-1,-2]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> h = H.HeegaardCode(lc, Dict("a"=>:alpha,"b"=>:beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing); r = first_homology(h);

julia> (r.free_rank, Tuple(r.torsion))
(0, (2,))
```

See also `Heegaard.hc_presentation`, [`linking_form`](@ref).

Role: invariant
"""
function first_homology end

"""
    linking_form(x::HeegaardCode; boundary_side=:alpha) -> Topology.HeegaardCodes.HCLinkingForm

Return the exact torsion-linking-form record of the supplied Heegaard code.

# Arguments
`x` and `boundary_side` have the domain of
`Heegaard.hc_linking_form`. The diagram, passage, meridian and
reference-lattice checks of that entry point are performed; their validation
exceptions propagate. Rejected linking data raise `HCLinkingInputError`
with the failing reason and detail. An unsupported boundary side raises
`ArgumentError`; malformed dimensions can raise `DimensionMismatch`.
The geometric data model is [Heegaard code](../specs/heegaard-code.md).

# Returns
`HCLinkingForm`: `free_rank`, torsion-generator orders `torsion`, rational
pairing `matrix`, integral generator lifts `generators`, and the integral
decomposition/Smith `certificate`, as in `Heegaard.HCLinkingForm`.

# Guarantees
1. The returned record equals `Heegaard.hc_linking_form(x; boundary_side)`;
   it is the exact coordinate pairing produced by that checked lattice
   construction. [proved/agent: by construction]
2. The input code is not modified. [proved/agent: by construction]

# Not guaranteed
This wrapper adds no geometric certification beyond its delegate. A complete
owning specification for the linking-form orientation and source theorem is
still missing; the delegated source records its literature identification
as unverified. Generator coordinates are presentation-dependent. No method
for `VirtualLink` or for an arbitrary abstract first-homology group is supplied.

# Sources
cf:lib/topology.jl/src/Topology.jl linking_form; lib/topology.jl/src/heegaard/hc_invariants_linking.jl hc_linking_form — the exact lattice construction and forwarding interface; geometric source identification remains unverified there.

# Examples
```jldoctest
julia> H = Topology.Heegaard; words = Dict("l1"=>[-9,-11,-13,-1,-10,-12], "m1"=>[14,-15,1,17,16], "alpha1"=>[-7,6,4,2,5,3,8,15], "beta1"=>[-8,10,-2,11,-16,-3,12,-4,13,-17,-5,9,-6,7,-14]);

julia> lc = H.LickorishCode(1, words, Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing); h = H.resolve_hc(lc, Dict("alpha1"=>:alpha,"beta1"=>:beta); disks=Dict("joining"=>["F13"]));

julia> f = linking_form(h); (f.free_rank, Tuple(f.torsion))
(0, (5,))
```

See also `first_homology`, `Heegaard.linking_pairing`.

Role: invariant
"""
function linking_form end

"""
    Z_Kuperberg(x, algebra::FiniteHopfAlgebra; kwargs...) -> ComplexF64

Return the normalized Kuperberg coordinate contraction of the supplied Heegaard presentation.

# Arguments
- `x` is a `HeegaardCode`, a path accepted by `Heegaard.read_hc`, a
  `HeegaardDiagram`, or a `PuncturedHeegaardDiagram`. The last form evaluates
  its stored closed diagram; it does not compute a relative invariant.
  The caller supplies coherent data in the
  [Heegaard specification](../specs/heegaard-code.md) and H5 of the
  [Hopf specification](../specs/hopf-invariants.md#numbered-conditions).
- The code/path route validates the code and its complete curve-family
  counts. Invalid counts raise `DomainError`; missing words or malformed
  passages raise `ArgumentError`; parser and Heegaard validation errors
  propagate. The diagram routes use the supplied certified record and do
  not repeat this code validation.
- `algebra` must supply valid finite Hopf tensors in H1 and H3. Those axioms
  and the literature identification are caller obligations. Numerical
  involutivity and a nonzero distinguished pairing are checked; failure
  raises `DomainError` identifying the failed requirement. Tensor-shape,
  allocation and arithmetic errors propagate.
- The code/path route accepts `strategy=:auto`, `:qa`, `:network`, `:none`
  or `:numeric`; unknown values raise `ArgumentError`. The diagram and
  punctured-diagram routes accept no keywords, including `strategy`.
  Other keywords raise `ArgumentError` about extra keywords.

# Returns
A `ComplexF64` scalar in the
[Kuperberg coordinate convention](../specs/hopf-invariants.md#kuperberg-coordinate-target).

# Guarantees
1. The returned scalar is that normalized coordinate contraction for the
   supplied data, using its entrance's current sign assignment; the intended
   HC/cellular relation is the [open sign point](../conventions/Algebra.md#kuperberg-sign-contradiction).
   [proved/agent: by construction]
2. The accepted code-route strategies give the same contraction value.
   Input diagram and algebra containers are not modified.
   [proved/agent: by construction]

# Conventions
[Tensor ports and integral fields](../conventions/Algebra.md).

# Not guaranteed
A scalar alone does not certify H1–H3, general manifold invariance, or
preservation under every `MoveAPI` move. The source crossing-sign and
normalization comparison recorded in the Hopf specification is unverified.
The [intended HC/cellular sign convention](../conventions/Algebra.md#kuperberg-sign-contradiction)
is undecided. Equality of `Z_Kuperberg(hc, algebra)` and
`Z_Kuperberg(forget_reference(realize_hc(hc)).diagram, algebra)` on convertible
inputs is not established by the label-set observation.
The calculation is numerical, with no certified error bound. In particular,
this is not a general exact-arithmetic or relative-boundary invariant.

# Sources
cf:[[Kup91#^kup91-sec-0105]] the Heegaard-circle contraction defining the literature target; the repository coordinate comparison remains as stated in the Hopf specification.
cf:[[Kup91#^kup91-blk-0138]] group-algebra interpretation used for the finite homomorphism-count checks.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc"));

julia> isapprox(Z_Kuperberg(path, cyclic_group_algebra(2)), 1; atol=1e-8)
true
```

See also [`invariant_expression`](@ref), `evaluate_invariant`.

Role: invariant
"""
function Z_Kuperberg end

"""
    mirror(x::Union{HeegaardCode,VirtualLink,KnotGaussDiagram}) -> Union{HeegaardCode,VirtualLink,KnotGaussDiagram}

Return the mirror presentation selected by the input's domain-specific convention.

# Arguments
- For `HeegaardCode`, the input requirements and face-mark transport checks
  are those of `Heegaard.hc_mirror`; ambiguous or absent face-token
  matches raise `KodaInputError`. This wrapper adds no validation.
- For `VirtualLink`, the checked domain is that of [`mst_sign_change`](@ref)
  in the [symmetry specification](../specs/ograph-symmetries.md#meaning-bearing-api).
  Rejected weighted, virtual or nonclosed input raises `ArgumentError`.
- For `KnotGaussDiagram`, see its specific `mirror` contract and
  [K1](../specs/knot-gauss.md#numbered-conditions).

# Returns
A new presentation of the same concrete domain type as `x`.

# Guarantees
1. The Heegaard and VirtualLink results are respectively `hc_mirror(x)`
   and `mst_sign_change(x)`; their contracts determine the represented
   structures and transported data. [proved/agent: by construction]
2. On the closed-normal VirtualLink domain, the geometric meaning is that
   of the named MST sign change in the symmetry specification.
   [proved/literature: [[MST25#^mst25-blk-0188]] Proposition A.1]
3. Input containers are not modified. [proved/agent: by construction]

# Not guaranteed
The three domains do not identify the same ancillary data under mirroring.
This wrapper supplies no framing lift, general Reidemeister-torsion
conjugation theorem, or convention-independent orientation translation.
The complete owning entry for Heegaard mirror face-mark transport is missing;
its delegate's geometric interpretation must be established separately.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1–§3.2 — domain-specific generic dispatch.
cf:[[MST25#^mst25-blk-0188]] Proposition A.1 — VirtualLink sign-change meaning.

# Examples
```jldoctest
julia> d = KnotGaussDiagram([1,2,-1,-2], [1,-1]);

julia> Topology.KnotGauss.canonical_hash(d) != Topology.KnotGauss.canonical_hash(mirror(d))
true
```

See also [`mst_sign_change`](@ref), `Heegaard.hc_mirror`.

Role: convert
"""
function mirror end

"""
    connected_sum(x, y; kwargs...) -> Union{VirtualLink,NamedTuple}

Return a connected-sum presentation in the supported Heegaard or o-graph domain.

# Arguments
- Two `HeegaardCode` inputs require the domain of `hc_connected_sum` in
  [the Heegaard construction specification](../specs/heegaard-snappy.md).
  Its `left_gap` and `right_gap` keywords are forwarded. Unsupported input,
  collar or splice data raise the delegate's `LCMoveError`; realization and
  other validation errors propagate.
- A `VirtualLink` input instead takes a component selection `y=(1,2)` and
  `transport_weights::Bool=true`, with the domain and errors of
  [`connected_sum!`](@ref). The caller supplies disjoint closed-normal
  summands under the [connected-sum preconditions](../specs/ograph-connected-sum.md#preconditions-and-postconditions).
  Selection, cut, weight and arithmetic errors propagate from that entry
  point. Here `y` is not a second `VirtualLink` summand.

# Returns
For Heegaard input, the named tuple `(lc,hc,certificate)` described by
`hc_connected_sum`, including its pending joining-disk status. For link
input, a new `VirtualLink` with the selected operation applied.

# Guarantees
1. Each result is the corresponding domain-specific construction; the
   Heegaard curve roles and o-graph local operation have the meanings
   stated in their linked specifications. [proved/agent: by construction]
2. For o-graph inputs satisfying the linked geometric preconditions, the
   local operation represents the oriented connected sum.
   [proved/literature: [[K10#^k10-sec-0033]] Figure 6 and Theorem 1.1]
3. The supplied inputs are not modified. [proved/agent: by construction]

# Not guaranteed
This is not an arbitrary-genus Heegaard constructor or a two-VirtualLink
binary sum. The Heegaard output does not acquire a certified joining disk.
The integer weight lift is not established as a canonical connected sum of
framing classes. No torsion multiplication rule is asserted by this wrapper.

# Sources
cf:lib/topology.jl/docs/specs/heegaard-snappy.md hc_connected_sum — supported Heegaard construction and certificate boundary.
cf:[[K10#^k10-sec-0033]] Figure 6 and Theorem 1.1 — oriented o-graph connected sum under the local specification.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); disjoint_sum!(v, deepcopy(v));

julia> w = connected_sum(v); is_closed(w)
(true, (true, true, true))
```

See also [`connected_sum!`](@ref), `Heegaard.hc_connected_sum`.

Role: convert
"""
function connected_sum end

export first_homology, linking_form, Z_Kuperberg, mirror, connected_sum,
    sawtooth, dedekind_sum, casson_walker_lens, casson_walker_knot_surgery,
    casson_walker_surgery_pup, casson_walker_brieskorn, casson_walker_seifert

end
