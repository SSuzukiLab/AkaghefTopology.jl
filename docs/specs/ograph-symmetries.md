# VirtualLink symmetries and the switched 0–2 move

Owner decision: 2026-09-24 conversation; implementation: 2026-09-27.
This specification uses the [existing passage/vertex sign convention](ograph-terminology.md).
It does not supply an integer or Z2 framing lift.

## Source and encoding

**Evidence.**
- literature: [[MST25#^mst25-blk-0186]] (three diagram involutions).
- literature: [[MST25#^mst25-blk-0188]] (Appendix A, Proposition A.1).
- literature: [[MST25#^mst25-fig-0189]] (crossing change and butterfly).
- literature: [[MST25#^mst25-blk-0062]] (Proposition 2.1, combed move equivalence).
- literature: [[MST25#^mst25-fig-0066]] (standard 0–2 figure).

The managed [MST25 Work](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/MST25.md)
owns the source text and figures; they are not copied into another bibliography here.

For real crossings, let `g` denote each Gauss component and `ε` the vertex-sign array.

| Diagram API (also available with `!`) | Source operation | Array action |
|---|---|---|
| `switch_crossings` | Γ× | `g → -g`, `ε → -ε` |
| `reflect_diagram` | Γ* | `g → g`, `ε → -ε` |
| `reverse_strands` | →Γ | reverse every component, `ε → ε` |

Reversing both incidences at each crossing leaves its vertex sign unchanged.
For cut strands the trailing zero remains last; reversal exchanges endpoints.
All primitives reject weighted inputs. `switch_crossings` leaves explicit virtual
visits unchanged; `reflect_diagram` rejects explicit virtual crossings because their
cyclic-order reflection is not represented by negating a stored zero sign.

`mirror_manifold!` remains a legacy diagram-level alias of `switch_crossings!`.
It is not plane reflection and does not certify normality. Unlike its old behavior,
it rejects weighted inputs instead of silently preserving weights.

## Meaning-bearing API

`mst_sign_change[!]` implements Γ×, representing `(M^op, -v^op)`.
`mst_reflect_reverse[!]` implements →(Γ*), representing `(M^op, v^op)`.
Both first require the existing `is_closed` check to accept the input and reject
weights and explicit virtual crossings. This guard is a computational validator,
not a formal proof of its own correctness. Cut diagrams remain available to the
diagram primitives but are outside these semantic entrypoints.

The superscript `op` means transport by an orientation-reversing diffeomorphism;
the separate minus negates the transported vector field. These are different
combed structures. Proposition A.1 supplies these identities, not a formula for
transforming integer framing weights. No such formula is claimed here.

## 0–2 parameters

Existing `param="dir"` / `"inv"` are unchanged. New `"dir-switched"` /
`"inv-switched"` exchange over/under at BOTH inserted crossings without reversing
strand directions. Both are unweighted-only, including exact-candidate application.
This extends the VirtualLink move API. The existing DB search selector still
enumerates standard 0–2 only; no new certificate/verifier support is claimed.

| Insertion | Local visits | Vertex signs |
|---|---|---|
| standard | `[-a,-b,b,a]` | `[+1,-1]` |
| switched | `[a,b,-b,-a]` | `[-1,+1]` |

```julia
v = lens_space_link(2,1)
w = move(v, "02"; param="dir-switched", vertices=[-2,-1])
v2 = move(w, "02"; param="inv-switched", vertices=[3,4])
```

Derivation (`proved/agent`, not an additional figure attributed to the source):
apply Γ× globally, perform the standard 0–2, then apply Γ× globally again.
Propositions 2.1 and A.1 imply preservation of the initial combed equivalence
class. Exterior crossing changes cancel; only the two new crossings are switched.
No framed lift follows from this argument.

There is deliberately no strand-direction-reversed parameter. The fragment
`[a,b,-b,-a]` with signs `[+1,-1]` is NOT the switched fragment above.
For the finite regression input `[-2,-1,1,2]`, signs `[1,1]`, inserting that
different fragment gives `[-2,3,4,-4,-3,-1,1,2]`, signs `[1,1,1,-1]`;
the current validator returns `(C1,C2,C3)=(true,true,false)`.
This is a counterexample to unconditional acceptance, not to all conditional
applications of a direction-reversed local picture.

## Verification boundary

[Executable tests](../../test/virtual_link_symmetries.jl) exercise the public APIs:
involutions, mutation/copy equivalence, cut handling, unsupported weighted input,
and both 0–2 inverses. The conjugation and closedness checks run at every directed
edge of the standard `lens_space_link(s,t)` presentations with `2 ≤ s ≤ 7`,
`1 ≤ t < s`, `gcd(s,t)=1`. Tests do not establish a complete invariant or certify
the entire library. Standard integer-weighted 0–2 transport is regression-tested.

Run from the research root with Julia 1.12.6 and the library project:

```powershell
julia --startup-file=no --pkgimages=no --project=lib/topology.jl -e 'include("lib/topology.jl/test/virtual_link_symmetries.jl")'
```

The source-resolution receipt is stored at
`threads/PJ260830/data/260927_ograph-symmetries-literature-receipt.json`.
Resolution checks source bindings; it does not prove the conjugation argument.
