"""
Generic names for presentation-independent 3-manifold invariants.

`first_homology`, `linking_form`, `Z_Kuperberg`, `mirror`, and `connected_sum`
have no methods in this module. `casson_walker.jl` gives Casson–Walker values
from surgery coefficients and Seifert data. Those functions do not read a Heegaard code.
"""
module Invariants

include("casson_walker.jl")

"""
    first_homology(x)

First homology of the closed 3-manifold presented by `x`, as an abelian group.
A `MoveAPI` move that is declared to preserve the manifold leaves the value unchanged.
"""
function first_homology end

"""
    linking_form(x)

Linking form on the torsion of the first homology of the closed 3-manifold presented by `x`.
A `MoveAPI` move that is declared to preserve the manifold leaves the value unchanged.
"""
function linking_form end

"""
    Z_Kuperberg(x, algebra)

Kuperberg invariant of the closed 3-manifold presented by `x`, with values determined by `algebra`.
A `MoveAPI` move that is declared to preserve the manifold leaves the value unchanged.
"""
function Z_Kuperberg end

"""
    mirror(x)

Orientation reversal of the closed 3-manifold presented by `x`. `x` is unchanged.

The first homology is the same. The linking form is multiplied by −1.
The Reidemeister torsion is sent to its complex conjugate.
"""
function mirror end

"""
    connected_sum(x, y)

Connected sum of the closed 3-manifolds presented by `x` and `y`.
The inputs are unchanged.

The first homology is the direct sum. The linking form is the orthogonal sum.
The Reidemeister torsion is Theorem XII.1.2 of Turaev, *Torsions of
3-dimensional manifolds*, Birkhäuser 2002, p. 176, not the product of the
two torsions. A character that is nontrivial on both summands makes the
twisted complex non-acyclic. For ``b_1(M_1) = b_1(M_2) = 0`` the theorem states

```
τ(M₁ # M₂, e₁ # e₂) = τ(M₁, e₁) Σ_{H₁(M₂)} + τ(M₂, e₂) Σ_{H₁(M₁)}.
```

Both signs on this line are plus. When both first Betti numbers are positive,
the same theorem gives ``τ = 0``.
"""
function connected_sum end

export first_homology, linking_form, Z_Kuperberg, mirror, connected_sum,
    sawtooth, dedekind_sum, casson_walker_lens, casson_walker_knot_surgery,
    casson_walker_surgery_pup, casson_walker_brieskorn, casson_walker_seifert

end
