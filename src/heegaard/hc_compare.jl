# cf: lib/topology.jl/docs/specs/lickorish-code.md §3 — passage sign; −k on b and +k on a means the frame (b, a) is positive
# cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — intersection rows are α, columns are β; H₁ is the abelianization of the Heegaard presentation
# cf: threads/PJ260922/code/hc_matrix.py — (α, β) entry is +1 when + is on α and − is on β
# cf: lib/topology.jl/src/heegaard/groups.jl — smith_normal_form, abelianization_matrix, first_homology(::GroupPresentation); α generates, β is the relator
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3.1 — face numbers come from the rotation system, not from stored label integers

"""計算と比較。入口の引数は HeegaardCode だけである。"""

"""
    HCClosedDiagramExcluded(condition, witness) -> HCClosedDiagramExcluded

Store a named obstruction to the closed-diagram domain of an HC computation.

# Arguments
`condition::Symbol` and `witness::Dict{String,Any}`, or values convertible to
those types. Conversion errors propagate. Construction does not throw the
result or evaluate the named condition.

# Returns
`HCClosedDiagramExcluded` with the supplied condition tag and diagnostic
data. `showerror` displays the condition tag.

# Guarantees
1. Fields preserve supplied values after conversion; a typed witness
   dictionary is retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. When thrown by `hc_first_homology`, the record identifies that operation's
   rejected domain check. The distinction between curve data and a valid
   manifold diagram is owned by
   [the Heegaard conditions](../specs/heegaard-koda-conversion.md#51-曲線系と多様体). [proved/agent: by construction]

# Not guaranteed
Arbitrary conditions and witnesses are accepted. A count mismatch is not
itself a classification of the represented surface or manifold, and the
absence of this exception is not proof of complete cut systems.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = HCClosedDiagramExcluded(:alpha_beta_count_ne_genus, Dict{String,Any}("genus" => 1, "n_alpha" => 2, "n_beta" => 1));

julia> (e.condition, e.witness["genus"], e.witness["n_alpha"], e.witness["n_beta"])
(:alpha_beta_count_ne_genus, 1, 2, 1)
```

See also [`hc_first_homology`](@ref).

Role: construct
"""
struct HCClosedDiagramExcluded <: Exception
    condition::Symbol
    witness::Dict{String,Any}
end

function Base.showerror(io::IO, err::HCClosedDiagramExcluded)
    print(io, "HCClosedDiagramExcluded(", err.condition, ")")
end

"""
    HCNotComputed(reason, detail) -> HCNotComputed

Store a reason why an HC comparison computation did not finish its decision.

# Arguments
`reason::Symbol` and `detail::Dict{String,Any}`, or convertible values.
Conversion errors propagate. No reason tag or diagnostic content is checked.

# Returns
`HCNotComputed` with the reason and diagnostic data. `showerror` prints the
reason tag. Constructing the value does not throw it.

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed
   dictionary is retained and not modified. [proved/agent: by construction]
2. When raised by the canonical comparison search, this is a computation
   failure rather than a verdict that the diagrams or manifolds differ;
   see [comparison scope](../specs/heegaard-terminology.md#9-操作比較検証).
   [proved/agent: by construction]

# Not guaranteed
The record does not decide any equivalence relation or certify that retrying
will finish. A caller can construct it without any failed computation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> HCNotComputed(:canonical_tie_group, Dict{String,Any}()).reason
:canonical_tie_group
```

See also [`compare_hc`](@ref), [`hc_canonical`](@ref).

Role: construct
"""
struct HCNotComputed <: Exception
    reason::Symbol
    detail::Dict{String,Any}
end

function Base.showerror(io::IO, err::HCNotComputed)
    print(io, "HCNotComputed(", err.reason, ")")
end

"""
    HCIntersection(rows, columns, matrix, row_order, column_order, sign_convention) -> HCIntersection

Store an alpha–beta intersection matrix with curve labels and convention metadata.

# Arguments
`rows::Vector{String}`, `columns::Vector{String}`, `matrix::Matrix{Int}`,
`row_order::Symbol`, `column_order::Symbol`, `sign_convention::String`, or
values convertible to those types. Conversion errors propagate. Dimensions,
labels and convention metadata are not checked.

# Returns
`HCIntersection` with supplied row and column curve names, integer matrix,
ordering tags and sign-convention text. The associated mathematical notion is
[the intersection pairing](../specs/heegaard-terminology.md#91-群と不変量で使う用語);
interpret the stored signs and ordering using
[comparison coordinates](../specs/heegaard-code.md#comparison-records).

# Guarantees
1. Fields equal the supplied values after conversion; typed vectors and
   matrix are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. For output of [`hc_intersection`](@ref), the matrix records the signed
   alpha–beta intersection data of the input curve arrangement under its
   recorded sign convention. A closed
   manifold diagram is not required for this matrix to be returned. [proved/agent: by construction]

# Not guaranteed
The constructor does not compute intersections or establish complete cut
systems. An arbitrary matrix and convention string need not describe any
surface arrangement. Agreement with the ordered oriented intersection pairing
is not guaranteed by this record. This record alone is not manifold homology.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("a" => [1], "b" => [-1]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> h = H.HeegaardCode(lc, Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.hc_intersection(h); s = H.HCIntersection(r.rows, r.columns, r.matrix, r.row_order, r.column_order, r.sign_convention);

julia> (Tuple(s.rows), Tuple(s.columns), s.matrix[1,1], s.row_order)
(("a",), ("b",), 1, :sorted_component_name)
```

See also [`hc_intersection`](@ref), [`hc_invariant_factors`](@ref).

Role: construct
"""
struct HCIntersection
    rows::Vector{String}
    columns::Vector{String}
    matrix::Matrix{Int}
    row_order::Symbol
    column_order::Symbol
    sign_convention::String
end

"""
    HCInvariantFactors(factors, source, manifold_homology) -> HCInvariantFactors

Store Smith invariant factors and their declared interpretation.

# Arguments
`factors::Vector{Int}`, `source::Symbol`, `manifold_homology::Bool`, or
values convertible to those types. Conversion errors propagate. Positivity,
divisibility and the interpretation flag are not validated.

# Returns
`HCInvariantFactors` with the supplied factor list, source tag and
manifold-homology flag.

# Guarantees
1. Fields preserve supplied values after conversion; a typed factor vector
   is retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. For output of [`hc_invariant_factors`](@ref), these are the nonzero Smith
   factors of the input HC's intersection matrix, and `manifold_homology`
   is false. The relevant matrix is the
   [intersection object](../specs/heegaard-terminology.md#91-群と不変量で使う用語). [proved/agent: by construction]

# Not guaranteed
Construction does not compute a Smith form. A manually supplied true flag
does not establish a manifold interpretation. These factors alone omit the
free rank needed to specify first homology.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("a" => [1,2], "b" => [-1,-2]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> h = H.HeegaardCode(lc, Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.hc_invariant_factors(h); s = H.HCInvariantFactors(r.factors, r.source, r.manifold_homology);

julia> Tuple(s.factors)
(2,)
```

See also [`hc_invariant_factors`](@ref), [`hc_first_homology`](@ref).

Role: construct
"""
struct HCInvariantFactors
    factors::Vector{Int}
    source::Symbol
    manifold_homology::Bool
end

"""
    HCPresentationResult(presentation, generator_family, relator_family, generator_names, relator_names) -> HCPresentationResult

Store an HC group presentation with its generator and relator curve names.

# Arguments
`presentation::GroupPresentation`, `generator_family::Symbol`,
`relator_family::Symbol`, and `generator_names, relator_names::Vector{String}`,
or convertible values. Conversion errors propagate. No correspondence of
names, families or words is checked.

# Returns
`HCPresentationResult` with the presentation, the two family tags, and curve
names in the respective generator and relator orders.

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are retained and not modified. [proved/agent: by construction]
2. For output of [`hc_presentation`](@ref), the record associates its group
   presentation with the input's alpha and beta curve data. Its manifold
   interpretation requires the Heegaard hypotheses of
   [the presentation specification](../specs/heegaard-terminology.md#91-群と不変量で使う用語).
   [proved/agent: by construction]

# Not guaranteed
Construction does not establish complete cut systems, check a presentation,
or identify an arbitrary supplied group with a 3-manifold group. Matching
alpha and beta counts alone does not establish those hypotheses.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("a" => [1], "b" => [-1]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing); h = H.HeegaardCode(lc, Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.hc_presentation(h); s = H.HCPresentationResult(r.presentation, r.generator_family, r.relator_family, r.generator_names, r.relator_names);

julia> (s.generator_family, s.relator_family, H.abelianization_matrix(s.presentation) == transpose(H.hc_intersection(h).matrix))
(:alpha, :beta, true)
```

See also [`hc_presentation`](@ref), [`hc_first_homology`](@ref).

Role: convert
"""
struct HCPresentationResult
    presentation::GroupPresentation
    generator_family::Symbol
    relator_family::Symbol
    generator_names::Vector{String}
    relator_names::Vector{String}
end

"""
    HCFirstHomology(free_rank, torsion, invariant_factors, generator_family, relator_family, condition) -> HCFirstHomology

Store the abelianization data of an HC presentation and its domain-check label.

# Arguments
`free_rank::Int`, `torsion::Vector{Int}`, `invariant_factors::Vector{Int}`,
`generator_family::Symbol`, `relator_family::Symbol`, `condition::Symbol`,
or values convertible to those types. Conversion errors propagate. No rank,
divisibility, family or condition check is performed.

# Returns
`HCFirstHomology` with the supplied free rank, torsion orders, nonzero Smith
factors, generator and relator families, and condition tag. The homology
interpretation is specified in
[groups and invariants](../specs/heegaard-terminology.md#91-群と不変量で使う用語).

# Guarantees
1. Fields preserve supplied values after conversion; typed vectors are
   retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. Data returned by `hc_first_homology` describe the abelianization of
   `hc_presentation` for that input. If the input `H` is a valid closed Heegaard
   diagram with the required complete cut systems, the intended manifold
   interpretation is the first integral homology of `Mhat(H)` from
   [the Heegaard manifold specification](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
   [unverified: general identification of the computed presentation with the manifold group]

# Not guaranteed
The constructor does not compute homology. `hc_first_homology` checks curve
counts but does not establish complete cut systems; its condition tag alone
does not certify the manifold interpretation. Stored data can be inconsistent.
The general identification of the computed presentation with the manifold
group is not certified here.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — the manifold Mhat(H) associated with the input diagram.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("a" => [1,2], "b" => [-1,-2]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> h = H.HeegaardCode(lc, Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.hc_first_homology(h); s = H.HCFirstHomology(r.free_rank, r.torsion, r.invariant_factors, r.generator_family, r.relator_family, r.condition);

julia> (s.free_rank, Tuple(s.torsion))
(0, (2,))
```

See also [`hc_first_homology`](@ref), [`hc_presentation`](@ref).

Role: construct
"""
struct HCFirstHomology
    free_rank::Int
    torsion::Vector{Int}
    invariant_factors::Vector{Int}
    generator_family::Symbol
    relator_family::Symbol
    condition::Symbol
end

"""
    HCCanonical(genus, k_words, alpha_words, beta_words, other_words, places, lc_disks, hc_disks, lc_points, hc_points, lc_weights, hc_weights) -> HCCanonical

Store the word and annotation data of a canonical HC representation.

# Arguments
`genus::Int`; `k_words` and `other_words` are
`Vector{Pair{String,Vector{Int}}}`; `alpha_words` and `beta_words` are
`Vector{Vector{Int}}`; `places`, `lc_points`, `hc_points` are
`Vector{Pair{String,String}}`; `lc_disks`, `hc_disks` are
`Vector{Pair{String,Vector{String}}}`; `lc_weights`, `hc_weights` are
`Vector{Pair{String,Int}}`. Convertible values are accepted; conversion
errors propagate. No canonicalization is performed by this constructor.

# Returns
`HCCanonical` with the supplied genus, reference and other named words,
alpha and beta word lists, placements and LC/HC annotation lists.
The underlying data and comparison order belong to
[HC comparison records](../specs/heegaard-code.md#comparison-records).

# Guarantees
1. Fields preserve supplied values after conversion; typed vectors are
   retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. `==` compares all fields by value. For output of [`hc_canonical`](@ref),
   the record is the representation used by that operation's identity test. [proved/agent: by construction]

# Not guaranteed
Direct construction does not establish validity, ordering or canonicality.
Equality of manually assembled records is not a checked correspondence of
their source diagrams, an isotopy certificate or a manifold classification.
Mutating stored containers can change equality and hashing.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1" => [-1,5], "m1" => [1], "alpha1" => [-5,2,-3], "alpha2" => [4,-6], "beta1" => [-2,3], "beta2" => [-4,6]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> h = H.HeegaardCode(lc, Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> c = H.hc_canonical(h); d = H.HCCanonical(c.genus, c.k_words, c.alpha_words, c.beta_words, c.other_words, c.places, c.lc_disks, c.hc_disks, c.lc_points, c.hc_points, c.lc_weights, c.hc_weights);

julia> d == H.hc_canonical(h)
true
```

See also [`hc_canonical`](@ref), [`hc_canonical_witness`](@ref).

Role: construct
"""
struct HCCanonical
    genus::Int
    k_words::Vector{Pair{String,Vector{Int}}}
    alpha_words::Vector{Vector{Int}}
    beta_words::Vector{Vector{Int}}
    other_words::Vector{Pair{String,Vector{Int}}}
    places::Vector{Pair{String,String}}
    lc_disks::Vector{Pair{String,Vector{String}}}
    hc_disks::Vector{Pair{String,Vector{String}}}
    lc_points::Vector{Pair{String,String}}
    hc_points::Vector{Pair{String,String}}
    lc_weights::Vector{Pair{String,Int}}
    hc_weights::Vector{Pair{String,Int}}
end

function Base.:(==)(a::HCCanonical, b::HCCanonical)
    all(getfield(a, i) == getfield(b, i) for i in 1:fieldcount(HCCanonical))
end

function Base.hash(canonical::HCCanonical, h::UInt)
    for i in 1:fieldcount(HCCanonical)
        h = hash(getfield(canonical, i), h)
    end
    hash(HCCanonical, h)
end

"""
    HCMoveStep(kind, detail) -> HCMoveStep

Store a caller-supplied move tag and accompanying HC comparison data.

# Arguments
`kind::Symbol` and `detail::Dict{String,Any}`, or convertible values.
Conversion errors propagate. The constructor accepts every symbol and does
not check the contents of `detail`.

# Returns
`HCMoveStep` with the move-kind tag and supplied detail dictionary.

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed
   dictionary is retained and not modified. [proved/agent: by construction]
2. The record can be passed in `HCMoveSearchResult` as an assertion about
   [comparison move assertions](../specs/heegaard-moves.md#comparison-and-move-witnesses). It carries no computed
   source or target diagram. [proved/agent: by construction]

# Not guaranteed
The constructor and `compare_hc` do not replay or validate this move.
Its presence does not prove isotopy, move equivalence or homeomorphism.
The caller must supply a valid move witness independently of the tag.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> make(w) = H.HeegaardCode(H.LickorishCode(1, Dict("a" => w, "b" => -w), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing), Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.HCMoveSearchResult([H.HCMoveStep(:isotopy, Dict{String,Any}())]);

julia> H.compare_hc(make([1,-2]), make([1,-2,3,-4]); move_result=r).isotopy.status
:equal
```

See also [`HCMoveSearchResult`](@ref), [`compare_hc`](@ref).

Role: construct
"""
struct HCMoveStep
    kind::Symbol
    detail::Dict{String,Any}
end

"""
    HCMoveSearchResult(steps) -> HCMoveSearchResult

Store a sequence of caller-supplied move assertions for an HC comparison.

# Arguments
`steps::Vector{HCMoveStep}`, or a value convertible to that type. Conversion
errors propagate. Empty sequences and arbitrary step tags are accepted.
The caller is responsible for any claimed geometric witness.

# Returns
`HCMoveSearchResult` whose `steps` field contains the supplied sequence.

# Guarantees
1. An already-typed step vector is retained in order and not modified.
   [proved/agent: by construction]
2. The record is an input to the move-equivalence part of
   [comparison move assertions](../specs/heegaard-moves.md#comparison-and-move-witnesses);
   constructing it performs no search or equivalence check.
   [proved/agent: by construction]

# Not guaranteed
No endpoint, consecutive-step compatibility or mathematical validity is
checked. `compare_hc` can issue a move-based equality verdict from the step
tags alone; that verdict is conditional on the caller's witness. An empty
sequence supplies no such evidence to `compare_hc`.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> make(w) = H.HeegaardCode(H.LickorishCode(1, Dict("a" => w, "b" => -w), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing), Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.HCMoveSearchResult([H.HCMoveStep(:handle_slide, Dict{String,Any}("source" => "caller"))]);

julia> c = H.compare_hc(make([1,-2]), make([1,-2,3,-4]); move_result=r); (c.identity.status, c.isotopy.status, c.move_equivalence.status, c.manifold.status, c.manifold.reason)
(:different, :undetermined, :equal, :equal, :move_sequence)
```

See also [`HCMoveStep`](@ref), [`compare_hc`](@ref).

Role: convert
"""
struct HCMoveSearchResult
    steps::Vector{HCMoveStep}
end

"""
    HCVerdict(level, status, reason, evidence) -> HCVerdict

Store a verdict at one HC comparison level with its reason and supporting data.

# Arguments
`level::Symbol`, `status::Symbol`, `reason::Symbol`, and
`evidence::Dict{String,Any}`, or convertible values. Conversion errors
propagate. The constructor does not check any tag or evidence.

# Returns
`HCVerdict` with the comparison level, status, reason and evidence dictionary.
The intended levels are distinguished by
[the comparison specification](../specs/heegaard-terminology.md#9-操作比較検証).

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed
   evidence dictionary is retained and not modified.
   [proved/agent: by construction]
2. For results of [`compare_hc`](@ref), each verdict addresses its indicated
   comparison level. Agreement of the selected group invariants alone yields
   an undetermined manifold verdict. [proved/agent: by construction]

# Not guaranteed
An arbitrary `:equal` record is not proof of any equivalence. A move-based
verdict from `compare_hc` relies on the caller's move assertions, which are
not replayed or validated. Interpretation as a manifold comparison also
requires valid Heegaard diagrams.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> make(w) = H.HeegaardCode(H.LickorishCode(1, Dict("a" => w, "b" => -w), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing), Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> v = H.compare_hc(make([1,-2]), make([1,-2,3,-4])).manifold; w = H.HCVerdict(v.level, v.status, v.reason, v.evidence);

julia> (w.status, w.reason)
(:undetermined, :invariants_agree_not_identity)
```

See also [`HCComparison`](@ref), [`compare_hc`](@ref).

Role: construct
"""
struct HCVerdict
    level::Symbol
    status::Symbol
    reason::Symbol
    evidence::Dict{String,Any}
end

"""
    HCComparison(identity, isotopy, move_equivalence, manifold) -> HCComparison

Store the separate identity, isotopy, move and manifold verdicts of an HC comparison.

# Arguments
Four `HCVerdict` values, or values convertible to that type. Conversion
errors propagate. The constructor does not check the verdicts or their
agreement with each other.

# Returns
`HCComparison` with fields `identity`, `isotopy`, `move_equivalence` and
`manifold` retaining the supplied verdicts.

# Guarantees
1. Fields preserve the supplied values after conversion. Already-typed arguments are not modified. [proved/agent: by construction]
2. A result of [`compare_hc`](@ref) separates these comparison questions
   rather than treating agreement of selected invariants as manifold
   identification; see [comparison scope](../specs/heegaard-terminology.md#9-操作比較検証). [proved/agent: by construction]

# Not guaranteed
Direct construction does not compare diagrams or verify any evidence.
In particular, `:equal` in a supplied verdict is not a new proof of isotopy
or homeomorphism. Caller-supplied move records used by `compare_hc` are not
replayed by this constructor.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> make(w) = H.HeegaardCode(H.LickorishCode(1, Dict("a" => w, "b" => -w), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing), Dict("a" => :alpha, "b" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> c = H.compare_hc(make([1,-2]), make([1,-2,3,-4]));

julia> d = H.HCComparison(c.identity, c.isotopy, c.move_equivalence, c.manifold);

julia> (d.manifold.status, d.manifold.reason, d.isotopy.status, d.move_equivalence.status)
(:undetermined, :invariants_agree_not_identity, :undetermined, :undetermined)
```

See also [`compare_hc`](@ref), [`HCVerdict`](@ref).

Role: convert
"""
struct HCComparison
    identity::HCVerdict
    isotopy::HCVerdict
    move_equivalence::HCVerdict
    manifold::HCVerdict
end

"""
    HCNamedRelations(identical, orientation_reversal, mirror) -> HCNamedRelations

Store separate identity, orientation-reversal and mirror relation results for HC data.

# Arguments
`identical::Bool`, `orientation_reversal::Symbol`, `mirror::Symbol`, or
convertible values. Conversion errors propagate. Relation tags and their
consistency are not checked.

# Returns
`HCNamedRelations` with the identity flag and the two named-relation tags.

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. For output of [`hc_named_relations`](@ref), the relations concern the
   input codes and their explicitly transformed curve data, in the scope of
   [diagram comparison](../specs/heegaard-terminology.md#9-操作比較検証).
   [proved/agent: by construction]

# Not guaranteed
These fields do not classify underlying 3-manifolds. Arbitrary symbolic
tags can be constructed; construction does not establish any relation.
An orientation-reversal search may return `:undetermined`.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> words = Dict("l1" => [-1,5], "m1" => [1], "alpha1" => [-5,2,-3], "alpha2" => [4,-6], "beta1" => [-2,3], "beta2" => [-4,6]);

julia> lc = H.LickorishCode(1, words, Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing); h = H.HeegaardCode(lc, Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> r = H.hc_named_relations(h, H.hc_reverse_orientation(h, ["alpha1"])); s = H.HCNamedRelations(r.identical, r.orientation_reversal, r.mirror);

julia> (s.identical, s.orientation_reversal)
(false, :yes)
```

See also [`hc_named_relations`](@ref), [`compare_hc`](@ref).

Role: construct
"""
struct HCNamedRelations
    identical::Bool
    orientation_reversal::Symbol
    mirror::Symbol
end

const _HC_SIGN = "+1 when + is on alpha and - is on beta; -1 for the opposite; labels that are not alpha-beta do not contribute"
const _HC_MOVE_KINDS = (:isotopy, :handle_slide, :stabilization, :destabilization)
const _HC_TIE_LIMIT = 2_000_000

function _hc_k_set(genus::Integer)
    Set(_lc_k_names(genus))
end

function _hc_owners(words::AbstractDict)
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in words
        for x in word
            x == 0 && throw(KodaInputError(:hc_compare, "label 0 is not a passage", Dict("curve" => String(name))))
            lab = abs(x)
            dest = x > 0 ? plus : minus
            haskey(dest, lab) && throw(KodaInputError(
                :hc_compare, "passage sign is not unique", Dict("label" => lab, "curve" => String(name))))
            dest[lab] = String(name)
        end
    end
    owners = Dict{Int,Tuple{String,String}}()
    for lab in keys(plus)
        haskey(minus, lab) || throw(KodaInputError(
            :hc_compare, "label is missing a minus passage", Dict("label" => lab)))
        owners[lab] = (plus[lab], minus[lab])
    end
    for lab in keys(minus)
        haskey(plus, lab) || throw(KodaInputError(
            :hc_compare, "label is missing a plus passage", Dict("label" => lab)))
    end
    owners
end

function _hc_curve_lists(hc::HeegaardCode)
    genus = hc.lc.genus
    kn = _hc_k_set(genus)
    alphas = String[]
    betas = String[]
    for (name, role) in hc.roles
        haskey(hc.lc.words, name) || throw(KodaInputError(
            :hc_compare, "role names a missing component", Dict("curve" => String(name))))
        if role == :alpha
            push!(alphas, String(name))
        elseif role == :beta
            push!(betas, String(name))
        else
            throw(KodaInputError(:hc_compare, "invalid role", Dict("curve" => String(name), "role" => String(role))))
        end
    end
    others = String[]
    for name in keys(hc.lc.words)
        sname = String(name)
        sname in kn && continue
        haskey(hc.roles, name) && continue
        push!(others, sname)
    end
    sort!(alphas)
    sort!(betas)
    sort!(others)
    alphas, betas, others
end

function _hc_index(names::Vector{String})
    Dict(name => i for (i, name) in enumerate(names))
end

"""
    hc_intersection(hc::HeegaardCode) -> HCIntersection

Return the signed alpha–beta intersection matrix in HC coordinates.

# Arguments
`hc` has valid roles and paired passages under [the HC input convention](../specs/heegaard-code.md).
Role, word and passage-ownership errors propagate as `KodaInputError`.
The caller supplies geometrically valid curve data for an intersection-pairing
interpretation; complete cut systems and count-equals-genus are not checked.

# Returns
`HCIntersection` with row and column names, an integer `matrix`, order tags
and the sign-convention string. Coordinates follow
[comparison records and sign translation](../specs/heegaard-code.md#comparison-records).

# Guarantees
1. The result is the algebraic alpha–beta intersection array of the encoded
   curve data in the linked HC sign convention.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The array alone is not manifold homology. Its sign convention must be
translated before treating it as the ordered tangent intersection pairing;
no silent identification of these conventions is made.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#intersection-sign-translation — HC sign versus ordered intersection pairing.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> a = H.hc_intersection(hc); (Tuple(a.rows), Tuple(a.columns), a.matrix[1,1])
(("m1",), ("l1",), 1)
```

Role: invariant
"""
function hc_intersection(hc::HeegaardCode)
    alphas, betas, _ = _hc_curve_lists(hc)
    owners = _hc_owners(hc.lc.words)
    ai = _hc_index(alphas)
    bi = _hc_index(betas)
    matrix = zeros(Int, length(alphas), length(betas))
    for (_, (plus_c, minus_c)) in owners
        if haskey(ai, plus_c) && haskey(bi, minus_c)
            matrix[ai[plus_c], bi[minus_c]] += 1
        elseif haskey(ai, minus_c) && haskey(bi, plus_c)
            matrix[ai[minus_c], bi[plus_c]] -= 1
        end
    end
    HCIntersection(alphas, betas, matrix, :sorted_component_name, :sorted_component_name, _HC_SIGN)
end

"""
    hc_invariant_factors(hc::HeegaardCode) -> HCInvariantFactors

Return the nonzero Smith invariant factors of the HC intersection matrix.

# Arguments
`hc` has the domain of `hc_intersection`; its input errors and Smith-form
arithmetic errors propagate. The
[comparison coordinate convention](../specs/heegaard-code.md#comparison-records)
applies. No closed-manifold or complete-cut-system check is performed.

# Returns
`HCInvariantFactors` with positive `factors` including units, the
`:intersection_matrix` source tag and a false manifold-homology flag.

# Guarantees
1. The factors are those of the integer intersection matrix; zero diagonal
   entries are absent. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This record omits free rank and is not a certified manifold-homology result.
Use `hc_first_homology` with its separate domain and interpretation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — matrix factors versus first homology.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> Tuple(H.hc_invariant_factors(hc).factors)
(1,)
```

Role: invariant
"""
function hc_invariant_factors(hc::HeegaardCode)
    matrix = hc_intersection(hc).matrix
    smith = smith_normal_form(matrix)
    diagonal = [abs(smith[i, i]) for i in 1:min(size(smith)...)]
    HCInvariantFactors(filter(!iszero, diagonal), :intersection_matrix, false)
end

function _hc_relator(word::Vector{Int}, owners, ai)
    letters = Tuple{Int,Int}[]
    for x in word
        plus_c, minus_c = owners[abs(x)]
        if haskey(ai, plus_c) && !haskey(ai, minus_c)
            push!(letters, (ai[plus_c], 1))
        elseif haskey(ai, minus_c) && !haskey(ai, plus_c)
            push!(letters, (ai[minus_c], -1))
        end
    end
    letters
end

"""
    hc_presentation(hc::HeegaardCode) -> HCPresentationResult

Return the group presentation associated with the alpha and beta words of an HC.

# Arguments
`hc` satisfies the [HC input domain](../specs/heegaard-code.md).
Invalid roles and passage ownership raise `KodaInputError`; malformed words
and missing role curves can also raise lookup errors. Complete cut systems
and a valid manifold diagram are the caller's responsibility under
[the Heegaard conditions](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
Curve counts equal to genus are not checked here.

# Returns
`HCPresentationResult` containing a `GroupPresentation`, family tags and
generator/relator names in [the comparison coordinates](../specs/heegaard-code.md#comparison-records).
Word encoding is owned by [legacy presentation words](../specs/heegaard-terminology.md#legacy-cellular-records).

# Guarantees
1. Abelianizing the result gives the transpose of `hc_intersection(hc).matrix`.
   [proved/agent: by construction]
2. The result is the presentation attached to the supplied diagram; for a
   valid complete Heegaard diagram its intended group is that of the
   represented closed 3-manifold. [unverified: general geometric identification is not certified by this entry]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Successful output does not certify the manifold hypotheses of guarantee 2.
The stored relator words depend on stored starting points, and no group
simplification or isomorphism decision is performed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — Heegaard presentation and applicability.
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — generator, relator and sign coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha)); p = H.hc_presentation(hc);

julia> (Tuple(p.presentation.generators), Tuple(only(p.presentation.relators)))
(("m1",), ((1, 1),))
```

Role: invariant
"""
function hc_presentation(hc::HeegaardCode)
    alphas, betas, _ = _hc_curve_lists(hc)
    owners = isempty(hc.lc.words) ? Dict{Int,Tuple{String,String}}() : _hc_owners(hc.lc.words)
    ai = _hc_index(alphas)
    relators = [_hc_relator(hc.lc.words[name], owners, ai) for name in betas]
    meta = Dict{String,Any}(
        "generator_family" => "alpha",
        "relator_family" => "beta",
        "generator_names" => alphas,
        "relator_names" => betas,
        "generator_order" => "sorted_component_name",
        "relator_order" => "sorted_component_name",
        "sign_convention" => _HC_SIGN,
        "beta_ids" => "positions in relator_names, not CombinatorialHeegaardCode curve ids",
    )
    presentation = GroupPresentation(alphas, relators, collect(1:length(betas)), meta)
    HCPresentationResult(presentation, :alpha, :beta, alphas, betas)
end

"""
    hc_first_homology(hc::HeegaardCode) -> HCFirstHomology

Return the abelianization of the presentation associated with a closed-count HC.

# Arguments
`hc` must satisfy the [HC input conditions](../specs/heegaard-code.md).
Checked: both curve counts equal genus, or `HCClosedDiagramExcluded` with
`:alpha_beta_count_ne_genus` is raised. Input, presentation and Smith-form
errors propagate. To interpret the result as a closed 3-manifold's H₁, the
caller supplies the complete cut systems of
[the Heegaard domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).

# Returns
`HCFirstHomology` with `free_rank`, torsion factors, nonzero
`invariant_factors`, generator/relator family tags and the checked-domain tag.

# Guarantees
1. The returned abelian group is the abelianization of `hc_presentation(hc)`.
   [proved/agent: by construction]
2. On the registered genus-one S³ fixture it is the trivial group.
   [verified-computationally: test/heegaard/hc_search.jl, genus-1 S3 and omitted selections]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The count check does not certify complete cut systems or identify an arbitrary
input with a closed 3-manifold. The supplied diagram's topological validity
is required for interpreting this presentation invariant as that manifold's H₁.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — presentation abelianization and manifold scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> h = H.hc_first_homology(hc); (h.free_rank, Tuple(h.torsion))
(0, ())
```

Role: invariant
"""
function hc_first_homology(hc::HeegaardCode)
    alphas, betas, _ = _hc_curve_lists(hc)
    genus = hc.lc.genus
    if length(alphas) != genus || length(betas) != genus
        throw(HCClosedDiagramExcluded(
            :alpha_beta_count_ne_genus,
            Dict{String,Any}(
                "genus" => genus,
                "n_alpha" => length(alphas),
                "n_beta" => length(betas),
                "alpha_names" => alphas,
                "beta_names" => betas,
            ),
        ))
    end
    presented = hc_presentation(hc)
    homology = first_homology(presented.presentation)
    HCFirstHomology(
        homology.free_rank,
        homology.torsion,
        homology.invariant_factors,
        :alpha,
        :beta,
        :alpha_beta_count_eq_genus,
    )
end

function _hc_s3_tables()
    perms = [(1, 2, 3), (2, 1, 3), (1, 3, 2), (3, 2, 1), (2, 3, 1), (3, 1, 2)]
    mul = zeros(Int, 6, 6)
    inv = zeros(Int, 6)
    for i in 1:6, j in 1:6
        p, q = perms[i], perms[j]
        composed = ntuple(k -> q[p[k]], 3)
        mul[i, j] = findfirst(==(composed), perms)
    end
    for i in 1:6
        p = perms[i]
        inverse = ntuple(k -> findfirst(==(k), p), 3)
        inv[i] = findfirst(==(inverse), perms)
    end
    mul, inv
end

const _HC_S3_MUL, _HC_S3_INV = _hc_s3_tables()

function _hc_eval_s3(relator, assign)
    acc = 1
    for (gen, exp) in relator
        step = exp >= 0 ? assign[gen] : _HC_S3_INV[assign[gen]]
        for _ in 1:abs(exp)
            acc = _HC_S3_MUL[acc, step]
        end
    end
    acc
end

function _hc_count_s3(relators, n::Int)
    n == 0 && return BigInt(1)
    used = falses(n)
    for relator in relators
        for (gen, _) in relator
            1 <= gen <= n || throw(KodaInputError(
                :hc_hom_count, "relator names an unknown generator", Dict("generator" => gen)))
            used[gen] = true
        end
    end
    active = [i for i in 1:n if used[i]]
    remap = Dict(old => i for (i, old) in enumerate(active))
    mapped = Vector{Tuple{Int,Int}}[]
    for relator in relators
        isempty(relator) && continue
        push!(mapped, [(remap[gen], exp) for (gen, exp) in relator])
    end
    free = n - length(active)
    length(active) == 0 && return BigInt(6)^free
    assign = zeros(Int, length(active))
    nodes = Ref(0)
    count = Ref(0)
    function rec(i::Int)
        if i > length(active)
            count[] += 1
            return
        end
        for g in 1:6
            nodes[] += 1
            nodes[] > _HC_TIE_LIMIT && throw(HCNotComputed(
                :hom_search_budget,
                Dict{String,Any}("group" => "S3", "generators" => n, "budget" => _HC_TIE_LIMIT),
            ))
            assign[i] = g
            ok = true
            for relator in mapped
                any(letter -> letter[1] > i, relator) && continue
                if _hc_eval_s3(relator, assign) != 1
                    ok = false
                    break
                end
            end
            ok && rec(i + 1)
        end
    end
    rec(1)
    count[] * BigInt(6)^free
end

"""
    hc_hom_count(presented::HCPresentationResult, group::Symbol) -> BigInt

Count homomorphisms from an HC presentation to a supported finite group.

# Arguments
`presented.presentation` must contain valid generator indices and integer
exponents in [the presentation format](../specs/heegaard-terminology.md#legacy-cellular-records).
The caller supplies a valid presentation, including when there are no generators.
`group` is `:Z2`, `:Z3`, `:Z5` or `:S3`; other values raise `ArgumentError`.
Unknown generator references in the S₃ search raise `KodaInputError`, and
search-budget exhaustion raises `HCNotComputed(:hom_search_budget)`.
Abelianization and arithmetic errors propagate.

# Returns
The exact nonnegative number of homomorphisms as a `BigInt`.

# Guarantees
1. For a valid presentation the result counts homomorphisms of its presented
   group to the selected target. [proved/agent: by construction]
2. The presentation is not modified. [proved/agent: by construction]

# Not guaranteed
The count is not a complete group or manifold invariant, and does not require
curve counts equal to genus or certify a closed-manifold presentation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — group presentation and finite-group comparison.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> p = H.hc_presentation(hc); (H.hc_hom_count(p,:Z2), H.hc_hom_count(p,:S3))
(1, 1)
```

Role: invariant
"""
function hc_hom_count(presented::HCPresentationResult, group::Symbol)
    presentation = presented.presentation
    if group in (:Z2, :Z3, :Z5)
        modulus = group == :Z2 ? 2 : group == :Z3 ? 3 : 5
        homology = first_homology(presentation)
        exponent = homology.free_rank + count(d -> d % modulus == 0, homology.torsion)
        return BigInt(modulus)^exponent
    elseif group == :S3
        return _hc_count_s3(presentation.relators, length(presentation.generators))
    else
        throw(ArgumentError("unsupported hom target: $group"))
    end
end

"""
    hc_hom_count(hc::HeegaardCode, group::Symbol) -> BigInt

Count homomorphisms from an HC's associated presentation to a supported finite group.

# Arguments
`hc` satisfies the [HC input conditions](../specs/heegaard-code.md).
`group` is `:Z2`, `:Z3`, `:Z5` or `:S3`, otherwise `ArgumentError` is raised.
Presentation errors propagate and bounded S₃ search can raise
`HCNotComputed(:hom_search_budget)`. Complete cut systems are the caller's
responsibility for a manifold interpretation.

# Returns
The homomorphism count as a `BigInt`.

# Guarantees
1. The result is the count for `hc_presentation(hc)` and the selected target.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Neither equal counts nor successful computation establish manifold or group
isomorphism. Curve counts need not equal genus.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — presentation-based comparison.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> (H.hc_hom_count(hc,:Z2), H.hc_hom_count(hc,:S3))
(1, 1)
```

Role: invariant
"""
hc_hom_count(hc::HeegaardCode, group::Symbol) = hc_hom_count(hc_presentation(hc), group)

function _hc_permutations(items::Vector{String})
    n = length(items)
    n == 0 && return [String[]]
    n > 8 && throw(HCNotComputed(:canonical_tie_group, Dict{String,Any}("n" => n)))
    out = Vector{String}[]
    used = falses(n)
    current = Vector{String}(undef, n)
    function rec(depth::Int)
        if depth > n
            push!(out, copy(current))
            return
        end
        for i in 1:n
            used[i] && continue
            used[i] = true
            current[depth] = items[i]
            rec(depth + 1)
            used[i] = false
        end
    end
    rec(1)
    out
end

function _hc_emit_rotation(word::Vector{Int}, start::Int, label_map::Dict{Int,Int}, next_id::Int)
    n = length(word)
    n == 0 && return Int[], label_map, next_id
    local_map = copy(label_map)
    nid = next_id
    emitted = Vector{Int}(undef, n)
    for i in 1:n
        x = word[mod1(start + i - 1, n)]
        lab = abs(x)
        if !haskey(local_map, lab)
            local_map[lab] = nid
            nid += 1
        end
        emitted[i] = sign(x) * local_map[lab]
    end
    emitted, local_map, nid
end

include(joinpath(@__DIR__, "hc_compare_form.jl"))
include(joinpath(@__DIR__, "hc_compare_canonical.jl"))


"""
    hc_identical(a::HeegaardCode, b::HeegaardCode) -> Bool

Compare two oriented, marked HC records by their canonical forms.

# Arguments
`a` and `b` must satisfy the comparison domain in
[HC comparison](../specs/heegaard-code.md#comparison-records).
The caller supplies valid arrangements and annotations; this is not full HC
validation. Canonicalization checks can raise `KodaInputError` for malformed
incidence or annotations and `HCNotComputed` for a tie group or search budget
that exceeds the supported limits. Such failures are not returned as `false`.

# Returns
`Bool`: equality of the two `HCCanonical` records.

# Guarantees
1. The result is exactly `hc_canonical(a) == hc_canonical(b)` under the
   named comparison conventions. [proved/agent: by construction]
2. Relabeling, cyclic rotation and permitted component renaming preserve
   the comparison in the registered examples.
   [verified-computationally: test/heegaard/hc_compare.jl, canonical identity examples]
3. Neither input is modified. [proved/agent: by construction]

# Conventions
[Canonical correspondence and ordering](../conventions/HeegaardCodes.md#comparison-conventions).

# Not guaranteed
This is not a test of homeomorphism of the represented 3-manifolds, nor does
it return a realized curve correspondence. Orientation reversal and mirror
equivalence are separate comparisons.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> h = HeegaardCode(lc, Dict("l1" => :beta, "m1" => :alpha), Dict(), Dict(), Dict(), nothing);

julia> hc_identical(h, h)
true
```

See also [`hc_identity_correspondence`](@ref), [`hc_named_relations`](@ref).

Role: condition
"""
hc_identical(a::HeegaardCode, b::HeegaardCode) = hc_canonical(a) == hc_canonical(b)

function _hc_replace_words(hc::HeegaardCode, words::Dict{String,Vector{Int}})
    lc = hc.lc
    new_lc = LickorishCode(lc.genus, words, copy(lc.places), copy(lc.disks), copy(lc.points), copy(lc.weights), nothing)
    HeegaardCode(new_lc, copy(hc.roles), copy(hc.disks), copy(hc.points), copy(hc.weights), nothing)
end

"""
    hc_reverse_orientation(hc::HeegaardCode, curves) -> HeegaardCode

Return an HC with selected curve traversals reversed and their incident passage signs changed.

# Arguments
`curves` is an iterable of names convertible to `String`; every name must
occur in `hc.lc.words`, checked with `KodaInputError` naming a missing curve.
The input must satisfy [LC passage conditions](../specs/lickorish-code.md);
conversion and copying errors propagate. Geometric annotation compatibility
is the caller's responsibility.

# Returns
A new HC with changed words, copied roles and annotations, and cleared
source paths.

# Guarantees
1. Each selected curve has reversed traversal and negated passages; an
   unselected curve has its passages negated at labels incident to a selected
   curve. [proved/agent: by construction]
2. Paired passage signs remain opposite and the input is not modified.
   [proved/agent: by construction]

# Not guaranteed
At an intersection of two selected curves, both stored passages are negated
once, not twice. The geometric interpretation as independent simultaneous
orientation reversals is therefore not established here; see [the simultaneous-reversal question](../specs/heegaard-code.md#simultaneous-curve-reversal). Annotation supports
are copied rather than geometrically revalidated.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#intersection-sign-translation — passage signs and ordered intersection pairing.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> r = H.hc_reverse_orientation(hc,["m1"]); (Tuple(r.lc.words["l1"]), Tuple(r.lc.words["m1"]))
((1,), (-1,))
```

Role: convert
"""
function hc_reverse_orientation(hc::HeegaardCode, curves)
    names = Set(String.(curves))
    for name in names
        haskey(hc.lc.words, name) || throw(KodaInputError(
            :hc_reverse_orientation, "missing curve", Dict("curve" => name)))
    end
    labels = Set{Int}()
    for name in names
        for x in hc.lc.words[name]
            push!(labels, abs(x))
        end
    end
    words = Dict{String,Vector{Int}}()
    for (name, word) in hc.lc.words
        if String(name) in names
            words[String(name)] = Int[-x for x in reverse(word)]
        else
            words[String(name)] = Int[abs(x) in labels ? -x : x for x in word]
        end
    end
    _hc_replace_words(hc, words)
end

function _hc_face_cycle_key(corners)
    n = length(corners)
    n == 0 && return NTuple{0,Int}()
    labels = [Int(corner[1]) for corner in corners]
    best = nothing
    for base in (labels, reverse(labels))
        for shift in 0:(n - 1)
            cyc = ntuple(j -> base[mod1(shift + j, n)], n)
            if best === nothing || cyc < best
                best = cyc
            end
        end
    end
    best
end

function _hc_unique_face_keys(words)
    faces, = _lc_faces(words)
    groups = Dict{Any,Vector{String}}()
    keys = Dict{String,Any}()
    for face in faces
        key = _hc_face_cycle_key(face.corners)
        push!(get!(groups, key, String[]), face.id)
        keys[face.id] = key
    end
    unique = Dict{String,Any}()
    ambiguous = Set{String}()
    for face in faces
        if length(groups[keys[face.id]]) == 1
            unique[face.id] = keys[face.id]
        else
            push!(ambiguous, face.id)
        end
    end
    unique, ambiguous
end

"""Face ids transported by the sign-invariant label cycle. A cycle shared by two faces is not a key."""
function _hc_mirror_face_map(before_words, after_words)
    before_unique, ambiguous = _hc_unique_face_keys(before_words)
    after_unique, = _hc_unique_face_keys(after_words)
    by_key = Dict(key => id for (id, key) in after_unique)
    mapping = Dict{String,String}()
    for (id, key) in before_unique
        haskey(by_key, key) && (mapping[id] = by_key[key])
    end
    mapping, ambiguous
end

function _hc_rewrite_disks(disks, mapping, ambiguous)
    out = Dict{String,Vector{String}}()
    for (name, refs) in disks
        rewritten = String[]
        for ref in refs
            token = String(ref)
            haskey(mapping, token) || throw(KodaInputError(:hc_mirror,
                "mirrored disk face has no unique sign-invariant",
                Dict{String,Any}("disk" => String(name), "face" => token, "ambiguous" => token in ambiguous)))
            push!(rewritten, mapping[token])
        end
        out[String(name)] = rewritten
    end
    out
end

function _hc_rewrite_face_token(value, mapping, ambiguous)
    token = String(value)
    haskey(mapping, token) && return mapping[token]
    token in ambiguous && throw(KodaInputError(:hc_mirror,
        "mirrored disk face has no unique sign-invariant",
        Dict{String,Any}("face" => token, "ambiguous" => true)))
    token
end

function _hc_rewrite_face_dict(items, mapping, ambiguous)
    Dict{String,String}(String(name) => _hc_rewrite_face_token(value, mapping, ambiguous) for (name, value) in items)
end

"""
    hc_mirror(hc::HeegaardCode) -> HeegaardCode

Return the passage-sign mirror of an HC with transported face annotations.

# Arguments
`hc` must have valid paired passages and transportable annotations under
[the HC surface convention](../specs/heegaard-code.md).
Ambiguous marked-face transport raises `KodaInputError(:hc_mirror)`;
malformed word and annotation errors propagate. This entry does not run a
complete manifold-diagram certification.

# Returns
A new HC with mirrored words, copied roles and weights, transported places,
disks and points, and cleared source paths.

# Guarantees
1. Applying the operation changes the encoded surface orientation while
   retaining curve traversal order in the LC passage convention.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
General annotations whose support is not a direct face token need not be
transported geometrically. The output is not a certificate that the input
represents any particular manifold or that it is amphichiral.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#hc-mirror-convention — passage mirror, annotation transport and invariant-coordinate translation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> m = H.hc_mirror(hc); (Tuple(m.lc.words["l1"]), Tuple(m.lc.words["m1"]))
((1,), (-1,))
```

Role: convert
"""
function hc_mirror(hc::HeegaardCode)
    words = Dict{String,Vector{Int}}(String(name) => Int[-x for x in word] for (name, word) in hc.lc.words)
    mapping, ambiguous = _hc_mirror_face_map(hc.lc.words, words)
    lc = hc.lc
    new_lc = LickorishCode(lc.genus, words,
        _hc_rewrite_face_dict(lc.places, mapping, ambiguous),
        _hc_rewrite_disks(lc.disks, mapping, ambiguous),
        _hc_rewrite_face_dict(lc.points, mapping, ambiguous),
        copy(lc.weights), nothing)
    HeegaardCode(new_lc, copy(hc.roles),
        _hc_rewrite_disks(hc.disks, mapping, ambiguous),
        _hc_rewrite_face_dict(hc.points, mapping, ambiguous),
        copy(hc.weights), nothing)
end

function _hc_orientation_related(a::HeegaardCode, b::HeegaardCode)
    names = sort!(collect(String.(keys(a.lc.words))))
    length(names) > 12 && return :undetermined
    n = length(names)
    for mask in 1:(1 << n) - 1
        chosen = String[names[i] for i in 1:n if (mask >> (i - 1)) & 1 == 1]
        hc_identical(hc_reverse_orientation(a, chosen), b) && return :yes
    end
    :no
end

"""
    hc_named_relations(a::HeegaardCode, b::HeegaardCode) -> HCNamedRelations

Compare HC identity, curve-orientation changes and the passage-sign mirror separately.

# Arguments
`a` and `b` must be accepted by the canonical, orientation-reversal and
mirror operations in [the comparison domain](../specs/heegaard-code.md#comparison-records).
Their input, annotation and search-budget errors propagate.

# Returns
`HCNamedRelations` with a Boolean identity result, an orientation-change
status and a mirror status.

# Guarantees
1. The identity field reports canonical diagram identity; the mirror field
   reports identity after mirroring `a`. [proved/agent: by construction]
2. The orientation field tests nonempty curve-orientation changes of `a`,
   and is `:undetermined` when `a` has more than twelve curves.
   [proved/agent: by construction]
3. Neither input is modified. [proved/agent: by construction]

# Not guaranteed
A `:no` orientation result does not negate identity under an empty change,
and none of these relations is a complete manifold-equivalence test.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — scope of diagram identity.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> H.hc_named_relations(hc,hc).identical
true
```

Role: condition
"""
function hc_named_relations(a::HeegaardCode, b::HeegaardCode)
    HCNamedRelations(hc_identical(a, b), _hc_orientation_related(a, b), hc_identical(hc_mirror(a), b) ? :yes : :no)
end

function _hc_move_levels(move_result::HCMoveSearchResult)
    isempty(move_result.steps) && return Symbol[]
    kinds = [step.kind for step in move_result.steps]
    all(kind -> kind in _HC_MOVE_KINDS, kinds) || return Symbol[]
    levels = Symbol[:move_equivalence, :manifold]
    all(kind -> kind == :isotopy, kinds) && push!(levels, :isotopy)
    levels
end

_hc_move_levels(::Nothing) = Symbol[]

function _hc_invariant_snapshot(hc::HeegaardCode)
    presented = hc_presentation(hc)
    hom = Dict{Symbol,BigInt}()
    hom_error = Dict{Symbol,Symbol}()
    for group in (:Z2, :Z3, :Z5, :S3)
        try
            hom[group] = hc_hom_count(presented, group)
        catch err
            err isa HCNotComputed || rethrow()
            hom_error[group] = err.reason
        end
    end
    homology = nothing
    homology_error = nothing
    try
        homology = hc_first_homology(hc)
    catch err
        err isa HCClosedDiagramExcluded || rethrow()
        homology_error = err.condition
    end
    (; hom, hom_error, homology, homology_error)
end

function _hc_invariants_differ(left, right)
    for group in (:Z2, :Z3, :Z5, :S3)
        haskey(left.hom, group) && haskey(right.hom, group) && left.hom[group] != right.hom[group] &&
            return true, Symbol(:hom_, group)
    end
    if left.homology !== nothing && right.homology !== nothing
        if left.homology.free_rank != right.homology.free_rank || left.homology.torsion != right.homology.torsion
            return true, :h1
        end
    end
    false, :none
end

function _hc_agreeing(left, right)
    agreeing = String[]
    if left.homology !== nothing && right.homology !== nothing &&
            left.homology.free_rank == right.homology.free_rank &&
            left.homology.torsion == right.homology.torsion
        push!(agreeing, "h1")
    end
    for group in (:Z2, :Z3, :Z5, :S3)
        haskey(left.hom, group) && haskey(right.hom, group) && left.hom[group] == right.hom[group] &&
            push!(agreeing, string("hom_", group))
    end
    agreeing
end

function _hc_invariant_evidence(left, right)
    evidence = Dict{String,Any}(
        "agreeing_invariants" => _hc_agreeing(left, right),
        "hom" => Dict(string(group) => haskey(left.hom, group) ? left.hom[group] : "not_computed" for group in (:Z2, :Z3, :Z5, :S3)),
    )
    if left.homology !== nothing
        evidence["h1"] = Dict("free_rank" => left.homology.free_rank, "torsion" => left.homology.torsion)
    else
        evidence["h1_condition"] = left.homology_error
    end
    evidence
end

function _hc_equal(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :equal, reason, evidence)
end

function _hc_different(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :different, reason, evidence)
end

function _hc_undetermined(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :undetermined, reason, evidence)
end

"""
    compare_hc(a::HeegaardCode, b::HeegaardCode; move_result=nothing) -> HCComparison

Return HC identity, isotopy, move-equivalence and manifold comparison verdicts with their recorded evidence.

# Arguments
`a` and `b` must satisfy `hc_identical` and, when distinct, the invariant
snapshot domain; canonicalization, validation and invariant errors propagate.
A canonicalization budget failure raises `HCNotComputed` rather than a
negative identity verdict. Optional `move_result::HCMoveSearchResult` is a
supplied assertion record under [the comparison-witness boundary](../specs/heegaard-moves.md#comparison-and-move-witnesses);
its endpoint compatibility and geometric validity are the caller's responsibility.

# Returns
An `HCComparison` with four `HCVerdict` records carrying status, reason and
evidence for identity, isotopy, move equivalence and manifold comparison.

# Guarantees
1. Canonical identity yields equal verdicts at all four levels. Otherwise,
   the implemented invariant discrepancies and supplied move claims determine
   the coarse verdicts; invariant agreement alone does not yield manifold
   equality. [proved/agent: by construction]
2. The inputs are not modified and no move search or replay is performed.
   [proved/agent: by construction]

# Not guaranteed
An unchecked `move_result` can cause an asserted equality without a validated
geometric witness. An undetermined verdict is not inequivalence. The invariant
producers' mathematical identification limits apply.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §comparison-and-move-witnesses — comparison assertions and validated moves are distinct.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc,Dict("m1"=>:alpha,"l1"=>:beta),Dict(),Dict(),Dict(),nothing); compare_hc(h,h).manifold.status
:equal
```

See also [`hc_identical`](@ref).

Role: condition
"""
function compare_hc(a::HeegaardCode, b::HeegaardCode; move_result::Union{Nothing,HCMoveSearchResult}=nothing)
    identical = hc_identical(a, b)
    if identical
        evidence = Dict{String,Any}("basis" => "canonical_form")
        return HCComparison(
            _hc_equal(:identity, :canonical_form, evidence),
            _hc_equal(:isotopy, :identical, evidence),
            _hc_equal(:move_equivalence, :identical, evidence),
            _hc_equal(:manifold, :identical, evidence),
        )
    end
    identity = _hc_different(:identity, :canonical_form)
    left = _hc_invariant_snapshot(a)
    right = _hc_invariant_snapshot(b)
    differ, differ_reason = _hc_invariants_differ(left, right)
    witnessed = _hc_move_levels(move_result)
    move_evidence = Dict{String,Any}(
        "steps" => move_result === nothing ? Any[] : [Dict("kind" => String(step.kind)) for step in move_result.steps],
    )
    function coarse(level::Symbol)
        differ && return _hc_different(level, differ_reason, _hc_invariant_evidence(left, right))
        level in witnessed && return _hc_equal(level, :move_sequence, move_evidence)
        _hc_undetermined(level, :no_move_search, _hc_invariant_evidence(left, right))
    end
    manifold = if differ
        _hc_different(:manifold, differ_reason, _hc_invariant_evidence(left, right))
    elseif :manifold in witnessed
        _hc_equal(:manifold, :move_sequence, move_evidence)
    else
        _hc_undetermined(:manifold, :invariants_agree_not_identity, _hc_invariant_evidence(left, right))
    end
    HCComparison(identity, coarse(:isotopy), coarse(:move_equivalence), manifold)
end

export HCClosedDiagramExcluded, HCNotComputed,
    HCIntersection, HCInvariantFactors, HCPresentationResult, HCFirstHomology, HCCanonical,
    HCCanonicalWitness, HCIdentityCorrespondence,
    HCMoveStep, HCMoveSearchResult, HCVerdict, HCComparison, HCNamedRelations,
    hc_intersection, hc_invariant_factors, hc_presentation, hc_first_homology, hc_hom_count,
    hc_canonical, hc_canonical_witness, hc_identity_correspondence, hc_apply_correspondence,
    hc_identical, hc_reverse_orientation, hc_mirror, hc_named_relations, compare_hc
