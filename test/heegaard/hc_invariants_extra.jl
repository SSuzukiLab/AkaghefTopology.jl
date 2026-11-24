# Reidemeister torsion from the Heegaard Fox matrix, Massuyeau arXiv:1003.2517 Lemma 5.3.
# Lens-space values are Lemma 5.7 with ε_p = +1. The quotient is F^× / ±φ(H₁).
using Test
using Topology
import AkaghefAlgebra

const HCT_H = Topology.HeegaardCodes
const HCT_TE = AkaghefAlgebra.TaftExact
const HCT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const HCT_CASES = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-usecases")

# Straight (p, q)-curve on the square torus. Same events as
# threads/PJ260922/data/261114-usecases/generate.jl::class_words.
function hct_class_words(p::Int, q::Int)
    p > 0 || error("p must be positive")
    gcd(p, abs(q)) == 1 || error("not a primitive class")
    events = NamedTuple[]
    for k in 0:(p - 1)
        push!(events, (t=k // p, kind=:m, idx=k))
    end
    js = q > 0 ? (1:q) : q < 0 ? ((q + 1):0) : (1:0)
    for j in js
        t = (2 * p * j - 1) // (2 * p * q)
        0 <= t < 1 || error("l-crossing parameter $t")
        push!(events, (t=t, kind=:l, idx=j))
    end
    sort!(events, by=e -> e.t)
    for i in 2:length(events)
        events[i].t == events[i - 1].t && error("curve meets a grid vertex")
    end
    mlab = Dict(k => 2 + k for k in 0:(p - 1))
    lev = [e for e in events if e.kind == :l]
    llab = Dict(e.idx => p + 1 + i for (i, e) in enumerate(lev))
    bparts = String[]
    for e in events
        if e.kind == :m
            push!(bparts, "-$(mlab[e.idx])")
        else
            push!(bparts, (q > 0 ? "+" : "-") * string(llab[e.idx]))
        end
    end
    lcross = [(x=mod(p * e.t, 1), idx=e.idx) for e in lev]
    sort!(lcross, by=c -> c.x)
    lparts = ["-1"]
    for c in lcross
        push!(lparts, (q > 0 ? "-" : "+") * string(llab[c.idx]))
    end
    mcross = [(y=mod(q * (k // p) + 1 // (2 * p), 1), idx=k) for k in 0:(p - 1)]
    sort!(mcross, by=c -> c.y)
    mparts = ["+1"]
    for c in mcross
        push!(mparts, "+$(mlab[c.idx])")
    end
    (l=join(lparts, " "), m=join(mparts, " "), b=join(bparts, " "))
end

function hct_lens_hc(p::Int, q::Int)
    words = hct_class_words(p, q)
    mktempdir() do dir
        write(joinpath(dir, "diagram.lc"), "genus: 1\nl1 : $(words.l)\nm1 : $(words.m)\nbeta1 : $(words.b)\n")
        write(joinpath(dir, "diagram.hc"), "format : hc/2\nlc : diagram.lc\nstandard\nrole beta1 : beta\n")
        HCT_H.read_hc(joinpath(dir, "diagram.hc"))
    end
end

function hct_classical(p::Int, q::Int, k::Int)
    z = HCT_TE.zeta(p)
    r = invmod(mod(q, p), p)
    inv((z^k - 1) * (z^(k * r) - 1))
end

function hct_in_class(order::Int, value, expected)
    z = HCT_TE.zeta(order)
    for s in 0:(order - 1), sgn in (1, -1)
        sgn * z^s * expected == value && return true
    end
    false
end

function hct_same_values(order::Int, left, right)
    for k in 1:(order - 1)
        hct_in_class(order, left(k), right(k)) || return false
        # One global sign and one global shift, not a fresh choice at each k.
    end
    z = HCT_TE.zeta(order)
    for s in 0:(order - 1), sgn in (1, -1)
        good = true
        for k in 1:(order - 1)
            if left(k) != sgn * z^(k * s) * right(k)
                good = false
                break
            end
        end
        good && return true
    end
    false
end

@testset "Reidemeister torsion of lens spaces p ≤ 7" begin
    for p in 1:7, q in 0:(p - 1)
        gcd(p, q) == 1 || continue
        hc = hct_lens_hc(p, q)
        torsion = HCT_H.hc_reidemeister_torsion(hc)
        if p == 1
            @test torsion.order == 1
            @test torsion.acyclic == false
            @test_throws HCT_H.HCTorsionError HCT_H.hc_torsion_character(torsion, 1)
        else
            @test torsion.order == p
            @test torsion.acyclic
            for k in 1:(p - 1)
                @test HCT_H.hc_torsion_character(torsion, k) == hct_classical(p, q, k)
            end
        end
    end

    t71 = HCT_H.hc_reidemeister_torsion(hct_lens_hc(7, 1))
    t72 = HCT_H.hc_reidemeister_torsion(hct_lens_hc(7, 2))
    t76 = HCT_H.hc_reidemeister_torsion(hct_lens_hc(7, 6))
    char71(k) = HCT_H.hc_torsion_character(t71, k)
    char72(k) = HCT_H.hc_torsion_character(t72, k)
    char76(k) = HCT_H.hc_torsion_character(t76, k)
    @test !hct_same_values(7, char71, char72)
    @test hct_same_values(7, char71, char76)
    p71 = HCT_H.hc_presentation(hct_lens_hc(7, 1)).presentation.relators
    p72 = HCT_H.hc_presentation(hct_lens_hc(7, 2)).presentation.relators
    @test p71 == p72
end

@testset "saved genus-one diagrams" begin
    for name in readdir(HCT_CASES)
        m = match(r"^g1-D-(\d+)-(m)?(\d+)$", name)
        m === nothing && continue
        p = parse(Int, m.captures[1])
        q = parse(Int, m.captures[3]) * (m.captures[2] === nothing ? 1 : -1)
        p <= 7 || continue
        hc = HCT_H.read_hc(joinpath(HCT_CASES, name, "diagram.hc"))
        torsion = HCT_H.hc_reidemeister_torsion(hc)
        if p == 1
            @test !torsion.acyclic
        else
            @test torsion.order == p
            for k in 1:(p - 1)
                @test HCT_H.hc_torsion_character(torsion, k) == hct_classical(p, mod(q, p), k)
            end
        end
    end
    s1s2 = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    @test_throws HCT_H.HCTorsionError HCT_H.hc_reidemeister_torsion(s1s2)
end

@testset "stabilization of L(5, 2)" begin
    genus_one = HCT_H.read_hc(joinpath(HCT_CASES, "g1-D-5-2", "diagram.hc"))
    stabilized = HCT_H.read_hc(joinpath(HCT_CASES, "stab-D-5-2", "diagram.hc"))
    @test HCT_H.hc_presentation(stabilized).presentation.generators == ["alpha2", "m1"]
    torsion_one = HCT_H.hc_reidemeister_torsion(genus_one)
    torsion_stab = HCT_H.hc_reidemeister_torsion(stabilized)
    @test torsion_stab.invariant_factors == [1, 5]
    for k in 1:4
        classical = hct_classical(5, 2, k)
        @test HCT_H.hc_torsion_character(torsion_one, k) == classical
        @test HCT_H.hc_torsion_character(torsion_stab, k) == -classical
    end
end

@testset "connected sum is not the product of the lens torsions" begin
    diagram = HCT_H.read_hc(joinpath(HCT_CASES, "sum-D-2-1__D-3-1", "diagram.hc"))
    relators = HCT_H.hc_presentation(diagram).presentation.relators
    @test relators == [[(1, 1), (1, 1)], [(2, 1), (2, 1), (2, 1)]]
    torsion = HCT_H.hc_reidemeister_torsion(diagram)
    @test torsion.order == 6
    @test length(torsion.keys) == 5
    acyclic = [value for value in torsion.values if value !== nothing]
    # ⟨x, y | x^2, y^3⟩ has Fox matrix diag(1+x, 1+y+y^2). At the two characters
    # nontrivial on both factors this matrix is zero, so the complex is not acyclic.
    @test length(acyclic) == 3
    z = HCT_TE.zeta(6)
    product_r1 = (1 // 4) * inv((z^2 - 1)^2)
    product_r2 = (1 // 4) * inv((z^4 - 1)^2)
    @test !iszero(product_r1)
    @test !iszero(product_r2)
    for value in acyclic
        @test !hct_in_class(6, value, product_r1)
        @test !hct_in_class(6, value, product_r2)
    end
end

# Lemma 5.7 inside the cyclotomic field of order `order`, which is a multiple of p.
function hct_lens_in(order::Int, p::Int, q::Int, k::Int)
    order % p == 0 || error("order $order is not a multiple of $p")
    z = HCT_TE.zeta(order)
    root = z^(order ÷ p)
    r = invmod(mod(q, p), p)
    inv((root^k - 1) * (root^(k * r) - 1))
end

@testset "connected sum, character trivial on one factor" begin
    # Turaev, Torsions of 3-dimensional Manifolds, Birkhäuser 2002,
    # Theorem XII.1.2, p. 176, the case b₁(M₁) = b₁(M₂) = 0:
    # τ(M₁#M₂, e₁#e₂) = τ(M₁, e₁) Σ_{H₁(M₂)} + τ(M₂, e₂) Σ_{H₁(M₁)}.
    # A character trivial on H₁(M₁) and nontrivial on H₁(M₂) kills the first
    # summand and sends Σ_{H₁(M₁)} to |H₁(M₁)|. The class is taken in
    # F^× / ±φ(H₁). The lens factor is Lemma 5.7 with ε_p = +1.
    diagram = HCT_H.read_hc(joinpath(HCT_CASES, "sum-D-2-1__D-3-1", "diagram.hc"))
    torsion = HCT_H.hc_reidemeister_torsion(diagram)
    acyclic = [value for value in torsion.values if value !== nothing]
    @test length(acyclic) == 3
    on_L3_k1 = 2 * hct_lens_in(6, 3, 1, 1)
    on_L3_k2 = 2 * hct_lens_in(6, 3, 1, 2)
    on_L2 = 3 * hct_lens_in(6, 2, 1, 1)
    @test hct_in_class(6, on_L3_k1, on_L3_k2)
    @test !hct_in_class(6, on_L2, on_L3_k1)
    @test count(v -> hct_in_class(6, v, on_L3_k1), acyclic) == 2
    @test count(v -> hct_in_class(6, v, on_L2), acyclic) == 1
end

# ρ(t) = ± t^k Δ(t^ε), ε ∈ {±1}. Both sides already have valuation 0 and a
# positive leading coefficient, so the relation is equality or reversal.
function hca_equivalent(got::Vector{Int}, expected::Vector{Int})
    reversed = reverse(expected)
    if !isempty(reversed) && reversed[end] < 0
        reversed = .-reversed
    end
    got == expected || got == reversed
end

@testset "Alexander polynomial of S1xS2" begin
    # Nicolaescu, Notes on the Reidemeister Torsion, January 2002.
    # Example 2.43 (1), printed p. 72: T_{S^1×S^2} = (1−t)^{−2}, H₁ = Z.
    # Printed pp. 22–23: aug is the identity when Tors H₁ = 0.
    # Theorem 2.37, printed p. 69: b₁ = 1 gives T^aug ∼ Δ/(t−1)^2.
    # Hence Δ ∼ 1.
    diagram = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    poly = HCT_H.hc_alexander_polynomial(diagram)
    @test poly.betti == 1
    @test hca_equivalent(poly.coefficients, [1])
end

@testset "Alexander polynomial of five 0-surgeries" begin
    # Expected coefficients are the body of SnapPy 3.3.2 alexander_polynomial
    # on fundamental_group(), computed by alexander_snappy.py. The decorated
    # method raises SageNotAvailable here. Each closed filling agrees with the
    # classical Alexander polynomial of the knot: trefoil 1−t+t^2, figure-eight
    # 1−3t+t^2, cinquefoil 1−t+t^2−t^3+t^4, 5_2 is 2−3t+2t^2, 6_1 is 2−5t+2t^2.
    root = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-invariants", "alexander-b1", "hc")
    expected = [
        ("00", [1, -1, 1]),
        ("01", [1, -3, 1]),
        ("02", [1, -1, 1, -1, 1]),
        ("03", [2, -3, 2]),
        ("04", [2, -5, 2]),
    ]
    for (stem, coefficients) in expected
        diagram = HCT_H.read_hc(joinpath(root, stem * ".hc"))
        homology = HCT_H.hc_first_homology(diagram)
        @test homology.free_rank == 1
        @test homology.torsion == Int[]
        poly = HCT_H.hc_alexander_polynomial(diagram)
        @test poly.betti == 1
        @test hca_equivalent(poly.coefficients, coefficients)
    end
end

# Complex conjugation is the automorphism ζ ↦ ζ^{−1} of Q(ζ_N).
function hct_conjugate(value)
    order = Int(only(typeof(value).parameters))
    root = HCT_TE.zeta(order)
    conjugated = zero(value)
    for (index, coefficient) in enumerate(HCT_TE.coefficients(value))
        iszero(coefficient) && continue
        conjugated += coefficient * root^(-(index - 1))
    end
    conjugated
end

function hct_signs_flipped(original, mirrored)
    left = Dict(String(name) => Int[-x for x in word] for (name, word) in original.lc.words)
    right = Dict(String(name) => Int.(word) for (name, word) in mirrored.lc.words)
    left == right
end

@testset "mirror sends lens torsion to the conjugate at the inverse character" begin
    # Lemma 5.7 with ε_p = +1. The automorphism ζ ↦ ζ^{−1} satisfies
    # conjugate(τ(k)) = ζ^{k(1+q^{−1})} τ(k) and conjugate(τ(−k)) = τ(k).
    # For L(5, 2), q^{−1} ≡ 3 (mod 5), so the same-index conjugate is ζ^4 τ(k).
    root = HCT_TE.zeta(5)
    classical = hct_classical(5, 2, 1)
    @test hct_conjugate(classical) == root^4 * classical
    @test hct_conjugate(classical) != classical
    @test hct_conjugate(hct_classical(5, 2, 4)) == classical
    for p in 1:7, q in 0:(p - 1)
        gcd(p, q) == 1 || continue
        hc = hct_lens_hc(p, q)
        mirrored = HCT_H.hc_mirror(hc)
        @test hct_signs_flipped(hc, mirrored)
        torsion = HCT_H.hc_reidemeister_torsion(mirrored)
        if p == 1
            @test torsion.order == 1
            @test torsion.acyclic == false
        else
            @test torsion.order == p
            @test torsion.acyclic
            for k in 1:(p - 1)
                expected = hct_conjugate(hct_classical(p, q, p - k))
                @test expected == hct_classical(p, q, k)
                @test HCT_H.hc_torsion_character(torsion, k) == expected
            end
        end
    end
end

@testset "Alexander polynomial of a mirror" begin
    # ρ(t) = ± t^k Δ_*(t^ε), ε ∈ {±1}. These polynomials are reciprocal, so the
    # mirror is equivalent to the same normalized coefficients. The coefficients
    # are the literature values already used for the unmirrored diagrams.
    s1s2 = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    mirrored_s1s2 = HCT_H.hc_mirror(s1s2)
    @test hct_signs_flipped(s1s2, mirrored_s1s2)
    poly = HCT_H.hc_alexander_polynomial(mirrored_s1s2)
    @test poly.betti == 1
    @test hca_equivalent(poly.coefficients, [1])
    root = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-invariants", "alexander-b1", "hc")
    expected = [
        ("00", [1, -1, 1]),
        ("01", [1, -3, 1]),
        ("02", [1, -1, 1, -1, 1]),
        ("03", [2, -3, 2]),
        ("04", [2, -5, 2]),
    ]
    for (stem, coefficients) in expected
        diagram = HCT_H.read_hc(joinpath(root, stem * ".hc"))
        mirrored = HCT_H.hc_mirror(diagram)
        @test hct_signs_flipped(diagram, mirrored)
        poly = HCT_H.hc_alexander_polynomial(mirrored)
        @test poly.betti == 1
        @test hca_equivalent(poly.coefficients, coefficients)
    end
end

# Characters trivial on exactly one factor, Theorem XII.1.2, both b₁ = 0.
# |H₁| of the trivial factor times Lemma 5.7 of the other factor.
function hct_b0_expected(order::Int, p::Int, q::Int, r::Int, s::Int)
    on_left = [r * hct_lens_in(order, p, q, k) for k in 1:(p - 1)]
    on_right = [p * hct_lens_in(order, r, s, k) for k in 1:(r - 1)]
    vcat(on_left, on_right)
end

function hct_match_classes(order::Int, values, expected)
    length(values) == length(expected) || return false
    used = falses(length(values))
    for target in expected
        hit = findfirst(index -> !used[index] && hct_in_class(order, values[index], target), eachindex(values))
        hit === nothing && return false
        used[hit] = true
    end
    true
end

@testset "connected sum built by hc_connected_sum, both b1 = 0" begin
    # Turaev, Torsions of 3-dimensional Manifolds, Birkhäuser 2002,
    # Theorem XII.1.2, p. 176. Both signs in the b₁ = b₁ = 0 line are plus.
    # A character trivial on one factor specializes that line to |H₁| times
    # the torsion of the other factor, in F^× / ±φ(H₁).
    pairs = [(2, 1, 3, 1), (3, 1, 2, 1), (5, 2, 2, 1), (2, 1, 5, 2)]
    for (p, q, r, s) in pairs
        built = HCT_H.hc_connected_sum(hct_lens_hc(p, q), hct_lens_hc(r, s))
        @test built.hc.lc.genus == 2
        homology = HCT_H.hc_first_homology(built.hc)
        order = p * r
        @test homology.free_rank == 0
        @test homology.torsion == [order]
        torsion = HCT_H.hc_reidemeister_torsion(built.hc)
        @test torsion.order == order
        acyclic = [value for value in torsion.values if value !== nothing]
        expected = hct_b0_expected(order, p, q, r, s)
        @test length(acyclic) == length(expected)
        @test hct_match_classes(order, acyclic, expected)
    end
end

@testset "connected sum built by hc_connected_sum, one b1 positive" begin
    # Turaev, Theorem XII.1.2, p. 176, the mixed line: τ = τ(S¹×S²) Σ_{H₁(L(p,q))}.
    # Nicolaescu, Example 2.43 (1), printed p. 72: τ(S¹×S²) = (1−t)^{−2}.
    # Printed pp. 22–23: aug sends the torsion summand onto 1, so Σ augments to p.
    # Theorem 2.37, printed p. 69: b₁ = 1 gives T^aug ∼ Δ/(t−1)^2. Hence Δ ∼ p.
    s1s2 = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    samples = [(2, 1), (3, 1), (5, 1), (5, 2), (7, 1), (7, 2)]
    for (p, q) in samples
        lens = hct_lens_hc(p, q)
        for (left, right) in ((s1s2, lens), (lens, s1s2))
            built = HCT_H.hc_connected_sum(left, right)
            homology = HCT_H.hc_first_homology(built.hc)
            @test homology.free_rank == 1
            @test homology.torsion == [p]
            poly = HCT_H.hc_alexander_polynomial(built.hc)
            @test poly.betti == 1
            @test hca_equivalent(poly.coefficients, [p])
        end
    end
end

@testset "reverse prism diagrams have a Fox boundary" begin
    # |H₁| = 4|q| is the order computed with the prism homology, before torsion.
    # The assertion is that no nontrivial character raises bad_differential.
    for (name, order) in (("p3q1", 4), ("p3qm1", 4), ("p5q1", 4), ("p3q2", 8), ("p5q2", 8))
        hc = HCT_H.read_hc(joinpath(HCT_CASES, name, "diagram.hc"))
        torsion = HCT_H.hc_reidemeister_torsion(hc)
        @test torsion.order == order
        @test torsion.acyclic
        for (key, value) in zip(torsion.keys, torsion.values)
            if value === nothing
                @test_throws HCT_H.HCTorsionError HCT_H.hc_torsion_character(torsion, key)
            else
                @test HCT_H.hc_torsion_character(torsion, key) == value
            end
        end
    end
end

function hct_exact_multiset(left, right)
    length(left) == length(right) || return false
    used = falses(length(left))
    for value in right
        hit = findfirst(index -> !used[index] && left[index] == value, eachindex(left))
        hit === nothing && return false
        used[hit] = true
    end
    true
end

@testset "a slide keeps Reidemeister torsion" begin
    # The expected values are the torsion of the same diagram before the slide,
    # recorded from the source HC. A slide does not change the manifold, so the
    # prepared diagram is compared as a set of acyclic values in F^× / ±φ(H₁).
    z6 = HCT_TE.zeta(6)
    cases = [
        (joinpath(HCT_CASES, "sum-D-2-1__D-2-1", "diagram.hc"), 2, 1,
            Any[-1//2, -1//2]),
        (joinpath(HCT_CASES, "sum-D-2-1__D-3-1", "diagram.hc"), 6, 2,
            Any[-(2//3) * z6, -3//4, -(2//3) + (2//3) * z6]),
        (joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261112-example-table", "E1-1", "turn10-public-sum-m", "diagram.hc"), 6, 2,
            Any[-(2//3) + (2//3) * z6, 3//4, -(2//3) * z6]),
    ]
    for (path, order, gaps, expected) in cases
        source = HCT_H.read_hc(path)
        before = HCT_H.hc_reidemeister_torsion(source)
        after = HCT_H.hc_reidemeister_torsion(HCT_H.hc_prepare_joining_disk(source).hc)
        @test before.order == order
        @test after.order == order
        @test count(isnothing, before.values) == gaps
        @test count(isnothing, after.values) == gaps
        before_acyclic = [value for value in before.values if value !== nothing]
        after_acyclic = [value for value in after.values if value !== nothing]
        @test hct_exact_multiset(before_acyclic, expected)
        @test hct_match_classes(order, after_acyclic, expected)
        for (key, value) in zip(after.keys, after.values)
            if value === nothing
                @test_throws HCT_H.HCTorsionError HCT_H.hc_torsion_character(after, key)
            else
                @test HCT_H.hc_torsion_character(after, key) == value
            end
        end
    end
end

@testset "connected sum of two copies of S1xS2" begin
    # b₁ = 1 + 1. Turaev, Theorem XII.1.2, says τ = 0 when both factors have
    # positive b₁. H₁ is Z⊕Z. The torsion and one-variable Alexander functions
    # do not return that 0; this test only checks that the sum is the manifold
    # with that homology.
    s1s2 = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    built = HCT_H.hc_connected_sum(s1s2, s1s2)
    @test built.hc.lc.genus == 2
    homology = HCT_H.hc_first_homology(built.hc)
    @test homology.free_rank == 2
    @test homology.torsion == Int[]
end

# Group-ring matrix whose +1 class and its inverse both square the Fox boundary,
# but the torsion ratio is not one character of H₁. H = (Z/3)^2, β♯ = (α₁, α₂).
# Columns lie along (α₂−1, 1−α₁) and rows do too, so (φ(α)−1) is a nullvector.
# The factor (α₂−α₁²) kills the two characters where the inverse nullvector is
# not parallel to it; there the matrix is zero and the character is not acyclic.
function hct_poly_add(polys...)
    acc = Dict{Vector{Int},Int}()
    for poly in polys, (exponents, coefficient) in poly
        acc[exponents] = get(acc, exponents, 0) + coefficient
    end
    filter!(pair -> !iszero(pair.second), acc)
    acc
end

function hct_poly_scale(poly, coefficient)
    scaled = Dict{Vector{Int},Int}()
    for (exponents, value) in poly
        product = value * coefficient
        iszero(product) && continue
        scaled[exponents] = product
    end
    scaled
end

function hct_poly_mul(left, right)
    acc = Dict{Vector{Int},Int}()
    for (left_exponents, left_coefficient) in left, (right_exponents, right_coefficient) in right
        acc[left_exponents .+ right_exponents] = get(acc, left_exponents .+ right_exponents, 0) +
            left_coefficient * right_coefficient
    end
    filter!(pair -> !iszero(pair.second), acc)
    acc
end

function hct_split_fox()
    x = Dict([1, 0] => 1)
    y = Dict([0, 1] => 1)
    one_poly = Dict([0, 0] => 1)
    x2 = hct_poly_mul(x, x)
    mu = hct_poly_add(y, hct_poly_scale(x2, -1))
    w1 = hct_poly_add(y, hct_poly_scale(one_poly, -1))
    w2 = hct_poly_add(one_poly, hct_poly_scale(x, -1))
    column = hct_poly_mul(mu, hct_poly_mul(w1, w2))
    [
        hct_poly_mul(mu, hct_poly_mul(w1, w1)) column;
        column hct_poly_mul(mu, hct_poly_mul(w2, w2));
    ]
end

@testset "two 2-cell orientations outside one torsion class stop" begin
    order = 3
    factors = [3, 3]
    zeta = HCT_TE.zeta(order)
    right = [1 0; 0 1]
    dual = [[1, 0], [0, 1]]
    err = try
        HCT_H._hct_select_orientation(hct_split_fox(), right, factors, order, zeta, dual)
        nothing
    catch exc
        exc
    end
    @test err isa HCT_H.HCTorsionError
    @test err.reason == :ambiguous_orientation

    onez = one(zeta)
    keys = [[1], [2]]
    @test_throws HCT_H.HCTorsionError HCT_H._hct_require_unique_class(
        order, [3], keys, [Any[onez, onez], Any[onez, 2 * onez]])
    @test HCT_H._hct_require_unique_class(
        order, [3], keys, [Any[onez, onez], Any[-onez, -onez]]) === nothing
    @test HCT_H._hct_require_unique_class(
        order, [3], keys, [Any[onez, onez], Any[zeta, zeta^2]]) === nothing
end

import AbstractAlgebra

function hcm_univariate(ring, coefficients::Vector{Int})
    variable = AbstractAlgebra.gen(ring, 1)
    polynomial = zero(ring)
    for (index, coefficient) in enumerate(coefficients)
        polynomial += coefficient * variable^(index - 1)
    end
    polynomial
end

@testset "fraction-free determinant matches Laplace" begin
    ring, variables = AbstractAlgebra.laurent_polynomial_ring(AbstractAlgebra.ZZ, [:p, :q])
    matrix = AbstractAlgebra.zero_matrix(ring, 3, 3)
    matrix[1, 1] = variables[1]
    matrix[1, 2] = one(ring)
    matrix[2, 1] = variables[2] - 1
    matrix[2, 3] = variables[1]^(-1)
    matrix[3, 2] = variables[2]
    matrix[3, 3] = one(ring) + variables[1]
    @test HCT_H._hcm_det(matrix) == HCT_H._hcm_det_laplace(matrix)
end

@testset "GL(2,Z) search swaps the variables" begin
    ring, variables = AbstractAlgebra.laurent_polynomial_ring(AbstractAlgebra.ZZ, [:u, :v])
    polynomial = one(ring) + variables[2]^2
    matrix = HCT_H.hc_laurent_gl_matrix(polynomial, [(Int[0, 0], 1), (Int[2, 0], 1)]; bound=1)
    @test matrix !== nothing
end

@testset "positive torsion of S1xS2" begin
    # Nicolaescu, Example 2.43 (1): T = (1−t)^{−2}. Theorem 2.37: b₁ = 1 gives T ∼ Δ/(t−1)^2, so Δ ∼ 1.
    # The cyclotomic function still refuses this diagram.
    diagram = HCT_H.read_hc(joinpath(HCT_CASES, "g1-S1xS2", "diagram.hc"))
    @test_throws HCT_H.HCTorsionError HCT_H.hc_reidemeister_torsion(diagram)
    torsion = HCT_H.hc_positive_torsion(diagram)
    @test torsion.betti == 1
    @test size(torsion.free_basis) == (diagram.lc.genus, 1)
    @test HCT_H.hc_turaev_matches(torsion, one(parent(torsion.numerator)))
    @test_throws HCT_H.HCTorsionError HCT_H.hc_multivariable_alexander(diagram)
end

@testset "positive torsion of five 0-surgeries" begin
    # Same SnapPy Fox coefficients as the Alexander test above. The relation is
    # Massuyeau Theorem 4.4: τ ∼ Δ/(t−1)^2 on the recorded basis, modulo ± a monomial.
    root = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-invariants", "alexander-b1", "hc")
    expected = [
        ("00", [1, -1, 1]),
        ("01", [1, -3, 1]),
        ("02", [1, -1, 1, -1, 1]),
        ("03", [2, -3, 2]),
        ("04", [2, -5, 2]),
    ]
    for (stem, coefficients) in expected
        diagram = HCT_H.read_hc(joinpath(root, stem * ".hc"))
        alexander = HCT_H.hc_alexander_polynomial(diagram)
        @test hca_equivalent(alexander.coefficients, coefficients)
        torsion = HCT_H.hc_positive_torsion(diagram)
        @test torsion.betti == 1
        @test size(torsion.free_basis) == (diagram.lc.genus, 1)
        polynomial = hcm_univariate(parent(torsion.numerator), alexander.coefficients)
        @test HCT_H.hc_turaev_matches(torsion, polynomial)
    end
end

@testset "multivariable Alexander of a 2-component 0-surgery and of T3" begin
    # Expected terms are the gcd of the (n−1) Fox minors of SnapPy 3.3.2
    # fundamental_group(), written in positive-b1/expected.json. HC is not used
    # on that side. Whitehead link L5a1(0,1)(0,1) and the Borromean rings
    # L6a4(0,1)(0,1)(0,1). Both closed triangulations have 6 tetrahedra.
    root = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-invariants", "positive-b1", "hc")
    expected = [
        ("whitehead", 2, [(Int[0, 0], 1)]),
        ("borromean", 3, [(Int[0, 0, 0], 1)]),
    ]
    for (stem, betti, terms) in expected
        diagram = HCT_H.read_hc(joinpath(root, stem, "diagram.hc"))
        homology = HCT_H.hc_first_homology(diagram)
        @test homology.free_rank == betti
        @test homology.torsion == Int[]
        alexander = HCT_H.hc_multivariable_alexander(diagram)
        torsion = HCT_H.hc_positive_torsion(diagram)
        @test alexander.betti == betti
        @test size(alexander.free_basis) == (diagram.lc.genus, betti)
        @test alexander.free_basis == torsion.free_basis
        @test HCT_H.hc_laurent_gl_matrix(alexander.polynomial, terms; bound=1) !== nothing
        @test HCT_H.hc_turaev_matches(torsion, alexander.polynomial)
    end
end

@testset "nonconstant multivariable Alexander for b1 >= 2" begin
    # Fox gcd of SnapPy 3.3.2 fundamental_group(), written to
    # positive-b1/expected2.json before these diagrams. L11n219(0,1)(0,1)
    # is 1 + y + y^2. L11n218(0,1)(0,1) is 1 + 2y + y^2. Bound 3 was fixed
    # in that file. The coefficient patterns are not a GL(2, Z) image of
    # each other: that action keeps coefficients and moves exponents, and the
    # patterns are 1,1,1 and 1,2,1.
    root = joinpath(HCT_ROOT, "threads", "PJ260922", "data", "261114-invariants", "positive-b1", "hc")
    expected = [
        ("l11n219", 2, [(Int[0, 0], 1), (Int[0, 1], 1), (Int[0, 2], 1)]),
        ("l11n218", 2, [(Int[0, 0], 1), (Int[0, 1], 2), (Int[0, 2], 1)]),
    ]
    for (stem, betti, terms) in expected
        diagram = HCT_H.read_hc(joinpath(root, stem, "diagram.hc"))
        homology = HCT_H.hc_first_homology(diagram)
        @test homology.free_rank == betti
        @test homology.torsion == Int[]
        alexander = HCT_H.hc_multivariable_alexander(diagram)
        torsion = HCT_H.hc_positive_torsion(diagram)
        @test alexander.betti == betti
        @test size(alexander.free_basis) == (diagram.lc.genus, betti)
        @test alexander.free_basis == torsion.free_basis
        @test any(pair -> any(!iszero, pair[1]), terms)
        @test HCT_H.hc_laurent_gl_matrix(alexander.polynomial, terms; bound=3) !== nothing
        @test HCT_H.hc_turaev_matches(torsion, alexander.polynomial)
    end
end
