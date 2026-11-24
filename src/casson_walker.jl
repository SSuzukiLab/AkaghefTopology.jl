# Casson–Walker λ from surgery data. These functions do not read a Heegaard code.
# Normalization: [[Les07#^anonnd70-blk-0138]]. Draft: docs/conventions/CassonWalker.md.

"""
    sawtooth(x)

((x)) from [[Sav12#^sav12-sec-0540]]: x − [x] − 1/2 off the integers, and 0 on the integers.
"""
function sawtooth(x::Rational)
    isinteger(x) && return zero(x)
    return x - floor(x) - 1//2
end

"""
    dedekind_sum(b, a)

Dedekind sum s(b, a) for a > 0. The sum from k = 1 to a − 1 is empty when a = 1.
"""
function dedekind_sum(b::Integer, a::Integer)
    a > 0 || throw(ArgumentError("the second argument of a Dedekind sum is positive"))
    bb = mod(b, a)
    s = 0//1
    for k in 1:(a - 1)
        s += sawtooth(k//a) * sawtooth((k * bb)//a)
    end
    return s
end

"""
    casson_walker_lens(p, q)

λ(L(p, q)) = −s(q, p)/2 for p > 0 and gcd(p, q) = 1.
L(p, q) is the p/q-surgery on the unknot in S³.
"""
function casson_walker_lens(p::Integer, q::Integer)
    p > 0 || throw(ArgumentError("lens order p is positive in this naming"))
    gcd(p, abs(q)) == 1 || throw(ArgumentError("p and q are coprime"))
    return -dedekind_sum(q, p) // 2
end

"""
    casson_walker_knot_surgery(p, q, delta2)

λ of p/q-surgery on a knot in S³ with Alexander second derivative `delta2` at 1.
The knot term is (q/p)·Δ''(1)/2 and the lens term is λ(L(p, q)) for p > 0.
A negative p is the orientation reverse of the positive surgery with the same |slope|.
"""
casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Integer) =
    casson_walker_knot_surgery(p, q, delta2//1)

function casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Rational)
    q == 0 && throw(ArgumentError("q = 0 is not a finite surgery coefficient"))
    if q < 0
        p = -p
        q = -q
    end
    gcd(abs(p), q) == 1 || throw(ArgumentError("p and q are coprime"))
    knot = (q//p) * (delta2 // 2)
    lens = p > 0 ? casson_walker_lens(p, q) : -casson_walker_lens(-p, q)
    return knot + lens
end

"""
    casson_walker_surgery_pup(p, q, delta2)

The q > 0 expansion in the proof of [[Les07#^anonnd70-blk-0314]]:
(q/p)(Δ''(1)/2 − 1/24 − (p²+1)/(24q²)) + sign(pq)/8 + s(p, q)/2.
For a knot in S³ the self-linking is 0.
"""
casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Integer) =
    casson_walker_surgery_pup(p, q, delta2//1)

function casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Rational)
    q > 0 || throw(ArgumentError("this expansion is stated for q > 0"))
    gcd(abs(p), q) == 1 || throw(ArgumentError("p and q are coprime"))
    sign = p > 0 ? 1 : -1
    return (q//p) * (delta2 // 2 - 1//24 - (p * p + 1)//(24 * q * q)) +
           sign//8 + dedekind_sum(p, q) // 2
end

"""
    casson_walker_brieskorn(p, q, r)

λ(Σ(p, q, r)) from [[Sav12#^sav12-sec-0544]], with d(a, b) = −4 s(b, a).
"""
function casson_walker_brieskorn(p::Integer, q::Integer, r::Integer)
    (p > 0 && q > 0 && r > 0) || throw(ArgumentError("Seifert invariants are positive"))
    term = (1 - (p * q * r)^2 + (p * q)^2 + (q * r)^2 + (p * r)^2)//(3 * p * q * r)
    d(a, b) = -4 * dedekind_sum(b, a)
    bracket = 1 - term - d(p, q * r) - d(q, p * r) - d(r, p * q)
    return -bracket // 8
end

"""
    casson_walker_seifert(b, fibers)

Closed form for a Seifert manifold with integer Euler contribution `b` and
exceptional fibers `(α, β)`. Returns `(lambda, h1, E)`.
`E = −b − Σ β/α`. The value is the expression that matches the circle bundle
L(p, 1) and Sav12's λ(Σ(2, 3, 5)) and λ(Σ(2, 3, 7)). Returns `nothing` when E = 0.
"""
function casson_walker_seifert(b::Integer, fibers::Vector{Tuple{Int,Int}})
    E = -b//1
    ssum = 0//1
    quad = (2 - length(fibers))//1
    order = 1
    for (alpha, beta) in fibers
        alpha == 0 && throw(ArgumentError("a fiber multiplicity is nonzero"))
        gcd(abs(alpha), abs(beta)) == 1 || throw(ArgumentError("a fiber pair is coprime"))
        E -= beta//alpha
        quad += 1//(alpha * alpha)
        order *= abs(alpha)
        s = dedekind_sum(beta, abs(alpha))
        ssum += 12 * (alpha > 0 ? 1 : -1) * s
    end
    E == 0 && return nothing
    inside = 3 * (E > 0 ? 1 : -1) - E - ssum - quad / E
    return (lambda = inside / 24, h1 = abs(E) * order, E = E)
end
