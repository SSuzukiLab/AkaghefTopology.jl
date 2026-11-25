"""
    LaurentPolynomial(variables::Vector{String}, terms::Dict{Tuple{Vararg{Int}},Rational{BigInt}}) -> LaurentPolynomial

Store exact rational coefficients for a Laurent polynomial in ordered variables.

# Arguments
- `variables` and `terms` supply the polynomial fields. The general constructor
  converts arguments to their declared field types; conversion errors propagate.
- The caller supplies L1 of the
  [specification](../specs/laurent-polynomials.md#numbered-conditions).
  Construction does not check L1 or normalization; malformed fields
  can invalidate later polynomial operations. Use [`laurent_variable`](@ref) or
  [`laurent_variables`](@ref) to obtain supported polynomial inputs.

# Returns
`LaurentPolynomial` with `variables::Vector{String}` holding the ordered names
and `terms::Dict{Tuple{Vararg{Int}},Rational{BigInt}}` holding exponent-indexed
coefficients, in the [Algebra API](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. For consistent fields, the result represents the supplied Laurent
   polynomial over the rational numbers, with no floating-point coefficient
   approximation. [proved/agent: by construction]
2. Construction does not modify the arguments. Containers already having the
   declared field types are retained, so subsequent mutations are shared.
   [proved/agent: by construction]

# Not guaranteed
The constructor does not remove zero coefficients or establish a canonical
representation. Equality compares the stored fields; a malformed or
unnormalized value need not behave like the corresponding normalized
polynomial. No knot-invariant interpretation follows from the type alone.

# Sources
cf:lib/topology.jl/src/laurent.jl, LaurentPolynomial and evaluate — exact coefficient storage and its polynomial evaluation semantics.

# Examples
```jldoctest
julia> t = laurent_variable(:t); p = LaurentPolynomial(t.variables,t.terms);

julia> evaluate(p^-2 + 3p - 1, [2//1])
21//4
```

See also [`LaurentRational`](@ref), [`evaluate`](@ref).

Role: construct
"""
struct LaurentPolynomial
    variables::Vector{String}
    terms::Dict{Tuple{Vararg{Int}},Rational{BigInt}}
end

"""
    LaurentRational(numerator::LaurentPolynomial, denominator::LaurentPolynomial) -> LaurentRational

Return a stored quotient of Laurent polynomials or embed one polynomial as a quotient.

# Arguments
- In the two-argument form, the arguments must be `LaurentPolynomial` values
  or convertible to that type. Failed conversions propagate. A meaningful
  quotient requires L1–L2 of the
  [specification](../specs/laurent-polynomials.md#numbered-conditions);
  this form does not check those conditions or cancel factors.
- The one-argument form accepts a well-formed `LaurentPolynomial` and embeds
  it with denominator one. Malformed polynomial fields are not validated.

# Returns
`LaurentRational` with `numerator` and `denominator` fields, both
`LaurentPolynomial`. This representation belongs to the
[Laurent-polynomial API](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. The two-argument form retains the supplied polynomial values, including
   their mutable coefficient and variable containers. [proved/agent: by construction]
2. The one-argument form represents the same rational function as its input
   polynomial and has denominator equal to `one(polynomial)`. [proved/agent: by construction]
3. For matching variables and a nonzero denominator, the stored pair
   represents the quotient of those Laurent polynomials over the rationals. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Not guaranteed
Direct two-argument construction does not establish L2 or the
[fraction normalization](../specs/laurent-polynomials.md#fraction-data-and-equality).
No general rational-function identity test is provided. No manifold or knot
interpretation follows from this type.

# Sources
cf:lib/topology.jl/src/laurent.jl, LaurentRational — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> x, y = laurent_variables(2; prefix="x");

julia> ((x+y)*(x-y))/(x-y) == LaurentRational(x+y)
true
```

See also [`evaluate`](@ref), [`limit_zero`](@ref).

Role: convert
"""
struct LaurentRational
    numerator::LaurentPolynomial
    denominator::LaurentPolynomial
end

function _laurent(variables,terms)
    filtered=Dict{Tuple{Vararg{Int}},Rational{BigInt}}()
    for (exponents,coefficient) in terms
        coefficient==0 && continue
        key=Tuple(Int.(exponents))
        filtered[key]=get(filtered,key,0//1)+coefficient
    end
    filter!(pair -> last(pair)!=0,filtered)
    LaurentPolynomial(String.(variables),filtered)
end

"""
    laurent_variable(name::Union{Symbol,String}) -> LaurentPolynomial

Return one indeterminate in an exact univariate Laurent-polynomial ring.

# Arguments
`name` is a `Symbol` or `String`. It is used as the variable name without
identifier or nonemptiness validation.

# Returns
A `LaurentPolynomial` with one variable and exact rational coefficients.
The polynomial API belongs to
[Algebra](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. The result represents the indeterminate named by `String(name)`. [proved/agent: by construction]
2. The coefficient of that indeterminate is one and all other coefficients
   are zero. Each call supplies new variable and coefficient containers. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Not guaranteed
The name is not checked against another polynomial's variable order. There
is no automatic identification of separately named rings or interpretation
as a knot invariant.

# Sources
cf:lib/topology.jl/src/laurent.jl, laurent_variable — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> t = laurent_variable(:t);

julia> evaluate(t^-2 + 3t - 1, [2//1])
21//4
```

See also [`laurent_variables`](@ref), [`evaluate`](@ref).

Role: construct
"""
laurent_variable(name::Union{Symbol,String})=_laurent([String(name)],Dict((1,)=>1//1))

"""
    laurent_variables(count::Integer; prefix="t") -> Vector{LaurentPolynomial}

Return the ordered indeterminates of one exact multivariate Laurent-polynomial ring.

# Arguments
`count` is the requested number of variables. Nonpositive counts return an
empty vector. `prefix` supplies the textual prefix of their names; no name
validation is performed. The count must fit available allocation resources.

# Returns
A vector of `LaurentPolynomial` values, one per variable, all with the same
ordered variable names and exact rational coefficients. See the
[Algebra scope](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. For positive `count`, the result has that many elements. Element `i`
   represents the `i`th indeterminate in the common ring. [proved/agent: by construction]
2. Its variable names are `string(prefix, i)` in increasing `i` order. [proved/agent: by construction]
3. Each polynomial has independently allocated variable and coefficient
   containers. No global variable bindings are created. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Not guaranteed
The result supplies no canonical renaming between different variable lists.
It does not certify a topological interpretation of any polynomial formed
from these variables.

# Sources
cf:lib/topology.jl/src/laurent.jl, laurent_variables — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> x, y = laurent_variables(2; prefix="x");

julia> x/(x+y) + y/(x+y) == LaurentRational(one(x))
true
```

See also [`laurent_variable`](@ref), [`LaurentRational`](@ref).

Role: construct
"""
laurent_variables(count::Integer;prefix="t")=[_laurent(["$prefix$index" for index in 1:count],
    Dict(Tuple(index==position ? 1 : 0 for index in 1:count)=>1//1)) for position in 1:count]

function _constant(variables,value)
    value==0 && return _laurent(variables,Dict{Tuple{Vararg{Int}},Rational{BigInt}}())
    _laurent(variables,Dict(Tuple(zeros(Int,length(variables)))=>Rational{BigInt}(value)))
end

Base.zero(p::LaurentPolynomial)=_constant(p.variables,0)
Base.one(p::LaurentPolynomial)=_constant(p.variables,1)
Base.iszero(p::LaurentPolynomial)=isempty(p.terms)
Base.:(==)(left::LaurentPolynomial,right::LaurentPolynomial)=left.variables==right.variables && left.terms==right.terms
Base.:(==)(left::LaurentRational,right::LaurentRational)=left.numerator==right.numerator && left.denominator==right.denominator
Base.show(io::IO,p::LaurentPolynomial)=print(io,"LaurentPolynomial(",p.variables,", ",p.terms,")")
Base.show(io::IO,p::LaurentRational)=print(io,"(",p.numerator,")/(",p.denominator,")")

function _same_variables(left::LaurentPolynomial,right::LaurentPolynomial)
    left.variables==right.variables || throw(ArgumentError("Laurent variables differ"))
end

function Base.:+(left::LaurentPolynomial,right::LaurentPolynomial)
    _same_variables(left,right); terms=copy(left.terms)
    for (exponents,coefficient) in right.terms
        terms[exponents]=get(terms,exponents,0//1)+coefficient
    end
    _laurent(left.variables,terms)
end
Base.:-(left::LaurentPolynomial)=_laurent(left.variables,Dict(exponents=>-coefficient for (exponents,coefficient) in left.terms))
Base.:-(left::LaurentPolynomial,right::LaurentPolynomial)=left+(-right)
function Base.:*(left::LaurentPolynomial,right::LaurentPolynomial)
    _same_variables(left,right); terms=Dict{Tuple{Vararg{Int}},Rational{BigInt}}()
    for (a,ac) in left.terms,(b,bc) in right.terms
        exponents=Tuple(a[index]+b[index] for index in eachindex(a))
        terms[exponents]=get(terms,exponents,0//1)+ac*bc
    end
    _laurent(left.variables,terms)
end
function Base.:^(polynomial::LaurentPolynomial,power::Integer)
    power>=0 || return LaurentRational(one(polynomial),polynomial^(-power))
    result=one(polynomial); base=polynomial
    while power>0
        isodd(power) && (result=result*base)
        power>>=1; power>0 && (base=base*base)
    end
    result
end

for operator in (:+,:-)
    @eval Base.$operator(left::LaurentPolynomial,right::Number)=$operator(left,_constant(left.variables,right))
    @eval Base.$operator(left::Number,right::LaurentPolynomial)=$operator(_constant(right.variables,left),right)
end
Base.:*(scalar::Number,polynomial::LaurentPolynomial)=_laurent(polynomial.variables,
    Dict(exponents=>Rational{BigInt}(scalar)*coefficient for (exponents,coefficient) in polynomial.terms))
Base.:*(polynomial::LaurentPolynomial,scalar::Number)=scalar*polynomial

function _shift_univariate(p::LaurentPolynomial,shift::Integer)
    _laurent(p.variables,Dict((exponents[1]+shift,)=>coefficient for (exponents,coefficient) in p.terms))
end

_degree_univariate(p::LaurentPolynomial)=isempty(p.terms) ? -1 : maximum(exponents[1] for exponents in keys(p.terms))
_minimum_univariate(p::LaurentPolynomial)=isempty(p.terms) ? 0 : minimum(exponents[1] for exponents in keys(p.terms))
_leading_univariate(p::LaurentPolynomial)=p.terms[(_degree_univariate(p),)]

function _monic_univariate(p::LaurentPolynomial)
    iszero(p) && return p
    (1/_leading_univariate(p))*p
end

function _divrem_univariate(dividend::LaurentPolynomial,divisor::LaurentPolynomial)
    iszero(divisor) && throw(DivideError())
    quotient=_constant(dividend.variables,0)
    remainder=dividend
    divisor_degree=_degree_univariate(divisor)
    divisor_leading=_leading_univariate(divisor)
    while !iszero(remainder) && _degree_univariate(remainder)>=divisor_degree
        degree=_degree_univariate(remainder)-divisor_degree
        term=_laurent(dividend.variables,Dict((degree,)=>_leading_univariate(remainder)/divisor_leading))
        quotient+=term
        remainder-=term*divisor
    end
    quotient,remainder
end

function _gcd_univariate(left::LaurentPolynomial,right::LaurentPolynomial)
    a,b=_monic_univariate(left),_monic_univariate(right)
    while !iszero(b)
        _,remainder=_divrem_univariate(a,b)
        a,b=b,_monic_univariate(remainder)
    end
    _monic_univariate(a)
end

"""Cancel exact common univariate polynomial factors without a CAS."""
function _cancel_univariate(numerator::LaurentPolynomial,denominator::LaurentPolynomial)
    isempty(numerator.terms) && return numerator,one(denominator)
    numerator_shift=_minimum_univariate(numerator)
    denominator_shift=_minimum_univariate(denominator)
    normalized_numerator=_shift_univariate(numerator,-numerator_shift)
    normalized_denominator=_shift_univariate(denominator,-denominator_shift)
    common=_gcd_univariate(normalized_numerator,normalized_denominator)
    reduced_numerator,remainder=_divrem_univariate(normalized_numerator,common)
    iszero(remainder) || error("internal Laurent GCD division failed")
    reduced_denominator,remainder=_divrem_univariate(normalized_denominator,common)
    iszero(remainder) || error("internal Laurent GCD division failed")
    # Fix the remaining scalar ambiguity so equivalent fractions have one
    # deterministic representation (denominator leading coefficient = 1).
    scale=_leading_univariate(reduced_denominator)
    reduced_numerator=(1/scale)*reduced_numerator
    reduced_denominator=(1/scale)*reduced_denominator
    _shift_univariate(reduced_numerator,numerator_shift-denominator_shift),reduced_denominator
end

_leading_multivariate(p::LaurentPolynomial)=maximum(keys(p.terms))
_minimum_multivariate(p::LaurentPolynomial)=Tuple(minimum(exponents[index] for exponents in keys(p.terms))
    for index in eachindex(p.variables))
function _shift_multivariate(p::LaurentPolynomial,shift)
    _laurent(p.variables,Dict(Tuple(exponents[index]+shift[index] for index in eachindex(exponents))=>coefficient
        for (exponents,coefficient) in p.terms))
end

"""Exact multivariate polynomial division, after removing Laurent monomials.

This is deliberately a *division* primitive rather than a full multivariate
GCD implementation.  It catches the dominant cancellation pattern created by
fraction arithmetic (a whole expanded denominator remaining in the numerator)
without adding a CAS dependency or guessing a monomial ordering.
"""
function _divrem_multivariate(dividend::LaurentPolynomial,divisor::LaurentPolynomial)
    iszero(divisor) && throw(DivideError())
    _same_variables(dividend,divisor)
    dividend_minimum=_minimum_multivariate(dividend)
    divisor_minimum=_minimum_multivariate(divisor)
    normalized_dividend=_shift_multivariate(dividend,Tuple(-value for value in dividend_minimum))
    normalized_divisor=_shift_multivariate(divisor,Tuple(-value for value in divisor_minimum))
    quotient=_constant(dividend.variables,0); remainder=normalized_dividend
    divisor_leading=_leading_multivariate(normalized_divisor)
    divisor_coefficient=normalized_divisor.terms[divisor_leading]
    while !iszero(remainder)
        leading=_leading_multivariate(remainder)
        all(leading[index]>=divisor_leading[index] for index in eachindex(leading)) || break
        exponent=Tuple(leading[index]-divisor_leading[index] for index in eachindex(leading))
        term=_laurent(dividend.variables,Dict(exponent=>remainder.terms[leading]/divisor_coefficient))
        quotient+=term; remainder-=term*normalized_divisor
    end
    # Restore the Laurent monomial removed independently from dividend/divisor.
    _shift_multivariate(quotient,Tuple(dividend_minimum[index]-divisor_minimum[index]
        for index in eachindex(dividend_minimum))),remainder
end

function _cancel_multivariate_exact(numerator::LaurentPolynomial,denominator::LaurentPolynomial)
    isempty(numerator.terms) && return numerator,one(denominator)
    quotient,remainder=_divrem_multivariate(numerator,denominator)
    iszero(remainder) && return quotient,one(denominator)
    numerator,denominator
end

"""Remove shared unit-coefficient binomial factors from a multivariate fraction.

`simplify2!` creates these factors when it eliminates a crossing.  They are
usually not the entire denominator, so whole-denominator division alone leaves
them to grow at every subsequent elimination.  We test each Laurent binomial
`1 - x^a` with exponents in `{-1,0,1}` exactly; this is a finite, deterministic
subset tailored to the Live Script's Alexander matrices, not a heuristic CAS
factorization.
"""
function _cancel_multivariate_binomials(numerator::LaurentPolynomial,denominator::LaurentPolynomial)
    # Exact division costs grow with both expanded supports.  Once an
    # intermediate has escaped the small Alexander-minor regime, defer to the
    # caller rather than turning every ring operation into 13 trial divisions.
    # This is a performance guard, not an approximation.
    (length(numerator.terms)>32 || length(denominator.terms)>32) && return numerator,denominator
    count=length(numerator.variables)
    deltas=Tuple{Vararg{Int}}[]
    for delta in Iterators.product(((-1):1 for _ in 1:count)...)
        all(==(0),delta) && continue
        first_nonzero=first(value for value in delta if value!=0)
        first_nonzero>0 || continue  # x^a-1 and x^-a-1 differ by a monomial.
        push!(deltas,Tuple(delta))
    end
    changed=true
    while changed
        changed=false
        for delta in deltas
            factor=_laurent(numerator.variables,Dict(Tuple(zeros(Int,count))=>1//1,delta=>-1//1))
            nq,nr=_divrem_multivariate(numerator,factor)
            iszero(nr) || continue
            dq,dr=_divrem_multivariate(denominator,factor)
            iszero(dr) || continue
            numerator,denominator=nq,dq
            changed=true
            break
        end
    end
    numerator,denominator
end

function _monic_multivariate_fraction(numerator::LaurentPolynomial,denominator::LaurentPolynomial)
    scale=denominator.terms[_leading_multivariate(denominator)]
    (1/scale)*numerator,(1/scale)*denominator
end

function _rational(numerator::LaurentPolynomial,denominator::LaurentPolynomial)
    _same_variables(numerator,denominator); iszero(denominator) && throw(DivideError())
    iszero(numerator) && return LaurentRational(numerator,one(denominator))
    if length(numerator.variables)==1
        numerator,denominator=_cancel_univariate(numerator,denominator)
    else
        numerator,denominator=_cancel_multivariate_exact(numerator,denominator)
        numerator,denominator=_cancel_multivariate_binomials(numerator,denominator)
        numerator,denominator=_cancel_multivariate_exact(numerator,denominator)
        numerator,denominator=_monic_multivariate_fraction(numerator,denominator)
    end
    LaurentRational(numerator,denominator)
end
"""
    LaurentRational(polynomial::LaurentPolynomial) -> LaurentRational

Embed a Laurent polynomial in the stored rational-function representation.

# Arguments
`polynomial` must have well-formed variable and exponent data. These fields
are not validated; conversion or arithmetic errors from malformed data
propagate. The representation belongs to
[Algebra](../specs/laurent-polynomials.md#definitions).

# Returns
`LaurentRational` with Laurent-polynomial `numerator` and `denominator` fields.

# Guarantees
1. The result represents the same rational function as the input polynomial. [proved/agent: by construction]
2. Its denominator equals `one(polynomial)`; its ordered variables agree
   with those of the input. The input is unchanged. [proved/agent: by construction]

# Not guaranteed
The result is not a knot or manifold invariant by virtue of its type.
Construction does not certify the well-formedness of directly assembled
polynomial fields or establish a general rational-function equality test.

# Sources
cf:lib/topology.jl/src/laurent.jl, LaurentRational — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> x, y = laurent_variables(2; prefix="x");

julia> ((x+y)*(x-y))/(x-y) == LaurentRational(x+y)
true
```

See also [`evaluate`](@ref).

Role: convert
"""
LaurentRational(polynomial::LaurentPolynomial)=_rational(polynomial,one(polynomial))
Base.:(==)(left::LaurentRational,right::LaurentPolynomial)=left==LaurentRational(right)
Base.:(==)(left::LaurentPolynomial,right::LaurentRational)=LaurentRational(left)==right
Base.zero(value::LaurentRational)=LaurentRational(zero(value.numerator))
Base.one(value::LaurentRational)=LaurentRational(one(value.numerator))
Base.iszero(value::LaurentRational)=iszero(value.numerator)
Base.:-(value::LaurentRational)=_rational(-value.numerator,value.denominator)
function Base.:+(left::LaurentRational,right::LaurentRational)
    left.denominator==right.denominator && return _rational(left.numerator+right.numerator,left.denominator)
    _rational(left.numerator*right.denominator+right.numerator*left.denominator,
              left.denominator*right.denominator)
end
Base.:-(left::LaurentRational,right::LaurentRational)=left+(-right)
Base.:*(left::LaurentRational,right::LaurentRational)=_rational(left.numerator*right.numerator,left.denominator*right.denominator)
Base.:/(left::LaurentRational,right::LaurentRational)=_rational(left.numerator*right.denominator,left.denominator*right.numerator)
Base.inv(value::LaurentRational)=_rational(value.denominator,value.numerator)
function Base.:^(value::LaurentRational,power::Integer)
    power>=0 ? _rational(value.numerator^power,value.denominator^power) :
               _rational(value.denominator^(-power),value.numerator^(-power))
end
Base.:/(left::LaurentPolynomial,right::LaurentPolynomial)=_rational(left,right)
Base.:/(left::LaurentPolynomial,right::Number)=_rational(left,_constant(left.variables,right))
Base.:/(left::Number,right::LaurentPolynomial)=_rational(_constant(right.variables,left),right)
Base.inv(value::LaurentPolynomial)=_rational(one(value),value)

for operator in (:+,:-,:*)
    @eval Base.$operator(left::LaurentRational,right::LaurentPolynomial)=$operator(left,LaurentRational(right))
    @eval Base.$operator(left::LaurentPolynomial,right::LaurentRational)=$operator(LaurentRational(left),right)
    @eval Base.$operator(left::LaurentRational,right::Number)=$operator(left,LaurentRational(_constant(left.numerator.variables,right)))
    @eval Base.$operator(left::Number,right::LaurentRational)=$operator(LaurentRational(_constant(right.numerator.variables,left)),right)
end
Base.:/(left::LaurentRational,right::LaurentPolynomial)=left/LaurentRational(right)
Base.:/(left::LaurentPolynomial,right::LaurentRational)=LaurentRational(left)/right
Base.:/(left::LaurentRational,right::Number)=left/LaurentRational(_constant(left.numerator.variables,right))
Base.:/(left::Number,right::LaurentRational)=LaurentRational(_constant(right.numerator.variables,left))/right

"""
    evaluate(polynomial::LaurentPolynomial, values::AbstractVector{<:Number}) -> Rational{BigInt}

Return the exact rational value of a Laurent polynomial or stored quotient at supplied coordinates.

# Arguments
- Supply one numeric coordinate per ordered variable; a length mismatch
  throws `DimensionMismatch` with the message `one value per Laurent variable`.
- Each coordinate must convert to `Rational{BigInt}`; conversion failures
  propagate. The stored polynomials must have well-formed exponent tuples.
- For a finite mathematical value, all negative-power coordinates and the
  evaluated denominator must be nonzero. Singular evaluations are not
  regularized: rational infinities or arithmetic exceptions can result.

# Returns
`Rational{BigInt}` in the variable order stored by the polynomial. For a
quotient, this is the quotient of the numerator and denominator evaluations.
See the [Laurent API](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. For valid nonsingular data, the result is substitution into the represented
   Laurent polynomial or rational function, with exact rational arithmetic. [proved/agent: by construction]
2. Machine-integer coordinates do not restrict intermediate integer sizes. [proved/agent: by construction]
3. The input polynomial, quotient and coordinate vector are unchanged. [proved/agent: by construction]

# Not guaranteed
This is not a limit at a singularity and does not cancel removable poles in
a directly constructed quotient. Floating-point input is converted to a
rational value rather than interpreted as a measurement with an error bound.

# Sources
cf:lib/topology.jl/src/laurent.jl, evaluate — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> t = laurent_variable(:t);

julia> evaluate(t^-2 + 3t - 1, [2//1])
21//4
```

See also [`limit_zero`](@ref), [`LaurentRational`](@ref).

Role: invariant
"""
function evaluate(polynomial::LaurentPolynomial,values::AbstractVector{<:Number})
    length(values)==length(polynomial.variables) || throw(DimensionMismatch("one value per Laurent variable"))
    # Do not let a caller-provided `Int` / `Rational{Int}` cap the exact
    # arithmetic: the intermediate Laurent powers can exceed machine words.
    promoted=Rational{BigInt}.(values)
    sum(coefficient*prod(promoted[index]^exponents[index] for index in eachindex(promoted))
        for (exponents,coefficient) in polynomial.terms;init=BigInt(0)//BigInt(1))
end
"""
    evaluate(value::LaurentRational, values::AbstractVector{<:Number}) -> Rational{BigInt}

Return the exact rational evaluation of a stored Laurent quotient.

# Arguments
`values` must satisfy the coordinate and nonsingularity requirements of
[`evaluate(::LaurentPolynomial, ::AbstractVector)`](@ref evaluate).
Each polynomial's variable count is checked independently; a mismatch
throws `DimensionMismatch` naming the required coordinate count.
Conversion and rational arithmetic exceptions propagate. For a finite value,
the denominator evaluation must be nonzero.

# Returns
`Rational{BigInt}`, the ratio of the two polynomial evaluations in their
stored variable orders. This belongs to the
[Laurent API](../specs/laurent-polynomials.md#definitions).

# Guarantees
1. On a well-formed quotient at a nonsingular point, the result is the
   represented rational function's value with exact rational arithmetic. [proved/agent: by construction]
2. Neither the input quotient nor the coordinate vector is changed. [proved/agent: by construction]

# Not guaranteed
The variable lists are not compared here. A directly constructed quotient
need not be well formed. Singular points are not treated as limits;
rational infinities or arithmetic exceptions can result, including for
a removable singularity in an uncancelled pair.

# Sources
cf:lib/topology.jl/src/laurent.jl, evaluate — exact rational polynomial storage and evaluation semantics.

# Examples
```jldoctest
julia> t = laurent_variable(:t);

julia> evaluate(t^-2 + 3t - 1, [2//1])
21//4
```

See also [`limit_zero`](@ref), [`LaurentRational`](@ref).

Role: invariant
"""
evaluate(value::LaurentRational,values)=evaluate(value.numerator,values)/evaluate(value.denominator,values)

"""
    limit_zero(polynomial::LaurentPolynomial) -> Rational

Return the exact finite limit at zero of a univariate Laurent polynomial.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/laurent-polynomials.md).

`polynomial` must be a well-formed normalized Laurent polynomial in one
variable. Variable count is checked: `ArgumentError` reports the one-variable
requirement. Other field consistency and normalization are the caller's
responsibility; malformed data can cause indexing errors or an incorrect
limit. A pole raises `DomainError` reporting divergence at zero.

# Returns
The finite rational limit. A zero result can have machine-integer numerator
and denominator; a stored constant coefficient has `BigInt` numerator and
denominator.

# Guarantees
1. For valid univariate input with a finite limit, the result is its limit
   at zero with exact rational arithmetic. [proved/agent: by construction]
2. A genuine pole is reported by an exception, not an infinite return value.
   [proved/agent: by construction]
3. The polynomial and its containers are not modified.
   [proved/agent: by construction]

# Not guaranteed
This is not multivariate, path-dependent or symbolic analytic limit
evaluation. Direct construction of unnormalized polynomial fields is not
validated by this call.

# Sources
cf:lib/topology.jl/src/laurent.jl, limit_zero — exact univariate Laurent evaluation at the zero limit.

# Examples
```jldoctest
julia> t = laurent_variable(:t);

julia> limit_zero(t^2+3t) == 0
true
```

See also [`evaluate`](@ref), [`LaurentPolynomial`](@ref).

Role: invariant
"""
function limit_zero(polynomial::LaurentPolynomial)
    length(polynomial.variables)==1 || throw(ArgumentError("limit_zero requires one Laurent variable"))
    iszero(polynomial) && return 0//1
    valuation=_minimum_univariate(polynomial)
    valuation<0 && throw(DomainError(polynomial,"Laurent polynomial diverges at zero"))
    valuation==0 ? polynomial.terms[(0,)] : 0//1
end
"""
    limit_zero(value::LaurentRational) -> Rational

Return the exact finite limit at zero of a univariate Laurent rational function.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/laurent-polynomials.md).

`value` must have consistent normalized univariate polynomial fields with
matching variables and nonzero denominator. Only the numerator's variable
count is checked explicitly: `ArgumentError` reports the one-variable
requirement. Remaining field validity is the caller's responsibility;
malformed fields can cause exceptions or a result without the claimed
interpretation. A pole raises `DomainError` reporting divergence at zero.

# Returns
The finite rational limit, using `Rational{Int}` for immediate zero returns
and `Rational{BigInt}` for a ratio of stored coefficients.

# Guarantees
1. For valid input with a finite limit, the result is the rational function's
   exact limit at zero, including removable singularities.
   [proved/agent: by construction]
2. A genuine pole is reported by an exception. The input quotient and its
   polynomial containers are not modified. [proved/agent: by construction]

# Not guaranteed
This does not validate a directly constructed quotient, compare variable
lists, or handle multivariate limits. In particular, a zero numerator does
not trigger validation of a directly supplied zero denominator.

# Sources
cf:lib/topology.jl/src/laurent.jl, limit_zero — exact univariate quotient limit, including removable singularities.

# Examples
```jldoctest
julia> t = laurent_variable(:t);

julia> (limit_zero((t^2+3t)/(one(t)-t)) == 0, limit_zero((t+t^2)/(t+t^2)) == 1)
(true, true)
```

See also [`evaluate`](@ref), [`LaurentRational`](@ref).

Role: invariant
"""
function limit_zero(value::LaurentRational)
    numerator,denominator=value.numerator,value.denominator
    length(numerator.variables)==1 || throw(ArgumentError("limit_zero requires one Laurent variable"))
    iszero(numerator) && return 0//1
    valuation=_minimum_univariate(numerator)-_minimum_univariate(denominator)
    valuation<0 && throw(DomainError(value,"Laurent rational diverges at zero"))
    valuation>0 && return 0//1
    numerator.terms[(_minimum_univariate(numerator),)] / denominator.terms[(_minimum_univariate(denominator),)]
end
