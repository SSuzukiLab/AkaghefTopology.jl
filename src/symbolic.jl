"""
    SymExpr(text::String) -> SymExpr

Store formal expression text for the lightweight symbolic operations.

# Arguments
`text` is a string, or a value convertible to `String` through the general
constructor. Conversion errors propagate. Expression syntax and mathematical
meaning are the caller's responsibility and are not checked.

# Returns
`SymExpr` with the supplied `text::String`, in the
[Algebra symbolic API](../specs/algebra.md#definitions).

# Guarantees
1. The text is preserved without evaluation or simplification. The value is a
   formal textual expression, not a certified element of a specified symbolic
   ring. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Arbitrary text need not denote a well-formed mathematical expression. Equality
compares text and does not decide mathematical identity; arithmetic provides
only the simplifications supported by this API.

# Sources
cf:lib/topology.jl/src/symbolic.jl, SymExpr and its equality operation — the textual expression model.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [SymExpr("x"),SymExpr("y"),SymExpr("z")]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

See also [`sym`](@ref), [`LaurentPolynomial`](@ref).

Role: construct
"""
struct SymExpr
    text::String
end

"""
    sym(x::SymExpr) -> SymExpr

Return a textual symbolic expression for the supplied value.

# Arguments
`x` may be a `SymExpr`, symbol, string or number. Strings are accepted
verbatim without syntax validation. Conversion to text can propagate the
input type's string-conversion errors.

# Returns
`SymExpr`, whose `text::String` is the supplied expression text or the
textual representation of the supplied symbol or number. This lightweight
symbolic API belongs to [Algebra](../specs/algebra.md#definitions).

# Guarantees
1. A `SymExpr` input is returned unchanged. String and symbol inputs preserve
   their text; a number uses its Julia string representation. [proved/agent: by construction]
2. The result represents formal expression text for subsequent symbolic
   operations; construction does not evaluate the text. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
This is not a parser, a computer algebra system or a mathematical identity
test. Arbitrary text need not denote a well-formed expression. Equality of
these values compares text, so equivalent formulas can compare unequal.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

See also [`SymExpr`](@ref), [`laurent_variable`](@ref).

Role: convert
"""
sym(x::SymExpr) = x
"""
    sym(x::Symbol) -> SymExpr

Return formal expression text for a symbol.

# Arguments
`x` is a Julia symbol; no mathematical identifier validation is performed.

# Returns
`SymExpr` with `text == String(x)`, in the
[Algebra symbolic API](../specs/algebra.md#definitions).

# Guarantees
1. The symbol's text is preserved, with no evaluation or simplification. [proved/agent: by construction]
2. The result can be used as a formal symbol by the `SymExpr` operations. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The text does not certify a domain, variable binding or mathematical identity.
Equivalent formulas need not have equal text.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

See also [`SymExpr`](@ref).

Role: construct
"""
sym(x::Symbol) = SymExpr(String(x))

"""
    sym(x::AbstractString) -> SymExpr

Return formal expression text from a string.

# Arguments
`x` is an `AbstractString` convertible to `String`. Conversion failures
propagate. Arbitrary text is accepted without mathematical syntax checks.

# Returns
`SymExpr` with `text == String(x)`, in the
[Algebra symbolic API](../specs/algebra.md#definitions).

# Guarantees
1. The result retains the supplied expression text without evaluation. [proved/agent: by construction]
2. Constructing it does not create variable bindings or execute the text. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The text need not denote a valid mathematical expression. Equality compares
text rather than mathematical meaning.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym("x"),sym("y"),sym("z")]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

See also [`SymExpr`](@ref).

Role: construct
"""
sym(x::AbstractString) = SymExpr(String(x))

"""
    sym(x::Number) -> SymExpr

Return the textual symbolic form of a number.

# Arguments
`x` is a `Number`. Its `string` conversion must be defined; conversion
errors propagate.

# Returns
`SymExpr` with `text == string(x)`, in the
[Algebra symbolic API](../specs/algebra.md#definitions).

# Guarantees
1. The result stores the number's Julia string representation. [proved/agent: by construction]
2. No numeric computation, parsing or variable binding is performed. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Printed text need not preserve an arbitrary number type's full mathematical
meaning, precision or evaluation rules. This is not exact scalar embedding
into a general computer algebra system.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

See also [`SymExpr`](@ref).

Role: construct
"""
sym(x::Number) = SymExpr(string(x))

Base.show(io::IO, x::SymExpr) = print(io, x.text)
Base.string(x::SymExpr) = x.text
Base.zero(::Type{SymExpr}) = sym(0)
Base.one(::Type{SymExpr}) = sym(1)
Base.iszero(x::SymExpr) = x.text == "0"
Base.:(==)(a::SymExpr, b::SymExpr) = a.text == b.text
Base.:(==)(a::SymExpr, b::Number) = a.text == string(b)
Base.:(==)(a::Number, b::SymExpr) = b == a

_paren(x::SymExpr) = occursin(r"[+\-*/]", x.text) ? "(" * x.text * ")" : x.text
_sx(x) = sym(x)
function _bin(op::String, a, b)
    x, y = _sx(a), _sx(b)
    op == "+" && iszero(x) && return y
    op == "+" && iszero(y) && return x
    op == "-" && iszero(y) && return x
    op == "*" && (iszero(x) || iszero(y)) && return sym(0)
    op == "*" && x == 1 && return y
    op == "*" && y == 1 && return x
    op == "/" && y == 1 && return x
    SymExpr(_paren(x) * op * _paren(y))
end

for op in (:+, :-, :*, :/)
    text = string(op)
    @eval Base.$op(a::SymExpr, b::SymExpr) = _bin($text, a, b)
    @eval Base.$op(a::SymExpr, b::Number) = _bin($text, a, b)
    @eval Base.$op(a::Number, b::SymExpr) = _bin($text, a, b)
end
Base.:-(a::SymExpr) = SymExpr("-" * _paren(a))
Base.:^(a::SymExpr, n::Integer) = n == 0 ? sym(1) : n == 1 ? a : SymExpr(_paren(a) * "^" * string(n))

function substitute(x::SymExpr, values::Dict{String,<:Number})
    text = x.text
    for (name, value) in values
        text = replace(text, Regex("\\b" * name * "\\b") => string(value))
    end
    SymExpr(text)
end
