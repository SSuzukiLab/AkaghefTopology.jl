"""A lightweight exact symbolic expression used by the migrated Live Scripts."""
struct SymExpr
    text::String
end

sym(x::SymExpr) = x
sym(x::Symbol) = SymExpr(String(x))
sym(x::AbstractString) = SymExpr(String(x))
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
