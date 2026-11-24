function invariant_expression(v::VirtualLink, ::URAlgebra)
    _require_closed_ograph_input(v, "U(R) MST")
    components, signs = real_gauss_code(v)
    length(components) == 1 && !isempty(only(components)) ||
        throw(ArgumentError("U(R) MST requires one nonempty real-Gauss component"))
    word = reverse(only(components))
    n = length(signs)
    all(s -> s in (-1, 1), signs) &&
        sort(word) == vcat(collect(-n:-1), collect(1:n)) ||
        throw(ArgumentError("U(R) MST requires one positive and negative visit per vertex"))
    word = circshift(word, 1 - findfirst(>(0), word))
    section = Dict(zip(word, cumsum(word .> 0)))
    matrix = Matrix{BigInt}(I, n, n)
    matrix -= circshift(Matrix{BigInt}(I, n, n), (1, 0))
    for vertex in 1:n
        matrix[section[vertex], section[-vertex]] += signs[vertex]
    end
    URMSTExpression(matrix)
end

Z_MST(v::VirtualLink, H::URAlgebra; strategy=:auto) =
    evaluate_invariant(simplify_invariant(invariant_expression(v,H); strategy=strategy))
