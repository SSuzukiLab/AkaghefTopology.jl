using Test, Topology, LinearAlgebra

const S3_TABLE = [
    1 2 3 4 5 6
    2 3 1 6 4 5
    3 1 2 5 6 4
    4 5 6 1 2 3
    5 6 4 3 1 2
    6 4 5 2 3 1
]

function _group_inverses(table::Matrix{Int})
    order = size(table, 1)
    identity = findfirst(i -> all(j -> table[i, j] == j, 1:order), 1:order)
    [findfirst(j -> table[i, j] == identity, 1:order) for i in 1:order]
end

function count_pi1_homomorphisms(presentation::Heegaard.GroupPresentation, table::Matrix{Int})
    order = size(table, 1)
    ngen = length(presentation.generators)
    ngen == 0 && return 1
    inv = _group_inverses(table)
    id = findfirst(i -> all(j -> table[i, j] == j, 1:order), 1:order)
    count = 0
    for tuple in Iterators.product((1:order for _ in 1:ngen)...)
        imgs = collect(tuple)
        all(rel -> begin
            s = id
            for (g, σ) in rel
                h = imgs[g]
                s = σ > 0 ? table[s, h] : table[s, inv[h]]
            end
            s == id
        end, presentation.relators) && (count += 1)
    end
    count
end

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data"))
const DRAW_B = normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260922", "data", "261107-ld-draw-b"))

@testset "Kuperberg domain and algebra guards" begin
    path = joinpath(DATA, "group-15", "case-01.hc")
    @test_throws DomainError Z_Kuperberg(path, sweedler_algebra())
    hc = Heegaard.read_hc(path)
    alg = group_algebra(S3_TABLE)
    expected = Z_Kuperberg(hc, alg)
    broken = deepcopy(alg)
    broken.right_integral.terms[(1,)] = 2
    @test isapprox(Z_Kuperberg(hc, broken), expected; atol=1e-8)
    zeroed = deepcopy(alg)
    empty!(zeroed.right_integral.terms)
    @test_throws DomainError Z_Kuperberg(hc, zeroed)
end

@testset "Kuperberg pi1 relator uses beta word order (OS04)" begin
    hc = Heegaard.read_hc(joinpath(DATA, "group-15", "case-01.hc"))
    presentation = heegaard_pi1_presentation(hc)
    @test presentation.relators == [[(1, -1), (1, 1), (1, -1)]]
end


const Z2_TABLE = [1 2; 2 1]
const Z3_TABLE = [mod(i + j - 2, 3) + 1 for i in 1:3, j in 1:3]
const Z5_TABLE = [mod(i + j - 2, 5) + 1 for i in 1:5, j in 1:5]
const KUPERBERG_TABLES = (
    ("C[Z2]", Z2_TABLE),
    ("C[Z3]", Z3_TABLE),
    ("C[Z5]", Z5_TABLE),
    ("C[S3]", S3_TABLE),
)
function _hom_signature(path)
    hc = Heegaard.read_hc(path)
    pres = heegaard_pi1_presentation(hc)
    ntuple(i -> count_pi1_homomorphisms(pres, KUPERBERG_TABLES[i][2]), length(KUPERBERG_TABLES))
end

function _assert_hom_and_dual(path)
    hc = Heegaard.read_hc(path)
    pres = heegaard_pi1_presentation(hc)
    for (_, table) in KUPERBERG_TABLES
        expected = count_pi1_homomorphisms(pres, table)
        alg = group_algebra(table)
        @test isapprox(real(Z_Kuperberg(hc, alg)), expected; atol=1e-8)
        @test isapprox(real(Z_Kuperberg(hc, dual_algebra(alg))), expected; atol=1e-8)
    end
end

const KUPERBERG_S3_CASES = [
    joinpath(DATA, "group-01", "case-01.hc"),
    joinpath(DATA, "group-01", "case-03.hc"),
    joinpath(DATA, "group-15", "case-01.hc"),
]
const KUPERBERG_S2XS1 = joinpath(DATA, "group-13", "case-01.hc")
const KUPERBERG_G13 = joinpath(DRAW_B, "g13-c01.hc")
const KUPERBERG_G31 = joinpath(DRAW_B, "g31-c01.hc")

@testset "Kuperberg matches |Hom| on several S3 diagrams" begin
    for path in KUPERBERG_S3_CASES
        @test _hom_signature(path) == (1, 1, 1, 1)
        _assert_hom_and_dual(path)
    end
end

@testset "Kuperberg matches |Hom| on S2xS1" begin
    @test isfile(KUPERBERG_S2XS1)
    @test _hom_signature(KUPERBERG_S2XS1) == (2, 3, 5, 6)
    _assert_hom_and_dual(KUPERBERG_S2XS1)
end

@testset "Kuperberg matches |Hom| on g13-c01 and g31-c01" begin
    @test isfile(KUPERBERG_G13)
    @test isfile(KUPERBERG_G31)
    _assert_hom_and_dual(KUPERBERG_G13)
    _assert_hom_and_dual(KUPERBERG_G31)
    @test _hom_signature(KUPERBERG_G13) == _hom_signature(KUPERBERG_S2XS1)
end

const KUPERBERG_T3 = [
    joinpath(DATA, "group-19", "case-01.hc"),
    joinpath(DATA, "group-32", "case-02.hc"),
]

function _commuting_triple_signature()
    ntuple(length(KUPERBERG_TABLES)) do i
        table = KUPERBERG_TABLES[i][2]
        order = size(table, 1)
        count = 0
        for a in 1:order, b in 1:order, c in 1:order
            table[a, b] == table[b, a] && table[a, c] == table[c, a] &&
                table[b, c] == table[c, b] && (count += 1)
        end
        count
    end
end

function _lens_power(path)
    hc = Heegaard.read_hc(path)
    pres = heegaard_pi1_presentation(hc)
    length(pres.generators) == 1 || return nothing
    length(pres.relators) == 1 || return nothing
    relator = only(pres.relators)
    length(relator) >= 2 || return nothing
    all(letter -> letter[1] == 1 && letter[2] == relator[1][2], relator) || return nothing
    length(relator)
end

@testset "Kuperberg later diagrams" begin
@testset "Kuperberg matches |Hom| on lens diagrams built from o-graphs" begin
    # The .hc search for a genus-1 pure power found nothing. Replaced, not deleted,
    # by the registered lens o-graphs the conversion lane can turn into diagrams.
    include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))
    for name in ("lens_2_1", "lens_3_1_mst2_example23", "lens_3_2")
        v = vl_example(name)
        gauss, orientation = Topology.real_gauss_code(v)
        converted = Topology.Heegaard.odata_to_heegaard_koda(
            Topology.Heegaard.OData(deepcopy(gauss), copy(orientation), Dict{String,Any}()))
        @test converted.value !== nothing
        diagram = converted.value.diagram
        pres = Heegaard.fundamental_group_presentation(diagram)
        for (_, table) in KUPERBERG_TABLES
            expected = count_pi1_homomorphisms(pres, table)
            alg = group_algebra(table)
            @test isapprox(Z_Kuperberg(diagram, alg), expected; atol=1e-6)
            @test isapprox(Z_Kuperberg(diagram, dual_algebra(alg)), expected; atol=1e-6)
        end
    end
end

@testset "Kuperberg matches |Hom| on T3 diagrams" begin
    wanted = _commuting_triple_signature()
    for path in KUPERBERG_T3
        @test isfile(path)
        @test _hom_signature(path) == wanted
        _assert_hom_and_dual(path)
    end
end

@testset "Kuperberg passage orientations are not identically +1 when beta has + labels" begin
    hc = Heegaard.read_hc(joinpath(DATA, "group-15", "case-01.hc"))
    expr = invariant_expression(hc, exterior_algebra(2))
    _, orientation = real_gauss_code(expr.link)
    @test any(<(0), orientation)
end

@testset "Kac-Paljutkin is 1 on trivial pi1 and agrees on another diagram of the same manifold" begin
    kp = kac_paljutkin_algebra()
    trivial = String[]
    for path in vcat(KUPERBERG_S3_CASES, [KUPERBERG_G13, KUPERBERG_G31])
        _hom_signature(path) == (1, 1, 1, 1) && push!(trivial, path)
    end
    @test length(trivial) >= 2
    values = ComplexF64[]
    for path in trivial
        value = Z_Kuperberg(path, kp)
        push!(values, value)
        @test isapprox(value, 1; atol=1e-6)
    end
    @test all(isapprox(value, values[1]; atol=1e-6) for value in values)
    t3_values = [Z_Kuperberg(path, kp) for path in KUPERBERG_T3]
    @test length(t3_values) == 2
    @test all(isapprox(value, t3_values[1]; atol=1e-6) for value in t3_values)
end
end

function _basis_changed(algebra, M)
    n = algebra.dimension
    Minv = inv(M)
    vec_of(factor) = begin
        v = zeros(ComplexF64, n)
        for ((i,), value) in factor.terms
            v[i] = value
        end
        v
    end
    mul = Dict{Tuple, ComplexF64}()
    for ((i, j, k), value) in algebra.multiplication.terms, a in 1:n, b in 1:n, c in 1:n
        coeff = M[i, a] * M[j, b] * Minv[c, k] * value
        abs(coeff) <= 1e-12 && continue
        mul[(a, b, c)] = get(mul, (a, b, c), 0.0 + 0im) + coeff
    end
    cop = Dict{Tuple, ComplexF64}()
    for ((i, j, k), value) in algebra.coproduct.terms, a in 1:n, b in 1:n, c in 1:n
        coeff = M[i, a] * Minv[b, j] * Minv[c, k] * value
        abs(coeff) <= 1e-12 && continue
        cop[(a, b, c)] = get(cop, (a, b, c), 0.0 + 0im) + coeff
    end
    Sold = zeros(ComplexF64, n, n)
    for ((output, input), value) in algebra.antipode.terms
        Sold[output, input] = value
    end
    Snew = Minv * Sold * M
    antipode = Dict{Tuple, ComplexF64}((o, i) => Snew[o, i] for o in 1:n, i in 1:n if abs(Snew[o, i]) > 1e-12)
    unit_new = Minv * vec_of(algebra.unit)
    counit_new = transpose(M) * vec_of(algebra.counit)
    lambda_new = transpose(M) * vec_of(algebra.right_integral)
    Lambda_new = Minv * vec_of(algebra.right_cointegral)
    sparsev(v) = Dict{Tuple, ComplexF64}((i,) => v[i] for i in 1:n if abs(v[i]) > 1e-12)
    empty = SparseFactor(Int[], Dict{Tuple, ComplexF64}())
    FiniteHopfAlgebra(n,
        SparseFactor(Int[], mul), SparseFactor(Int[], sparsev(unit_new)),
        SparseFactor(Int[], cop), SparseFactor(Int[], sparsev(counit_new)),
        empty, empty, SparseFactor(Int[], antipode),
        SparseFactor(Int[], sparsev(lambda_new)), SparseFactor(Int[], sparsev(Lambda_new)),
        SparseFactor(Int[], sparsev(transpose(Snew) * lambda_new)),
        SparseFactor(Int[], sparsev(inv(Snew) * Lambda_new))), Snew
end

@testset "Kuperberg is unchanged by a non-orthogonal basis change" begin
    paths = [
        joinpath(DRAW_B, "g12-c01.hc"),
        joinpath(DRAW_B, "g15-c01.hc"),
        normpath(joinpath(DRAW_B, "..", "261107-ld-draw-batch", "src", "g19-c01.hc")),
    ]
    algebras = [
        group_algebra(S3_TABLE),
        dual_algebra(group_algebra(S3_TABLE)),
        kac_paljutkin_algebra(),
    ]
    for algebra in algebras
        n = algebra.dimension
        M = Matrix{ComplexF64}(I, n, n)
        for i in 2:n
            M[i, 1] = i
            M[1, i] = 1 / (i + 1)
        end
        changed, Snew = _basis_changed(algebra, M)
        @test !isapprox(Snew, transpose(Snew); atol=1e-8)
        for path in paths
            @test isapprox(Z_Kuperberg(path, algebra), Z_Kuperberg(path, changed); atol=1e-5)
        end
    end
end

@testset "HC and diagram entrances mark the same crossings negative" begin
    samples = [
        joinpath(DATA, "group-15", "case-01.hc"),
        joinpath(DRAW_B, "g12-c01.hc"),
        joinpath(DRAW_B, "g13-c01.hc"),
        joinpath(DRAW_B, "g15-c01.hc"),
        joinpath(DRAW_B, "g31-c01.hc"),
    ]
    checked = 0
    for path in samples
        isfile(path) || continue
        hc = Heegaard.read_hc(path)
        forgot = Heegaard.forget_reference(Heegaard.realize_hc(hc))
        forgot.diagram === nothing && continue
        checked += 1
        hc_neg = Set(string(label) for label in kuperberg_hc_negative_labels(hc))
        diag_neg = Set(point.source_label for point in forgot.diagram.code.code.intersections
            if point.sign < 0 && point.source_label !== nothing)
        @test hc_neg == diag_neg
    end
    @test checked >= 3
end

@testset "Z_Kuperberg on converted diagrams agrees with Z_MST" begin
    include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))
    names = ("lens_2_1", "lens_3_1_mst2_example23", "lens_3_2", "odata_variant",
        "poincare", "quaternion", "s2xs1", "s2xs1_2vertex", "s3_abalone")
    # All nine examples, but only C[Z2], C[S3], their duals, and Kac–Paljutkin.
    # Z3, Z5 and the Kac–Paljutkin dual are in the scratch comparison, not here.
    algebras = [
        group_algebra(Z2_TABLE),
        dual_algebra(group_algebra(Z2_TABLE)),
        group_algebra(S3_TABLE),
        dual_algebra(group_algebra(S3_TABLE)),
        kac_paljutkin_algebra(),
    ]
    for name in names, algebra in algebras
        v = vl_example(name)
        gauss, orientation = Topology.real_gauss_code(v)
        converted = Topology.Heegaard.odata_to_heegaard_koda(
            Topology.Heegaard.OData(deepcopy(gauss), copy(orientation), Dict{String,Any}()))
        @test converted.value !== nothing
        @test isapprox(Z_Kuperberg(converted.value.diagram, algebra), Z_MST(v, algebra); atol=1e-5)
    end
end


