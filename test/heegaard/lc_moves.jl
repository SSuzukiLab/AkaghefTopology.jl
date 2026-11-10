# cf: lib/topology.jl/docs/specs/heegaard-moves.md §8 — inverse, checker, genus, cyclic wrap
using Test
using Topology
using SHA

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const CHECKER = joinpath(ROOT, "threads", "PJ260922", "code", "lc_check.py")

if !isdefined(H, :load_heegaard_moves)
    Base.include(H, joinpath(ROOT, "lib", "topology.jl", "src", "heegaard", "lc_moves.jl"))
end

function parse_lc(text)
    H._lc_parse_text(text)
end

function checker(lc)
    dir = mktempdir()
    path = joinpath(dir, "probe.lc")
    H.write_lc(path, lc)
    out = read(Cmd(`python $CHECKER $path`; ignorestatus=true), String)
    occursin("\"status\": \"pass\"", out) || occursin("\"status\":\"pass\"", out), out
end

function cyclic_equal(left::Vector{Int}, right::Vector{Int})
    length(left) == length(right) || return false
    isempty(left) && return true
    any(i -> left == vcat(right[i:end], right[1:i - 1]), 1:length(right))
end

function relabel(word, mapping)
    [sign(entry) * mapping[abs(entry)] for entry in word]
end

function same_up_to(left, right)
    H._component_names(left) == H._component_names(right) || return false
    labels = sort(unique(abs.(reduce(vcat, values(left.words); init=Int[]))))
    right_labels = sort(unique(abs.(reduce(vcat, values(right.words); init=Int[]))))
    length(labels) == length(right_labels) || return false
    # K labels are pinned. Remaining labels are matched by trying the identity first,
    # then the bijection that sends left's non-K labels onto right's in sorted order
    # only when both sides use the same crossing pattern. The composed move maps are
    # checked by the caller; this predicate searches permutations of the non-K labels.
    limit = left.genus <= 0 ? 0 : 3 * left.genus - 2
    fixed = [label for label in labels if label <= limit]
    all(label in right_labels for label in fixed) || return false
    free_left = [label for label in labels if label > limit]
    free_right = [label for label in right_labels if label > limit]
    length(free_left) == length(free_right) || return false
    length(free_left) > 6 && return false
    for perm in permutations(free_right)
        mapping = Dict(label => label for label in fixed)
        for (src, dst) in zip(free_left, perm)
            mapping[src] = dst
        end
        good = true
        for name in keys(left.words)
            cyclic_equal(relabel(left.words[name], mapping), right.words[name]) || (good = false; break)
        end
        good && return true
    end
    false
end

function permutations(values)
    isempty(values) && return [Int[]]
    out = Vector{Vector{Int}}()
    for i in eachindex(values)
        rest = vcat(values[1:i - 1], values[i + 1:end])
        for tail in permutations(rest)
            push!(out, vcat(values[i], tail))
        end
    end
    out
end

function algebraic_intersections(lc, alphas, betas)
    words = lc.words
    owners = Dict{Int,Tuple{String,Int}}()
    for (name, word) in words
        for entry in word
            owners[sign(entry) * abs(entry) == entry ? abs(entry) + 0 : abs(entry)] = (name, sign(entry))
        end
    end
    # rebuild properly
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in words
        for entry in word
            entry > 0 ? (plus[abs(entry)] = name) : (minus[abs(entry)] = name)
        end
    end
    matrix = zeros(Int, length(alphas), length(betas))
    for (i, alpha) in enumerate(alphas), (j, beta) in enumerate(betas)
        for label in keys(plus)
            pair = (plus[label], minus[label])
            if pair == (alpha, beta)
                matrix[i, j] += 1
            elseif pair == (beta, alpha)
                matrix[i, j] -= 1
            end
        end
    end
    matrix
end

const DIGON = """
genus: 1
l1 : -1 2 -3
m1 : 1
c : -2 3
"""

const DIGON_WRAP = """
genus: 1
l1 : -1 2 -3
m1 : 1 5
c : 3 -5 -2
"""

const TRIANGLE = """
genus: 1
l1 : 3 -1
m1 : 2 1
c : -3 -2
"""

@testset "lc moves" begin
@testset "catalog and isotopy" begin
    ids = H.heegaard_move_ids()
    @test "lc.isotopy.digon.remove.00" in ids
    @test "lc.twist.annulus" in ids
    @test "hc.slide.pants" in ids
    @test "hc.stabilization.once" in ids

    digon = parse_lc(DIGON)
    hits = H.movable(digon, "lc.isotopy.digon.remove.00")
    @test length(hits) == 1
    forward = H.move(digon, only(hits))
    ok, report = checker(forward.output)
    @test ok
    @test forward.output.genus == digon.genus
    back = H.move(forward.output, forward.inverse)
    @test back.output.genus == digon.genus
    @test same_up_to(digon, back.output)

    wrapped = parse_lc(DIGON_WRAP)
    whits = H.movable(wrapped, "lc.isotopy.digon.remove.00")
    @test length(whits) == 1
    wforward = H.move(wrapped, only(whits))
    wok, wreport = checker(wforward.output)
    @test wok
    wback = H.move(wforward.output, wforward.inverse)
    @test same_up_to(wrapped, wback.output)

    tri = parse_lc(TRIANGLE)
    thits = H.movable(tri, "lc.isotopy.vertex.000")
    @test length(thits) == 1
    tforward = H.move(tri, only(thits))
    tok, treport = checker(tforward.output)
    @test tok
    tback = H.move(tforward.output, tforward.inverse)
    @test same_up_to(tri, tback.output)
    @test tforward.output.genus == 1
end

@testset "marks and unresolved" begin
    digon = parse_lc(DIGON)
    marked = Topology.Heegaard.LickorishConfiguration(
        digon.genus, digon.words, digon.places, Dict("D" => ["F1"]), digon.points, digon.weights, nothing)
    @test_throws H.LCMoveError H.movable(marked, "lc.isotopy.digon.remove.00")
    @test_throws H.LCMoveError H.movable(digon, "lc.fixed.unresolved")
end

# Spec §8: a partial twist with β′ = a+2b changes the intersection determinant from 1 to 2.
# This test states that property. It is not adjusted to the current expander.
@testset "partial twist determinant" begin
    torus = parse_lc("""
    genus: 1
    l1 : -2 -1
    m1 : 1 3
    a : -4 -3
    b : 2 4
    """)
    before = algebraic_intersections(torus, ["a"], ["b"])[1, 1]
    @test abs(before) == 1
    hits = H.movable(torus, "lc.twist.annulus"; q=2, support="a", acted_components=["b"])
    @test !isempty(hits)
    result = H.move(torus, only(hits))
    ok, report = checker(result.output)
    @test ok
    after = algebraic_intersections(result.output, ["a"], ["b"])[1, 1]
    @test abs(after) == 2
    @test result.output.genus == torus.genus
end

@testset "handle slide empty band" begin
    lc = H.read_lc(joinpath(ROOT, "threads", "PJ260922", "code", "lc-fixtures", "p28-after.lc"))
    hc = H.HeegaardLayer(lc, Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta),
        Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
    slides = H.movable(hc, "hc.slide.pants"; curves=["alpha1", "alpha2"], band=Any[])
    @test !isempty(slides)
    slid = H.move(hc, only(slides))
    sok, sreport = checker(slid.output.lc)
    @test sok
    @test slid.output.lc.genus == lc.genus
    back = H.move(slid.output, slid.inverse)
    @test back.output.lc.genus == lc.genus
end

@testset "stabilization genus" begin
    base = H.read_lc(joinpath(ROOT, "threads", "PJ260922", "code", "lc-fixtures", "operations", "g1-base.lc"))
    g1 = H.HeegaardLayer(base, Dict("alpha1" => :alpha, "beta1" => :beta),
        Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
    stabs = H.movable(g1, "hc.stabilization.once"; disk="F1", k_identification="ladder_end")
    @test !isempty(stabs)
    stabilized = H.move(g1, only(stabs))
    @test stabilized.output.lc.genus == base.genus + 1
    tok, treport = checker(stabilized.output.lc)
    @test tok
    undone = H.move(stabilized.output, stabilized.inverse)
    @test undone.output.lc.genus == base.genus
end
end

println("LC_MOVES_TESTS_DONE")
