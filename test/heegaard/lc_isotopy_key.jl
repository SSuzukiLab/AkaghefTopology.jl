# Relations from the accepted comparison: commute, braid, torus order 6,
# the genus-2 5-chain, the lantern, and the four pairs that must differ.
# The expected value is the relation (keys agree, or the pair has two keys),
# not a stored word.

using Test
using Topology

const H = Topology.HeegaardCodes

const DIGON_ROWS = (
    "lc.isotopy.digon.remove.00",
    "lc.isotopy.digon.remove.01",
    "lc.isotopy.digon.remove.10",
    "lc.isotopy.digon.remove.11",
)

function blank_lc(genus, words)
    H.LickorishCode(
        genus, words,
        Dict{String,String}(),
        Dict{String,Vector{String}}(),
        Dict{String,String}(),
        Dict{String,Int}(),
        nothing,
    )
end

function k_surface(g)
    words = Dict{String,Vector{Int}}(
        String(k) => Int.(copy(w)) for (k, w) in H._lc_k_constants(g))
    blank_lc(Int(g), words)
end

function owners(lc)
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in lc.words
        for x in word
            if x > 0
                plus[abs(x)] = String(name)
            else
                minus[abs(x)] = String(name)
            end
        end
    end
    plus, minus
end

function algebraic(lc, a, b)
    plus, minus = owners(lc)
    n = 0
    for lab in keys(plus)
        pair = (plus[lab], minus[lab])
        if pair == (a, b)
            n += 1
        elseif pair == (b, a)
            n -= 1
        end
    end
    n
end

function k_algebra(lc)
    kn = String.(H._lc_k_names(lc.genus))
    rows = Tuple{String,String,Int}[]
    for a in sort(collect(keys(lc.words))), b in kn
        String(a) == b && continue
        push!(rows, (String(a), b, algebraic(lc, String(a), b)))
    end
    rows
end

function reduce_digons(lc)
    state = lc
    n = 0
    while true
        hit = nothing
        row = ""
        for candidate in DIGON_ROWS
            hits = H.movable(state, candidate)
            isempty(hits) && continue
            hit = first(hits)
            row = candidate
            break
        end
        hit === nothing && return state
        before = k_algebra(state)
        state = H.move(state, hit).output
        after = k_algebra(state)
        before == after || error("algebra changed on $row")
        n += 1
        n > 4000 && error("digon removal did not finish")
        maximum(length(w) for w in values(state.words)) > 800 &&
            error("word stayed long after digon removal")
    end
end

function apply_twist(lc, support, acted, q)
    hits = H.movable(
        lc, "lc.twist.annulus";
        q=q, support=String(support), acted_components=String.(acted),
    )
    length(hits) == 1 || error("witness count $(length(hits)) for support=$support")
    H.move(lc, only(hits)).output
end

function twist_all(lc, steps, acted)
    state = lc
    for (support, q) in steps
        state = apply_twist(state, support, acted, q)
        state = reduce_digons(state)
    end
    state
end

function add_parallel(lc, curve; side="left")
    hits = H.movable(lc, "lc.isotopy.parallel"; curve=String(curve), side=side)
    w = only(hits)
    H.move(lc, w).output, String(only(w.extra["created_components"]))
end

function with_parallels(lc, names; side="left")
    out = lc
    mapping = Dict{String,String}()
    for curve in names
        out, copy_name = add_parallel(out, curve; side=side)
        mapping[String(curve)] = copy_name
    end
    out, mapping
end

function torus_lc()
    H._lc_parse_text("""
    genus: 1
    l1 : -2 -1
    m1 : 1 3
    a : -4 -3
    b : 2 4
    """)
end

function keyed(lc, curves)
    rows = H.lc_isotopy_key(lc, collect(String, curves))
    for row in rows
        @test row.algebra_ok
        @test row.failure === nothing
        @test row.key isa String
        @test !isempty(row.words)
    end
    rows
end

function agree(left, right, curves)
    L = keyed(left, curves)
    R = keyed(right, curves)
    @test [row.key for row in L] == [row.key for row in R]
end

function differ(left, right, curves)
    L = keyed(left, curves)
    R = keyed(right, curves)
    @test any(i -> L[i].key != R[i].key, eachindex(L))
end

function geometric(lc, a, b)
    plus, minus = owners(lc)
    n = 0
    for lab in keys(plus)
        pair = (plus[lab], minus[lab])
        (pair == (a, b) || pair == (b, a)) && (n += 1)
    end
    n
end

function cycle_start(seq, head)
    i = findfirst(==(head), seq)
    i === nothing && return seq
    vcat(seq[i:end], seq[1:i - 1])
end

function face_curve_sequences(lc)
    faces, _ = H._lc_faces(lc.words)
    plus, minus = owners(lc)
    seqs = Vector{String}[]
    for face in faces
        seq = String[]
        for (lab, sg, _) in face.corners
            name = sg == "+" ? plus[lab] : minus[lab]
            (isempty(seq) || seq[end] != name) && push!(seq, name)
        end
        if length(seq) > 1 && seq[1] == seq[end]
            pop!(seq)
        end
        push!(seqs, seq)
    end
    seqs
end

function lantern_lc()
    blank_lc(2, Dict{String,Vector{Int}}(
        "l1" => [-1, -3],
        "l2" => [-2, 4],
        "m1" => [1, -10, -14],
        "m2" => [2, -11, 17],
        "n1" => [3, -20, -8, -4, 9, -15],
        "x" => [8, -13, 16, -9, -18, 12],
        "y" => [10, -19, -12, 11, 13, 21],
        "z" => [14, 15, -16, -17, 18, 19, 20, -21],
    ))
end

@testset "lc_isotopy_key" begin
    torus = torus_lc()
    saved = Dict(name => copy(word) for (name, word) in torus.words)
    one = H.lc_isotopy_key(torus, "a")
    two = H.lc_isotopy_key(torus, ["a", "b"])
    @test torus.words == saved
    @test one.key == two[1].key
    @test one.algebra_ok
    @test two[2].key isa String
    @test_throws H.KodaInputError H.lc_isotopy_key(torus, "missing")

    g2, parallel_of = with_parallels(k_surface(2), H._lc_k_names(2))
    acted = [parallel_of[c] for c in H._lc_k_names(2)]

    @testset "commute" begin
        println("commute $(time())")
        flush(stdout)
        left = twist_all(g2, [("m2", 1), ("m1", 1)], acted)
        right = twist_all(g2, [("m1", 1), ("m2", 1)], acted)
        agree(left, right, acted)
    end

    @testset "braid" begin
        println("braid $(time())")
        flush(stdout)
        left = twist_all(g2, [("l1", 1), ("m1", 1), ("l1", 1)], acted)
        right = twist_all(g2, [("m1", 1), ("l1", 1), ("m1", 1)], acted)
        agree(left, right, acted)
    end

    @testset "torus order 6" begin
        println("torus $(time())")
        flush(stdout)
        steps = Tuple{String,Int}[]
        for _ in 1:6
            append!(steps, [("m1", 1), ("l1", 1)])
        end
        image = twist_all(torus, steps, ["a", "b"])
        agree(image, torus, ["a", "b"])
    end

    @testset "negative controls" begin
        println("negative $(time())")
        flush(stdout)
        differ(
            twist_all(g2, [("m1", 1), ("l1", 1)], acted),
            twist_all(g2, [("l1", 1), ("m1", 1)], acted),
            acted,
        )
        cubed = Tuple{String,Int}[]
        for _ in 1:3
            append!(cubed, [("b", 1), ("a", 1)])
        end
        differ(twist_all(torus, cubed, ["a", "b"]), torus, ["a", "b"])
        differ(twist_all(torus, [("a", 1)], ["a", "b"]), torus, ["a", "b"])
        differ(
            twist_all(torus, [("a", 1)], ["a", "b"]),
            twist_all(torus, [("a", -1)], ["a", "b"]),
            ["a", "b"],
        )
    end

    @testset "5-chain order 6" begin
        println("chain $(time())")
        flush(stdout)
        steps = Tuple{String,Int}[]
        for _ in 1:6
            append!(steps, [("m2", 1), ("l2", 1), ("n1", 1), ("l1", 1), ("m1", 1)])
        end
        image = twist_all(g2, steps, acted)
        agree(image, g2, acted)
    end

    @testset "lantern" begin
        println("lantern $(time())")
        flush(stdout)
        base = lantern_lc()
        H.validate_lc_input(base)
        built = H.realize_lc(base)
        @test built.correspondence["surface_chi"] == -2
        for a in ("x", "y", "z"), b in ("x", "y", "z")
            a >= b && continue
            @test geometric(base, a, b) == 2
            @test algebraic(base, a, b) == 0
        end
        triangles = [
            seq for seq in face_curve_sequences(base) if
            length(seq) == 3 && sort(seq) == ["x", "y", "z"]
        ]
        @test length(triangles) == 1
        positive = cycle_start(triangles[1], "x")
        clockwise = cycle_start(reverse(positive), "x")
        @test positive == ["x", "z", "y"]
        @test clockwise == ["x", "y", "z"]
        lc = base
        boundary = String[]
        for core in ("l1", "l2"), side in ("left", "right")
            lc, name = add_parallel(lc, core; side=side)
            push!(boundary, name)
        end
        acted_lantern = String[]
        for core in ("m1", "m2", "n1"), side in ("left", "right")
            lc, name = add_parallel(lc, core; side=side)
            push!(acted_lantern, name)
        end
        left_steps = [(clockwise[3], 1), (clockwise[2], 1), (clockwise[1], 1)]
        right_steps = [(name, 1) for name in boundary]
        agree(
            twist_all(lc, left_steps, acted_lantern),
            twist_all(lc, right_steps, acted_lantern),
            acted_lantern,
        )
    end
end
