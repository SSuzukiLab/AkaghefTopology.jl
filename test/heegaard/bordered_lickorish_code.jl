# Bordered LC, docs/specs/bordered-lickorish-code.md §2–§6.
# Counts are the ones written in threads/PJ260922/data/261116-bordered/euler-check.md.
using Test
using Topology

const H = Topology.Heegaard

function blc_roundtrip(text::AbstractString)
    piece = H.read_bordered_lc(text)
    again = H.read_bordered_lc(H.bordered_lc_text(piece))
    @test again == piece
    piece
end

function blc_counts(piece; reference::Bool=false)
    d = H.diagnose_bordered_lc(piece; reference=reference)
    (d.V, d.E, d.F, d.boundary_cycles, d.chi, d.connected)
end

function blc_section(f)
    try
        f()
        nothing
    catch e
        e isa H.BorderedLCError || rethrow()
        e
    end
end

function blc_rotations(xs)
    n = length(xs)
    [circshift(xs, -i) for i in 0:(n - 1)]
end

const ANN0 = """
piece ann0 : genus=0 boundaries=2
marking : reference
boundary B0 base=o0 : o0 s0 q0
boundary B1 base=o1 : o1 q1 s1
Q arc q q0 -> q1 : []
arc gamma s0 -> s1 : []
"""

const ANN1 = """
piece ann1 : genus=0 boundaries=2
marking : reference
boundary B0 base=o0 : o0 s0 q0
boundary B1 base=o1 : o1 q1 s1
Q arc q q0 -> q1 : [+1]
arc gamma s0 -> s1 : [-1]
"""

const BARE = """
piece bare-annulus : genus=0 boundaries=2
marking : reference
boundary B0 base=p0 : p0 s0
boundary B1 base=p1 : p1 s1
arc gamma s0 -> s1 : []
"""

const VERTEX_L = """
piece vertex-L : genus=0 boundaries=1
marking : reference
boundary d base=Bm : Bm Am Cp Bp Ap Cm
K arc A Am -> Ap : [+2 +1]
K arc B Bm -> Bp : [+3 -1]
arc C Cm -> Cp : [-3 -2]
"""

const VERTEX_R = """
piece vertex-R : genus=0 boundaries=1
marking : reference
boundary d base=Bm : Bm Am Cp Bp Ap Cm
K arc A Am -> Ap : [+1 +2]
K arc B Bm -> Bp : [-1 +3]
arc C Cm -> Cp : [-2 -3]
"""

const STAB_L = """
piece stab-L : genus=0 boundaries=1
marking : reference
boundary d base=u0 : u0 v0 u1 v1
Q arc q u0 -> u1 : [-4]
Q arc r v0 -> v1 : [+4]
"""

const STAB_R = """
piece stab-R : genus=1 boundaries=1
marking : reference
boundary d base=u0 : u0 v0 u1 v1
Q arc q u0 -> u1 : [+1]
Q arc r v0 -> v1 : [-2]
loop a : (-1 -3)
loop b : (+2 +3)
"""

const SLIDE_L = """
piece slide-L : genus=0 boundaries=3
marking : reference
boundary A base=oA : oA qA rA
boundary B base=oB : oB qB
boundary C base=oC : oC rC
Q arc q qA -> qB : [+1]
Q arc r rA -> rC : [+2]
loop a : (-1)
loop b : (-2)
"""

const SLIDE_R = """
piece slide-R : genus=0 boundaries=3
marking : reference
boundary A base=oA : oA qA rA
boundary B base=oB : oB qB
boundary C base=oC : oC rC
Q arc q qA -> qB : [+1]
Q arc r rA -> rC : [+3 +2]
loop a_new : (-1 -3)
loop b : (-2)
"""

const SLIDE_BAD = """
piece slide-bad : genus=0 boundaries=3
marking : reference
boundary A base=oA : oA qA rA
boundary B base=oB : oB qB
boundary C base=oC : oC rC
Q arc q qA -> qB : [+1]
Q arc r rA -> rC : [+2 +3]
loop a_new : (-1 -3)
loop b : (-2)
"""

const CUT_G2 = """
piece cut-g2 : genus=1 boundaries=2
marking : reference
boundary BL base=oL : oL 7L 12L
boundary BR base=oR : oR 12R 7R
K arc  m1 7L -> 7R origin=P.m1@+1:+7 : [+1]
K arc  n1 12L -> 12R origin=P.n1@+3:+12 : [-4 +3]
K loop l1 origin=P.l1@-3 : (-3 -1)
K loop l2 origin=P.l2@+4 : (4 -2)
K loop m2 origin=P.m2@+2 : (2)
"""

const DISK = """
piece disk-ps : genus=0 boundaries=1
marking : reference
boundary B base=p : p s
"""

function g2_parent()
    H.LickorishCode(
        2,
        Dict{String, Vector{Int}}(
            "l1" => [-3, -1],
            "l2" => [4, -2],
            "m1" => [1, 7],
            "m2" => [2],
            "n1" => [3, 12, -4],
            "c" => [-7, -12],
        ),
        Dict{String, String}(),
        Dict{String, Vector{String}}(),
        Dict{String, String}(),
        Dict{String, Int}(),
        nothing,
    )
end

function blc_orders_json(values)
    "[" * join(string.(values), ",") * "]"
end

function blc_matrix_json(matrix)
    lines = String[]
    for i in 1:size(matrix, 1)
        cells = String[]
        for j in 1:size(matrix, 2)
            entry = matrix[i, j]
            push!(cells, "[" * string(numerator(entry)) * "," * string(denominator(entry)) * "]")
        end
        push!(lines, "[" * join(cells, ",") * "]")
    end
    "[" * join(lines, ",") * "]"
end

function alpha_names(path::AbstractString)
    names = String[]
    for line in eachline(path)
        m = match(r"^role\s+(\S+)\s*:\s*alpha\s*$", line)
        m === nothing && continue
        push!(names, String(m.captures[1]))
    end
    names
end

@testset "bordered lickorish code" begin
    src = joinpath(@__DIR__, "..", "..", "src", "heegaard", "bordered_lickorish_code.jl")
    src_text = read(src, String)
    @test !occursin("HeegaardCode", src_text)
    @test !occursin("hc_", src_text)
    @test !occursin("role", src_text)

    @testset "examples and euler counts" begin
        ann0 = blc_roundtrip(ANN0)
        @test H.validate_bordered_lc(ann0)
        @test blc_counts(ann0) == (6, 8, 2, 2, 0, true)
        @test blc_counts(ann0; reference=true) == (6, 7, 1, 2, 0, true)

        ann1 = blc_roundtrip(ANN1)
        @test H.validate_bordered_lc(ann1)
        @test blc_counts(ann1) == (7, 10, 3, 2, 0, true)
        @test blc_counts(ann1; reference=true) == (6, 7, 1, 2, 0, true)

        unsigned = replace(ANN1, "[+1]" => "[1]", "[-1]" => "[-1]")
        @test H.read_bordered_lc(unsigned) == ann1

        for (text, counts) in ((VERTEX_L, (9, 15, 7, 1, 1, true)), (VERTEX_R, (9, 15, 7, 1, 1, true)))
            piece = blc_roundtrip(text)
            @test H.validate_bordered_lc(piece)
            @test blc_counts(piece) == counts
            @test blc_counts(piece; reference=true) == (7, 10, 4, 1, 1, true)
        end

        stab_l = blc_roundtrip(STAB_L)
        @test H.validate_bordered_lc(stab_l)
        @test blc_counts(stab_l) == (5, 8, 4, 1, 1, true)

        stab_r = blc_roundtrip(STAB_R)
        @test H.validate_bordered_lc(stab_r)
        @test blc_counts(stab_r) == (7, 12, 4, 1, -1, true)
        @test blc_counts(stab_r; reference=true) == (4, 6, 1, 1, -1, true)

        slide_l = blc_roundtrip(SLIDE_L)
        @test H.validate_bordered_lc(slide_l)
        @test blc_counts(slide_l) == (9, 13, 3, 3, -1, true)
        @test blc_counts(slide_l; reference=true) == (7, 9, 1, 3, -1, true)

        slide_r = blc_roundtrip(SLIDE_R)
        @test H.validate_bordered_lc(slide_r)
        @test blc_counts(slide_r) == (10, 15, 4, 3, -1, true)
        @test blc_counts(slide_r; reference=true) == (7, 9, 1, 3, -1, true)

        g2 = blc_roundtrip(CUT_G2)
        @test H.validate_bordered_lc(g2)
        @test blc_counts(g2) == (10, 16, 4, 2, -2, true)
        diag = H.diagnose_bordered_lc(g2)
        left = ["oL", "7L", "12L"]
        right = ["oR", "12R", "7R"]
        @test diag.boundary_ports[1] in blc_rotations(left)
        @test !(diag.boundary_ports[1] in blc_rotations(reverse(left)))
        @test diag.boundary_ports[2] in blc_rotations(right)
        @test !(diag.boundary_ports[2] in blc_rotations(reverse(right)))

        disk = blc_roundtrip(DISK)
        @test H.validate_bordered_lc(disk)
        @test blc_counts(disk) == (2, 2, 1, 1, 1, true)
    end

    @testset "rejected inputs" begin
        bare = blc_roundtrip(BARE)
        bare_ref = H.diagnose_bordered_lc(bare; reference=true)
        @test (bare_ref.V, bare_ref.E, bare_ref.F, bare_ref.boundary_cycles, bare_ref.chi, bare_ref.connected) ==
            (4, 4, 0, 4, 0, false)
        err = blc_section(() -> H.validate_bordered_lc(bare))
        @test err isa H.BorderedLCError
        @test err.section == "4.3"
        @test occursin("disconnected", err.detail)
        @test occursin("boundary_cycles=4", err.detail)
        @test !occursin("chi=", err.detail)

        bad = blc_roundtrip(SLIDE_BAD)
        @test blc_counts(bad) == (10, 15, 2, 3, -3, true)
        @test blc_counts(bad; reference=true) == (7, 9, 1, 3, -1, true)
        err = blc_section(() -> H.validate_bordered_lc(bad))
        @test err isa H.BorderedLCError
        @test err.section == "4.4"
        @test err.detail == "full chi -3, reference chi -1"
    end

    @testset "genus 2 cut and glue" begin
        parent = g2_parent()
        cut = H.cut_bordered_lc(parent, "c"; parent="P")
        hand = H.read_bordered_lc(CUT_G2)
        @test cut.piece.name == "cut-c"
        @test cut.piece.genus == hand.genus
        @test cut.piece.nboundaries == hand.nboundaries
        @test cut.piece.boundaries == hand.boundaries
        @test cut.piece.k_arcs == hand.k_arcs
        @test cut.piece.k_loops == hand.k_loops
        @test cut.piece.q_arcs == hand.q_arcs
        @test cut.piece.q_loops == hand.q_loops
        @test cut.piece.arcs == hand.arcs
        @test cut.piece.loops == hand.loops
        @test cut.piece.places == hand.places
        @test H.validate_bordered_lc(cut.piece)
        back = H.glue_bordered_lc(cut)
        @test back.words == parent.words
        @test back.genus == parent.genus
    end

    @testset "census 00-19" begin
        dir = joinpath(@__DIR__, "..", "..", "..", "..", "threads", "PJ260922", "data", "261114-census", "hc")
        for i in 0:19
            id = lpad(string(i), 2, '0')
            lc = H.read_lc(joinpath(dir, "$id.lc"))
            names = alpha_names(joinpath(dir, "$id.hc"))
            @test !isempty(names)
            cut = nothing
            used = ""
            for name in names
                try
                    cut = H.cut_bordered_lc(lc, name)
                    used = name
                    break
                catch e
                    e isa H.BorderedLCError && e.reason in (:separating, :not_simple) && continue
                    rethrow()
                end
            end
            @test cut !== nothing
            cut === nothing && continue
            back = H.glue_bordered_lc(cut)
            @test H.read_bordered_lc(H.bordered_lc_text(cut.piece)) == cut.piece
            @test back.words == lc.words
            if back.words != lc.words
                for key in sort!(collect(union(keys(back.words), keys(lc.words))))
                    got = get(back.words, key, Int[])
                    expected = get(lc.words, key, Int[])
                    got == expected && continue
                    error("$id cut along $used: $key got $got expected $expected")
                end
            end
            @test used != ""
        end
    end

    ISOLATED_BAD = """
piece isolated-bad : genus=0 boundaries=1
marking : reference
boundary B base=p : p s
loop z : (+1 +2 -1 -2)
place {z} : F1 ; outside=face(z/+1/out)
"""

    ISOLATED_OK = """
piece isolated-ok : genus=0 boundaries=1
marking : reference
boundary B base=p : p s
loop z : (+1 -1)
place {z} : F1 ; outside=face(z/+1/out)
"""

    QA = """
piece qa : genus=0 boundaries=1
marking : reference
boundary d base=q0 : q0 a0 a1 q1
Q arc q q0 -> q1 : []
arc a a0 -> a1 : []
loop z : ()
place {z} : F1 ; outside=z-
"""

    TRIANGLE = """
genus: 1
l1 : 3 -1
m1 : 2 1
c : -3 -2
"""

    DIGON = """
genus: 1
l1 : -1 2 -3
m1 : 1
c : -2 3
"""

    DIGON_WRAP = """
genus: 1
l1 : -1 2 -3
m1 : 1 5
c : 3 -5 -2
"""

    @testset "per face and visit order" begin
        bad = blc_roundtrip(ISOLATED_BAD)
        @test blc_counts(bad) == (2, 2, 1, 1, 1, true)
        @test blc_counts(bad; reference=true) == (2, 2, 1, 1, 1, true)
        err = blc_section(() -> H.validate_bordered_lc(bad))
        @test err isa H.BorderedLCError
        @test err.section == "4.4"
        @test occursin("isolated", err.detail)
        @test occursin("chi 0", err.detail)

        ok = blc_roundtrip(ISOLATED_OK)
        @test blc_counts(ok) == (2, 2, 1, 1, 1, true)
        @test H.validate_bordered_lc(ok)

        qa = blc_roundtrip(QA)
        faces = H.bordered_lc_faces(qa; drop=Set(["z"]))
        @test [face.id for face in faces] == [1, 2, 3]
        @test [face.ports for face in faces] == [
            ["a0", "a1"],
            ["a0", "q0", "q1", "a1"],
            ["q0", "q1"],
        ]
        @test all(!face.by_label for face in faces)
        @test H.validate_bordered_lc(qa)
        @test H.validate_bordered_lc(H.read_bordered_lc(replace(QA, "F1" => "F2")))

        labeled = H.bordered_lc_faces(H.read_bordered_lc(VERTEX_L); drop=Set{String}())
        @test any(face.by_label for face in labeled)
    end

    function blc_same(a, b)
        a.genus == b.genus && a.words == b.words && a.places == b.places
    end

    function blc_hit_for(hits, witness)
        curves = [witness.names[loc[1]] for loc in witness.locations]
        starts = [loc[2] for loc in witness.locations]
        for hit in hits
            [binding.second for binding in hit.bindings] == curves || continue
            hit.starts == starts || continue
            hit.labels == witness.vertices || continue
            return hit
        end
        nothing
    end

    function blc_mismatch(id, row, witness, got, expected)
        locations = [(witness.names[loc[1]], loc[2]) for loc in witness.locations]
        if got === nothing
            return "$id $row locations=$locations labels=$(witness.vertices) no embedding"
        end
        if got.genus != expected.genus
            return "$id $row locations=$locations genus $(got.genus) vs $(expected.genus)"
        end
        if got.places != expected.places
            return "$id $row locations=$locations places $(got.places) vs $(expected.places)"
        end
        keys = sort!(collect(union(keys(got.words), keys(expected.words))))
        for key in keys
            left = get(got.words, key, Int[])
            right = get(expected.words, key, Int[])
            left == right && continue
            return "$id $row locations=$locations $key pair=$left move=$right"
        end
        "$id $row locations=$locations words differ"
    end

    @testset "pair application" begin
        vertex = H.vertex_bordered_pair()
        @test vertex.left == H.read_bordered_lc(VERTEX_L)
        @test vertex.right == H.read_bordered_lc(VERTEX_R)
        @test vertex.arc_order == ["A", "B", "C"]
        @test vertex.corner_labels == [1, 2, 3]
        @test vertex.pin_below == [1]
        @test !vertex.above_k
        swapped = H.swap_bordered_pair(vertex)
        @test swapped.left == vertex.right
        @test swapped.right == vertex.left
        @test swapped.corner_labels == vertex.corner_labels
        @test swapped.pin_below == vertex.pin_below
        @test swapped.arc_order == vertex.arc_order

        digon_pair = H.digon_bordered_pair()
        @test digon_pair.arc_order == ["C", "K"]
        @test digon_pair.corner_labels == [1, 2]
        @test digon_pair.above_k
        @test isempty(digon_pair.pin_below)
        @test isempty(H.apply_bordered_pair(H._lc_parse_text(DIGON), H.swap_bordered_pair(digon_pair)))

        tri = H._lc_parse_text(TRIANGLE)
        tri_hits = H.movable(tri, "lc.isotopy.vertex.000")
        @test length(tri_hits) == 1
        tri_pair = H.apply_bordered_pair(tri, vertex)
        tri_hit = blc_hit_for(tri_pair, only(tri_hits))
        @test tri_hit !== nothing
        tri_move = H.move(tri, only(tri_hits)).output
        @test tri_hit.output.words == Dict(
            "l1" => [-1, 3],
            "m1" => [1, 2],
            "c" => [-2, -3],
        )
        @test blc_same(tri_hit.output, tri_move)
        back = H.apply_bordered_pair(tri_hit.output, swapped)
        back_hits = H.movable(tri_hit.output, "lc.isotopy.vertex.inverse.000")
        @test length(back_hits) == 1
        back_hit = blc_hit_for(back, only(back_hits))
        @test back_hit !== nothing
        @test blc_same(back_hit.output, tri)

        dig = H._lc_parse_text(DIGON)
        dig_hits = H.movable(dig, "lc.isotopy.digon.remove.00")
        @test length(dig_hits) == 1
        dig_pair = H.apply_bordered_pair(dig, digon_pair)
        dig_hit = blc_hit_for(dig_pair, only(dig_hits))
        @test dig_hit !== nothing
        dig_move = H.move(dig, only(dig_hits)).output
        @test dig_hit.output.words["l1"] == [-1]
        @test dig_hit.output.words["m1"] == [1]
        @test dig_hit.output.words["c"] == Int[]
        @test blc_same(dig_hit.output, dig_move)

        wrapped = H._lc_parse_text(DIGON_WRAP)
        wrap_hits = H.movable(wrapped, "lc.isotopy.digon.remove.00")
        @test length(wrap_hits) == 1
        wrap_pair = H.apply_bordered_pair(wrapped, digon_pair)
        wrap_hit = blc_hit_for(wrap_pair, only(wrap_hits))
        @test wrap_hit !== nothing
        wrap_move = H.move(wrapped, only(wrap_hits)).output
        @test wrap_hit.output.words["l1"] == [-1]
        @test wrap_hit.output.words["m1"] == [1, 3]
        @test wrap_hit.output.words["c"] == [-3]
        @test blc_same(wrap_hit.output, wrap_move)
    end

    @testset "census pair regression" begin
        # Agreement with move is a regression. A wrong map and its inverse can agree too.
        dir = joinpath(@__DIR__, "..", "..", "..", "..", "threads", "PJ260922", "data", "261114-census", "hc")
        vertex = H.vertex_bordered_pair()
        digon_pair = H.digon_bordered_pair()
        compared = 0
        mismatches = String[]
        for i in 0:19
            id = lpad(string(i), 2, '0')
            lc = H.read_lc(joinpath(dir, "$id.lc"))
            for (row, pair) in (
                ("lc.isotopy.vertex.000", vertex),
                ("lc.isotopy.digon.remove.00", digon_pair),
            )
                witnesses = H.movable(lc, row)
                hits = H.apply_bordered_pair(lc, pair)
                for witness in witnesses
                    compared += 1
                    hit = blc_hit_for(hits, witness)
                    moved = H.move(lc, witness).output
                    if hit === nothing || !blc_same(hit.output, moved)
                        push!(mismatches, blc_mismatch(id, row, witness,
                            hit === nothing ? nothing : hit.output, moved))
                    end
                end
            end
        end
        println("CENSUS compared=$compared mismatches=$(length(mismatches))")
        for line in mismatches
            println(line)
        end
        @test isempty(mismatches)
    end

    function blc_arc_entries(piece, name)
        for arc in Iterators.flatten((piece.k_arcs, piece.q_arcs, piece.arcs))
            arc.name == name || continue
            return [lab.sign * lab.n for lab in arc.word]
        end
        Int[]
    end

    function blc_digon_eps(push_word, cross_word)
        push_word == [-1, 2] && cross_word == [-2, 1] && return (1, 1)
        push_word == [1, -2] && cross_word == [-1, 2] && return (1, -1)
        push_word == [-1, 2] && cross_word == [1, -2] && return (-1, 1)
        push_word == [1, -2] && cross_word == [2, -1] && return (-1, -1)
        nothing
    end

    function blc_host_kind(lc, name)
        H._is_k(lc, name) ? :K : :input
    end

    function blc_band_from_route(route)
        point(p) = Dict{String,Any}("dart" => collect(p.dart), "u" => p.u)
        Dict{String,Any}(
            "start" => point(route.start),
            "finish" => point(route.finish),
            "crossings" => Any[point(p) for p in route.crossings],
        )
    end

    function blc_word_line(id, got, expected)
        got === nothing && return "$id no output"
        got.genus == expected.genus || return "$id genus $(got.genus) vs $(expected.genus)"
        got.places == expected.places || return "$id places $(got.places) vs $(expected.places)"
        names = sort!(collect(union(keys(got.words), keys(expected.words))))
        for name in names
            left = get(got.words, name, Int[])
            right = get(expected.words, name, Int[])
            left == right || return "$id $name pair=$left move=$right"
        end
        "$id words differ"
    end

    function blc_rehouse(lc, hc)
        H.HeegaardCode(lc, copy(hc.roles), copy(hc.disks), copy(hc.points), copy(hc.weights), nothing)
    end

    function blc_pack_form(hc)
        form = H.hc_linking_form(hc; boundary_side=:alpha)
        pres = nothing
        pres_error = ""
        try
            pres = H.first_homology(hc)
        catch err
            pres_error = sprint(showerror, err)
        end
        (
            free_rank=form.free_rank,
            torsion=copy(form.torsion),
            matrix=copy(form.matrix),
            pres=pres,
            pres_error=pres_error,
        )
    end

    function blc_h1(hc)
        try
            h = H.hc_first_homology(hc)
            (ok=true, free_rank=h.free_rank, torsion=Int.(h.torsion), error="")
        catch err
            (ok=false, free_rank=0, torsion=Int[], error=sprint(showerror, err))
        end
    end

    function blc_link_torsion(pack)
        sort!(Int[Int(x) for x in pack.torsion])
    end

    @testset "argument families" begin
        rows = Dict(
            (1, 1) => (points=["Cm", "Cp", "Dm", "Dp"], c=[-1, 2], d=[-2, 1]),
            (1, -1) => (points=["Cm", "Cp", "Dp", "Dm"], c=[1, -2], d=[-1, 2]),
            (-1, 1) => (points=["Cm", "Dm", "Dp", "Cp"], c=[-1, 2], d=[1, -2]),
            (-1, -1) => (points=["Cm", "Dp", "Dm", "Cp"], c=[1, -2], d=[2, -1]),
        )
        for (signs, row) in rows
            pair = H.digon_bordered_family(signs[1], signs[2]; create=false, c_kind=:input, d_kind=:K)
            @test pair.arc_order == ["C", "D"]
            @test pair.corner_labels == [1, 2]
            @test pair.above_k
            @test isempty(pair.pin_below)
            @test pair.left.boundaries[1].points == row.points
            @test pair.right.boundaries[1].points == row.points
            @test blc_arc_entries(pair.left, "C") == row.c
            @test blc_arc_entries(pair.left, "D") == row.d
            @test isempty(blc_arc_entries(pair.right, "C"))
            @test isempty(blc_arc_entries(pair.right, "D"))
            created = H.digon_bordered_family(signs[1], signs[2]; create=true, c_kind=:input, d_kind=:K)
            @test created == H.swap_bordered_pair(pair)
            qpair = H.digon_bordered_family(signs[1], signs[2]; create=false, c_kind=:input, d_kind=:Q)
            @test blc_arc_entries(qpair.left, "D") == row.d
            # Reference curve present: the §4.4 seed still runs, and χ stays 1.
            @test blc_counts(pair.left)[5] == 1
            @test blc_counts(pair.left; reference=true)[5] == 1
            @test blc_counts(qpair.left)[5] == 1
            @test blc_counts(qpair.left; reference=true)[5] == 1
        end
        # Both arcs are input, so there is no reference curve. All four sign rows
        # used to stop in §4.4 with "reference face 1 is empty".
        for signs in ((1, 1), (1, -1), (-1, 1), (-1, -1))
            both_input = H.digon_bordered_family(signs[1], signs[2]; create=false, c_kind=:input, d_kind=:input)
            @test blc_counts(both_input.left)[5] == 1
            @test blc_counts(both_input.left; reference=true)[5] == 1
            @test blc_counts(both_input.right)[5] == 1
            @test blc_counts(both_input.right; reference=true)[5] == 1
        end
        @test occursin("chi = _blc_region_chi(full, faces) - carrier + 1", src_text)
        removal = H.digon_bordered_family(-1, 1; create=false, c_kind=:input, d_kind=:K)
        @test blc_counts(removal.left) == (6, 10, 5, 1, 1, true)
        @test blc_counts(removal.right) == (4, 6, 3, 1, 1, true)
        sign_err = blc_section(() -> H.digon_bordered_family(2, 1))
        @test sign_err isa H.BorderedLCError
        @test sign_err.section == "7.7.1"
        @test sign_err.reason == :sign
        kind_err = blc_section(() -> H.digon_bordered_family(1, 1; c_kind=:other))
        @test kind_err isa H.BorderedLCError
        @test kind_err.section == "7.7.1"
        @test kind_err.reason == :kind

        slide = H.slide_bordered_family(1, 1)
        @test slide.left == H.read_bordered_lc(SLIDE_L)
        @test slide.right == H.read_bordered_lc(SLIDE_R)
        @test slide.arc_order == ["q", "r"]
        @test H.swap_bordered_pair(slide).left == slide.right
        @test H.swap_bordered_pair(slide).right == slide.left
        for ea in (1, -1), eb in (1, -1)
            built = H.slide_bordered_family(ea, eb)
            @test blc_counts(built.left)[5] == -1
            @test blc_counts(built.right)[5] == -1
            @test blc_arc_entries(built.right, "r")[1] == ea * 3
            @test blc_arc_entries(built.right, "r")[2] == eb * 2
        end
        slide_sign = blc_section(() -> H.slide_bordered_family(0, 1))
        @test slide_sign isa H.BorderedLCError
        @test slide_sign.section == "7.7.2"
        @test slide_sign.reason == :sign

        dig = H._lc_parse_text(DIGON)
        dig_remove = only(H.movable(dig, "lc.isotopy.digon.remove.00"))
        dig_hits = H.apply_bordered_family(dig, removal)
        dig_hit = blc_hit_for(dig_hits, dig_remove)
        @test dig_hit !== nothing
        @test blc_same(dig_hit.output, H.move(dig, dig_remove).output)

        host = H._lc_parse_text(DIGON)
        made = H.movable(host, "lc.isotopy.digon.create.ab")
        @test !isempty(made)
        restored = 0
        for witness in made
            eps = blc_digon_eps(witness.after_frags[1], witness.after_frags[2])
            @test eps !== nothing
            eps === nothing && continue
            c_name = witness.names[witness.locations[1][1]]
            d_name = witness.names[witness.locations[2][1]]
            pair = H.digon_bordered_family(eps[1], eps[2]; create=true,
                c_kind=blc_host_kind(host, c_name), d_kind=H._is_k(host, d_name) ? :K : :Q)
            hit = blc_hit_for(H.apply_bordered_family(host, pair), witness)
            @test hit !== nothing
            hit === nothing && continue
            @test blc_same(hit.output, H.move(host, witness).output)
            backs = H.apply_bordered_family(hit.output, H.swap_bordered_pair(pair))
            if any(item -> blc_same(item.output, host), backs)
                restored += 1
            end
            @test any(item -> blc_same(item.output, host), backs)
        end
        println("HAND digon create roundtrip=$restored")
        @test restored == length(made)

        p28_path = joinpath(@__DIR__, "..", "..", "..", "..", "threads", "PJ260922", "code", "lc-fixtures", "p28-after.lc")
        p28 = H.read_lc(p28_path)
        stops = Dict{Symbol,H.BorderedLCError}()
        function keep_stop(reason, err)
            @test err isa H.BorderedLCError
            err isa H.BorderedLCError || return
            @test err.section == "7.7.2"
            @test err.reason == reason
            stops[reason] = err
        end
        keep_stop(:curve, blc_section(() -> H.slide_bordered_family(
            p28, "missing", "alpha2", nothing; family_a="alpha", family_b="alpha")))
        keep_stop(:distinct, blc_section(() -> H.slide_bordered_family(
            p28, "alpha1", "alpha1", nothing; family_a="alpha", family_b="alpha")))
        keep_stop(:family, blc_section(() -> H.slide_bordered_family(
            p28, "alpha1", "alpha2", nothing; family_a="alpha", family_b="beta")))
        empty_lc = H._lc_parse_text("genus: 1\nl1 : +1\nm1 : -1\nb : +2\nc : -2\n")
        empty_lc.words["a"] = Int[]
        keep_stop(:empty, blc_section(() -> H.slide_bordered_family(
            empty_lc, "a", "b", nothing; family_a="alpha", family_b="alpha")))
        kink = H._lc_parse_text("genus: 1\nl1 : +2\nm1 : -2\na : +1 -1\nb : +3\nc : -3\n")
        keep_stop(:simple, blc_section(() -> H.slide_bordered_family(
            kink, "a", "b", nothing; family_a="alpha", family_b="alpha")))
        crossed = H._lc_parse_text("genus: 1\nl1 : +1\nm1 : -1\na : +2\nb : -2 +3\nc : -3\n")
        keep_stop(:intersect, blc_section(() -> H.slide_bordered_family(
            crossed, "a", "b", nothing; family_a="alpha", family_b="alpha")))
        keep_stop(:band, blc_section(() -> H.slide_bordered_family(
            p28, "alpha1", "alpha2", "not-a-band"; family_a="alpha", family_b="alpha")))
        absent = Dict{String,Any}(
            "start" => Dict{String,Any}("dart" => Any[999, "+", "out"], "u" => "1/2"),
            "finish" => Dict{String,Any}("dart" => Any[998, "-", "in"], "u" => "1/2"),
        )
        keep_stop(:band, blc_section(() -> H.slide_bordered_family(
            p28, "alpha1", "alpha2", absent; family_a="alpha", family_b="alpha")))

        faces, right = H._lc_faces(H._lc_nonempty_words(p28.words))
        occ = H._lc_rotation_occ(p28.words)
        function curve_darts(name)
            found = Tuple{String,Vector{Any}}[]
            for face in faces
                for corner in face.corners
                    occ[corner[1]][corner[2]][1] == name || continue
                    push!(found, (face.id, Any[corner[1], corner[2], corner[3]]))
                end
            end
            found
        end
        rectangle = nothing
        for (face_a, dart_a) in curve_darts("alpha1"), (face_b, dart_b) in curve_darts("alpha2")
            rectangle === nothing || break
            face_a == face_b && continue
            band = Dict{String,Any}(
                "start" => Dict{String,Any}("dart" => dart_a, "u" => "1/2"),
                "finish" => Dict{String,Any}("dart" => dart_b, "u" => "1/2"),
                "crossings" => Any[],
            )
            err = blc_section(() -> H.slide_bordered_family(
                p28, "alpha1", "alpha2", band; family_a="alpha", family_b="alpha"))
            err isa H.BorderedLCError && err.reason == :rectangle || continue
            rectangle = err
        end
        keep_stop(:rectangle, rectangle)

        winding = nothing
        for band0 in H._move_face_bands(p28, "alpha1", "alpha2")
            face_id = right[H._move_corner(band0["start"]["dart"])]
            face = only(item for item in faces if item.id == face_id)
            for corner in face.corners
                crossed_name = String(occ[corner[1]][corner[2]][1])
                crossed_name in ("alpha1", "alpha2") && continue
                trial = deepcopy(band0)
                trial["crossings"] = Any[
                    Dict{String,Any}("dart" => collect(corner), "u" => "1/3"),
                    Dict{String,Any}("dart" => collect(H._lc_tau(p28.words, occ, corner)), "u" => "2/3"),
                ]
                try
                    H._move_certify_route(p28, trial; avoid=String["alpha1", "alpha2"])
                catch
                    continue
                end
                winding = (trial, crossed_name)
                break
            end
            winding === nothing || break
        end
        @test winding !== nothing
        if winding !== nothing
            trial, crossed_name = winding
            keep_stop(:winding, blc_section(() -> H.slide_bordered_family(
                p28, "alpha1", "alpha2", trial; family_a="alpha", family_b="alpha",
                avoid=String["alpha1", "alpha2", crossed_name])))
        end
        @test sort!(collect(keys(stops))) == [
            :band, :curve, :distinct, :empty, :family, :intersect, :rectangle, :simple, :winding,
        ]

        p28_hc = H.HeegaardCode(p28, Dict(
            "alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta,
        ), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
        p28_bands = H._move_face_bands(p28, "alpha1", "alpha2")
        @test !isempty(p28_bands)
        hand_slide = 0
        hand_rows = NamedTuple[]
        for band in p28_bands
            result = H.apply_slide_bordered(p28, "alpha1", "alpha2", band;
                family_a="alpha", family_b="alpha", avoid=String["alpha1", "alpha2"])
            witness = only(H.movable(p28_hc, "hc.slide.pants"; curves=["alpha1", "alpha2"], band=band))
            moved = H.move(p28_hc, witness).output
            @test blc_same(result.output, moved.lc)
            back = H.apply_slide_inverse(result.output, result.inverse)
            @test blc_same(back, p28)
            geo = H.apply_reversed_slide(result.output, p28, "alpha1", "alpha2", band;
                avoid=String["alpha1", "alpha2"])
            @test blc_same(geo, back)
            housed = blc_rehouse(result.output, p28_hc)
            before = blc_pack_form(p28_hc)
            after = blc_pack_form(housed)
            @test before.free_rank == after.free_rank
            @test before.torsion == after.torsion
            hb = blc_h1(p28_hc)
            ha = blc_h1(housed)
            @test hb.ok
            @test ha.ok
            @test hb.free_rank == ha.free_rank
            @test hb.torsion == ha.torsion
            @test hb.free_rank == before.free_rank
            @test sort(hb.torsion) == blc_link_torsion(before)
            if before.pres !== nothing && after.pres !== nothing
                @test before.pres.free_rank == after.pres.free_rank
                @test before.pres.torsion == after.pres.torsion
            end
            hand_slide += 1
            push!(hand_rows, (tag="p28-$hand_slide-before", pack=before))
            push!(hand_rows, (tag="p28-$hand_slide-after", pack=after))
            q_word = blc_arc_entries(result.pair.left, "q")
            r_word = blc_arc_entries(result.pair.left, "r")
            sd = band["start"]["dart"]
            fd = band["finish"]["dart"]
            println("HAND slide $hand_slide q=$q_word r=$r_word start=$(sd[3]) finish=$(fd[3]) eps=$(H._move_dir(sd)),$(H._move_dir(fd))")
        end
        println("HAND slide compared=$hand_slide")
        @test hand_slide > 0

        global BLC_HAND_ROWS = hand_rows   # used by bordered_lickorish_code_corpus.jl
    end

    # Expectations were written in inverse-expected.json before any library run
    # (2026-11-18 10:28:59). Signs stay part of the cyclic match.
    @testset "geometric slide inverse" begin
        function closed_lc(words)
            H.LickorishCode(
                1,
                Dict(name => copy(word) for (name, word) in words),
                Dict{String,String}(),
                Dict{String,Vector{String}}(),
                Dict{String,String}(),
                Dict{String,Int}(),
                nothing,
            )
        end
        function band_of(sl, ss, sw, fl, fs, fw)
            Dict{String,Any}(
                "start" => Dict{String,Any}("dart" => Any[sl, ss, sw], "u" => "1/2"),
                "finish" => Dict{String,Any}("dart" => Any[fl, fs, fw], "u" => "1/2"),
                "crossings" => Any[],
            )
        end
        function cyclic_same(got::Vector{Int}, expect::Vector{Int})
            length(got) == length(expect) || return false
            n = length(got)
            n == 0 && return true
            for shift in 0:(n - 1)
                all(j -> got[mod1(j + shift, n)] == expect[j], 1:n) && return true
            end
            false
        end
        function same_words(got, expect)
            Set(keys(got)) == Set(keys(expect)) &&
                all(cyclic_same(Int[got[name]...], expect[name]) for name in keys(expect))
        end
        function surface_counts(lc)
            v, e, f = H._lc_surface_stats(lc.words)
            (v, e, f, v - e + f)
        end
        function digon_faces(lc)
            faces, _ = H._lc_faces(lc.words)
            count(face -> length(face.corners) == 2, faces)
        end
        function check_round(source, band, back_counts, start_dart, start_edge, finish_dart, finish_edge,
                attachment, forward_edge, minted)
            lc = closed_lc(source)
            slid = H.apply_slide_bordered(lc, "m1", "m2", band; family_a="slide", family_b="slide")
            back = H.apply_geometric_slide_inverse(slid.output, lc, "m1", "m2", band)
            @test same_words(back.diagram.words, source)
            @test surface_counts(back.diagram) == back_counts
            rec = back.record
            @test rec["inverse_start"][1] == start_dart[1]
            @test String(rec["inverse_start"][2]) == start_dart[2]
            @test String(rec["inverse_start"][3]) == start_dart[3]
            @test Int(rec["inverse_start_edge"]) == start_edge
            @test rec["inverse_finish"][1] == finish_dart[1]
            @test String(rec["inverse_finish"][2]) == finish_dart[2]
            @test String(rec["inverse_finish"][3]) == finish_dart[3]
            @test Int(rec["inverse_finish_edge"]) == finish_edge
            @test rec["attachment_index"] == attachment
            @test rec["forward_edge"] == forward_edge
            @test Int[rec["minted"]...] == minted
            @test rec["replaced"] == "m1"
            @test rec["stationary"] == "m2"
            back.diagram
        end

        slide_source = Dict("l1" => [-1, -2], "m1" => [1], "m2" => [2])
        slide_back = check_round(
            slide_source,
            band_of(1, "+", "out", 2, "+", "in"),
            (2, 4, 2, 0),
            [1, "+", "in"], 2,
            [2, "+", "out"], 1,
            2, 1, [3],
        )
        @test digon_faces(slide_back) == 0

        digon_source = Dict("l1" => [-1, 3, -2], "m1" => [1, -3], "m2" => [2])
        digon_back = check_round(
            digon_source,
            band_of(1, "+", "in", 2, "+", "in"),
            (3, 6, 3, 0),
            [1, "+", "in"], 3,
            [2, "+", "out"], 1,
            3, 2, [4],
        )
        @test digon_faces(closed_lc(digon_source)) == 2
        @test digon_faces(digon_back) == 2

        hand = closed_lc(Dict("l1" => [-1, 3, -4, -2], "m1" => [1, -3, 4], "m2" => [2]))
        flipped = band_of(1, "+", "in", 2, "+", "out")
        route = H._blc_slide_route(hand, "m1", "m2", flipped, String["m1", "m2"])
        @test route.start.edge == 3
        @test route.finish.edge == 1
        unflipped = blc_section(() -> H._blc_slide_route(
            hand, "m1", "m2", band_of(1, "+", "in", 2, "+", "in"), String["m1", "m2"]))
        @test unflipped isa H.BorderedLCError
        @test unflipped.section == "7.7.2"
        @test unflipped.reason == :rectangle
    end
end

