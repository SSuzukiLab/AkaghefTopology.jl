@testset "Lescop's decorated Heegaard formula" begin
    function diagram(alpha, beta, intersections)
        code = H.heegaard_code(alpha, beta, intersections)
        H.certify_heegaard(H.decorate(code; cellular=true))
    end

    # Lescop, https://arxiv.org/abs/1402.2261, Figures 1 and 3,
    # Examples 1.2 and 1.3. The complementary region containing w1
    # has oriented boundary turn -1 after omitting the outer boundary.
    # The two oriented beta arcs turn 1/2 and 3/2 in Fig. 3, so this
    # region follows the c-to-d beta edge forward and the d-to-c edge
    # backward, with two clockwise alpha half-turns and four corners.
    D1 = diagram([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2])
    @test H.lescop_intersection_terms(D1, [1]) == (ell2=0, sl=0)
    @test H.lescop_theta(D1, [1], 2) == (ell2=0, sl=0, e=-1//4, theta=1//4)
    @test H.lescop_theta(D1, [1], 1) == (ell2=0, sl=0, e=1//4, theta=-1//4)
    @test_throws ArgumentError H.lescop_theta(D1, [1], 3)
    @test_throws ArgumentError H.lescop_intersection_terms(D1, [1,2])

    # Fig. 3 D2: crossings c,d,e,f correspond to IDs 1,2,3,4.
    # The beta turns are 1/2, 0, -1/2 on beta1 and 0 on beta2.
    # With these turns the cycle beginning with alpha1:c->d has total
    # turn -1, identifying the region at infinity (w2).
    D2 = diagram([[1,2], [3,4]], [[1,2,4], [3]],
                 [H.Intersection(1,1,1,1), H.Intersection(2,1,1,1),
                  H.Intersection(3,2,2,1), H.Intersection(4,2,1,1)])
    @test H.lescop_theta(D2, [1,3], 1) == (ell2=0, sl=0, e=1//4, theta=-1//4)
end

@testset "Lescop rotation on the Koda lens-family ribbon" begin
    for p in 2:12
        x, y = p + 1, p + 2
        signs = [fill(-1, p); 1; -1]
        code = H.heegaard_code([[x; collect(p:-1:1); y]], [[y; collect(1:p); x]],
                               [H.Intersection(i,1,1,s) for (i,s) in enumerate(signs)])
        diagram = H.certify_heegaard(H.decorate(code; cellular=true))
        # The joining bigon follows alpha y->x and beta x->y.
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        dart = only(d for d in graph.darts if
            graph.edges[d.edge].family == :alpha && d.origin == y && d.target == x && isodd(d.id))
        region = only(r for (r, cycle) in enumerate(cycles) if dart.id in cycle)
        @test length(cycles[region]) == 2
        # Koda (2008), Fig. 25 and Theorem 5.14. Compare the full cyclic
        # word and vertex signs to the registered all-positive lens family.
        ograph = H.to_odata(diagram; joining_disk=H.joining_disk(diagram, region)).ograph
        expected = [-p:-1; 1:p]
        @test length(ograph.gauss) == 1
        @test any(circshift(only(ograph.gauss), k) == expected for k in eachindex(expected))
        @test ograph.orientation == ones(Int, p)
        @test H.lescop_theta(diagram, [x], region).theta == (p + 1)//(2p)
        @test H.lescop_theta(diagram, [x], region).e == 0
        @test H.lescop_theta(diagram, [y], region).theta == (1 - p)//(2p)
        @test H.lescop_theta(diagram, [y], region).e == 1//2
    end
end

@testset "Lescop Θ from a resolved Heegaard code" begin
    # The exterior region is the complement region whose boundary contains
    # the named arcs. It is not chosen by comparing Θ, and an LC face id
    # is not used as a region index.
    function resolved(text, roles)
        H.resolve_hc(H._lc_parse_text(text), roles)
    end

    function code_of(view)
        view.diagram.code.code
    end

    function visits_of(view, family, curve_label)
        code = code_of(view)
        curves = family == :alpha ? code.alpha : code.beta
        curve = only(c for c in curves if c.label == curve_label)
        names = Dict(id => string(label) for (id, label) in view.label_of_id)
        [names[id] for id in curve.visits]
    end

    function sigma_of(view, label)
        only(p.sign for p in code_of(view).intersections if p.source_label == string(label))
    end

    function ids_for(view, labels)
        wanted = Dict(label => id for (id, label) in view.label_of_id)
        [wanted[label] for label in labels]
    end

    function boundary_signatures(diagram)
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        blocks = diagram.code.regions.boundary_cycle_blocks
        [Set((graph.edges[graph.darts[id].edge].family,
              graph.edges[graph.darts[id].edge].curve_id,
              graph.darts[id].origin, graph.darts[id].target, isodd(id))
             for c in block for id in cycles[c]) for block in blocks]
    end

    function forward_arc(view, family, curve_label, from_label, to_label)
        code = code_of(view)
        curves = family == :alpha ? code.alpha : code.beta
        curve = only(c for c in curves if c.label == curve_label)
        labels = Dict(label => id for (id, label) in view.label_of_id)
        origin, target = labels[from_label], labels[to_label]
        (dart, edge) -> edge.family == family && edge.curve_id == curve.id &&
            dart.origin == origin && dart.target == target && isodd(dart.id)
    end

    function region_with(view, preds)
        graph = H.ribbon_graph(view.diagram)
        cycles = H.boundary_cycles(graph)
        blocks = view.diagram.code.regions.boundary_cycle_blocks
        found = Int[]
        for (r, block) in enumerate(blocks)
            ok = true
            for pred in preds
                hit = false
                for c in block, id in cycles[c]
                    dart = graph.darts[id]
                    if pred(dart, graph.edges[dart.edge])
                        hit = true
                        break
                    end
                end
                if !hit
                    ok = false
                    break
                end
            end
            ok && push!(found, r)
        end
        only(found)
    end

    # Les14 D1. m1 carries the α role, so this K curve is kept. l1 is
    # roleless and is deleted. The forgotten visits are the genus-one
    # diagram already compared with Examples 1.2–1.3. Label 2 is the
    # matching crossing. The α word carries a minus there, and σ is +1.
    D1 = """
    genus: 1
    l1 : -1 -4
    m1 : 1 -2 -3
    b : 4 2 3
    """
    hc1 = resolved(D1, Dict("m1" => :alpha, "b" => :beta))
    view1 = H.hc_to_cellular(hc1)
    @test visits_of(view1, :alpha, "m1") == ["2", "3"]
    @test visits_of(view1, :beta, "b") == ["2", "3"]
    @test sigma_of(view1, 2) == 1
    @test sigma_of(view1, 3) == 1
    @test -2 in hc1.lc.words["m1"]
    @test sort(parse.(Int, [p.source_label for p in code_of(view1).intersections])) == [2, 3]
    old1 = let
        code = H.heegaard_code([[1, 2]], [[1, 2]], [H.Intersection(i, 1, 1, 1) for i in 1:2])
        H.certify_heegaard(H.decorate(code; cellular=true))
    end
    w1 = only(r for (r, sig) in enumerate(boundary_signatures(view1.diagram))
              if sig == boundary_signatures(old1)[2])
    @test view1.face_region["F3"] == w1
    @test view1.face_region["F4"] == w1
    @test view1.face_region["F3"] != 3
    ids1 = ids_for(view1, [2])
    @test H.lescop_theta(view1.diagram, ids1, view1.face_region["F4"]) ==
          (ell2=0, sl=0, e=-1//4, theta=1//4)
    @test H.hc_lescop_theta(hc1, [2], "F3") == (ell2=0, sl=0, e=-1//4, theta=1//4)
    @test_throws ArgumentError H.hc_lescop_theta(hc1, [2], 3)

    # Les14 D2. β2 is parallel to ℓ2: m2 carries −5 and n1 carries +6.
    # β1 then meets m1 twice and m2 once, in the order c, d, f.
    # Labels: 13 = c, 14 = d, 12 = f, 5 = e. l1, l2, and n1 are roleless.
    # m1 and m2 stay because they carry the α role.
    D2 = """
    genus: 2
    l1 : -1 10 11 -3
    l2 : -2 4
    m1 : 9 -11 -13 -14 1
    m2 : 8 7 -10 -12 -5 2
    n1 : 6 3 -9 -7 -8 -4
    b1 : 13 14 12
    b2 : 5 -6
    """
    hc2 = resolved(D2, Dict("m1" => :alpha, "m2" => :alpha, "b1" => :beta, "b2" => :beta))
    view2 = H.hc_to_cellular(hc2)
    @test visits_of(view2, :alpha, "m1") == ["13", "14"]
    @test visits_of(view2, :alpha, "m2") == ["12", "5"]
    @test visits_of(view2, :beta, "b1") == ["13", "14", "12"]
    @test visits_of(view2, :beta, "b2") == ["5"]
    @test all(label -> sigma_of(view2, label) == 1, (13, 14, 12, 5))
    @test -13 in hc2.lc.words["m1"]
    w2 = region_with(view2, [
        forward_arc(view2, :alpha, "m1", 13, 14),
        forward_arc(view2, :beta, "b1", 14, 12),
        forward_arc(view2, :beta, "b1", 12, 13),
    ])
    @test view2.face_region["F1"] == w2
    @test view2.face_region["F11"] == w2
    @test view2.face_region["F11"] != 11
    ids2 = ids_for(view2, [13, 5])
    @test H.lescop_theta(view2.diagram, ids2, view2.face_region["F11"]) ==
          (ell2=0, sl=0, e=1//4, theta=-1//4)
    @test H.hc_lescop_theta(hc2, [13, 5], "F1") == (ell2=0, sl=0, e=1//4, theta=-1//4)
    @test_throws ArgumentError H.hc_lescop_theta(hc2, [13, 5], 1)

    # Les12a Figure 6. D3 is the D2 word plus g = 15 and h = 16.
    # σ(g) = +1 because m1 carries −15. σ(h) = −1 because m1 carries +16.
    # α1 reads c, d, g, h. β2 reads h, g, e, a rotation of e, h, g.
    # The exterior region is the unique region that still contains the
    # D2 arcs c→d, d→f, and f→c. Those arcs do not use the old β2 loop.
    D3 = """
    genus: 2
    l1 : -1 10 11 -3
    l2 : -2 4
    m1 : 9 -11 -13 -14 -15 16 1
    m2 : 8 7 -10 -12 -5 2
    n1 : 6 3 -9 -7 -8 -4
    b1 : 13 14 12
    b2 : -16 15 5 -6
    """
    hc3 = resolved(D3, Dict("m1" => :alpha, "m2" => :alpha, "b1" => :beta, "b2" => :beta))
    view3 = H.hc_to_cellular(hc3)
    @test visits_of(view3, :alpha, "m1") == ["13", "14", "15", "16"]
    @test visits_of(view3, :alpha, "m2") == ["12", "5"]
    @test visits_of(view3, :beta, "b1") == ["13", "14", "12"]
    @test visits_of(view3, :beta, "b2") == ["16", "15", "5"]
    @test all(label -> sigma_of(view3, label) == 1, (13, 14, 12, 5, 15))
    @test sigma_of(view3, 16) == -1
    @test -15 in hc3.lc.words["m1"]
    @test 16 in hc3.lc.words["m1"]
    w3 = region_with(view3, [
        forward_arc(view3, :alpha, "m1", 13, 14),
        forward_arc(view3, :beta, "b1", 14, 12),
        forward_arc(view3, :beta, "b1", 12, 13),
    ])
    @test view3.face_region["F1"] == w3
    @test view3.face_region["F12"] == w3
    @test view3.face_region["F12"] != 12
    ids3 = ids_for(view3, [13, 5])
    expected3 = (ell2=1//4, sl=0, e=1//2, theta=-1//4)
    @test H.lescop_theta(view3.diagram, ids3, view3.face_region["F12"]) == expected3
    @test H.hc_lescop_theta(hc3, [13, 5], "F1") == expected3
    other3 = first(face for (face, region) in view3.face_region if region != w3)
    other_theta3 = H.lescop_theta(view3.diagram, ids3, view3.face_region[other3])
    @test other_theta3.ell2 == 1//4
    @test other_theta3.sl == 0
    @test other_theta3.e != 1//2
    @test_throws ArgumentError H.hc_lescop_theta(hc3, [13, 5], 1)
end

using LinearAlgebra

@testset "Lescop variation, rotation sign, and singular or multi-boundary regions" begin
    # The right-hand side is the intersection matrix and the half-weight arcs.
    # It does not call lescop_theta, _lescop_rotation, or _lescop_arcs.
    function l2_rat(text)
        s = string(text)
        if occursin('/', s)
            a, b = split(s, '/')
            return parse(BigInt, a) // parse(BigInt, b)
        end
        parse(BigInt, s) // 1
    end

    function arc_on(visits, base, target)
        weights = Dict{Int,Rational{BigInt}}()
        edges = Tuple{Int,Int}[]
        target == base && return weights, edges
        start = findfirst(==(base), visits)
        start === nothing && error("base $base missing")
        n = length(visits)
        prev = base
        for step in 1:n
            nxt = visits[mod1(start + step, n)]
            push!(edges, (prev, nxt))
            if nxt == target
                weights[base] = 1//2
                weights[target] = get(weights, target, 0//1) + 1//2
                return weights, edges
            end
            weights[nxt] = get(weights, nxt, 0//1) + 1
            prev = nxt
        end
        error("target $target not reached")
    end

    function diagram_pack(diagram, matching)
        code = diagram.code.code
        points = Dict(p.id => p for p in code.intersections)
        alpha_index = Dict(c.id => i for (i, c) in enumerate(code.alpha))
        beta_index = Dict(c.id => j for (j, c) in enumerate(code.beta))
        J = inv(Rational{BigInt}.(H.intersection_matrix(diagram)))
        u = Dict{Int,Rational{BigInt}}()
        k = Dict{Int,Rational{BigInt}}()
        match = Set(matching)
        for p in code.intersections
            uc = J[beta_index[p.beta], alpha_index[p.alpha]] * p.sign
            u[p.id] = uc
            k[p.id] = (p.id in match ? 1 : 0) - uc
        end
        beta_match = Dict{Int,Int}()
        alpha_match = Dict{Int,Int}()
        for id in matching
            p = points[id]
            beta_match[p.beta] = id
            alpha_match[p.alpha] = id
        end
        (; code, points, alpha_index, beta_index, J, u, k, beta_match, alpha_match)
    end

    # Lemma 5.2 sums over every crossing. A crossing on another beta still
    # contributes through ⟨path, β_s⟩. all_crossings=false keeps the sum that
    # produced the frozen census-L(3,1) value, for comparison with that file.
    function exterior_change(diagram, matching, from_id, to_id; all_crossings::Bool=true)
        pack = diagram_pack(diagram, matching)
        g = length(pack.code.alpha)
        P = Dict(id => 0//1 for id in keys(pack.points))
        Pbeta = zeros(Rational{BigInt}, g)
        A = [Dict(id => 0//1 for id in keys(pack.points)) for _ in 1:g]
        owner = only(c for c in pack.code.beta if c.id == pack.points[from_id].beta)
        Pbeta[pack.beta_index[owner.id]] = 1
        for id in keys(pack.points)
            curve = only(c for c in pack.code.beta if c.id == pack.points[id].beta)
            on_owner = curve.id == owner.id
            all_crossings || on_owner || continue
            weights, edges = arc_on(curve.visits, pack.beta_match[curve.id], id)
            on_owner && (from_id, to_id) in edges && (P[id] = 1)
            for (x, w) in weights
                i = pack.alpha_index[pack.points[x].alpha]
                A[i][id] += w * pack.points[x].sign
            end
        end
        lk = sum(pack.k[id] * P[id] for id in keys(P); init=0//1)
        for j in 1:g, i in 1:g
            acc = sum(pack.k[id] * A[i][id] for id in keys(P); init=0//1)
            lk -= pack.J[j, i] * acc * Pbeta[j]
        end
        delta = 0//1
        for id in keys(pack.u)
            inner = -P[id]
            for s in 1:g, r in 1:g
                inner += pack.J[s, r] * A[r][id] * Pbeta[s]
            end
            delta += 2 * pack.u[id] * inner
        end
        (; lk, delta, two_lk=2lk)
    end

    function arc_product(pack, alpha_weights, beta_weights)
        acc = 0//1
        for (id, wa) in alpha_weights
            wb = get(beta_weights, id, 0//1)
            iszero(wb) && continue
            acc += wa * wb * pack.points[id].sign
        end
        acc
    end

    function arc_against_curve(pack, weights, family, curve_id)
        acc = 0//1
        for (id, w) in weights
            p = pack.points[id]
            hit = family == :beta ? p.beta == curve_id : p.alpha == curve_id
            hit || continue
            acc += w * p.sign
        end
        acc
    end

    function ell_tilde(pack, x, y)
        ax = only(c for c in pack.code.alpha if c.id == pack.points[x].alpha)
        by = only(c for c in pack.code.beta if c.id == pack.points[y].beta)
        aw, _ = arc_on(ax.visits, pack.alpha_match[ax.id], x)
        bw, _ = arc_on(by.visits, pack.beta_match[by.id], y)
        raw = arc_product(pack, aw, bw)
        corr = 0//1
        for (j, bcurve) in enumerate(pack.code.beta), (i, acurve) in enumerate(pack.code.alpha)
            corr += pack.J[j, i] * arc_against_curve(pack, aw, :beta, bcurve.id) *
                    arc_against_curve(pack, bw, :alpha, acurve.id)
        end
        raw - corr
    end

    function parallel_lk(diagram, matching, coeffs)
        pack = diagram_pack(diagram, matching)
        acc = 0//1
        for x in keys(pack.points), y in keys(pack.points)
            iszero(coeffs[x]) && continue
            iszero(coeffs[y]) && continue
            acc += coeffs[x] * coeffs[y] * ell_tilde(pack, x, y)
        end
        acc
    end

    function link_coeffs(diagram, plus_id, minus_id)
        Dict(p.id => (p.id == plus_id ? 1//1 : p.id == minus_id ? -1//1 : 0//1)
             for p in diagram.code.code.intersections)
    end

    function side_of_edge(diagram, from_id, to_id)
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        blocks = diagram.code.regions.boundary_cycle_blocks
        owner = Dict{Int,Int}()
        for (r, block) in enumerate(blocks), c in block, id in cycles[c]
            owner[id] = r
        end
        odd = even = nothing
        for e in graph.edges
            e.family == :beta || continue
            e.from == from_id && e.to == to_id || continue
            odd = owner[2e.id - 1]
            even = owner[2e.id]
        end
        (; odd, even)
    end

    function resolved(text, roles)
        H.resolve_hc(H._lc_parse_text(text), roles)
    end

    function id_of(view, label)
        only(id for (id, lab) in view.label_of_id if lab == label)
    end

    function positive_lens(p)
        code = H.heegaard_code([collect(1:p)], [collect(1:p)],
                               [H.Intersection(i, 1, 1, 1) for i in 1:p])
        H.certify_heegaard(H.decorate(code; cellular=true))
    end

    function tangent_degree(pts)
        acc = 0.0
        ang(a, b) = atan(b[2] - a[2], b[1] - a[1])
        prev = ang(pts[1], pts[2])
        for i in 2:length(pts)-1
            nxt = ang(pts[i], pts[i + 1])
            d = nxt - prev
            while d > π
                d -= 2π
            end
            while d < -π
                d += 2π
            end
            acc += d
            prev = nxt
        end
        acc / (2π)
    end

    function cubic(p0, c1, c2, p3, n)
        [(let t = k / n, s = 1 - t
              (s^3*p0[1] + 3s^2*t*c1[1] + 3s*t^2*c2[1] + t^3*p3[1],
               s^3*p0[2] + 3s^2*t*c1[2] + 3s*t^2*c2[2] + t^3*p3[2])
          end) for k in 0:n]
    end

    function chain(segments; n=800)
        pts = Tuple{Float64,Float64}[]
        for seg in segments
            chunk = cubic(seg..., n)
            isempty(pts) || pop!(pts)
            append!(pts, chunk)
        end
        pts
    end

    function check_exterior(diagram, matching, from_id, to_id, spec)
        got = exterior_change(diagram, matching, from_id, to_id)
        @test got.lk == l2_rat(spec["lk"])
        @test got.delta == l2_rat(spec["delta_theta"])
        @test got.two_lk == got.delta
        sides = side_of_edge(diagram, from_id, to_id)
        @test sides.odd !== nothing && sides.even !== nothing && sides.odd != sides.even
        right = H.lescop_theta(diagram, matching, sides.even)
        left = H.lescop_theta(diagram, matching, sides.odd)
        @test right.ell2 == left.ell2
        @test right.sl == left.sl
        @test right.theta - left.theta == l2_rat(spec["delta_theta"])
        sides
    end

    repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    expected = JSON.parsefile(joinpath(repo, "threads", "PJ260922", "data",
                                        "261117-lescop", "l2-expected.json"))
    corrected = JSON.parsefile(joinpath(repo, "threads", "PJ260922", "data",
                                         "261117-lescop", "l3-rhs-corrected.json"))
    spec_of(group, id) = only(row for row in expected[group]["cases"] if row["id"] == id)
    rotation = expected["rotation"]

    # Clockwise shores of one alpha hole. Both have planar degree -1/2.
    # East through south to west, then west through north to east: both clockwise.
    cw_east = [(cos(-π * k / 1024), sin(-π * k / 1024)) for k in 0:1024]
    cw_west = [(cos(-π - π * k / 1024), sin(-π - π * k / 1024)) for k in 0:1024]
    planar = l2_rat(rotation["backward_alpha_planar_degree"])
    @test abs(tangent_degree(cw_east) - Float64(planar)) < 1e-3
    @test abs(tangent_degree(cw_west) - Float64(planar)) < 1e-3

    # Les14 Figure 3, left panel, tex 426–430. The last leg writes one control.
    short = chain([
        ((2.4, 0.9), (1.9, 0.9), (1.4, 1.4), (1.0, 1.4)),
        ((1.0, 1.4), (0.6, 1.4), (0.3, 1.25), (0.3, 1.05)),
        ((0.3, 1.05), (0.3, 0.9), (0.3, 0.9), (0.6, 0.9)),
    ])
    long = chain([
        ((3.2, 0.9), (3.5, 0.9), (3.5, 1.1), (3.5, 1.3)),
        ((3.5, 1.3), (3.5, 1.45), (3.0, 1.55), (1.0, 1.55)),
        ((1.0, 1.55), (0.6, 1.55), (0.1, 1.4), (0.1, 0.5)),
        ((0.1, 0.5), (0.1, 0.3), (0.5, 0.2), (1.0, 0.2)),
        ((1.0, 0.2), (1.5, 0.2), (1.9, 0.9), (1.4, 0.9)),
    ])
    d1rot = rotation["D1_w1_matching_c"]
    @test abs(tangent_degree(short) - Float64(l2_rat(d1rot["d_rot_arc_c_to_d"]))) < 1e-3
    @test abs(tangent_degree(short) + tangent_degree(long) - Float64(l2_rat(d1rot["d_rot_beta1"]))) < 1e-3

    D1 = positive_lens(2)
    data1 = H._lescop_data(D1, [1])
    e1 = H._lescop_rotation(D1, data1, H._lescop_intersections(data1),
                            side_of_edge(D1, 1, 2).odd)
    @test e1 == l2_rat(d1rot["e"])

    hc1 = resolved("""
    genus: 1
    l1 : -1 -4
    m1 : 1 -2 -3
    b : 4 2 3
    """, Dict("m1" => :alpha, "b" => :beta))
    view1 = H.hc_to_cellular(hc1)
    c1, d1 = id_of(view1, 2), id_of(view1, 3)
    check_exterior(view1.diagram, [c1], c1, d1, spec_of("exterior", "D1-edge-c-to-d"))

    hc2 = resolved("""
    genus: 2
    l1 : -1 10 11 -3
    l2 : -2 4
    m1 : 9 -11 -13 -14 1
    m2 : 8 7 -10 -12 -5 2
    n1 : 6 3 -9 -7 -8 -4
    b1 : 13 14 12
    b2 : 5 -6
    """, Dict("m1" => :alpha, "m2" => :alpha, "b1" => :beta, "b2" => :beta))
    view2 = H.hc_to_cellular(hc2)
    c2, d2, e2 = id_of(view2, 13), id_of(view2, 14), id_of(view2, 5)
    check_exterior(view2.diagram, [c2, e2], c2, d2, spec_of("exterior", "D2-edge-c-to-d"))

    hc3 = resolved("""
    genus: 2
    l1 : -1 10 11 -3
    l2 : -2 4
    m1 : 9 -11 -13 -14 -15 16 1
    m2 : 8 7 -10 -12 -5 2
    n1 : 6 3 -9 -7 -8 -4
    b1 : 13 14 12
    b2 : -16 15 5 -6
    """, Dict("m1" => :alpha, "m2" => :alpha, "b1" => :beta, "b2" => :beta))
    view3 = H.hc_to_cellular(hc3)
    check_exterior(view3.diagram, [id_of(view3, 13), id_of(view3, 5)],
                   id_of(view3, 13), id_of(view3, 14),
                   spec_of("exterior", "D3-edge-c-to-d"))

    lens = spec_of("exterior", "standard-positive-L(p,1)")
    for sample in lens["samples"]
        p = Int(sample["p"])
        check_exterior(positive_lens(p), [1], 1, 2, sample)
    end

    for row in expected["exterior"]["cases"]
        haskey(row, "file") || continue
        hc = H.read_hc(joinpath(repo, split(row["file"], "/")...))
        view = H.hc_to_cellular(hc)
        labels = row["matching_labels"]
        edge = row["edge_labels"]
        matching = [id_of(view, Int(label)) for label in labels]
        if row["id"] == "census-L(3,1)"
            # l2-expected records the sum over the crossed beta only. Lemma 5.2
            # sums over every crossing. The corrected file is that sum.
            from_id = id_of(view, Int(edge[1]))
            to_id = id_of(view, Int(edge[2]))
            partial = exterior_change(view.diagram, matching, from_id, to_id;
                                      all_crossings=false)
            @test partial.lk == l2_rat(row["lk"])
            @test partial.delta == l2_rat(row["delta_theta"])
            named = corrected["named_edge"]
            @test Int.(named["labels"]) == Int.(edge)
            full_named = exterior_change(view.diagram, matching, from_id, to_id)
            @test full_named.lk == l2_rat(named["lk"])
            @test full_named.delta == l2_rat(named["delta_theta"])
            @test full_named.two_lk == full_named.delta
            @test full_named.delta - partial.delta ==
                  l2_rat(corrected["omitted_on_every_beta_edge"])
            for curve in view.diagram.code.code.beta
                visits = curve.visits
                for k in eachindex(visits)
                    fr = visits[k]
                    to = visits[mod1(k + 1, length(visits))]
                    labels = [view.label_of_id[fr], view.label_of_id[to]]
                    spec = only(item for item in corrected["edges"]
                                if Int.(item["labels"]) == labels)
                    partial_edge = exterior_change(view.diagram, matching, fr, to;
                                                   all_crossings=false)
                    change = exterior_change(view.diagram, matching, fr, to)
                    @test partial_edge.lk == l2_rat(spec["truncated_lk"])
                    @test partial_edge.delta == l2_rat(spec["truncated_delta_theta"])
                    @test change.lk == l2_rat(spec["lk"])
                    @test change.delta == l2_rat(spec["delta_theta"])
                    @test change.two_lk == change.delta
                    @test change.delta - partial_edge.delta ==
                          l2_rat(corrected["omitted_on_every_beta_edge"])
                    sides = side_of_edge(view.diagram, fr, to)
                    right = H.lescop_theta(view.diagram, matching, sides.even)
                    left = H.lescop_theta(view.diagram, matching, sides.odd)
                    @test right.ell2 - left.ell2 == 0
                    @test right.sl - left.sl == 0
                    @test right.theta - left.theta == -(right.e - left.e)
                    @test right.theta - left.theta == change.delta
                end
            end
        else
            check_exterior(view.diagram, matching, id_of(view, Int(edge[1])),
                           id_of(view, Int(edge[2])), row)
        end
    end

    # Proposition 5.5 at the old exterior region. The parallel term is the arc sum.
    match1 = spec_of("matching", "D1-matching-c-to-d")
    @test parallel_lk(view1.diagram, [c1], link_coeffs(view1.diagram, d1, c1)) ==
          l2_rat(match1["lk_delta_parallel"])
    w1 = side_of_edge(view1.diagram, c1, d1).odd
    theta_c = H.lescop_theta(view1.diagram, [c1], w1)
    theta_d = H.lescop_theta(view1.diagram, [d1], w1)
    @test theta_d.theta - theta_c.theta == l2_rat(match1["delta_theta"])
    @test l2_rat(match1["sum_d_e"]) - l2_rat(match1["lk_delta_parallel"]) ==
          l2_rat(match1["delta_theta"])

    match2 = spec_of("matching", "D2-matching-ce-to-de")
    @test parallel_lk(view2.diagram, [c2, e2], link_coeffs(view2.diagram, d2, c2)) ==
          l2_rat(match2["lk_delta_parallel"])
    graph2 = H.ribbon_graph(view2.diagram)
    cycles2 = H.boundary_cycles(graph2)
    blocks2 = view2.diagram.code.regions.boundary_cycle_blocks
    function has_forward(block, family, curve_label, from_label, to_label)
        curve = only(c for c in (family == :alpha ? view2.diagram.code.code.alpha :
                                  view2.diagram.code.code.beta) if c.label == curve_label)
        origin, target = id_of(view2, from_label), id_of(view2, to_label)
        any(block) do c
            any(cycles2[c]) do dart_id
                dart = graph2.darts[dart_id]
                edge = graph2.edges[dart.edge]
                edge.family == family && edge.curve_id == curve.id &&
                    dart.origin == origin && dart.target == target && isodd(dart_id)
            end
        end
    end
    w2 = only(r for (r, block) in enumerate(blocks2)
              if has_forward(block, :alpha, "m1", 13, 14) &&
                 has_forward(block, :beta, "b1", 14, 12) &&
                 has_forward(block, :beta, "b1", 12, 13))
    data2 = H._lescop_data(view2.diagram, [c2, e2])
    @test H._lescop_rotation(view2.diagram, data2, H._lescop_intersections(data2), w2) ==
          l2_rat(rotation["D2_w2_matching_c_e"]["e"])
    theta_ce = H.lescop_theta(view2.diagram, [c2, e2], w2)
    theta_de = H.lescop_theta(view2.diagram, [d2, e2], w2)
    @test theta_de.theta - theta_ce.theta == l2_rat(match2["delta_theta"])

    ann = expected["multi_boundary"]
    hc_ann = resolved(ann["hc"], Dict(name => Symbol(role) for (name, role) in ann["roles"]))
    view_ann = H.hc_to_cellular(hc_ann)
    blocks_ann = view_ann.diagram.code.regions.boundary_cycle_blocks
    cycles_ann = H.boundary_cycles(H.ribbon_graph(view_ann.diagram))
    @test any(length(block) > 1 for block in blocks_ann)
    @test sum(length, blocks_ann) == length(cycles_ann)
    @test H.intersection_matrix(view_ann.diagram) == [1 0; 0 1]
    @test length(unique(values(view_ann.face_region))) == 1
    face_ann = first(face for (face, region) in view_ann.face_region if region == 1)
    got_ann = H.hc_lescop_theta(hc_ann, Int.(ann["matching_labels"]), face_ann)
    want_ann = Tuple(l2_rat(x) for x in ann["theta"])
    @test (got_ann.ell2, got_ann.sl, got_ann.e, got_ann.theta) == want_ann
    loops = H.heegaard_code([[1], [2]], [[1], [2]],
                            [H.Intersection(1, 1, 1, 1), H.Intersection(2, 2, 2, 1)])
    loop_cycles = H.boundary_cycles(H.ribbon_graph(loops))
    loop_diagram = H.certify_heegaard(H.decorate(loops, [collect(eachindex(loop_cycles))]))
    direct = H.lescop_theta(loop_diagram, [1, 2], 1)
    @test (direct.ell2, direct.sl, direct.e, direct.theta) == want_ann
    @test direct == got_ann

    singular = expected["singular"]
    hc_sing = H.read_hc(joinpath(repo, split(singular["file"], "/")...))
    view_sing = H.hc_to_cellular(hc_sing)
    @test only(H.intersection_matrix(view_sing.diagram)) == 0
    @test_throws ArgumentError H.hc_lescop_theta(hc_sing, [first(values(view_sing.label_of_id))],
                                                 first(keys(view_sing.face_region)))
end

@testset "Lescop p1 difference across one beta edge" begin
    # Theorem 3.4 and Proposition 5.1, read from l6-p1-expected.json.
    # The linking side does not call lescop_theta. Corollary 5.3 is the
    # separate comparison 4 (Θ(right) − Θ(left)) = delta_p1.
    function l6_rat(text)
        s = string(text)
        if occursin('/', s)
            a, b = split(s, '/')
            return parse(BigInt, a) // parse(BigInt, b)
        end
        parse(BigInt, s) // 1
    end

    function id_of(view, label)
        only(id for (id, lab) in view.label_of_id if lab == label)
    end

    function side_of_edge(diagram, from_id, to_id)
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        blocks = diagram.code.regions.boundary_cycle_blocks
        owner = Dict{Int,Int}()
        for (r, block) in enumerate(blocks), c in block, id in cycles[c]
            owner[id] = r
        end
        odd = even = nothing
        for e in graph.edges
            e.family == :beta || continue
            e.from == from_id && e.to == to_id || continue
            odd = owner[2e.id - 1]
            even = owner[2e.id]
        end
        (; odd, even)
    end

    function check_edge(hc, view, matching_labels, from_label, to_label, lk, delta_p1)
        got = H.hc_lescop_p1_difference(hc, matching_labels, from_label, to_label)
        @test got.lk == l6_rat(lk)
        @test got.delta_p1 == l6_rat(delta_p1)
        @test got.delta_p1 == 8 * got.lk
        from_id = id_of(view, Int(from_label))
        to_id = id_of(view, Int(to_label))
        matching = [id_of(view, Int(label)) for label in matching_labels]
        direct = H.lescop_p1_difference(view.diagram, matching, from_id, to_id)
        @test direct == got
        sides = side_of_edge(view.diagram, from_id, to_id)
        @test sides.odd !== nothing && sides.even !== nothing && sides.odd != sides.even
        right = H.lescop_theta(view.diagram, matching, sides.even)
        left = H.lescop_theta(view.diagram, matching, sides.odd)
        @test right.ell2 == left.ell2
        @test right.sl == left.sl
        @test 4 * (right.theta - left.theta) == got.delta_p1
    end

    repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    spec = JSON.parsefile(joinpath(repo, "threads", "PJ260922", "data",
                                    "261117-lescop", "l6-p1-expected.json"))
    @test spec["frozen_before_evaluator"] == true
    @test spec["absolute_p1_of_each_combing"]["assigned"] == false
    @test spec["adopted_formula"]["does_not_read_theta"] == true

    for row in spec["edges"]
        if row["diagram"] == "inline"
            roles = Dict(name => Symbol(role) for (name, role) in row["roles"])
            hc = H.resolve_hc(H._lc_parse_text(row["lc"]), roles)
            view = H.hc_to_cellular(hc)
            edge = row["edge_labels"]
            check_edge(hc, view, row["matching_labels"], edge[1], edge[2],
                       row["lk"], row["delta_p1"])
            # The stored number equals 2 lk. Corollary 5.3 calls that ΔΘ.
            @test 2 * l6_rat(row["lk"]) ==
                  l6_rat(row["four_delta_theta_predicted_by_the_corollary"])
        else
            hc = H.read_hc(joinpath(repo, split(row["diagram"], "/")...))
            view = H.hc_to_cellular(hc)
            for edge in row["edges"]
                labels = edge["labels"]
                check_edge(hc, view, row["matching_labels"], labels[1], labels[2],
                           edge["lk"], edge["delta_p1"])
            end
        end
    end
end

@testset "Spin^c labels from matchings and exterior regions" begin
    # OS04a subsec:Disks and lemma:VarySpinC. The hand values are the
    # genus-one count and the coboundary residue. Theta and p1 are not read.
    repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    spec = JSON.parsefile(joinpath(repo, "threads", "PJ260922", "data",
                                    "261117-lescop", "l8-spinc-expected.json"))
    @test spec["frozen_before_evaluator"] == true
    @test spec["uses_theta"] == false
    @test spec["uses_p1"] == false

    function region_count(diagram)
        length(diagram.code.regions.boundary_cycle_blocks)
    end

    function homology_order(diagram)
        homology = H.first_homology(diagram)
        homology.free_rank == 0 || throw(ArgumentError("H1 has positive free rank"))
        isempty(homology.torsion) ? 1 : prod(homology.torsion)
    end

    function visit_step(diagram)
        code = diagram.code.code
        alpha = only(code.alpha)
        beta = only(code.beta)
        modulus = length(beta.visits)
        beta_index = Dict(id => position - 1 for (position, id) in enumerate(beta.visits))
        alpha_position = Dict(beta_index[id] => position - 1
                              for (position, id) in enumerate(alpha.visits))
        mod(alpha_position[1] - alpha_position[0], modulus)
    end

    function modular_inverse(step, modulus)
        for value in 0:modulus-1
            mod(step * value, modulus) == 1 && return value
        end
        throw(ArgumentError("the beta step has no inverse"))
    end

    function hand_class(beta_visits, first_id, second_id)
        modulus = length(beta_visits)
        first = findfirst(==(first_id), beta_visits) - 1
        second = findfirst(==(second_id), beta_visits) - 1
        Int[mod(-(second - first), modulus)]
    end

    function owners_of(diagram)
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        blocks = diagram.code.regions.boundary_cycle_blocks
        owner = Dict{Int,Int}()
        for (region, block) in enumerate(blocks), cycle in block, dart in cycles[cycle]
            owner[dart] = region
        end
        graph, owner
    end

    function labels_of(diagram)
        code = diagram.code.code
        matchings = H._lescop_matchings(code)
        found = Dict{Vector{Int},Int}()
        for matching in matchings, region in 1:region_count(diagram)
            label = H.lescop_spinc_label(diagram, matching, region)
            found[label] = get(found, label, 0) + 1
        end
        matchings, found
    end

    agreements = 0
    disagreements = 0
    changed_edges = Dict{String,Int}()
    unchanged_edges = Dict{String,Int}()

    for row in spec["diagrams"]
        row["id"] in ("L(5,2)", "L(7,2)") || continue
        roles = Dict(name => Symbol(role) for (name, role) in row["roles"])
        hc = H.resolve_hc(H._lc_parse_text(row["lc"]), roles)
        view = H.hc_to_cellular(hc)
        diagram = view.diagram
        code = diagram.code.code
        order = homology_order(diagram)
        @test order == row["h1_order"]
        @test order == abs(det(BigInt.(H.intersection_matrix(diagram))))
        @test all(point.sign == 1 for point in code.intersections)
        @test visit_step(diagram) == row["q"]
        @test modular_inverse(row["q"], order) == row["q_inverse_mod_h1"]
        residue = row["q_inverse_mod_h1"]
        beta = only(code.beta)
        matchings, found = labels_of(diagram)
        @test length(matchings) == order
        @test length(found) == order
        println("SPIN_C $(row["id"]) classes=$(length(found)) combings=$(sum(values(found))) regions=$(region_count(diagram))")
        pairs = 0
        for first in matchings, second in matchings
            first == second && continue
            pairs += 1
            got = H.lescop_epsilon(diagram, first, second)
            hand = hand_class(beta.visits, only(first), only(second))
            if got == hand
                agreements += 1
            else
                disagreements += 1
                println("EPSILON_MISMATCH $(row["id"]) $first $second got=$got hand=$hand")
            end
            left = H.lescop_spinc_label(diagram, first, 1)
            right = H.lescop_spinc_label(diagram, second, 1)
            mod(right[1] - left[1], order) == hand[1] ||
                println("LABEL_MISMATCH $(row["id"]) $first $second left=$left right=$right hand=$hand")
            @test mod(right[1] - left[1], order) == hand[1]
            @test H.lescop_epsilon(diagram, first, first) == Int[0]
        end
        @test pairs == row["ordered_distinct_pairs"]
        graph, owner = owners_of(diagram)
        changed = 0
        unchanged = 0
        for edge in graph.edges
            edge.family == :beta || continue
            left_region = owner[2edge.id - 1]
            right_region = owner[2edge.id]
            @test left_region != right_region
            label_left = H.lescop_spinc_label(diagram, matchings[1], left_region)
            label_right = H.lescop_spinc_label(diagram, matchings[1], right_region)
            @test mod(label_right[1] - label_left[1], order) == residue
            label_right == label_left ? (unchanged += 1) : (changed += 1)
        end
        @test changed == row["beta_edges"]["changed"]
        @test unchanged == row["beta_edges"]["unchanged"]
        changed_edges[row["id"]] = changed
        unchanged_edges[row["id"]] = unchanged
        face = first(keys(view.face_region))
        direct = H.lescop_spinc_label(diagram, matchings[1], view.face_region[face])
        @test H.hc_lescop_spinc_label(hc, [view.label_of_id[only(matchings[1])]], face) == direct
    end

    census = only(row for row in spec["diagrams"] if row["id"] == "census-L(3,1)")
    hc = H.read_hc(joinpath(repo, split(census["diagram"], "/")...))
    view = H.hc_to_cellular(hc)
    diagram = view.diagram
    order = homology_order(diagram)
    @test order == census["h1_order"]
    @test order == abs(det(BigInt.(H.intersection_matrix(diagram))))
    matchings, found = labels_of(diagram)
    @test length(found) == order
    census_pairs = count(pair -> pair[1] != pair[2], Iterators.product(matchings, matchings))
    matrix = H.intersection_matrix(diagram)
    for first in matchings, second in matchings
        epsilon = H.lescop_epsilon(diagram, first, second)
        reverse = H.lescop_epsilon(diagram, second, first)
        @test H._lescop_in_lattice(epsilon .+ reverse, matrix)
        left = H.lescop_spinc_label(diagram, first, 1)
        right = H.lescop_spinc_label(diagram, second, 1)
        @test H._lescop_in_lattice(right .- left .- epsilon, matrix)
    end
    graph, owner = owners_of(diagram)
    changed = 0
    unchanged = 0
    reference = first(sort(matchings))
    for edge in graph.edges
        edge.family == :beta || continue
        left_region = owner[2edge.id - 1]
        right_region = owner[2edge.id]
        label_left = H.lescop_spinc_label(diagram, reference, left_region)
        label_right = H.lescop_spinc_label(diagram, reference, right_region)
        label_left == label_right ? (unchanged += 1) : (changed += 1)
    end
    changed_edges["census-L(3,1)"] = changed
    unchanged_edges["census-L(3,1)"] = unchanged
    println("SPIN_C census-L(3,1) classes=$(length(found)) combings=$(sum(values(found))) regions=$(region_count(diagram)) ordered_pairs=$census_pairs changed=$changed unchanged=$unchanged")
    id_of = Dict(label => id for (id, label) in view.label_of_id)
    illumination_matching = [id_of[5], id_of[6]]
    for edge in spec["illumination_against_turn6"]["census_edges_from_l6"]
        from_label, to_label = edge["labels"]
        from_id = id_of[from_label]
        to_id = id_of[to_label]
        graph, owner = owners_of(diagram)
        ribbon = only(item for item in graph.edges
                      if item.family == :beta && item.from == from_id && item.to == to_id)
        left_region = owner[2ribbon.id - 1]
        right_region = owner[2ribbon.id]
        label_left = H.lescop_spinc_label(diagram, illumination_matching, left_region)
        label_right = H.lescop_spinc_label(diagram, illumination_matching, right_region)
        println("ILLUMINATION census $(from_label)->$(to_label) quarter=$(edge["quarter"]) left=$label_left right=$label_right changed=$(label_left != label_right)")
    end

    @test fieldtype(H.CombingFraming, :spin_c) == Nothing
    @test spec["koda"]["expected_correspondence"] == "missing"
    println("SPIN_C agreements=$agreements disagreements=$disagreements changed=$changed_edges unchanged=$unchanged_edges")
    @test disagreements == 0
    @test agreements == 62
end

@testset "S3 residual along the lemma anchor" begin
    # Les14 lemlambst and propaddconnsum. This is internal to Lescop.
    repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    spec = JSON.parsefile(joinpath(repo, "threads", "PJ260922", "data",
                                    "261117-lescop", "l8-s3-expected.json"))
    @test spec["frozen_before_evaluator"] == true
    @test spec["not_an_external_foundation"] == true
    @test spec["does_not_assign_lens_space_p1"] == true

    function positive_lens(p)
        code = H.heegaard_code([collect(1:p)], [collect(1:p)],
                               [H.Intersection(i, 1, 1, 1) for i in 1:p])
        H.certify_heegaard(H.decorate(code; cellular=true))
    end

    function region_count(diagram)
        length(diagram.code.regions.boundary_cycle_blocks)
    end

    function once_pairs(matrix)
        genus = size(matrix, 1)
        size(matrix, 2) == genus || return false
        for row in 1:genus
            count(!iszero, matrix[row, :]) == 1 || return false
            abs(sum(matrix[row, :])) == 1 || return false
        end
        for column in 1:genus
            count(!iszero, matrix[:, column]) == 1 || return false
        end
        true
    end

    reached = 0
    genus_one = positive_lens(1)
    @test region_count(genus_one) == 1
    @test only(point.sign for point in genus_one.code.code.intersections) == 1
    theta = H.lescop_theta(genus_one, [1], 1)
    residual = 4 * theta.theta
    println("S3 genus1 ell2=$(theta.ell2) sl=$(theta.sl) e=$(theta.e) theta=$(theta.theta) p1=0 lambda=0 R=$residual")
    @test theta.ell2 == 0
    @test theta.sl == 0
    @test theta.e == 0
    @test theta.theta == 0
    @test residual == 0
    reached += 1

    stabilized, certificate = H.stabilize(genus_one)
    @test certificate.attachment_data["intersection_id"] == 2
    @test region_count(stabilized) == 1
    matching = [1, certificate.attachment_data["intersection_id"]]
    theta = H.lescop_theta(stabilized, matching, 1)
    residual = 4 * theta.theta
    println("S3 stabilized ell2=$(theta.ell2) sl=$(theta.sl) e=$(theta.e) theta=$(theta.theta) p1=0 lambda=0 R=$residual")
    @test theta.ell2 == 0
    @test theta.sl == 0
    @test theta.e == 0
    @test theta.theta == 0
    @test residual == 0
    reached += 1
    graph = H.ribbon_graph(stabilized)
    cycles = H.boundary_cycles(graph)
    blocks = stabilized.code.regions.boundary_cycle_blocks
    owner = Dict{Int,Int}()
    for (region, block) in enumerate(blocks), cycle in block, dart in cycles[cycle]
        owner[dart] = region
    end
    for edge in graph.edges
        edge.family == :beta || continue
        difference = H.lescop_p1_difference(stabilized, matching, edge.from, edge.to)
        println("S3 stabilized edge $(edge.from)->$(edge.to) same_region=$(owner[2edge.id - 1] == owner[2edge.id]) delta_p1=$(difference.delta_p1)")
        @test owner[2edge.id - 1] == owner[2edge.id]
    end

    library = spec["genus_two_library_diagram"]
    hc = H.read_hc(joinpath(repo, split(library["path"], "/")...))
    view = H.hc_to_cellular(hc)
    diagram = view.diagram
    matrix = H.intersection_matrix(diagram)
    regions = region_count(diagram)
    println("S3 g2-S3-pairs matrix=$matrix regions=$regions")
    if once_pairs(matrix) && regions == 1
        oriented = diagram
        signs = copy(matrix)
        for column in 1:size(signs, 2)
            entry = sum(signs[:, column])
            entry < 0 || continue
            oriented, _ = H.reverse_curve_orientation(oriented,
                (:beta, oriented.code.code.beta[column].id))
        end
        oriented_matrix = H.intersection_matrix(oriented)
        oriented_regions = region_count(oriented)
        println("S3 g2-S3-pairs oriented matrix=$oriented_matrix regions=$oriented_regions")
        if once_pairs(oriented_matrix) && all(>=(0), oriented_matrix) && oriented_regions == 1
            ids = [point.id for point in oriented.code.code.intersections]
            theta = H.lescop_theta(oriented, ids, 1)
            residual = 4 * theta.theta
            println("S3 g2 ell2=$(theta.ell2) sl=$(theta.sl) e=$(theta.e) theta=$(theta.theta) p1=0 lambda=0 R=$residual")
            @test theta.ell2 == 0
            @test theta.sl == 0
            @test theta.e == 0
            @test theta.theta == 0
            @test residual == 0
            reached += 1
        else
            println("S3 g2-S3-pairs unassigned after orientation")
        end
    else
        println("S3 g2-S3-pairs unassigned before orientation")
    end
    println("S3 R_zero=$reached")
    @test reached >= 1
end
