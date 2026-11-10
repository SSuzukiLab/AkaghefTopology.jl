module KodaConversionTests
using Test, Topology
const H = Topology.Heegaard
const ROOT = abspath(joinpath(@__DIR__, "..", "..", "..", ".."))

@testset "koda local tables load" begin
    tabs = H.load_koda_local_tables()
    @test tabs.catalog["schema"] == "akaghef.heegaard.koda-local/2"
    @test tabs.catalog["implementation_ready"] == true
    model = tabs.vertex["models"]["1"]
    @test model["pd_ports_ccw"] == ["-i", "+o", "-o", "+i"]
    pending = H.koda_pending_table_fields(model)
    @test pending == String[]
end

@testset "odata_to_spine forward decode" begin
    G = H.OData([[1, -1]], [1], Dict{String, Any}())
    spine = H.odata_to_spine(G)
    @test spine.gauss == [[1, -1]]
    @test spine.vertex_signs[1] == 1
    punctured = H.spine_to_heegaard_koda(spine).value
    @test punctured isa H.PuncturedHeegaardDiagram
end

@testset "§12.1 p28-after stops without joining disk" begin
    lc_path = joinpath(ROOT, "threads", "PJ260922", "code", "lc-fixtures", "p28-after.lc")
    isfile(lc_path) || @test_skip "fixture missing"
    lc = H.read_lc(lc_path)
    roles = Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta)
    hc = H.HeegaardLayer(lc, roles, Dict(), Dict(), Dict(), nothing)
    err = nothing
    try
        H.hc_to_odata(hc)
    catch e
        err = e
    end
    @test err isa H.KodaInputError || err isa H.KodaImplementationError
    if err isa H.KodaInputError
        @test err.stage == :hc_to_odata
    end
end

@testset "HeegaardAbalone forward o-data matches reference" begin
    include(joinpath(@__DIR__, "..", "..", "examples", "HeegaardAbalone.jl"))
    ex = HeegaardAbaloneExample.fixture()
    b = H.to_odata(H.make_punctured(ex.diagram; joining_disk=ex.disk))
    G = H.OData([[1, -1]], [1], Dict{String, Any}())
    @test H.compare_odata_words(G, b.ograph)["equal_up_to_metadata"]
    ph = H.odata_to_heegaard_koda(G).value
    @test ph isa H.PuncturedHeegaardDiagram
    G2 = H.spine_to_odata(H.heegaard_to_spine(ph).value)
    @test H.compare_odata_words(G, G2)["equal_up_to_metadata"]
    hb = H.first_homology(ex.diagram)
    hr = H.first_homology(ph.diagram)
    @test hb.free_rank == hr.free_rank
    @test hb.torsion == hr.torsion
end

function _closed(G)
    link = Topology.VirtualLink()
    Topology.set_data!(link; rgauss=deepcopy(G.gauss), orientation=copy(G.orientation))
    Topology.is_closed(link)
end

function _odata(diagram, region)
    jd = H.joining_disk(diagram, region)
    ph = H.make_punctured_combinatorial(diagram, jd)
    H.spine_to_odata(H.heegaard_to_spine(ph).value)
end

@testset "short words do not invent a direction" begin
    @test H._koda_lc_arc_direction([1], 1, 1) === nothing
    @test H._koda_lc_arc_direction([1, 2], 1, 2) === nothing
    @test H._koda_lc_arc_direction([1, 2, 3], 1, 2) == 1
    @test H._koda_lc_arc_direction([1, 3, 2], 1, 2) == -1
    err = try
        H._koda_strand_oriented_cycles([[1, 2]], [[1]], Dict((:alpha, 1) => 1))
        nothing
    catch e
        e
    end
    @test err isa H.KodaInputError
end

@testset "vertex-free encoding is not computed, not undefined" begin
    result = H._koda_combinatorial_spine_core(
        [[1]], [[1]],
        [Dict("id" => 1, "sign" => 1, "alpha" => 1, "beta" => 1)],
        Set([1]),
        Dict((:alpha, 1) => 1, (:beta, 1) => 1),
        nothing,
        1,
    )
    @test result.value isa H.FlowSpine
    @test H.spine_topology(result.value)["vertex_free"] == true
    err = try
        H.spine_to_odata(result.value)
        nothing
    catch e
        e
    end
    @test err isa H.KodaNotComputed
    @test !(err isa H.KodaDomainExcluded)
    @test err.reason == :vertex_free_odata_not_encoded
end

@testset "explicit h is not chosen" begin
    G = H.OData([[1, -1]], [1], Dict{String, Any}())
    missing = try
        H.odata_to_hc(G)
        nothing
    catch e
        e
    end
    @test missing isa H.KodaSpecificationPending
    @test missing.topic == :surface_identification
    bad = try
        H.odata_to_hc(G; h_map=Dict("orientation_preserving" => false))
        nothing
    catch e
        e
    end
    @test bad isa H.KodaInputError
    accepted = try
        H.odata_to_hc(G; h_map=Dict("orientation_preserving" => true))
        nothing
    catch e
        e
    end
    @test accepted isa H.KodaSpecificationPending
    @test accepted.topic == :lc_extraction
end

@testset "abalone properties, not a stored o-data" begin
    include(joinpath(@__DIR__, "..", "..", "examples", "HeegaardAbalone.jl"))
    ex = HeegaardAbaloneExample.fixture()
    produced = Tuple{Int, Any}[]
    for r in 1:length(ex.diagram.code.regions.boundary_cycle_blocks)
        G = try
            _odata(ex.diagram, r)
        catch
            continue
        end
        push!(produced, (r, G))
        closed, = _closed(G)
        @test closed
    end
    @test !isempty(produced)
    for i in 2:length(produced)
        @test H.compare_odata_words(produced[1][2], produced[i][2])["equal_up_to_metadata"]
    end
    rev, _ = H.reverse_curve_orientation(ex.diagram, (:alpha, 1))
    rev_ok = false
    for r in 1:length(rev.code.regions.boundary_cycle_blocks)
        G = try
            _odata(rev, r)
        catch
            continue
        end
        closed, = _closed(G)
        @test closed
        if H.compare_odata_words(produced[1][2], G)["equal_up_to_metadata"]
            rev_ok = true
        end
    end
    @test rev_ok
    code = ex.diagram.code.code
    points = [H.Intersection(p.id, p.alpha, p.beta, -p.sign) for p in code.intersections]
    mirror = H.certify_heegaard(H.decorate(H.heegaard_code(code.alpha, code.beta, points); cellular=true))
    mirror_closed = false
    for r in 1:length(mirror.code.regions.boundary_cycle_blocks)
        G = try
            _odata(mirror, r)
        catch
            continue
        end
        closed, = _closed(G)
        closed && (mirror_closed = true)
        @test closed
    end
    @test mirror_closed
    back = H.spine_to_odata(H.heegaard_to_spine(H.odata_to_heegaard_koda(produced[1][2]).value).value)
    @test H.compare_odata_words(produced[1][2], back)["equal_up_to_metadata"]
end

@testset "one- and two-visit curves keep o-data under orientation reversal" begin
    function core_odata(α, β, crossings, joining, dγ)
        sp = H._koda_combinatorial_spine_core(α, β, crossings, joining, dγ, nothing, 1)
        H.spine_to_odata(sp.value)
    end
    # A 1-visit curve meets J at its only crossing, so that crossing is not a true
    # vertex. Reversing the word does not change it. Flipping the contact step and
    # the incident sign must still return the same o-data. The pure 1-crossing
    # certified diagram is rejected by joining_disk (repeated boundary vertex);
    # this is the words-and-ports calculation the label test could not do.
    α = [[1], [2, 3]]
    β = [[1, 2], [3]]
    cross = [
        Dict("id" => 1, "sign" => 1, "alpha" => 1, "beta" => 1),
        Dict("id" => 2, "sign" => 1, "alpha" => 2, "beta" => 1),
        Dict("id" => 3, "sign" => -1, "alpha" => 2, "beta" => 2),
    ]
    joining = Set([1, 2])
    d_plus = Dict((:alpha, 1) => 1, (:alpha, 2) => 1, (:beta, 1) => 1, (:beta, 2) => 1)
    d_flip = Dict((:alpha, 1) => -1, (:alpha, 2) => 1, (:beta, 1) => 1, (:beta, 2) => 1)
    Gplus = core_odata(α, β, cross, joining, d_plus)
    cross_rev = [
        Dict("id" => x["id"], "sign" => x["alpha"] == 1 ? -x["sign"] : x["sign"], "alpha" => x["alpha"], "beta" => x["beta"])
        for x in cross
    ]
    Grev = core_odata(α, β, cross_rev, joining, d_flip)
    @test H.compare_odata_words(Gplus, Grev)["equal_up_to_metadata"]
    # 2-visit words: both cyclic neighbors carry the same label, so the old
    # label test returned +1. The contact step is a position. Forcing +1 changes
    # the vertex sign. Reversing the word, flipping its signs, and flipping its
    # step restores the o-data.
    α2 = [[1, 2], [3, 4]]
    β2 = [[1, 3], [2, 4]]
    cross2 = [
        Dict("id" => 1, "sign" => 1, "alpha" => 1, "beta" => 1),
        Dict("id" => 2, "sign" => -1, "alpha" => 1, "beta" => 2),
        Dict("id" => 3, "sign" => 1, "alpha" => 2, "beta" => 1),
        Dict("id" => 4, "sign" => -1, "alpha" => 2, "beta" => 2),
    ]
    join2 = Set([1, 2, 3])
    d2 = Dict((:alpha, 1) => -1, (:alpha, 2) => 1, (:beta, 1) => 1, (:beta, 2) => -1)
    d2old = Dict((:alpha, 1) => 1, (:alpha, 2) => 1, (:beta, 1) => 1, (:beta, 2) => 1)
    G2 = core_odata(α2, β2, cross2, join2, d2)
    Gold = core_odata(α2, β2, cross2, join2, d2old)
    @test !H.compare_odata_words(G2, Gold)["equal_up_to_metadata"]
    closed, = _closed(G2)
    @test closed
    α2r = [[2, 1], [3, 4]]
    cross2r = [
        Dict("id" => x["id"], "sign" => x["alpha"] == 1 ? -x["sign"] : x["sign"], "alpha" => x["alpha"], "beta" => x["beta"])
        for x in cross2
    ]
    d2r = Dict((:alpha, 1) => 1, (:alpha, 2) => 1, (:beta, 1) => 1, (:beta, 2) => -1)
    G2r = core_odata(α2r, β2, cross2r, join2, d2r)
    @test H.compare_odata_words(G2, G2r)["equal_up_to_metadata"]
end

function _test_word_presentation(hc)
    owners = Dict{Int, Dict{Symbol, String}}()
    for (name, word) in hc.lc.words
        for x in word
            get!(owners, abs(x), Dict{Symbol, String}())[x > 0 ? :+ : :-] = name
        end
    end
    alphas = sort([name for (name, family) in hc.roles if family == :alpha])
    betas = sort([name for (name, family) in hc.roles if family == :beta])
    index = Dict(name => i for (i, name) in enumerate(alphas))
    relators = Vector{Vector{Tuple{Int, Int}}}()
    for beta_name in betas
        relator = Tuple{Int, Int}[]
        for x in get(hc.lc.words, beta_name, Int[])
            label = abs(x)
            sides = get(owners, label, nothing)
            sides === nothing && continue
            haskey(sides, :+) && haskey(sides, :-) || continue
            ra, rb = get(hc.roles, sides[:+], :none), get(hc.roles, sides[:-], :none)
            ((ra == :alpha && rb == :beta) || (ra == :beta && rb == :alpha)) || continue
            if ra == :beta
                push!(relator, (index[sides[:-]], 1))
            elseif rb == :beta
                push!(relator, (index[sides[:+]], -1))
            end
        end
        push!(relators, relator)
    end
    H.GroupPresentation(["x$name" for name in alphas], relators, Int[], Dict{String, Any}())
end

function _test_count_hom(pres, table)
    order = size(table, 1)
    id = findfirst(i -> all(j -> table[i, j] == j, 1:order), 1:order)
    inv = [findfirst(j -> table[i, j] == id, 1:order) for i in 1:order]
    gens = length(pres.generators)
    gens == 0 && return 1
    total = 0
    for assign in Iterators.product(ntuple(_ -> 1:order, gens)...)
        imgs = collect(assign)
        ok = true
        for rel in pres.relators
            s = id
            for (g, sigma) in rel
                h = imgs[g]
                s = sigma > 0 ? table[s, h] : table[s, inv[h]]
            end
            s == id || (ok = false; break)
        end
        total += ok
    end
    total
end

const _TEST_TABLES = (
    ("Z2", [1 2; 2 1]),
    ("Z3", [mod(i + j - 2, 3) + 1 for i in 1:3, j in 1:3]),
    ("Z5", [mod(i + j - 2, 5) + 1 for i in 1:5, j in 1:5]),
    ("S3", [1 2 3 4 5 6; 2 3 1 6 4 5; 3 1 2 5 6 4; 4 5 6 1 2 3; 5 6 4 3 1 2; 6 4 5 2 3 1]),
)

@testset "§5.2 zero-intersection contact is not a joining disk" begin
    path = joinpath(ROOT, "threads", "PJ260922", "data", "261107-lickorish-data", "group-13", "case-01.hc")
    hc = H.read_hc(path)
    forgot = H.forget_reference(H.realize_hc(hc))
    region = only(r for r in forgot.regions if r.id == 3)
    @test region.joining_disk == false
    reason = H.joining_disk_failure(region, merge(forgot.alpha_words, forgot.beta_words), hc.roles)
    @test reason == :whole_curve || reason == :no_ab_corner
    err = try
        H.hc_to_odata(hc; joining_disk=3)
    catch e
        e
    end
    @test err isa H.KodaInputError
    @test err.detail["condition"] in ("whole_curve", "no_ab_corner")
end

@testset "§6.3 vertex-free spine, o-data still not computed" begin
    path = joinpath(ROOT, "threads", "PJ260922", "data", "261107-lickorish-data", "group-12", "case-01.hc")
    hc = H.read_hc(path)
    forgot = H.forget_reference(H.realize_hc(hc))
    region = only(r for r in forgot.regions if r.joining_disk)
    err = try
        H.hc_to_odata(hc; joining_disk=region.id)
    catch e
        e
    end
    @test err isa H.KodaNotComputed
    @test err.reason == :vertex_free_odata_not_encoded
    spine = H.heegaard_to_spine(forgot; joining_region=region.id).value
    topo = H.spine_topology(spine)
    @test topo["vertex_free"] == true
    @test haskey(topo, "vertices") && haskey(topo, "edges")
    @test all(haskey(s, "genus") && haskey(s, "boundary_count") && haskey(s, "chi") && haskey(s, "boundaries") && haskey(s, "sources") for s in topo["sectors"])
    counts = Dict{Int, Int}()
    for e in topo["edges"]
        counts[e["id"]] = 0
    end
    for sec in topo["sectors"]
        @test sec["chi"] == 2 - 2 * sec["genus"] - sec["boundary_count"]
        for word in sec["boundaries"]
            for x in word
                counts[abs(x)] = get(counts, abs(x), 0) + 1
            end
        end
    end
    @test all(n == 3 for n in values(counts))
    @test length(topo["vertices"]) - length(topo["edges"]) + sum(s["chi"] for s in topo["sectors"]) == 1
    pres_spine = H.spine_pi1_presentation(spine)
    pres_word = _test_word_presentation(hc)
    for (name, table) in _TEST_TABLES
        @test _test_count_hom(pres_spine, table) == _test_count_hom(pres_word, table)
    end
end

@testset "§6.3 true-vertex spine π1 matches HC words without o-data" begin
    path = joinpath(ROOT, "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc")
    hc = H.read_hc(path)
    forgot = H.forget_reference(H.realize_hc(hc))
    region = first(r for r in forgot.regions if r.joining_disk)
    spine = H.heegaard_to_spine(forgot; joining_region=region.id).value
    topo = H.spine_topology(spine)
    @test topo["vertex_free"] == false
    @test haskey(topo, "vertices") && haskey(topo, "edges")
    counts = Dict{Int, Int}()
    for e in topo["edges"]
        counts[e["id"]] = 0
    end
    for sec in topo["sectors"]
        @test sec["chi"] == 2 - 2 * sec["genus"] - sec["boundary_count"]
        for word in sec["boundaries"]
            for x in word
                counts[abs(x)] = get(counts, abs(x), 0) + 1
            end
        end
    end
    @test all(n == 3 for n in values(counts))
    @test length(topo["vertices"]) - length(topo["edges"]) + sum(s["chi"] for s in topo["sectors"]) == 1
    pres_spine = H.spine_pi1_presentation(spine)
    pres_word = _test_word_presentation(hc)
    for (name, table) in _TEST_TABLES
        @test _test_count_hom(pres_spine, table) == _test_count_hom(pres_word, table)
    end
end
end
