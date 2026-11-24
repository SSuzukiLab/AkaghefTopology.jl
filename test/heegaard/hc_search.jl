# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.2–§6.4 — omitted selections enumerate candidates
# cf: lib/topology.jl/docs/specs/module-structure.md §3.1 — hc_search_* on MoveGraphStoreDB
using Test
using Topology

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const USE = joinpath(ROOT, "threads", "PJ260922", "data", "261114-usecases")

function genus_one_s3()
    lc = H._lc_parse_text("genus: 1\nl1 : -1\nm1 : +1\n")
    H.resolve_hc(lc, Dict("l1" => :beta); standard=true)
end

function copy_hc(hc; words=hc.lc.words, roles=hc.roles)
    copied = Dict(name => copy(word) for (name, word) in words)
    lc = H.LickorishCode(hc.lc.genus, copied, copy(hc.lc.places), deepcopy(hc.lc.disks),
        copy(hc.lc.points), copy(hc.lc.weights), nothing)
    H.HeegaardCode(lc, copy(roles), deepcopy(hc.disks), copy(hc.points), copy(hc.weights), nothing)
end

@testset "genus-1 S3 and omitted selections" begin
    hc = genus_one_s3()
    @test H.is_standard(hc)
    @test H.hc_curve_names(hc, :alpha) == ["m1"]
    @test H.hc_curve_names(hc, :beta) == ["l1"]
    homology = H.hc_first_homology(hc)
    @test homology.free_rank == 0 && isempty(homology.torsion)
    linking = H.hc_linking_form(hc)
    @test linking.free_rank == 0 && isempty(linking.torsion)
    witnesses = H.movable(hc)
    @test length(Topology.movable(hc)) == length(witnesses)
    @test length(H.hc_search_candidates(hc)) == length(witnesses)
    counts = Dict{String,Int}()
    for witness in witnesses
        counts[witness.row_id] = get(counts, witness.row_id, 0) + 1
    end
    @test counts["lc.twist.annulus"] == 4
    @test counts["hc.stabilization.once"] == 1
    @test counts["hc.role-support.separate"] == 4
    @test !haskey(counts, "hc.slide.pants")
    @test isempty(H.movable(hc, "hc.slide.pants"))
    @test isempty(H.movable(hc, "hc.stabilization.inverse"))
    for witness in H.movable(hc, "lc.twist.annulus")
        steps = witness.extra["sequence"]
        action = steps[end]
        @test action.extra["q"] in (1, -1)
        @test witness.extra["preparation"]["acted_before"] == ["l1", "m1"]
    end
end

@testset "explicit partial twist still reaches the applicator" begin
    lc = H._lc_parse_text("""
        genus: 1
        l1 : -1 -3 -2
        m1 : 1 4 5
        a : -4 -7 -6
        b : 2 6 8
        c : 7 -8 -5 3
        """)
    witness = only(H.movable(lc, "lc.twist.annulus"; q=1, support="c", acted_components=["b"]))
    @test witness.row_id == "lc.twist.annulus"
end

@testset "canonical key" begin
    lens5 = H.read_hc(joinpath(USE, "g1-D-5-2", "diagram.hc"))
    lens7 = H.read_hc(joinpath(USE, "g1-D-7-2", "diagram.hc"))
    base = H.hc_canonical(lens5)
    @test [pair.first for pair in base.k_words] == ["l1", "m1"]
    rotated = copy_hc(lens5)
    word = rotated.lc.words["beta1"]
    rotated.lc.words["beta1"] = vcat(word[2:end], word[1:1])
    @test H.hc_canonical(rotated) == base
    labels = sort!(unique!(abs.([entry for word in values(lens5.lc.words) for entry in word])))
    perm = Dict(label => label == labels[1] ? labels[2] : label == labels[2] ? labels[1] : label for label in labels)
    relabeled_words = Dict(name => [sign(entry) * perm[abs(entry)] for entry in word] for (name, word) in lens5.lc.words)
    @test H.hc_canonical(copy_hc(lens5; words=relabeled_words)) == base
    renamed_words = Dict(name => copy(word) for (name, word) in lens5.lc.words)
    renamed_words["b"] = pop!(renamed_words, "beta1")
    renamed_roles = copy(lens5.roles)
    renamed_roles["b"] = pop!(renamed_roles, "beta1")
    @test H.hc_canonical(copy_hc(lens5; words=renamed_words, roles=renamed_roles)) == base
    swapped = Dict(name => copy(word) for (name, word) in lens5.lc.words)
    swapped["l1"], swapped["m1"] = copy(swapped["m1"]), copy(swapped["l1"])
    @test H.hc_canonical(copy_hc(lens5; words=swapped)) != base
    @test H.hc_canonical(lens7) != base
end

@testset "search db round trip and one inverse" begin
    hc = genus_one_s3()
    store = Topology.open_store(":memory:")
    H.migrate_hc_search_store!(store)
    state = H.put_hc_search_state!(store, hc; identity=:canonical)
    loaded = H.load_hc_search_state(store, state.state_id)
    @test H.hc_search_state(loaded; identity=:canonical).state_id == state.state_id
    witness = first(H.movable(hc, "lc.twist.annulus"))
    result = H.move(hc, witness)
    back = H.move(result.output, result.inverse)
    @test H.hc_search_state(back.output; identity=:canonical).state_id == state.state_id
    @test hc.lc.words["l1"] == [-1] && hc.lc.words["m1"] == [1]
    script = joinpath(@__DIR__, "hc_search.jl")
    run_id = H.begin_hc_search_run!(store, script; label="round-trip", scope=Dict{String,Any}("identity" => "canonical"))
    expansion = H.expand_hc_search_state!(store, run_id, hc; depth=0, identity=:canonical, room=20)
    @test expansion.candidate_count == 9
    @test all(record.inverse_status == "ok" for record in expansion.records)
    H.finish_hc_search_run!(store, run_id; status="completed", counts=Dict{String,Any}("states" => 1 + expansion.created))
    @test !isempty(H.hc_search_neighbors(store, state.state_id))
    Topology.close_store!(store)
end
