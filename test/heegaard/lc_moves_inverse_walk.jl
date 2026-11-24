# Genus-2 S³, walk 3 of the 20261115 random walk, first four steps.
# The stored inverses, applied from the end, must return the start key.
using Test
using Topology
using Random

const WalkH = Topology.Heegaard

function walk_child_rng(label)
    parent = Random.Xoshiro(20261115)
    drawn = UInt64(0)
    for name in ("s3", "L52", "L72", "g2s3", "stab52", "sum23")
        drawn = rand(parent, UInt64)
        name == label && return Random.Xoshiro(drawn), drawn
    end
    error("unknown walk label " * label)
end

function walk_choose_pool(hc, witnesses)
    genus = hc.lc.genus
    if genus > 4
        return [witness for witness in witnesses if witness.genus_delta < 0]
    elseif genus >= 4
        return [witness for witness in witnesses if witness.genus_delta <= 0]
    end
    witnesses
end

function walk_replay!(rng, root, steps)
    hc = root
    rows = String[]
    inverses = WalkH.LCMoveWitness[]
    states = Any[hc]
    for _ in 1:steps
        witnesses = WalkH.movable(hc)
        sort!(witnesses, by=witness -> (witness.row_id, witness.digest, witness.fingerprint))
        pool = walk_choose_pool(hc, witnesses)
        isempty(pool) && error("empty walk pool at genus $(hc.lc.genus)")
        witness = pool[rand(rng, eachindex(pool))]
        result = WalkH.move(hc, witness)
        push!(rows, witness.row_id)
        push!(inverses, result.inverse)
        hc = result.output
        push!(states, hc)
    end
    (rows=rows, inverses=inverses, states=states)
end

if get(ENV, "G2_ONLY", "") != "1"
@testset "genus-2 walk returns through a stabilization mark and a twist" begin
    use = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "threads", "PJ260922", "data", "261114-usecases"))
    root = WalkH.read_hc(joinpath(use, "g2-S3-pairs", "diagram.hc"))
    rng, drawn = walk_child_rng("g2s3")
    @test string(drawn) == "10663057169867408135"
    start_key = WalkH.hc_search_state(root; identity=:canonical).state_id
    for _ in 1:2
        walk_replay!(rng, root, 40)
    end
    replay = walk_replay!(rng, root, 4)
    @test replay.rows == [
        "hc.role-support.separate",
        "hc.stabilization.once",
        "hc.mark.disk.move",
        "lc.twist.annulus",
    ]
    back = replay.states[end]
    err = nothing
    for index in length(replay.inverses):-1:1
        try
            back = WalkH.move(back, replay.inverses[index]).output
        catch thrown
            err = thrown
            break
        end
    end
    @test err === nothing
    if err === nothing
        @test WalkH.hc_search_state(back; identity=:canonical).state_id == start_key
        @test WalkH.lc_text(back.lc) == WalkH.lc_text(root.lc)
        @test back.disks == root.disks
        @test back.roles == root.roles
    end
end
end # G2_ONLY skips the four-step replay

function genus_one_s3()
    lc = WalkH._lc_parse_text("genus: 1\nl1 : -1\nm1 : +1\n")
    WalkH.resolve_hc(lc, Dict("l1" => :beta); standard=true)
end

function genus_two_s3()
    use = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "threads", "PJ260922", "data", "261114-usecases"))
    WalkH.read_hc(joinpath(use, "g2-S3-pairs", "diagram.hc"))
end

function assert_prefix_returns(replay, root, start_key, last)
    back = replay.states[last + 1]
    err = nothing
    for index in last:-1:1
        try
            back = WalkH.move(back, replay.inverses[index]).output
        catch thrown
            err = thrown
            break
        end
    end
    @test err === nothing
    if err === nothing
        @test WalkH.hc_search_state(back; identity=:canonical).state_id == start_key
        @test WalkH.lc_text(back.lc) == WalkH.lc_text(root.lc)
        @test back.disks == root.disks
        @test back.roles == root.roles
    end
end

# The shortest failing prefixes of the 20261115 genus-1 walks.
# Each inverse that changes crossings must restore the linear words.
if get(ENV, "G2_ONLY", "") != "1"
@testset "genus-1 isotopy inverses return the start" begin
    root = genus_one_s3()
    rng, drawn = walk_child_rng("s3")
    @test string(drawn) == "18354032652717575510"
    start_key = WalkH.hc_search_state(root; identity=:canonical).state_id
    checks = Dict(
        1 => (19, "lc.isotopy.digon.remove.10"),
        3 => (11, "lc.isotopy.digon.remove.11"),
        5 => (34, "lc.isotopy.vertex.100"),
        10 => (14, "lc.isotopy.crossing.inverse.010"),
        12 => (38, "lc.isotopy.digon.remove.11"),
        13 => (4, "lc.isotopy.crossing.inverse.001"),
    )
    for walk_index in 1:13
        println("INVERSE_WALK s3 ", walk_index)
        flush(stdout)
        if haskey(checks, walk_index)
            last, row = checks[walk_index]
            replay = walk_replay!(rng, root, last)
            @test replay.rows[end] == row
            assert_prefix_returns(replay, root, start_key, last)
            walk_index < 13 && walk_replay!(rng, replay.states[end], 40 - last)
        else
            walk_replay!(rng, root, 40)
        end
    end
end
end # G2_ONLY skips the genus-1 replay

@testset "genus-2 isotopy inverses return the start" begin
    root = genus_two_s3()
    rng, drawn = walk_child_rng("g2s3")
    @test string(drawn) == "10663057169867408135"
    start_key = WalkH.hc_search_state(root; identity=:canonical).state_id
    checks = Dict(
        3 => (30, "lc.isotopy.digon.remove.00"),
        4 => (36, "lc.isotopy.digon.remove.10"),
        6 => (32, "lc.isotopy.digon.remove.10"),
    )
    for walk_index in 1:6
        println("INVERSE_WALK g2s3 ", walk_index)
        flush(stdout)
        if haskey(checks, walk_index)
            last, row = checks[walk_index]
            replay = walk_replay!(rng, root, last)
            @test replay.rows[end] == row
            assert_prefix_returns(replay, root, start_key, last)
            walk_replay!(rng, replay.states[end], 40 - last)
        else
            walk_replay!(rng, root, 40)
        end
    end
end
