# Bordered LC corpus sweeps: every digon-create witness and every pants slide of the
# census 00-19 compared against the bordered families (acceptance tier).
# Split out of bordered_lickorish_code.jl on 2026-11-18; that file defines the helpers
# (blc_*) and the hand-checked rows BLC_HAND_ROWS this sweep extends, and the manifest
# lists it under `requires`.
using Test
using Topology

const H = Topology.Heegaard

@testset "bordered lickorish code corpus" begin
        hand_rows = BLC_HAND_ROWS
        repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
        census = joinpath(repo, "threads", "PJ260922", "data", "261114-census", "hc")
        create_compared = 0
        create_roundtrip = 0
        create_miss = String[]
        for i in 0:19
            id = lpad(string(i), 2, '0')
            println("CENSUS create start $id")
            flush(stdout)
            lc = H.read_lc(joinpath(census, "$id.lc"))
            witnesses = H.movable(lc, "lc.isotopy.digon.create.ab")
            println("CENSUS create $id witnesses=$(length(witnesses))")
            flush(stdout)
            groups = Dict{Tuple{Int,Int,Symbol,Symbol},Vector{Any}}()
            for witness in witnesses
                eps = blc_digon_eps(witness.after_frags[1], witness.after_frags[2])
                if eps === nothing
                    create_compared += 1
                    push!(create_miss, "$id pattern $(witness.after_frags) not in the digon table")
                    continue
                end
                c_name = witness.names[witness.locations[1][1]]
                d_name = witness.names[witness.locations[2][1]]
                key = (eps[1], eps[2], blc_host_kind(lc, c_name), H._is_k(lc, d_name) ? :K : :Q)
                push!(get!(groups, key, Any[]), witness)
            end
            for (key, group) in groups
                pair = H.digon_bordered_family(key[1], key[2]; create=true, c_kind=key[3], d_kind=key[4])
                hits = H.apply_bordered_family(lc, pair)
                turned = false
                for witness in group
                    create_compared += 1
                    hit = blc_hit_for(hits, witness)
                    moved = H.move(lc, witness).output
                    if hit === nothing || !blc_same(hit.output, moved)
                        push!(create_miss, blc_mismatch(id, "lc.isotopy.digon.create.ab", witness,
                            hit === nothing ? nothing : hit.output, moved))
                        continue
                    end
                    turned && continue
                    turned = true
                    backs = H.apply_bordered_family(hit.output, H.swap_bordered_pair(pair))
                    if any(item -> blc_same(item.output, lc), backs)
                        create_roundtrip += 1
                    else
                        push!(create_miss, "$id create roundtrip locations=$([witness.names[loc[1]] for loc in witness.locations]) did not restore")
                    end
                end
            end
            println("CENSUS create $id compared=$create_compared roundtrip=$create_roundtrip miss=$(length(create_miss))")
            flush(stdout)
        end
        println("CENSUS create compared=$create_compared roundtrip=$create_roundtrip mismatches=$(length(create_miss))")
        for line in create_miss
            println(line)
        end
        @test create_compared > 0
        @test isempty(create_miss)

        slide_compared = 0
        slide_restored = 0
        slide_direct = 0
        slide_prepared = 0
        slide_miss = String[]
        class_rows = NamedTuple[]
        pres_threw = 0
        pres_both = 0
        h1_ok = 0
        h1_threw = 0
        h1_miss = String[]
        link_miss = String[]
        h1_reasons = Dict{String,Int}()
        geo_match = 0
        geo_miss = String[]
        sign_seen = Dict{Tuple{String,String,Int,Int},String}()
        for i in 0:19
            id = lpad(string(i), 2, '0')
            println("CENSUS slide start $id")
            flush(stdout)
            hc = H.read_hc(joinpath(census, "$id.hc"))
            witnesses = H.movable(hc, "hc.slide.pants")
            println("CENSUS slide $id witnesses=$(length(witnesses))")
            flush(stdout)
            for (nwatch, witness) in enumerate(witnesses)
                if nwatch % 25 == 0
                    println("CENSUS slide $id $nwatch / $(length(witnesses))")
                    flush(stdout)
                end
                prepared = haskey(witness.extra, "sequence")
                local state, action
                counted = false
                try
                    if prepared
                        state = hc
                        steps = witness.extra["sequence"]
                        for step in steps[1:end-1]
                            state = H.move(state, step).output
                        end
                        action = last(steps)
                    else
                        state = hc
                        action = witness
                    end
                    support = action.extra["support"]
                    a = String(action.extra["boundary_target"])
                    b = String(support["stationary"])
                    avoid = String[String(name) for (name, role) in state.roles if role == state.roles[a]]
                    result = H.apply_slide_bordered(state.lc, a, b, blc_band_from_route(support["route"]);
                        family_a=String(state.roles[a]), family_b=String(state.roles[b]), avoid=avoid)
                    moved = H.move(hc, witness).output
                    slide_compared += 1
                    prepared ? (slide_prepared += 1) : (slide_direct += 1)
                    counted = true
                    if !blc_same(result.output, moved.lc)
                        push!(slide_miss, blc_word_line("$id slide $a $b prepared=$prepared", result.output, moved.lc))
                    end
                    back = H.apply_slide_inverse(result.output, result.inverse)
                    if blc_same(back, state.lc)
                        slide_restored += 1
                    else
                        push!(slide_miss, blc_word_line("$id slide inverse $a $b prepared=$prepared", back, state.lc))
                    end
                    band = blc_band_from_route(support["route"])
                    sd = band["start"]["dart"]
                    fd = band["finish"]["dart"]
                    skey = (String(sd[3]), String(fd[3]), H._move_dir(sd), H._move_dir(fd))
                    if !haskey(sign_seen, skey)
                        q_word = blc_arc_entries(result.pair.left, "q")
                        r_word = blc_arc_entries(result.pair.left, "r")
                        sign_seen[skey] = "$id $a $b q=$q_word r=$r_word"
                    end
                    try
                        geo = H.apply_reversed_slide(result.output, state.lc, a, b, band; avoid=avoid)
                        if blc_same(geo, back)
                            geo_match += 1
                        else
                            push!(geo_miss, blc_word_line("$id geo $a $b prepared=$prepared", geo, back))
                        end
                    catch err
                        push!(geo_miss, "$id geo $a $b prepared=$prepared $(sprint(showerror, err))")
                    end
                    before = blc_pack_form(state)
                    after = blc_pack_form(blc_rehouse(result.output, state))
                    if before.free_rank != after.free_rank || before.torsion != after.torsion
                        push!(slide_miss, "$id slide H1 $a $b prepared=$prepared before=$(before.free_rank),$(before.torsion) after=$(after.free_rank),$(after.torsion)")
                    else
                        tag = "$id-$nwatch"
                        push!(class_rows, (tag=tag * "-before", pack=before))
                        push!(class_rows, (tag=tag * "-after", pack=after))
                    end
                    if before.pres === nothing || after.pres === nothing
                        pres_threw += 1
                    else
                        pres_both += 1
                        if before.pres.free_rank != after.pres.free_rank || before.pres.torsion != after.pres.torsion
                            push!(slide_miss, "$id slide presentation $a $b prepared=$prepared")
                        end
                    end
                    housed = blc_rehouse(result.output, state)
                    hb = blc_h1(state)
                    ha = blc_h1(housed)
                    if hb.ok && ha.ok
                        h1_ok += 1
                        if hb.free_rank != ha.free_rank || hb.torsion != ha.torsion
                            push!(h1_miss, "$id h1 $a $b prepared=$prepared before=$(hb.free_rank),$(hb.torsion) after=$(ha.free_rank),$(ha.torsion)")
                        end
                        linked = blc_link_torsion(before)
                        if hb.free_rank != before.free_rank || hb.torsion != linked
                            push!(link_miss, "$id link $a $b prepared=$prepared presentation=$(hb.free_rank),$(hb.torsion) linking=$(before.free_rank),$linked")
                        end
                    else
                        h1_threw += 1
                        reason = hb.ok ? ha.error : hb.error
                        h1_reasons[reason] = get(h1_reasons, reason, 0) + 1
                    end
                catch err
                    if !counted
                        slide_compared += 1
                        prepared ? (slide_prepared += 1) : (slide_direct += 1)
                    end
                    push!(slide_miss, "$id slide witness=$nwatch prepared=$prepared $(sprint(showerror, err))")
                end
            end
            println("CENSUS slide $id compared=$slide_compared restored=$slide_restored miss=$(length(slide_miss))")
            flush(stdout)
        end
        println("CENSUS slide compared=$slide_compared direct=$slide_direct prepared=$slide_prepared restored=$slide_restored mismatches=$(length(slide_miss)) presentation_both=$pres_both presentation_threw=$pres_threw")
        println("CENSUS h1 ok=$h1_ok threw=$h1_threw disagree=$(length(h1_miss)) link_disagree=$(length(link_miss))")
        println("CENSUS geo match=$geo_match miss=$(length(geo_miss))")
        for (key, where) in sort!(collect(pairs(sign_seen)); by=first)
            println("SIGN $(key[1]) $(key[2]) eps=$(key[3]),$(key[4]) $where")
        end
        for (reason, count) in sort!(collect(h1_reasons); by=last, rev=true)
            println("CENSUS h1 reason $count $reason")
        end
        for line in slide_miss
            println(line)
        end

        append!(class_rows, hand_rows)
        class_miss = String[]
        kk_py = joinpath(tempdir(), "blc007_kk_form.py")
        kk_in = joinpath(tempdir(), "blc007_kk_in.json")
        kk_out = joinpath(tempdir(), "blc007_kk_out.tsv")
        try
            open(kk_in, "w") do io
                println(io, "[")
                for (index, row) in enumerate(class_rows)
                    index > 1 && println(io, ",")
                    orders = blc_orders_json(row.pack.torsion)
                    matrix = blc_matrix_json(row.pack.matrix)
                    println(io, "{\"id\":\"$(row.tag)\",\"orders\":$orders,\"matrix\":$matrix}")
                end
                println(io, "]")
            end
            write(kk_py, raw"""
import importlib.util, json, sys
from fractions import Fraction
kk_path, src, dest = sys.argv[1:4]
spec = importlib.util.spec_from_file_location("kk_form", kk_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
rows = json.loads(open(src, encoding="utf-8").read())
with open(dest, "w", encoding="utf-8") as out:
    for row in rows:
        orders = [int(x) for x in row["orders"]]
        matrix = [[Fraction(cell[0], cell[1]) for cell in line] for line in row["matrix"]]
        form = module.normal_form(orders, matrix)
        out.write("{}\t{}\t{}\t{}\t{}\n".format(
            row["id"], form["rank"], form["sigma"], form["legendre"], form["sigma_stability"]))
""")
            kk = joinpath(repo, "threads", "PJ260922", "data", "261114-census", "linking-0-19", "kk_form.py")
            run(`py -3 $kk_py $kk $kk_in $kk_out`)
            classes = Dict{String,NTuple{4,String}}()
            for line in eachline(kk_out)
                parts = split(line, '\t')
                length(parts) == 5 || continue
                classes[parts[1]] = (parts[2], parts[3], parts[4], parts[5])
            end
            seen = Set{String}()
            for row in class_rows
                endswith(row.tag, "-before") || continue
                stem = row.tag[1:end-7]
                stem in seen && continue
                push!(seen, stem)
                left = get(classes, stem * "-before", nothing)
                right = get(classes, stem * "-after", nothing)
                if left === nothing || right === nothing
                    push!(class_miss, "$stem class missing")
                    continue
                end
                if left[4] != "exact-threshold" || right[4] != "exact-threshold" || left[1:3] != right[1:3]
                    push!(class_miss, "$stem before=$(left) after=$(right)")
                end
            end
        finally
            for path in (kk_py, kk_in, kk_out)
                isfile(path) && rm(path)
            end
        end
        println("CENSUS class compared=$(length(class_rows) ÷ 2) mismatches=$(length(class_miss))")
        for line in class_miss
            println(line)
        end
        @test slide_compared > 0
        @test slide_restored == slide_compared
        for line in geo_miss[1:min(end, 8)]
            println(line)
        end
        for line in h1_miss
            println(line)
        end
        for line in link_miss[1:min(end, 8)]
            println(line)
        end
        @test isempty(slide_miss)
        @test isempty(class_miss)
        @test isempty(h1_miss)
        @test isempty(link_miss)
        @test isempty(geo_miss)
        @test geo_match == slide_compared
end
