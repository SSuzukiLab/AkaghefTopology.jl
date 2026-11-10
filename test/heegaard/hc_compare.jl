using Test
using Topology

function _hc_compare_source()
    specified = joinpath("lib", "topology.jl", "src", "heegaard", "hc_compare.jl")
    dir = @__DIR__
    for _ in 1:10
        candidate = joinpath(dir, specified)
        isfile(candidate) && return abspath(candidate)
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    abspath(specified)
end

Base.include(Topology.Heegaard, _hc_compare_source())

function layer(genus, words, roles;
        places=Dict{String,String}(),
        disks=Dict{String,Vector{String}}(),
        points=Dict{String,String}(),
        weights=Dict{String,Int}())
    lc = Heegaard.LickorishConfiguration(
        genus, words, places,
        Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
    Heegaard.HeegaardLayer(lc, roles, disks, points, weights, nothing)
end

function copied_words(h)
    Dict{String,Vector{Int}}(String(name) => copy(word) for (name, word) in h.lc.words)
end

function copied_roles(h)
    Dict{String,Symbol}(String(name) => role for (name, role) in h.roles)
end

function with_words(h, words, roles=copied_roles(h); places=copy(h.lc.places))
    layer(h.lc.genus, words, roles; places=places, disks=copy(h.disks), points=copy(h.points), weights=copy(h.weights))
end

function relabel_words(words, mapping)
    Dict(name => Int[sign(x) * mapping[abs(x)] for x in word] for (name, word) in words)
end

function rotate_word(word, start)
    n = length(word)
    n == 0 && return Int[]
    Int[word[mod1(start + i - 1, n)] for i in 1:n]
end

const SAMPLE = let
    words = Dict(
        "l1" => [-1, 5],
        "m1" => [1],
        "alpha1" => [-5, 2, -3],
        "alpha2" => [4, -6],
        "beta1" => [-2, 3],
        "beta2" => [-4, 6],
    )
    roles = Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta)
    layer(1, words, roles)
end

function status_of(comparison, level)
    getfield(comparison, level).status
end

ts = @testset "hc_compare" begin
    @testset "sign convention and presentation orientation" begin
        positive = layer(1, Dict("a" => [1], "b" => [-1]), Dict("a" => :alpha, "b" => :beta))
        negative = layer(1, Dict("a" => [-1], "b" => [1]), Dict("a" => :alpha, "b" => :beta))
        plus = Heegaard.hc_intersection(positive)
        minus = Heegaard.hc_intersection(negative)
        @test plus.rows == ["a"] && plus.columns == ["b"]
        @test plus.matrix == ones(Int, 1, 1)
        @test minus.matrix == -ones(Int, 1, 1)
        @test plus.row_order == :sorted_component_name
        presented = Heegaard.hc_presentation(positive)
        @test presented.generator_family == :alpha
        @test presented.relator_family == :beta
        @test presented.presentation.metadata["generator_family"] == "alpha"
        @test presented.presentation.metadata["relator_family"] == "beta"
        @test Heegaard.abelianization_matrix(presented.presentation) == transpose(plus.matrix)
    end

    @testset "closed diagram condition is named and does not block the matrix" begin
        extra = layer(1, Dict("a" => [1], "c" => [-2], "b" => [-1, 2]),
            Dict("a" => :alpha, "c" => :alpha, "b" => :beta))
        @test Heegaard.hc_intersection(extra).matrix == reshape([1, -1], 2, 1)
        @test Heegaard.hc_invariant_factors(extra).manifold_homology == false
        err = try
            Heegaard.hc_first_homology(extra)
            nothing
        catch e
            e
        end
        @test err isa Heegaard.HCClosedDiagramExcluded
        @test err.condition == :alpha_beta_count_ne_genus
        @test err.witness["genus"] == 1
        @test err.witness["n_alpha"] == 2
        @test err.witness["n_beta"] == 1
    end

    @testset "empty words and genus 0 are presentations, not missing values" begin
        empty_curves = layer(0, Dict("a" => Int[], "b" => Int[]), Dict("a" => :alpha, "b" => :beta);
            places=Dict("a" => "outer", "b" => "a+"))
        presented = Heegaard.hc_presentation(empty_curves)
        @test presented.generator_names == ["a"]
        @test presented.relator_names == ["b"]
        @test presented.presentation.relators == [Tuple{Int,Int}[]]
        err = try
            Heegaard.hc_first_homology(empty_curves)
            nothing
        catch e
            e
        end
        @test err isa Heegaard.HCClosedDiagramExcluded
        @test err.condition == :alpha_beta_count_ne_genus
        trivial = layer(0, Dict{String,Vector{Int}}(), Dict{String,Symbol}())
        trivial_h = Heegaard.hc_first_homology(trivial)
        @test trivial_h.free_rank == 0 && isempty(trivial_h.torsion)
        @test Heegaard.hc_hom_count(trivial, :Z2) == 1
        @test Heegaard.hc_hom_count(trivial, :S3) == 1
        no_crossing_alpha = layer(1, Dict("a" => Int[], "b" => Int[]), Dict("a" => :alpha, "b" => :beta))
        @test length(Heegaard.hc_presentation(no_crossing_alpha).presentation.relators) == 1
        @test isempty(only(Heegaard.hc_presentation(no_crossing_alpha).presentation.relators))
    end

    @testset "homology and hom counts follow the presentation" begin
        meridian = layer(1, Dict("a" => [1], "b" => [-1]), Dict("a" => :alpha, "b" => :beta))
        h = Heegaard.hc_first_homology(meridian)
        @test h.free_rank == 0 && isempty(h.torsion)
        @test Heegaard.hc_hom_count(meridian, :Z2) == 1
        @test Heegaard.hc_hom_count(meridian, :Z3) == 1
        @test Heegaard.hc_hom_count(meridian, :Z5) == 1
        @test Heegaard.hc_hom_count(meridian, :S3) == 1
        torsion = layer(1, Dict("a" => [1, 2], "b" => [-1, -2]), Dict("a" => :alpha, "b" => :beta))
        ht = Heegaard.hc_first_homology(torsion)
        @test ht.free_rank == 0 && ht.torsion == [2]
        @test Heegaard.hc_invariant_factors(torsion).factors == [2]
        @test Heegaard.hc_hom_count(torsion, :Z2) == 2
        @test Heegaard.hc_hom_count(torsion, :Z3) == 1
        @test Heegaard.hc_hom_count(torsion, :Z5) == 1
        @test Heegaard.hc_hom_count(torsion, :S3) == 4
        free = layer(1, Dict("a" => [1, -2], "b" => [-1, 2]), Dict("a" => :alpha, "b" => :beta))
        hf = Heegaard.hc_first_homology(free)
        @test hf.free_rank == 1 && isempty(hf.torsion)
        @test Heegaard.hc_hom_count(free, :Z2) == 2
        @test Heegaard.hc_hom_count(free, :Z5) == 5
        @test Heegaard.hc_hom_count(free, :S3) == 6
    end

    @testset "identity allows relabeling, rotation, and within-role renaming only" begin
        mapping = Dict(1 => 9, 2 => 4, 3 => 8, 4 => 2, 5 => 7, 6 => 3)
        renamed_labels = with_words(SAMPLE, relabel_words(copied_words(SAMPLE), mapping))
        rotated = with_words(SAMPLE, Dict(name => rotate_word(word, 2) for (name, word) in copied_words(SAMPLE)))
        swapped_words = copied_words(SAMPLE)
        swapped_words["alpha2"] = copy(SAMPLE.lc.words["alpha1"])
        swapped_words["alpha1"] = copy(SAMPLE.lc.words["alpha2"])
        swapped_roles = Dict("alpha2" => :alpha, "alpha1" => :alpha, "beta1" => :beta, "beta2" => :beta)
        swapped = with_words(SAMPLE, swapped_words, swapped_roles)
        @test Heegaard.hc_identical(SAMPLE, renamed_labels)
        @test Heegaard.hc_identical(SAMPLE, rotated)
        @test Heegaard.hc_identical(SAMPLE, swapped)
        @test status_of(Heegaard.compare_hc(SAMPLE, swapped), :identity) == :equal
        changed = copied_words(SAMPLE)
        changed["alpha1"] = [-5, -3, 2]
        @test !Heegaard.hc_identical(SAMPLE, with_words(SAMPLE, changed))
        @test status_of(Heegaard.compare_hc(SAMPLE, with_words(SAMPLE, changed)), :identity) == :different
        import Random
        rng = Random.MersenneTwister(20261110)
        labels = sort!(unique(abs(x) for word in values(SAMPLE.lc.words) for x in word))
        for _ in 1:8
            perm = Random.randperm(rng, length(labels))
            label_map = Dict(labels[i] => 20 + perm[i] for i in eachindex(labels))
            words = relabel_words(copied_words(SAMPLE), label_map)
            for name in keys(words)
                n = length(words[name])
                n == 0 && continue
                words[name] = rotate_word(words[name], rand(rng, 1:n))
            end
            roles = copied_roles(SAMPLE)
            alpha_names = sort!(collect(keys(filter(pair -> pair.second == :alpha, roles))))
            shuffle_alpha = alpha_names[Random.randperm(rng, length(alpha_names))]
            rename = Dict(alpha_names[i] => shuffle_alpha[i] for i in eachindex(alpha_names))
            new_words = Dict(get(rename, name, name) => word for (name, word) in words)
            new_roles = Dict(get(rename, name, name) => role for (name, role) in roles)
            @test Heegaard.hc_identical(SAMPLE, with_words(SAMPLE, new_words, new_roles))
        end
    end

    @testset "orientation reversal and mirror are not identity" begin
        reversed = Heegaard.hc_reverse_orientation(SAMPLE, ["alpha1"])
        mirrored = Heegaard.hc_mirror(SAMPLE)
        @test !Heegaard.hc_identical(SAMPLE, reversed)
        @test !Heegaard.hc_identical(SAMPLE, mirrored)
        relations = Heegaard.hc_named_relations(SAMPLE, reversed)
        @test relations.identical == false
        @test relations.orientation_reversal == :yes
        mirror_relations = Heegaard.hc_named_relations(SAMPLE, mirrored)
        @test mirror_relations.identical == false
        @test mirror_relations.mirror == :yes
        @test status_of(Heegaard.compare_hc(SAMPLE, reversed), :identity) == :different
        @test status_of(Heegaard.compare_hc(SAMPLE, mirrored), :identity) == :different
    end

    @testset "agreeing invariants are not manifold equality" begin
        left = layer(1, Dict("a" => [1, -2], "b" => [-1, 2]), Dict("a" => :alpha, "b" => :beta))
        right = layer(1, Dict("a" => [1, -2, 3, -4], "b" => [-1, 2, -3, 4]), Dict("a" => :alpha, "b" => :beta))
        @test !Heegaard.hc_identical(left, right)
        @test Heegaard.hc_first_homology(left).free_rank == Heegaard.hc_first_homology(right).free_rank == 1
        for group in (:Z2, :Z3, :Z5, :S3)
            @test Heegaard.hc_hom_count(left, group) == Heegaard.hc_hom_count(right, group)
        end
        comparison = Heegaard.compare_hc(left, right)
        @test comparison.manifold.status == :undetermined
        @test comparison.manifold.reason == :invariants_agree_not_identity
        @test "h1" in comparison.manifold.evidence["agreeing_invariants"]
        @test "hom_S3" in comparison.manifold.evidence["agreeing_invariants"]
        @test comparison.isotopy.status == :undetermined
        @test comparison.move_equivalence.status == :undetermined
        different = layer(1, Dict("a" => [1], "b" => [-1]), Dict("a" => :alpha, "b" => :beta))
        separated = Heegaard.compare_hc(left, different)
        @test separated.isotopy.status == :different
        @test separated.move_equivalence.status == :different
        @test separated.manifold.status == :different
        witnessed = Heegaard.HCMoveSearchResult([
            Heegaard.HCMoveStep(:handle_slide, Dict{String,Any}("source" => "caller")),
        ])
        moved = Heegaard.compare_hc(left, right; move_result=witnessed)
        @test moved.identity.status == :different
        @test moved.isotopy.status == :undetermined
        @test moved.move_equivalence.status == :equal
        @test moved.manifold.status == :equal
        @test moved.manifold.reason == :move_sequence
        isotopic = Heegaard.HCMoveSearchResult([Heegaard.HCMoveStep(:isotopy, Dict{String,Any}())])
        by_isotopy = Heegaard.compare_hc(left, right; move_result=isotopic)
        @test by_isotopy.isotopy.status == :equal
        contradicted = Heegaard.compare_hc(left, different; move_result=witnessed)
        @test contradicted.manifold.status == :different
        @test contradicted.move_equivalence.status == :different
    end

    @testset "old HeegaardCode is not an entry type" begin
        @test !hasmethod(Heegaard.compare_hc, (Heegaard.HeegaardCode, Heegaard.HeegaardCode))
        @test !hasmethod(Heegaard.hc_intersection, (Heegaard.HeegaardCode,))
        @test !hasmethod(Heegaard.hc_first_homology, (Heegaard.HeegaardCode,))
        @test !hasmethod(Heegaard.hc_presentation, (Heegaard.HeegaardCode,))
        @test hasmethod(Heegaard.compare_hc, (Heegaard.HeegaardLayer, Heegaard.HeegaardLayer))
    end

    @testset "a stored file is identical to its own rotation and relabeling" begin
        loaded = Heegaard.read_hc(joinpath("threads", "PJ260922", "code", "lc-fixtures", "p28-before.hc"))
        words = copied_words(loaded)
        name = "beta1"
        words[name] = rotate_word(words[name], 2)
        labels = sort!(unique(abs(x) for word in values(words) for x in word))
        mapping = Dict(labels[i] => 100 + i for i in eachindex(labels))
        @test Heegaard.hc_identical(loaded, with_words(loaded, relabel_words(words, mapping)))
        changed = copied_words(loaded)
        word = changed["alpha1"]
        @test length(word) >= 3
        changed["alpha1"] = [word[1], word[3], word[2]]
        @test !Heegaard.hc_identical(loaded, with_words(loaded, changed))
    end

    @testset "an unfinished tie search throws and is not a difference" begin
        words = Dict{String,Vector{Int}}()
        roles = Dict{String,Symbol}()
        for i in 1:9
            words["a$i"] = [i]
            words["b$i"] = [-i]
            roles["a$i"] = :alpha
            roles["b$i"] = :beta
        end
        tied = layer(0, words, roles)
        err = try
            Heegaard.compare_hc(tied, tied)
            nothing
        catch e
            e
        end
        @test err isa Heegaard.HCNotComputed
        @test err.reason == :canonical_tie_group
    end

    @testset "allowed transforms of every stored diagram stay identical" begin
        root = let
            dir = @__DIR__
            found = pwd()
            for _ in 1:10
                if isfile(joinpath(dir, "AGENTS.md"))
                    found = dir
                    break
                end
                parent = dirname(dir)
                parent == dir && break
                dir = parent
            end
            found
        end
        files = String[]
        data = joinpath(root, "threads", "PJ260922", "data")
        for name in readdir(data)
            startswith(name, "261107-") || continue
            for (dirpath, _, filenames) in walkdir(joinpath(data, name))
                for filename in filenames
                    endswith(filename, ".hc") && push!(files, joinpath(dirpath, filename))
                end
            end
        end
        fixtures = joinpath(root, "threads", "PJ260922", "code", "lc-fixtures")
        for (dirpath, _, filenames) in walkdir(fixtures)
            for filename in filenames
                endswith(filename, ".hc") && push!(files, joinpath(dirpath, filename))
            end
        end
        sort!(files)
        @test length(files) == 138
        import Random
        rng = Random.MersenneTwister(20261110)
        checked = 0
        function transport_disk(dict, transport)
            Dict{String,Vector{String}}(String(name) => String[transport(ref) for ref in refs] for (name, refs) in dict)
        end
        function transport_point(dict, transport)
            Dict{String,String}(String(name) => transport(ref) for (name, ref) in dict)
        end
        function transport_weight(dict, transport)
            Dict{String,Int}(transport(String(name)) => Int(value) for (name, value) in dict)
        end
        function allowed_copy(h; start_of=nothing)
            words0 = copied_words(h)
            roles0 = copied_roles(h)
            labels = sort!(unique(abs(x) for word in values(words0) for x in word))
            label_map = Dict{Int,Int}()
            if !isempty(labels)
                perm = Random.randperm(rng, length(labels))
                for i in eachindex(labels)
                    label_map[labels[i]] = 1000 + perm[i]
                end
            end
            name_map = Dict{String,String}()
            for role in (:alpha, :beta)
                names = sort!(collect(keys(filter(pair -> pair.second == role, roles0))))
                isempty(names) && continue
                shuffled = names[Random.randperm(rng, length(names))]
                for (old, new) in zip(names, shuffled)
                    name_map[old] = new
                end
            end
            for name in keys(words0)
                haskey(name_map, name) || (name_map[name] = name)
            end
            for name in keys(h.lc.places)
                sname = String(name)
                haskey(name_map, sname) || (name_map[sname] = sname)
            end
            new_words = Dict{String,Vector{Int}}()
            for (name, word) in words0
                n = length(word)
                start = n == 0 ? 1 : (start_of === nothing ? rand(rng, 1:n) : start_of(n))
                rotated = rotate_word(word, start)
                new_words[name_map[name]] = isempty(label_map) ? Int[] : Int[sign(x) * label_map[abs(x)] for x in rotated]
            end
            new_roles = Dict{String,Symbol}(get(name_map, name, name) => role for (name, role) in roles0)
            face_map = Heegaard._hc_face_map(h.lc.words, label_map)
            transport = ref -> Heegaard._hc_transport_ref(ref, name_map, face_map, label_map)
            places = Dict{String,String}(name_map[String(name)] => transport(ref) for (name, ref) in h.lc.places)
            lc = Heegaard.LickorishConfiguration(
                h.lc.genus, new_words, places,
                transport_disk(h.lc.disks, transport),
                transport_point(h.lc.points, transport),
                transport_weight(h.lc.weights, transport),
                nothing,
            )
            Heegaard.HeegaardLayer(
                lc, new_roles,
                transport_disk(h.disks, transport),
                transport_point(h.points, transport),
                transport_weight(h.weights, transport),
                nothing,
            )
        end
        function carried(original, copy)
            status_of(Heegaard.compare_hc(original, copy), :identity) == :equal || return false
            corr = Heegaard.hc_identity_correspondence(original, copy)
            applied = Heegaard.hc_apply_correspondence(original, corr)
            copied_words(applied) == copied_words(copy) && Heegaard.hc_identical(applied, copy)
        end
        for path in files
            original = Heegaard.read_hc(path)
            @test carried(original, allowed_copy(original; start_of=n -> n == 0 ? 1 : min(2, n)))
            checked += 1
            for _ in 1:5
                @test carried(original, allowed_copy(original))
                checked += 1
            end
        end
        println("HC_COMPARE_PROPERTY checked=$checked")
        @test checked == 138 * 6
    end

    @testset "group-21 beta1 reversal is a real symmetry, recovered by the correspondence" begin
        root = let
            dir = @__DIR__
            found = pwd()
            for _ in 1:10
                if isfile(joinpath(dir, "AGENTS.md"))
                    found = dir
                    break
                end
                parent = dirname(dir)
                parent == dir && break
                dir = parent
            end
            found
        end
        original = Heegaard.read_hc(joinpath(
            root, "threads", "PJ260922", "data", "261107-lickorish-data", "group-21", "case-01.hc"))
        reversed = Heegaard.hc_reverse_orientation(original, ["beta1"])
        @test status_of(Heegaard.compare_hc(original, reversed), :identity) == :equal
        corr = Heegaard.hc_identity_correspondence(reversed, original)
        applied = Heegaard.hc_apply_correspondence(reversed, corr)
        @test copied_words(applied) == copied_words(original)
        @test Heegaard.hc_identical(applied, original)
        println("HC_COMPARE_GROUP21 beta1_label_map=", sort!(collect(pairs(corr.label_map))))
        println("HC_COMPARE_GROUP21 name_map=", sort!(collect(pairs(corr.name_map))))
        println("HC_COMPARE_GROUP21 starts=", sort!(collect(pairs(corr.rotation_starts))))
    end
end

function _count_tests(testset)
    passed = failed = errored = 0
    for result in testset.results
        if result isa Test.DefaultTestSet
            p, f, e = _count_tests(result)
            passed += p
            failed += f
            errored += e
        elseif result isa Test.Pass
            passed += 1
        elseif result isa Test.Fail
            failed += 1
        elseif result isa Test.Error
            errored += 1
        end
    end
    passed, failed, errored
end

passed, failed, errored = _count_tests(ts)
println("HC_COMPARE_TESTS failed=$failed errored=$errored")
println("HC_COMPARE_TESTS note=Julia does not store passing results; use the Test Summary line")
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit((failed + errored) == 0 ? 0 : 1)
end
