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
import JSON

const ROOT = pwd()
const OUT = joinpath(ROOT, "threads", "PJ260922", "data", "261110-hc-compare")
const LOG = joinpath(ROOT, ".research-scratch", "PJ260922", "grok", "hc-compare", "table.log")

function collect_hc(root)
    files = String[]
    data = joinpath(root, "threads", "PJ260922", "data")
    for name in readdir(data)
        startswith(name, "261107-") || continue
        dir = joinpath(data, name)
        isdir(dir) || continue
        for (dirpath, _, filenames) in walkdir(dir)
            for filename in filenames
                endswith(filename, ".hc") || continue
                push!(files, joinpath(dirpath, filename))
            end
        end
    end
    fixtures = joinpath(root, "threads", "PJ260922", "code", "lc-fixtures")
    for (dirpath, _, filenames) in walkdir(fixtures)
        for filename in filenames
            endswith(filename, ".hc") || continue
            push!(files, joinpath(dirpath, filename))
        end
    end
    sort!(files)
end

function rel_path(path)
    replace(relpath(path, ROOT), "\\" => "/")
end

function exception_record(err)
    err isa InterruptException && rethrow()
    record = Dict{String,Any}(
        "exception" => string(nameof(typeof(err))),
        "message" => sprint(showerror, err),
    )
    if err isa Heegaard.HCClosedDiagramExcluded
        record["condition"] = String(err.condition)
    elseif err isa Heegaard.HCNotComputed
        record["reason"] = String(err.reason)
    elseif err isa Heegaard.KodaInputError
        record["stage"] = String(err.stage)
    end
    record
end

function json_number(value)
    value isa BigInt && value <= typemax(Int) && value >= typemin(Int) && return Int(value)
    value isa BigInt && return string(value)
    value
end

function canonical_key(canonical)
    [
        canonical.genus,
        [[pair.first, pair.second] for pair in canonical.k_words],
        canonical.alpha_words,
        canonical.beta_words,
        [[pair.first, pair.second] for pair in canonical.other_words],
        [[pair.first, pair.second] for pair in canonical.places],
        [[pair.first, pair.second] for pair in canonical.lc_disks],
        [[pair.first, pair.second] for pair in canonical.hc_disks],
        [[pair.first, pair.second] for pair in canonical.lc_points],
        [[pair.first, pair.second] for pair in canonical.hc_points],
        [[pair.first, pair.second] for pair in canonical.lc_weights],
        [[pair.first, pair.second] for pair in canonical.hc_weights],
    ]
end

function process_file(path)
    record = Dict{String,Any}("path" => rel_path(path))
    local hc
    try
        hc = Heegaard.read_hc(path)
        Heegaard.validate_hc_input(hc)
        record["input_check"] = Dict{String,Any}("ok" => true)
    catch err
        record["input_check"] = merge(Dict{String,Any}("ok" => false), exception_record(err))
        return record
    end
    try
        intersection = Heegaard.hc_intersection(hc)
        record["intersection"] = Dict{String,Any}(
            "rows" => intersection.rows,
            "columns" => intersection.columns,
            "row_order" => String(intersection.row_order),
            "column_order" => String(intersection.column_order),
            "matrix" => [collect(intersection.matrix[i, :]) for i in axes(intersection.matrix, 1)],
            "sign_convention" => intersection.sign_convention,
        )
    catch err
        record["intersection"] = exception_record(err)
    end
    try
        factors = Heegaard.hc_invariant_factors(hc)
        record["invariant_factors"] = factors.factors
        record["invariant_factors_are_manifold_homology"] = factors.manifold_homology
    catch err
        record["invariant_factors"] = exception_record(err)
    end
    try
        presented = Heegaard.hc_presentation(hc)
        letters = sum(length(relator) for relator in presented.presentation.relators; init=0)
        record["presentation"] = Dict{String,Any}(
            "generator_family" => String(presented.generator_family),
            "relator_family" => String(presented.relator_family),
            "n_generators" => length(presented.generator_names),
            "n_relators" => length(presented.relator_names),
            "n_letters" => letters,
            "generator_names" => presented.generator_names,
            "relator_names" => presented.relator_names,
        )
    catch err
        record["presentation"] = exception_record(err)
    end
    try
        homology = Heegaard.hc_first_homology(hc)
        record["h1"] = Dict{String,Any}(
            "free_rank" => homology.free_rank,
            "torsion" => homology.torsion,
            "invariant_factors" => homology.invariant_factors,
            "condition" => String(homology.condition),
        )
    catch err
        record["h1"] = exception_record(err)
    end
    hom = Dict{String,Any}()
    hom_ok = true
    for group in (:Z2, :Z3, :Z5, :S3)
        try
            hom[String(group)] = json_number(Heegaard.hc_hom_count(hc, group))
        catch err
            hom[String(group)] = exception_record(err)
            hom_ok = false
        end
    end
    record["hom"] = hom
    record["hom_complete"] = hom_ok
    try
        record["canonical_key"] = canonical_key(Heegaard.hc_canonical(hc))
    catch err
        record["canonical"] = exception_record(err)
    end
    record
end

function class_groups(records, key_of)
    groups = Dict{String,Vector{String}}()
    for record in records
        key = key_of(record)
        key === nothing && continue
        text = JSON.json(key)
        push!(get!(() -> String[], groups, text), record["path"])
    end
    classes = [Dict("members" => sort!(members), "key" => JSON.parse(text)) for (text, members) in groups]
    sort!(classes, by=class -> class["members"][1])
    classes
end

mkpath(OUT)
mkpath(dirname(LOG))
open(LOG, "w") do io
    println(io, "start")
    flush(io)
end

files = collect_hc(ROOT)
records = Dict{String,Any}[]
open(LOG, "a") do io
    println(io, "files ", length(files))
    flush(io)
    for (index, path) in enumerate(files)
        push!(records, process_file(path))
        if index % 10 == 0 || index == length(files)
            println(io, "processed ", index, "/", length(files))
            flush(io)
        end
    end
end

identity_classes = class_groups(records, record -> get(record, "canonical_key", nothing))
function manifold_class_key(record)
    get(get(record, "input_check", Dict()), "ok", false) || return nothing
    haskey(record, "hom") || return nothing
    h1 = get(record, "h1", nothing)
    h1_key = if h1 isa AbstractDict && haskey(h1, "free_rank") && haskey(h1, "torsion")
        Dict("status" => "returned", "free_rank" => h1["free_rank"], "torsion" => h1["torsion"])
    else
        Dict(
            "status" => "not_returned",
            "mark" => h1 isa AbstractDict ? get(h1, "condition", get(h1, "exception", "missing")) : "missing",
        )
    end
    Dict("h1" => h1_key, "hom" => record["hom"])
end
invariant_classes = class_groups(records, manifold_class_key)

input_ok = count(record -> record["input_check"]["ok"], records)
matrix_ok = count(record -> haskey(record, "intersection") && haskey(record["intersection"], "matrix"), records)
h1_ok = count(record -> haskey(record, "h1") && haskey(record["h1"], "free_rank"), records)
presentation_ok = count(record -> haskey(record, "presentation") && haskey(record["presentation"], "n_generators"), records)
hom4_ok = count(record -> get(record, "hom_complete", false), records)

summary = Dict{String,Any}(
    "schema" => "hc-compare-summary/1",
    "date" => "2026-11-10",
    "hc_total" => length(records),
    "input_ok" => input_ok,
    "matrix_ok" => matrix_ok,
    "h1_ok" => h1_ok,
    "presentation_ok" => presentation_ok,
    "hom4_ok" => hom4_ok,
    "identity_classes" => length(identity_classes),
    "invariant_classes" => length(invariant_classes),
    "conventions" => Dict(
        "row_order" => "sorted_component_name",
        "generator_family" => "alpha",
        "relator_family" => "beta",
        "sign" => "+1 when + is on alpha and - is on beta",
        "h1_condition" => "alpha_beta_count_eq_genus",
        "h1_exception" => "HCClosedDiagramExcluded",
        "h1_exception_condition" => "alpha_beta_count_ne_genus",
        "invariant_class_key" => "h1 free_rank and torsion, or not_returned mark; hom Z2 Z3 Z5 S3",
    ),
)
document = Dict{String,Any}(
    "schema" => "hc-compare-table/1",
    "date" => "2026-11-10",
    "summary" => summary,
    "records" => records,
    "identity_classes" => identity_classes,
    "invariant_classes" => invariant_classes,
)
open(joinpath(OUT, "summary.json"), "w") do io
    JSON.print(io, summary, 2)
    print(io, "\n")
end
open(joinpath(OUT, "table.json"), "w") do io
    JSON.print(io, document)
    print(io, "\n")
end
println("HC_COMPARE_SUMMARY total=$(summary["hc_total"]) input_ok=$(input_ok) matrix=$(matrix_ok) h1=$(h1_ok) presentation=$(presentation_ok) hom4=$(hom4_ok) identity_classes=$(length(identity_classes))")
open(LOG, "a") do io
    println(io, "HC_COMPARE_SUMMARY total=$(summary["hc_total"]) input_ok=$(input_ok) matrix=$(matrix_ok) h1=$(h1_ok) presentation=$(presentation_ok) hom4=$(hom4_ok) identity_classes=$(length(identity_classes))")
end
