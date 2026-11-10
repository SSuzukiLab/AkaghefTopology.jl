# 全部を値ごとの別プロセスで回すと、1 プロセスの起動がおよそ 20–40 秒なので、1040 個でおよそ 6–12 時間かかる。
using Pkg, JSON3, Printf, Dates
const PKG = dirname(dirname(@__FILE__))
const REPO = dirname(dirname(PKG))
const SCRIPT = @__FILE__

if get(ENV, "KUPERBERG_CHILD", "") == "1"
    Pkg.activate(PKG)
    using Topology
    child_z2 = [1 2; 2 1]
    child_z3 = [mod(i + j - 2, 3) + 1 for i in 1:3, j in 1:3]
    child_z5 = [mod(i + j - 2, 5) + 1 for i in 1:5, j in 1:5]
    child_s3 = [1 2 3 4 5 6; 2 3 1 6 4 5; 3 1 2 5 6 4; 4 5 6 1 2 3; 5 6 4 3 1 2; 6 4 5 2 3 1]
    name = ENV["KUPERBERG_CHILD_ALGEBRA"]
    path = ENV["KUPERBERG_CHILD_PATH"]
    table = Dict("C[Z2]" => child_z2, "C[Z3]" => child_z3, "C[Z5]" => child_z5, "C[S3]" => child_s3)
    base = replace(name, "*" => "")
    algebra = base == "KP" ? kac_paljutkin_algebra() : group_algebra(table[base])
    endswith(name, "*") && (algebra = dual_algebra(algebra))
    println(string(Z_Kuperberg(path, algebra)))
    exit(0)
end

Pkg.activate(PKG)
using Topology
const VALUE_LIMIT_SECONDS = let raw = get(ENV, "KUPERBERG_VALUE_SECONDS", "300")
    seconds = tryparse(Float64, raw)
    seconds === nothing && error("KUPERBERG_VALUE_SECONDS is not a number: $raw")
    seconds
end
const Z2 = [1 2; 2 1]
const Z3 = [mod(i + j - 2, 3) + 1 for i in 1:3, j in 1:3]
const Z5 = [mod(i + j - 2, 5) + 1 for i in 1:5, j in 1:5]
const S3_TABLE = [1 2 3 4 5 6; 2 3 1 6 4 5; 3 1 2 5 6 4; 4 5 6 1 2 3; 5 6 4 3 1 2; 6 4 5 2 3 1]
const TABLES = (("C[Z2]", Z2), ("C[Z3]", Z3), ("C[Z5]", Z5), ("C[S3]", S3_TABLE))

function group_inverses(table)
    order = size(table, 1)
    id = findfirst(i -> all(j -> table[i, j] == j, 1:order), 1:order)
    [findfirst(j -> table[i, j] == id, 1:order) for i in 1:order]
end

function count_homomorphisms(pres::Heegaard.GroupPresentation, table)
    order, ngen = size(table, 1), length(pres.generators)
    ngen == 0 && return 1
    ngen > 7 && return nothing
    inv = group_inverses(table)
    id = findfirst(i -> all(j -> table[i, j] == j, 1:order), 1:order)
    count = 0
    for tuple in Iterators.product((1:order for _ in 1:ngen)...)
        imgs = collect(tuple)
        all(rel -> begin
            s = id
            for (g, sigma) in rel
                h = imgs[g]
                s = sigma > 0 ? table[s, h] : table[s, inv[h]]
            end
            s == id
        end, pres.relators) && (count += 1)
    end
    count
end

function stop_process(pid::Integer)
    run(pipeline(`powershell.exe -NoProfile -Command "Stop-Process -Id $pid -Force"`; stdout=devnull, stderr=devnull))
end

function evaluate_bounded(path, algebra_name)
    out_path = tempname()
    err_path = tempname()
    cmd = setenv(`$(joinpath(Sys.BINDIR, Base.julia_exename())) --startup-file=no --project=$PKG $SCRIPT`,
        merge(ENV, Dict(
            "KUPERBERG_CHILD" => "1",
            "KUPERBERG_CHILD_PATH" => path,
            "KUPERBERG_CHILD_ALGEBRA" => algebra_name,
        )))
    t0 = time()
    proc = run(pipeline(cmd; stdout=out_path, stderr=err_path); wait=false)
    while process_running(proc) && time() - t0 < VALUE_LIMIT_SECONDS
        sleep(0.25)
    end
    seconds = time() - t0
    if process_running(proc)
        stop_process(getpid(proc))
        process_running(proc) && kill(proc)
        return (status="undecided", value=nothing, seconds=seconds, error="timeout")
    end
    success(proc) || return (status="error", value=nothing, seconds=seconds, error=read(err_path, String))
    text = strip(read(out_path, String))
    value = tryparse(ComplexF64, text)
    value === nothing && return (status="error", value=nothing, seconds=seconds, error="unparsed: " * text)
    (status="ok", value=value, seconds=seconds)
end

input_list = get(ENV, "KUPERBERG_INPUTS",
    joinpath(REPO, ".research-scratch", "PJ260922", "grok", "hc-kuperberg", "inputs.txt"))
scratch = get(ENV, "KUPERBERG_SCRATCH", joinpath(REPO, ".research-scratch", "PJ260922", "grok", "hc-kuperberg"))
mkpath(scratch)
stamp = Dates.format(Dates.now(), "yymmdd-HHMMSS")
log_path = joinpath(scratch, stamp * "-acceptance.log")
logio = open(log_path, "w")
paths = filter(!isempty, strip.(readlines(input_list)))
problems = Dict{String, Any}[]
judged = 0
matched = 0
undecided = 0
errored = 0
mismatched = 0
slowest = 0.0
slowest_label = ""
dual_judged = 0
dual_matched = 0

algebras = Tuple{String, Any}[]
for (name, table) in TABLES
    alg = group_algebra(table)
    push!(algebras, (name, alg))
    push!(algebras, (name * "*", dual_algebra(alg)))
end
push!(algebras, ("KP", kac_paljutkin_algebra()))
push!(algebras, ("KP*", dual_algebra(kac_paljutkin_algebra())))

println("value_limit_seconds=", VALUE_LIMIT_SECONDS, " inputs=", length(paths))
flush(stdout)

for (index, rel) in enumerate(paths)
    global judged, matched, undecided, errored, mismatched, slowest, slowest_label, dual_judged, dual_matched
    path = joinpath(REPO, rel)
    if index % 10 == 1
        println("progress ", index, "/", length(paths), " matched=", matched,
                " dual_matched=", dual_matched, " undecided=", undecided, " slowest=", round(slowest; digits=3))
        flush(stdout)
    end
    if !isfile(path)
        errored += 1
        push!(problems, Dict("input" => rel, "status" => "error", "error" => "missing file"))
        continue
    end
    local hc, pres
    try
        hc = Heegaard.read_hc(path)
        pres = heegaard_pi1_presentation(hc)
    catch err
        errored += 1
        push!(problems, Dict("input" => rel, "status" => "error", "error" => sprint(showerror, err)))
        continue
    end
    hom = Dict(name => count_homomorphisms(pres, table) for (name, table) in TABLES)
    for (name, alg) in algebras
        outcome = evaluate_bounded(path, name)
        seconds = outcome.seconds
        if seconds > slowest
            slowest = seconds
            slowest_label = rel * " " * name
        end
        is_kp = startswith(name, "KP")
        is_dual = endswith(name, "*") && !is_kp
        base = replace(name, "*" => "")
        expected = is_kp ? nothing : hom[base]
        row = Dict{String, Any}(
            "input" => rel,
            "algebra" => name,
            "seconds" => seconds,
            "status" => outcome.status,
            "expected" => expected,
        )
        if outcome.status == "undecided"
            undecided += 1
            row["error"] = get(outcome, :error, "timeout")
            push!(problems, row)
        elseif outcome.status == "error"
            errored += 1
            row["error"] = outcome.error
            push!(problems, row)
        else
            row["got"] = string(outcome.value)
            if is_kp || expected === nothing
                row["match"] = is_kp ? "not-applicable" : "hom-skipped"
            else
            ok = isapprox(outcome.value, expected; atol=1e-6)
            row["match"] = ok
            if is_dual
                dual_judged += 1
                ok && (dual_matched += 1)
            else
                judged += 1
                ok && (matched += 1)
            end
            ok || (mismatched += 1; push!(problems, row))
            end
        end
        println(logio, rel, "\t", name, "\t", outcome.status, "\t", seconds, "\t",
            something(outcome.value, ""), "\t", expected)
        flush(logio)
    end
end
close(logio)

summary = Dict(
    "value_limit_seconds" => VALUE_LIMIT_SECONDS,
    "inputs" => length(paths),
    "judged" => judged,
    "matched" => matched,
    "dual_judged" => dual_judged,
    "dual_matched" => dual_matched,
    "mismatched" => mismatched,
    "undecided" => undecided,
    "errored" => errored,
    "slowest_seconds" => slowest,
    "slowest_label" => slowest_label,
    "problems" => problems,
)
scratch = get(ENV, "KUPERBERG_SCRATCH", joinpath(REPO, ".research-scratch", "PJ260922", "grok", "hc-kuperberg"))
mkpath(scratch)
out = joinpath(scratch, stamp * "-acceptance_summary.json")
println("judged=", judged, " matched=", matched, " dual_judged=", dual_judged, " dual_matched=", dual_matched)
println("mismatched=", mismatched, " undecided=", undecided, " errored=", errored)
println("value_limit_seconds=", VALUE_LIMIT_SECONDS, " slowest=", slowest, " ", slowest_label)
println("log: ", log_path)
open(out, "w") do io
    JSON3.pretty(io, summary)
end
println("written: ", out)
