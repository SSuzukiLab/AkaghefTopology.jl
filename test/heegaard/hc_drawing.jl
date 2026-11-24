# cf: lib/topology.jl/docs/specs/drawing-objects.md §9
using Test, Topology

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

function load_lc(path)
    lc = H.read_lc(path)
    roles = Dict{String,Symbol}("alpha$i" => :alpha for i in 1:lc.genus)
    for i in 1:lc.genus
        roles["beta$i"] = :beta
    end
    H.resolve_hc(lc, roles)
end

function load_input(spec)
    path = normpath(joinpath(ROOT, spec.path))
    spec.kind == :hc ? H.read_hc(path) : load_lc(path)
end

const CASES = [
    (name="Weeks inner m", kind=:lc, path=".research-scratch/PJ260922/261111-weeks-manifold-verification/PJ260922_source_transcriptions/cases/Weeks_fig1_inner_m/Weeks_fig1_inner_m.lc", placing=false),
    (name="Koda genus 3", kind=:lc, path=".research-scratch/PJ260922/ld-draw-astra/261109_koda.ld", placing=false),
    (name="P28", kind=:lc, path="threads/PJ260922/code/lc-fixtures/p28-before.lc", placing=false),
    (name="Ozsvath-Szabo figure 2", kind=:hc, path=".research-scratch/PJ260922/transcribe/os04-fig2/figure2.hc", placing=false),
    (name="L(2,1) standard", kind=:hc, path="threads/PJ260922/data/261111-plot-acceptance/package-run-5/inputs/lens21-m/diagram.hc", placing=true),
    (name="Weeks standard", kind=:hc, path="threads/PJ260922/data/261111-plot-acceptance/package-run-5/inputs/weeks-m/diagram.hc", placing=true),
    (name="point z", kind=:hc, path="threads/PJ260922/data/261111-plot-acceptance/package-run-6/point-z/diagram.hc", placing=true),
]

function expected_counts(hc, placing)
    g = hc.lc.genus
    lset = Set{Int}(abs(x) for i in 1:g for x in get(hc.lc.words, "l$i", Int[]))
    alphas = H.hc_curve_names(hc, :alpha)
    betas = H.hc_curve_names(hc, :beta)
    aset = Set{Int}(abs(x) for n in alphas for x in hc.lc.words[n])
    bset = Set{Int}(abs(x) for n in betas for x in hc.lc.words[n])
    keys = Tuple{String,Int}[]
    if placing
        for n in vcat(alphas, betas), x in hc.lc.words[n]
            push!(keys, (n, x))
        end
    else
        for (names, other) in ((alphas, bset), (betas, aset))
            for n in names
                events = [x for x in hc.lc.words[n] if abs(x) in other || abs(x) in lset]
                if isempty(events)
                    error("case needs the unspecified key of a curve with no drawn event: $n")
                end
                append!(keys, [(n, x) for x in events])
            end
        end
    end
    face_ids = String[]
    if placing
        faces, _ = H._lc_faces(H._lc_nonempty_words(hc.lc.words))
        face_ids = [f.id for f in faces]
    end
    (curves=vcat(alphas, betas), crossings=sort!(collect(intersect(aset, bset))), edges=keys, faces=face_ids)
end

function on_segment(p, a, b; tol=1e-4)
    d = H._hp_sub(b, a)
    den = H._hp_dot(d, d)
    t = den == 0 ? 0.0 : clamp(H._hp_dot(H._hp_sub(p, a), d) / den, 0.0, 1.0)
    H._hp_dist(p, H._hp_lerp(a, b, t)) <= tol
end

function kind(d, k)
    [el for el in d.elements if el.kind == k]
end

println("makie_loaded=", any(k -> occursin("Makie", string(k)), keys(Base.loaded_modules)))
@test isdefined(H, :drawing)
@test :drawing in names(H)

const SELECTED = isempty(ARGS) ? CASES : [spec for spec in CASES if any(occursin(a, spec.name) for a in ARGS)]

@testset "Heegaard drawing objects" begin
    for spec in SELECTED
        @testset "$(spec.name)" begin
            println("CASE ", spec.name)
            flush(stdout)
            hc = load_input(spec)
            want = expected_counts(hc, spec.placing)
            seconds = @elapsed d = H.drawing(hc)
            println("  drawing $(round(seconds; digits=1))s path=", d.path)
            flush(stdout)
            @test d.path == (spec.placing ? :reference : :ordinary)
            @test sort!(String[el.key for el in kind(d, :curve)]) == sort!(collect(want.curves))
            @test sort!(Int[el.key for el in kind(d, :crossing)]) == want.crossings
            @test sort!(Tuple{String,Int}[el.key for el in kind(d, :edge)]) == sort!(want.edges)
            @test sort!(String[el.key for el in kind(d, :face)]) == sort!(collect(want.faces))
            spec.placing || @test isempty(kind(d, :face))
            result = H.plot(hc)
            for r in result.placement.runs
                pts = Tuple(r.points)
                @test any(prim.points == pts for prim in d.primitives)
            end
            for el in kind(d, :opening)
                if result.placement.square
                    y = el.key[2] == :+ ? result.placement.height : 0.0
                    @test all(abs(p[2] - y) <= 1e-6 && 0 <= p[1] <= result.placement.width for p in d.primitives[el.primitives[1]].points)
                else
                    hole = only(h for h in result.placement.holes if begin
                        m = match(r"^(top|bot)(\d+)$", h.first)
                        ("l$(m.captures[2])", startswith(h.first, "top") ? (:+) : (:-)) == el.key
                    end)
                    center, radius = hole.second
                    @test all(abs(H._hp_dist(p, center) - radius) <= 1e-4 for p in d.primitives[el.primitives[1]].points)
                end
            end
            for el in kind(d, :edge)
                runs = [r for r in result.placement.runs if r.component == el.input.curve]
                for pid in el.primitives, p in d.primitives[pid].points
                    on_run = any(on_segment(p, r.points[i], r.points[mod1(i + 1, length(r.points))]) for r in runs for i in 1:(r.closed ? length(r.points) : length(r.points) - 1))
                    on_cut = result.placement.square ?
                        (abs(p[1]) <= 1e-4 || abs(p[1] - result.placement.width) <= 1e-4 || abs(p[2]) <= 1e-4 || abs(p[2] - result.placement.height) <= 1e-4) :
                        any(abs(H._hp_dist(p, c) - r) <= 1e-3 for (_, (c, r)) in result.placement.holes)
                    @test on_run || on_cut
                end
            end
            for el in kind(d, :crossing)
                curves = (el.input.alpha, el.input.beta)
                pts = H.HCPlotPoint[]
                for (i, c) in enumerate(result.audit.crossings)
                    get(result.audit.labels, i, 0) == el.key && push!(pts, c.point)
                end
                for curve in curves, r in result.placement.runs
                    r.component == curve && !r.closed || continue
                    for side in (:start, :end)
                        get(result.contacts, (r.name, side), 0) == el.key || continue
                        push!(pts, side == :start ? first(r.points) : last(r.points))
                    end
                end
                got = H.HCPlotPoint[]
                for p in pts
                    any(H._hp_dist(p, q) <= 1e-5 for q in got) || push!(got, p)
                end
                @test length(el.primitives) == length(got)
                for pid in el.primitives
                    q = only(d.primitives[pid].points)
                    @test any(H._hp_dist(q, p) <= 1e-5 for p in got)
                end
            end
            edge = first(kind(d, :edge))
            @test H.element_at(d, edge.primitives[1], 1).key == edge.key
            cross = first(kind(d, :crossing))
            @test H.element_at(d, cross.primitives[1], 1).kind == :crossing
            @test H.crossing_curves(d, cross.key) == (alpha=cross.input.alpha, beta=cross.input.beta)
            @test occursin("V$(cross.key)", H.drawing_inspector_label(d, cross))
            if spec.placing
                face = first(kind(d, :face))
                @test H.element_at(d, face.primitives[1], 1).kind == :face
                @test !isempty(H.face_boundary(d, face.key))
            else
                @test_throws ArgumentError H.face_boundary(d, "F1")
            end
            sample = spec.name == "L(2,1) standard" ? "b" : want.curves[findfirst(n -> hc.roles[n] == :beta, want.curves)]
            n = H._hd_layout_calls[]
            prims = d.primitives
            H.style!(d, :curve; name=sample, color="#ff00aa")
            @test H._hd_layout_calls[] == n
            @test d.primitives === prims
            styled = only(el for el in d.elements if el.kind == :curve && el.key == sample)
            @test H.resolved_style(d, styled)[:color] == "#ff00aa"
            child = first(el for el in kind(d, :edge) if el.input.curve == sample)
            @test H.resolved_style(d, child)[:color] == "#ff00aa"
            H.style!(d, :curve; role=:alpha, visible=false)
            @test H.resolved_style(d, first(el for el in kind(d, :edge) if el.input.role == :alpha))[:visible] == false
            @test H.resolved_style(d, child)[:visible] == true
        end
    end
end

# cf: lib/topology.jl/docs/specs/drawing-objects.md §3.3.1 — regions are the components of Σ ∖ (α ∪ β)
@testset "regions glue faces across K without a role" begin
    for spec in CASES
        spec.placing || continue
        hc = load_input(spec)
        d = H.drawing(hc)
        faces = [el for el in d.elements if el.kind == :face]
        regions = [el for el in d.elements if el.kind == :region]
        region_of = Dict(f => r.key for r in regions for f in r.input.faces)
        # every face lies in exactly one region
        @test sort!(collect(keys(region_of))) == sort!([f.key for f in faces])
        @test sum(length(r.input.faces) for r in regions) == length(faces)
        # the two sides of an edge of a K curve without a role lie in one region; the face boundaries come from the
        # placement (_hp_mark_regions), the regions from forget_reference, so this compares two computations
        sides = Dict{Any,Vector{String}}()
        for f in faces, b in f.input.boundary
            push!(get!(sides, b.edge, String[]), f.key)
        end
        free = Set(c for c in keys(hc.lc.words) if !haskey(hc.roles, c))
        for (edge, fs) in sides
            edge[1] in free || continue
            @test length(unique(region_of[f] for f in fs)) == 1
        end
        # a region of more than one face has a polygon for each face
        @test all(length(r.primitives) == length(r.input.faces) for r in regions)
        # selecting a region by one of its faces
        f = first(faces).key
        H.style!(d, :region; face=f, fill=("#ffd700", 0.45))
        hit = [r for r in regions if H.resolved_style(d, r)[:fill] == ("#ffd700", 0.45)]
        @test length(hit) == 1 && f in only(hit).input.faces
    end
end

if any(k -> occursin("CairoMakie", string(k)), keys(Base.loaded_modules))
    @testset "Makie Weeks standard" begin
        hc = H.read_hc(normpath(joinpath(ROOT, "threads/PJ260922/data/261111-plot-acceptance/package-run-5/inputs/weeks-m/diagram.hc")))
        d = H.drawing(hc)
        beta = first(el.key for el in d.elements if el.kind == :curve && el.input.role == :beta)
        face = first(el.key for el in d.elements if el.kind == :face)
        H.style!(d, :curve; name=beta, color="#00aa46")
        H.style!(d, :region; face=face, fill=("#ffd700", 0.45))
        fig = H.drawingplot(d)
        @test only(block for block in fig.content if block isa CairoMakie.Axis).yreversed[]
        out = normpath(joinpath(ROOT, ".research-scratch/codex-lanes/pj922-drawing/out"))
        mkpath(out)
        svg = joinpath(out, "weeks-standard.svg")
        pdf = joinpath(out, "weeks-standard.pdf")
        CairoMakie.save(svg, fig)
        CairoMakie.save(pdf, fig)
        svg_text = read(svg, String)
        pdf_text = String(read(pdf))
        @test occursin("00aa46", lowercase(svg_text)) || occursin("0,170,70", svg_text) || occursin("0 0.666", lowercase(svg_text))
        @test occursin("ffd700", lowercase(svg_text)) || occursin("255,215,0", svg_text) || occursin("1,0.843", lowercase(svg_text))
        @test filesize(pdf) > 1000
        @test occursin("00aa46", lowercase(pdf_text)) || occursin("ffd700", lowercase(pdf_text)) || occursin("0 0.666", pdf_text)
        println("MAKIE saved ", svg)
    end
else
    println("CairoMakie is not loaded; Makie save check is skipped in this process")
end
