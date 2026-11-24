# Draft of lib/topology.jl/test/ograph_vertex_convention.jl.
# Unit tier: lib fixtures only. The n≤4 corpus is not this file.
#
# Definition 6 (docs/conventions/OGraphVertex.md):
#   ε(v) = sign det(t_o, t_u, X) in a right-handed frame of T_v M.
# Reading (convention page, definitions 2–6, and
# docs/specs/heegaard-koda-conversion.md §6.1 and §7.3.1):
#   X = (0,0,1) is the positive normal of S, β side toward the reader.
#   Over is β and under is α, so t_o = d_β t_β and t_u = d_α t_α.
#   t_α, t_β are the outgoing parameter tangents of the vertex chart.
#   d_γ = +1 iff the joining-disk boundary dart is the forward (odd) ribbon dart.
#
# Where the two sides come from
#   Input sign. OData.orientation, taken from vl_example, lens_space_link,
#   or connected_sum! / disjoint_sum!. It is not computed by Koda.
#   Chart tangents. This file reads load_koda_local_tables and selects the
#   model with that input sign. The port rule is the same sentence as
#   koda_hc_from_boundary (stratum "e", C_side, family, port name ending in
#   "o"). For these two models, det(t_β, t_α, X) equals the input sign when
#   d_α = d_β = +1. That equality is a rewrite of the table.
#   hc epsilon. koda_hc_from_boundary stores det(t_α, t_β) from those ports.
#   Comparing it with this file's chart det checks that the copied port rule
#   still matches that function. It is not a second geometric source.
#   d_γ. Parity of the joining-disk boundary darts on the diagram built by
#   koda_heegaard_from_hc and _koda_select_joining_region. The parity rule is
#   the same sentence as _koda_curve_directions. This file does not call
#   _koda_curve_directions, koda_vertex_sign_from_frame, or
#   koda_forward_gauss_from_hc. Those last two multiply -ε d_α d_β.
#   Ribbon dependence. ribbon_graph rotates ports by the intersection sign,
#   and that sign is hc epsilon, hence a function of the input sign. d_γ is
#   downstream of the sign. It is not assigned the value that makes the
#   product +1.
#
# For these two models the chart arithmetic is
#   det(t_o, t_u, X) = (input sign) * d_α * d_β.
# Agreement with the input sign is exactly d_α d_β = +1. The product is the
# part that is not a rewrite of the table. Reversing both strand tangents
# leaves the determinant unchanged, so definition 6 sees the product, not the
# absolute dart orientation beyond the copied parity rule.
#
# θ_G = +2 is the other planned check on the convention page. It is not here.
# koda_heegaard_from_hc and _koda_select_joining_region are unexported. They
# are the two calls spine_to_heegaard_koda makes after koda_hc_from_boundary.

module OGraphVertexConventionTests

using Test
using Topology

const H = Topology.Heegaard
const T_AFTER_USING = time()
const X_FRAME = (0, 0, 1)

function _include_examples()
    isdefined(@__MODULE__, :vl_example) && return
    candidates = [
        # Path after the lead moves this file to lib/topology.jl/test/.
        joinpath(@__DIR__, "..", "examples", "VLExample.jl"),
        # Path while the draft still sits in the lane directory.
        joinpath(@__DIR__, "..", "..", "..", "..", "lib", "topology.jl", "examples", "VLExample.jl"),
    ]
    for candidate in candidates
        isfile(candidate) || continue
        include(candidate)
        return
    end
    error("VLExample.jl not found from $(@__DIR__)")
end

_include_examples()

function det3(a, b, c)
    a[1] * (b[2] * c[3] - b[3] * c[2]) -
    a[2] * (b[1] * c[3] - b[3] * c[1]) +
    a[3] * (b[1] * c[2] - b[2] * c[1])
end

function det_otu(to, tu)
    det3((to[1], to[2], 0), (tu[1], tu[2], 0), X_FRAME)
end

function chart_positions(tables)
    Dict(parse(Int, string(k)) => (Int(v[1]), Int(v[2]))
         for (k, v) in tables.vertex["chart_port_vectors"])
end

"""Outgoing parameter tangents of one vertex model.

The model key is the input sign. This is a table lookup, not a measurement
on a diagram.
"""
function outgoing_parameter_tangents(tables, sign::Int)
    sign in (-1, 1) || throw(ArgumentError("vertex sign must be ±1"))
    model = tables.vertex["models"][string(sign)]
    C = tables.c_selection
    rule = Dict(string(r["spur_side"]) => r for r in C["vertex_rules"])
    positions = chart_positions(tables)
    tangents = Dict{String, Tuple{Int, Int}}()
    ports = Dict{String, Int}()
    for (_, row) in model["boundary_vertices"]
        string(row["stratum"]) == "e" || continue
        side = string(rule[string(row["spur_side"])]["C_side"])
        family = string(C["disk_types"]["C_in_" * side]["family"])
        chosen = Int[]
        for p in row["e_ports"]
            name = string(model["ports"][string(Int(p))]["name"])
            endswith(name, "o") && push!(chosen, Int(p))
        end
        length(chosen) == 1 || throw(ArgumentError("expected one outgoing e-port for $family"))
        port = only(chosen)
        tangents[family] = positions[port]
        ports[family] = port
    end
    haskey(tangents, "alpha") && haskey(tangents, "beta") ||
        throw(ArgumentError("model $sign has no alpha/beta outgoing tangent"))
    (tangents=tangents, ports=ports)
end

"""Dart parity on the joining disk. Odd id = forward = +1.

Same rule as H._koda_curve_directions. Not a call to that function, and not
the product -ε d_α d_β.
"""
function curve_directions(diagram, disk)
    graph = H.ribbon_graph(diagram)
    directions = Dict{Tuple{Symbol, Int}, Int}()
    for d in disk.boundary_darts
        dart = graph.darts[d]
        dart.id == d || error("ribbon dart $d is not stored at that index")
        edge = graph.edges[dart.edge]
        key = (edge.family, edge.curve_id)
        val = isodd(d) ? 1 : -1
        if haskey(directions, key) && directions[key] != val
            error("incoherent joining boundary direction for $key")
        end
        directions[key] = val
    end
    directions
end

function curve_index(words, vertex::Int)
    hits = findall(w -> string(vertex) in w, words)
    length(hits) == 1 || error("vertex $vertex meets $(length(hits)) curves")
    hits[1]
end

function vertex_frame(gauss, signs, tables)
    signs = Int[Int(x) for x in signs]
    gauss = [Int[Int(x) for x in w] for w in gauss]
    spine = H.odata_to_spine(H.OData(gauss, signs, Dict{String, Any}()); tables=tables)
    sigma = Dict(i => signs[i] for i in eachindex(signs))
    hc = H.koda_hc_from_boundary(spine.metadata["table_spine"], sigma; tables=tables)
    diagram = H.koda_heegaard_from_hc(hc)
    region = H._koda_select_joining_region(diagram, hc)
    directions = curve_directions(diagram, H.joining_disk(diagram, region))
    models = Dict(s => outgoing_parameter_tangents(tables, s) for s in (1, -1))
    rows = NamedTuple[]
    for v in eachindex(signs)
        model = models[signs[v]]
        dα = directions[(:alpha, curve_index(hc["alpha"], v))]
        dβ = directions[(:beta, curve_index(hc["beta"], v))]
        tα = model.tangents["alpha"]
        tβ = model.tangents["beta"]
        to = (dβ * tβ[1], dβ * tβ[2])
        tu = (dα * tα[1], dα * tα[2])
        chart_eps = tα[1] * tβ[2] - tα[2] * tβ[1]
        push!(rows, (
            v=v,
            sign=signs[v],
            dalpha=dα,
            dbeta=dβ,
            d_product=dα * dβ,
            det=det_otu(to, tu),
            chart_epsilon=chart_eps,
            hc_epsilon=Int(hc["epsilon"][string(v)]),
        ))
    end
    rows
end

function flipped(code, signs)
    code, Int[Int(-s) for s in signs]
end

function abalone_connected_sum(sign1::Int, sign2::Int)
    code, _ = real_gauss_code(vl_example("s3_abalone"))
    first = VirtualLink()
    second = VirtualLink()
    set_data!(first; rgauss=code, orientation=[sign1])
    set_data!(second; rgauss=code, orientation=[sign2])
    disjoint_sum!(first, second; preserve_weights=false)
    connected_sum!(first, (1, 2); transport_weights=false)
    real_gauss_code(first)
end

function lib_diagrams()
    # koda_fig13 fails C2 (the underlying graph is not connected), so
    # odata_to_spine rejects it. It is checked separately and is not a
    # definition 6 input.
    diagrams = Pair{String, Tuple{Vector{Vector{Int}}, Vector{Int}}}[]
    for name in vl_example_names()
        name == "koda_fig13" && continue
        code, signs = real_gauss_code(vl_example(name))
        push!(diagrams, "vl_example:$(name)" => (code, Int.(signs)))
    end
    code, signs = real_gauss_code(vl_example("s3_abalone"))
    push!(diagrams, "abalone_negative" => flipped(code, signs))
    push!(diagrams, "connected_sum:abalone++" => abalone_connected_sum(1, 1))
    push!(diagrams, "connected_sum:abalone+-" => abalone_connected_sum(1, -1))
    for (s, t) in ((3, 1), (4, 1), (5, 1))
        lens_code, lens_signs = real_gauss_code(lens_space_link(s, t))
        push!(diagrams, "lens_space_link:($s,$t)" => (lens_code, Int.(lens_signs)))
        push!(diagrams, "lens_space_link:($s,$t)_flipped" => flipped(lens_code, lens_signs))
    end
    diagrams
end

const TABLES = H.load_koda_local_tables()
const PLUS_MODEL = outgoing_parameter_tangents(TABLES, 1)
const MINUS_MODEL = outgoing_parameter_tangents(TABLES, -1)

function build_frames()
    frames = Dict{String, Vector{NamedTuple}}()
    seconds = Dict{String, Float64}()
    for (name, (gauss, signs)) in lib_diagrams()
        t0 = time()
        frames[name] = vertex_frame(gauss, signs, TABLES)
        seconds[name] = time() - t0
    end
    frames, seconds
end

const FRAMES, DIAGRAM_SECONDS = build_frames()

@testset "tautological: chart table equals the sign when d = +1" begin
    @test det3((-1, 0, 0), (0, -1, 0), (0, 0, 1)) == 1
    @test det_otu((-1, 0), (0, -1)) == 1
    @test PLUS_MODEL.ports == Dict("alpha" => 3, "beta" => 2)
    @test PLUS_MODEL.tangents["alpha"] == (-1, 0)
    @test PLUS_MODEL.tangents["beta"] == (0, 1)
    @test det_otu(PLUS_MODEL.tangents["beta"], PLUS_MODEL.tangents["alpha"]) == 1
    @test MINUS_MODEL.ports == Dict("alpha" => 3, "beta" => 4)
    @test MINUS_MODEL.tangents["alpha"] == (-1, 0)
    @test MINUS_MODEL.tangents["beta"] == (0, -1)
    @test det_otu(MINUS_MODEL.tangents["beta"], MINUS_MODEL.tangents["alpha"]) == -1
end

@testset "koda_fig13 is outside the closed normal o-graphs" begin
    code, signs = real_gauss_code(vl_example("koda_fig13"))
    link = VirtualLink()
    set_data!(link; rgauss=code, orientation=Int.(signs))
    closed, flags = is_closed(link)
    @test closed == false
    @test flags[2] == false
end

@testset "tautological: hc epsilon uses the same outgoing ports" begin
    @test length(FRAMES) == 18
    for (_, rows) in FRAMES
        for row in rows
            @test row.chart_epsilon == row.hc_epsilon
        end
    end
end

@testset "joining disk: d_alpha * d_beta = +1, so definition 6 matches the input sign" begin
    signs = Int[]
    for (_, rows) in FRAMES
        @test !isempty(rows)
        for row in rows
            push!(signs, row.sign)
            # Rewrite of the two chart models, on the vectors this diagram used.
            @test row.det == row.sign * row.d_product
            # Not that rewrite: the joining disk produced this product.
            @test row.d_product == 1
            @test row.det == row.sign
        end
    end
    @test any(==(1), signs)
    @test any(==(-1), signs)
end

const VERTEX_COUNT = sum(length(rows) for rows in values(FRAMES))
println("TIMING after_using_s=$(round(time() - T_AFTER_USING; digits=3)) diagrams=$(length(FRAMES)) vertices=$(VERTEX_COUNT) slowest_s=$(round(maximum(values(DIAGRAM_SECONDS)); digits=3))")

end
