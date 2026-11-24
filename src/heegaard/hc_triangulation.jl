# cf: lib/topology.jl/src/ograph_tetrahedra.jl::disk_tetrahedra — consume its checked face maps, without a second port convention
# cf: lib/topology.jl/docs/specs/heegaard-snappy.md — bridge contract and indexing

struct HCTriangulationError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::HCTriangulationError) = print(io, "HCTriangulationError(", e.reason, "): ", e.detail)

struct SnapPyUnavailable <: Exception
    detail::String
end
Base.showerror(io::IO, e::SnapPyUnavailable) = print(io, "SnapPyUnavailable: ", e.detail)

struct SnapPyExecutionError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::SnapPyExecutionError) = print(io, "SnapPyExecutionError(", e.reason, "): ", e.detail)

"""An unoriented isoSig and an explicitly oriented, finite closed triangulation.

Tetrahedron indices in `gluings` are zero based; `permutation[i+1]` is the
image of vertex i, including the omitted vertex. Each paired face occurs once.
`isosig === nothing` only when `signature=false` was explicitly requested.
`branch_orders` is the GKS vertex order carried into the positive vertex frame,
one permutation of `0,1,2,3` per tetrahedron. It is a branching, not an
orientation. Empty means this triangulation has no recorded branching.
"""
struct HCClosedTriangulation
    tetrahedra::Int
    gluings::Vector{NamedTuple}
    signs::Vector{Int}
    isosig::Union{Nothing,String}
    joining_disk::Any
    correspondence::Dict{String,Any}
    topology::NamedTuple
    branch_orders::Vector{Vector{Int}}
    function HCClosedTriangulation(tetrahedra, gluings, signs, isosig, joining_disk,
                                   correspondence, topology,
                                   branch_orders::Vector{Vector{Int}}=Vector{Vector{Int}}())
        new(tetrahedra, gluings, signs, isosig, joining_disk, correspondence, topology, branch_orders)
    end
end

"""Count quotient cells from explicit simplicial pairings; no o-graph recipe here."""
function _hc_triangulation_topology(n::Int, gluings)
    n > 0 || throw(HCTriangulationError(:no_tetrahedra, "a true vertex is required"))
    vertices = collect(1:4n)
    edges = collect(1:6n)
    edge_pairs = [(a,b) for a in 0:3 for b in a+1:3]
    edge_id = Dict(e => i for (i,e) in enumerate(edge_pairs))
    root(p, i) = (while p[i] != i; i = p[i]; end; i)
    unite!(p, i, j) = (p[root(p,j)] = root(p,i))
    seen = Set{Tuple{Int,Int}}()
    for g in gluings
        s, t, a, b = g.tet, g.to_tet, g.omitted, g.to_omitted
        0 <= s < n && 0 <= t < n && 0 <= a <= 3 && 0 <= b <= 3 ||
            throw(HCTriangulationError(:invalid_face, "face index out of range"))
        p = collect(g.permutation)
        sort(p) == collect(0:3) && p[a+1] == b ||
            throw(HCTriangulationError(:invalid_permutation, "face pairing is not a compatible permutation"))
        for face in ((s,a),(t,b))
            face in seen && throw(HCTriangulationError(:reused_face, string(face)))
            push!(seen, face)
        end
        for i in 0:3
            i == a && continue
            unite!(vertices, 4s+i+1, 4t+p[i+1]+1)
        end
        for (i,j) in edge_pairs
            a in (i,j) && continue
            target = minmax(p[i+1], p[j+1])
            unite!(edges, 6s+edge_id[(i,j)], 6t+edge_id[target])
        end
    end
    length(seen) == 4n || throw(HCTriangulationError(:unpaired_face, "not all 4n faces are paired"))
    nv = length(Set(root(vertices,i) for i in eachindex(vertices)))
    ne = length(Set(root(edges,i) for i in eachindex(edges)))
    nf = length(gluings)
    (; vertices=nv, edges=ne, faces=nf, euler_characteristic=nv-ne+nf-n)
end

"""Serialize the explicit triangulation, keeping the selected disk and source maps."""
function hc_triangulation_data(t::HCClosedTriangulation)
    Dict{String,Any}(
        "schema" => "heegaard.closed-triangulation/1", "tetrahedra" => t.tetrahedra,
        "gluings" => [Dict(string(k)=>v for (k,v) in pairs(g)) for g in t.gluings],
        "signs" => t.signs, "isosig" => t.isosig, "joining_disk" => t.joining_disk,
        "branch_orders" => t.branch_orders, "correspondence" => t.correspondence,
        "topology" => Dict(string(k)=>v for (k,v) in pairs(t.topology)))
end

function _hc_closed_triangulation(v; joining_disk=nothing, correspondence=Dict{String,Any}(),
                                  signature::Bool=true, runner=nothing, timeout::Real=120)
    data = parentmodule(@__MODULE__).disk_tetrahedra(v)
    hasproperty(data, :face_gluings) || throw(HCTriangulationError(
        :face_gluings_unavailable,
        "disk_tetrahedra must expose its already checked face_gluings; parent integration patch is required"))
    gluings = NamedTuple[g for g in data.face_gluings]
    n = length(data.tetrahedra)
    topology = _hc_triangulation_topology(n, gluings)
    topology.euler_characteristic == 0 || throw(HCTriangulationError(:nonzero_euler, string(topology)))
    topology.vertices == 1 || throw(HCTriangulationError(:not_one_vertex,
        "a certified closed normal o-graph must give one capped sphere-link vertex; got $(topology.vertices)"))
    corr = merge(Dict{String,Any}(correspondence), Dict{String,Any}(
        "tetrahedron_to_ograph_vertex" => [t.vertex for t in data.tetrahedra],
        "tetrahedron_edges_to_disks" => [collect(t.edges) for t in data.tetrahedra],
        "node_to_disk" => data.node_to_disk,
        "ograph_gauss" => deepcopy(v.gauss), "ograph_orientation" => copy(v.orientation)))
    result = HCClosedTriangulation(n, gluings, [t.sign for t in data.tetrahedra],
                                  nothing, joining_disk, corr, topology)
    if signature
        answer = hc_snappy(result; signature_only=true, runner=runner, timeout=timeout)
        result = HCClosedTriangulation(n, gluings, result.signs, String(answer["isosig"]),
                                      joining_disk, corr, topology)
    end
    result
end

# Odd relabeling that reverses one tetrahedron in the positive vertex frame.
const _HC_VERTEX_SWAP = (1, 0, 2, 3)
const _HC_ORDERED_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

function _hc_swapped_vertex(vertex::Int, flip::Bool)
    flip ? _HC_VERTEX_SWAP[vertex + 1] : vertex
end

"""GKS vertex order after the positive-frame swap. Identity stays a branching only for sign +1."""
function _hc_gks_branch_order(sign::Int)
    sign == -1 ? [1, 0, 2, 3] : [0, 1, 2, 3]
end

"""Disk ids of one tetrahedron, reindexed by the odd swap (0 1)."""
function _hc_relabel_edge_disks(edges)
    by_pair = Dict{Tuple{Int,Int},Int}(
        minmax(_HC_VERTEX_SWAP[a + 1], _HC_VERTEX_SWAP[b + 1]) => Int(edges[i])
        for (i, (a, b)) in enumerate(_HC_ORDERED_EDGES))
    [by_pair[pair] for pair in _HC_ORDERED_EDGES]
end

"""Edge-disk map in the same vertex numbers as the framed gluings."""
function _hc_framed_correspondence(correspondence, pre_signs::Vector{Int})
    corr = deepcopy(correspondence)
    edges = get(corr, "tetrahedron_edges_to_disks", nothing)
    if edges !== nothing
        length(edges) == length(pre_signs) || throw(HCTriangulationError(
            :correspondence_failure, "edge disks do not match tetrahedron signs"))
        corr["tetrahedron_edges_to_disks"] = [
            pre_signs[i] == -1 ? _hc_relabel_edge_disks(edges[i]) : collect(Int, edges[i])
            for i in eachindex(pre_signs)]
    end
    corr
end

"""Move every negative tetrahedron into the positive vertex frame.

SnapPy and Regina read vertices 0,1,2,3 as positive and do not read `signs`.
Conjugating each negative tetrahedron by the odd swap (0 1) keeps the signed
orientation and writes it into the permutations. The returned signs are all +1.
"""
function _hc_positive_vertex_frame(gluings, signs::Vector{Int})
    all(s -> s == 1 || s == -1, signs) || throw(HCTriangulationError(
        :invalid_sign, "tetrahedron sign must be ±1"))
    if all(==(1), signs)
        return gluings, signs
    end
    flips = [sign == -1 for sign in signs]
    framed = map(gluings) do g
        0 <= g.tet < length(flips) && 0 <= g.to_tet < length(flips) || throw(HCTriangulationError(
            :invalid_face, "tetrahedron index out of range"))
        flip_source = flips[g.tet + 1]
        flip_target = flips[g.to_tet + 1]
        perm = g.permutation
        new_perm = ntuple(i -> _hc_swapped_vertex(
            perm[_hc_swapped_vertex(i - 1, flip_source) + 1], flip_target), 4)
        omitted = _hc_swapped_vertex(g.omitted, flip_source)
        to_omitted = _hc_swapped_vertex(g.to_omitted, flip_target)
        new_perm[omitted + 1] == to_omitted || throw(HCTriangulationError(
            :invalid_permutation, "positive vertex frame broke the omitted-vertex image"))
        merge(g, (; omitted=omitted, to_omitted=to_omitted, permutation=new_perm))
    end
    framed, fill(1, length(signs))
end

function _hc_finish_frame(result::HCClosedTriangulation, pre_signs::Vector{Int};
                          signature::Bool=true, runner=nothing, timeout::Real=120)
    # F2. Orientation is the positive vertex frame. branch_orders is only the
    # transported GKS branching, so hc_triangulation_to_odata can reuse it.
    corr = _hc_framed_correspondence(result.correspondence, pre_signs)
    orders = [_hc_gks_branch_order(sign) for sign in pre_signs]
    gluings, signs = _hc_positive_vertex_frame(result.gluings, pre_signs)
    isosig = nothing
    if signature
        framed = HCClosedTriangulation(result.tetrahedra, gluings, signs, nothing,
            result.joining_disk, corr, result.topology, orders)
        answer = hc_snappy(framed; signature_only=true, runner=runner, timeout=timeout)
        isosig = String(answer["isosig"])
    end
    HCClosedTriangulation(result.tetrahedra, gluings, signs, isosig,
        result.joining_disk, corr, result.topology, orders)
end

"""Convert a resolved HC via the existing Koda conversion and checked dual tetrahedra.

An explicit disk is passed through unchanged. With no disk and no HC disk line,
`hc_to_odata` selects the candidate of minimum region id. An empty candidate
list and the other Koda domain exclusions keep their existing exception types.
No diagram is altered.
Use `signature=false` for the Julia-only combinatorial part.
"""
function hc_triangulation(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing,
                          disk_id=nothing, signature::Bool=true, runner=nothing, timeout::Real=120)
    conversion = hc_to_odata(hc; joining_disk=joining_disk, disk_name=disk_name, disk_id=disk_id)
    v = conversion.value
    v isa parentmodule(@__MODULE__).VirtualLink || throw(HCTriangulationError(
        :unexpected_koda_output, "hc_to_odata must return VirtualLink"))
    selected = get(conversion.choices, "joining_disk_selection", nothing)
    result = _hc_closed_triangulation(v; joining_disk=selected,
        correspondence=Dict{String,Any}("koda"=>conversion.correspondence,
                                        "koda_choices"=>conversion.choices),
        signature=false, runner=runner, timeout=timeout)
    # F2. Census 19 dual signs, written in the positive vertex frame, match the
    # source triangulation. The K08 left-hand minus stays at R1, not here.
    _hc_finish_frame(result, result.signs; signature=signature, runner=runner, timeout=timeout)
end

# The public HC type is HeegaardCode; the legacy local HeegaardDiagram alias
# is deliberately not used for dispatch.
@eval function hc_triangulation(v::$(parentmodule(@__MODULE__).VirtualLink);
                               signature::Bool=true, runner=nothing, timeout::Real=120)
    result = _hc_closed_triangulation(v; signature=false, runner=runner, timeout=timeout)
    _hc_finish_frame(result, result.signs; signature=signature, runner=runner, timeout=timeout)
end

"""The pinned optional backend command; construction does not install anything."""
function hc_snappy_runner()
    uv = Sys.which("uv")
    uv === nothing && throw(SnapPyUnavailable("uv executable not found"))
    `$uv run --offline --no-project --python 3.12 --with snappy==3.3.2 python`
end

"""Run the Python bridge once, with a deadline and JSON errors.

`runner` may be an explicit `uv run ... python` Cmd for an already provisioned
environment. The default is offline and never installs from the network.
"""
function _hc_snappy_call(payload; signature_only=false, references=String[], runner=nothing,
                         timeout::Real=120, extra_args=String[])
    timeout > 0 || throw(ArgumentError("timeout must be positive"))
    command = runner === nothing ? hc_snappy_runner() : runner
    command isa Cmd || throw(ArgumentError("runner must be a Cmd ending in a Python executable"))
    script = normpath(joinpath(@__DIR__, "..", "..", "scripts", "heegaard", "hc_snappy.py"))
    mktempdir() do dir
        input, output, diagnostic = (joinpath(dir,x) for x in ("input.json","output.json","stderr.txt"))
        write(input, JSON3.write(payload))
        args = String["--input", input, "--timeout", string(timeout)]
        signature_only && push!(args, "--isosig-only")
        for reference in references; append!(args, ["--reference", String(reference)]); end
        append!(args, extra_args)
        process = try
            run(pipeline(`$command $script $args`, stdout=output, stderr=diagnostic); wait=false)
        catch e
            throw(SnapPyUnavailable(sprint(showerror,e)))
        end
        if timedwait(() -> process_exited(process), timeout+15; pollint=0.1) == :timed_out
            kill(process) # only the process created by this invocation
            wait(process)
            throw(SnapPyExecutionError(:timeout, "deadline $(timeout)s exceeded"))
        end
        wait(process)
        raw = read(output,String)
        err = read(diagnostic,String)
        result = try
            JSON3.read(raw, Dict{String,Any})
        catch
            if !success(process)
                throw(SnapPyUnavailable("runner failed before a JSON reply: " * err))
            end
            throw(SnapPyExecutionError(:invalid_json, raw * "\n" * err))
        end
        if !get(result,"ok",false)
            error_data = get(result,"error",Dict{String,Any}())
            get(error_data,"type","") == "SnapPyUnavailable" &&
                throw(SnapPyUnavailable(string(get(error_data,"message","import failed"))))
            throw(SnapPyExecutionError(Symbol(get(error_data,"type","backend_error")), string(error_data)))
        end
        success(process) || throw(SnapPyExecutionError(:exit_status, err))
        result
    end
end

hc_snappy(t::HCClosedTriangulation; kwargs...) = _hc_snappy_call(hc_triangulation_data(t); kwargs...)

"""One Python invocation from HC: triangulate in Julia, then analyze and obtain isoSig."""
function hc_snappy(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, kwargs...)
    t = hc_triangulation(hc; joining_disk=joining_disk, disk_name=disk_name, disk_id=disk_id,
                         signature=false)
    hc_snappy(t; kwargs...)
end

"""Export a catalog manifold's finite triangulation and dual spine, without inventing branching."""
function hc_snappy_catalog(name::AbstractString; kwargs...)
    _hc_snappy_call(Dict("catalog_name"=>String(name)); extra_args=["--export-catalog"], kwargs...)
end

function hc_snappy_compare(left::HCClosedTriangulation, right::HCClosedTriangulation; kwargs...)
    _hc_snappy_call(Dict("left"=>hc_triangulation_data(left), "right"=>hc_triangulation_data(right));
                    extra_args=["--compare"], kwargs...)
end

include("hc_triangulation_reverse.jl")
include("hc_triangulation_lc_cut.jl")
include("hc_triangulation_standard.jl")
export hc_find_branching, hc_triangulation_to_odata, hc_triangulation_to_heegaard_koda
