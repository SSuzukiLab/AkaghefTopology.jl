# cf: lib/topology.jl/src/ograph_tetrahedra.jl::disk_tetrahedra — consume its checked face maps, without a second port convention
# cf: lib/topology.jl/docs/specs/heegaard-snappy.md — bridge contract and indexing

"""
    HCTriangulationError(reason, detail) -> HCTriangulationError

Store a rejected check in the HC finite-triangulation conversion API.

# Arguments
`reason::Symbol` and `detail::String`, or convertible values. Conversion
errors propagate. Construction neither validates the tag nor throws the
resulting exception.

# Returns
`HCTriangulationError` with the reason tag and diagnostic text. `showerror`
displays both fields.

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. When emitted by a triangulation operation, the record concerns a failed
   check in [the finite-triangulation bridge](../specs/heegaard-snappy.md).
   It is diagnostic data, not a triangulation or a manifold invariant.
   [proved/agent: by construction]

# Not guaranteed
The constructor does not verify that the diagnostic is accurate. This type
alone does not distinguish a mathematical domain failure from a consistency
failure; inspect the reason and detail supplied by the operation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = HCTriangulationError(:branching, "invalid order"); (e.reason, e.detail)
(:branching, "invalid order")
```

See also [`hc_triangulation`](@ref), [`hc_find_branching`](@ref).

Role: construct
"""
struct HCTriangulationError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::HCTriangulationError) = print(io, "HCTriangulationError(", e.reason, "): ", e.detail)

"""
    SnapPyUnavailable(detail) -> SnapPyUnavailable

Store a diagnostic that the requested SnapPy backend is unavailable.

# Arguments
`detail::String`, or a value convertible to `String`. Conversion errors
propagate. No backend check is performed by this constructor.

# Returns
`SnapPyUnavailable` with diagnostic text in `detail`; `showerror` displays
that text with the exception name.

# Guarantees
1. The supplied text is retained after conversion and is not modified.
   [proved/agent: by construction]
2. The result is a diagnostic for the [SnapPy bridge](../specs/heegaard-snappy.md),
   with no mathematical invariant attached. Constructing it neither launches
   a backend nor throws the resulting value. [proved/agent: by construction]

# Not guaranteed
A manually constructed value does not establish the actual availability of
SnapPy, and does not assert that any input manifold is invalid.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = SnapPyUnavailable("backend unavailable"); e.detail
"backend unavailable"
```

See also [`SnapPyExecutionError`](@ref).

Role: construct
"""
struct SnapPyUnavailable <: Exception
    detail::String
end
Base.showerror(io::IO, e::SnapPyUnavailable) = print(io, "SnapPyUnavailable: ", e.detail)

"""
    SnapPyExecutionError(reason, detail) -> SnapPyExecutionError

Store a reason and diagnostic text for a failed SnapPy bridge operation.

# Arguments
`reason::Symbol` and `detail::String`, or values convertible to those types.
Conversion errors propagate. Arbitrary reason tags are accepted.

# Returns
`SnapPyExecutionError` with the reason tag and diagnostic text; `showerror`
displays both fields.

# Guarantees
1. Fields preserve the supplied values after conversion; arguments are not
   modified. [proved/agent: by construction]
2. The result records a failure in the scope of the
   [SnapPy bridge](../specs/heegaard-snappy.md). Construction invokes no
   backend and neither computes an invariant nor throws the resulting value.
   [proved/agent: by construction]

# Not guaranteed
Construction does not verify the diagnostic. A backend failure is not a
proof of any mathematical property of the input manifold.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = SnapPyExecutionError(:process, "nonzero exit"); (e.reason, e.detail)
(:process, "nonzero exit")
```

See also [`SnapPyUnavailable`](@ref).

Role: construct
"""
struct SnapPyExecutionError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::SnapPyExecutionError) = print(io, "SnapPyExecutionError(", e.reason, "): ", e.detail)

"""
    HCClosedTriangulation(tetrahedra, gluings, signs, isosig, joining_disk, correspondence, topology, branch_orders=Vector{Vector{Int}}()) -> HCClosedTriangulation

Store a finite triangulation with orientation, optional branching and HC correspondence data.

# Arguments
`tetrahedra::Int`, `gluings::Vector{NamedTuple}`, `signs::Vector{Int}`,
`isosig::Union{Nothing,String}`, arbitrary `joining_disk`,
`correspondence::Dict{String,Any}`, `topology::NamedTuple`, and
`branch_orders::Vector{Vector{Int}}`, or convertible values where accepted.
Conversion errors propagate. The constructor does not check face pairings,
orientation, branching, closedness or a claimed signature.

# Returns
An `HCClosedTriangulation` holding the supplied finite gluing record and
annotations in [the triangulation interface](../specs/heegaard-snappy.md).

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers and
   the supplied joining-disk value are retained without modification.
   [proved/agent: by construction]
2. The optional signature and branching fields can be absent; construction
   performs no backend call or triangulation validation.
   [proved/agent: by construction]

# Not guaranteed
A hand-built value need not triangulate a closed 3-manifold. Its type and
stored topology fields do not certify the HC correspondence or its isoSig.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — finite closed triangulation data and optional backend boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> t = HCClosedTriangulation(0, NamedTuple[], Int[], nothing, nothing, Dict{String,Any}(), NamedTuple()); (t.tetrahedra, t.isosig, isempty(t.branch_orders))
(0, nothing, true)
```

See also [`hc_triangulation`](@ref).

Role: convert
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

"""
    hc_triangulation_data(t::HCClosedTriangulation) -> Dict{String,Any}

Return the bridge serialization of a closed triangulation and its source correspondences.

# Arguments
`t` is a record in [the bridge format](../specs/heegaard-snappy.md).
No mathematical validity check is performed; a caller-created record must
have consistent fields for consumers that interpret it as a triangulation.

# Returns
A string-keyed dictionary with schema, tetrahedra, gluing dictionaries, signs,
isoSig, joining disk, branching orders, correspondence and topology data.

# Guarantees
1. The dictionary preserves the record's declared triangulation and
   correspondence fields in bridge form. [proved/agent: by construction]
2. The record is not modified. Nested signs, branching orders, joining disk
   and correspondence data can be shared. [proved/agent: by construction]

# Not guaranteed
This is not a deep copy or a manifold-validity certificate. Mutating shared
nested values can change the original record.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — explicit finite-triangulation bridge data.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> t = H.HCClosedTriangulation(0,NamedTuple[],Int[],nothing,nothing,Dict{String,Any}(),(;)); x = H.hc_triangulation_data(t);

julia> (x["tetrahedra"], x["isosig"], x["signs"] === t.signs)
(0, nothing, true)
```

Role: convert
"""
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

"""
    hc_triangulation(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, signature=true, runner=nothing, timeout=120) -> HCClosedTriangulation

Return the oriented closed triangulation constructed from an HC's Koda o-graph.

# Arguments
`hc` and disk-selection keywords have the
[forward-conversion domain](../specs/heegaard-snappy.md#入力と定義域) of `hc_to_odata`;
those input and domain exceptions propagate. Dual gluing, orientation and
capped-link checks raise `HCTriangulationError` with their stage.
`signature=true` requests the optional backend with the runner/timeout domain
of `hc_snappy`; `signature=false` requires no backend.

# Returns
`HCClosedTriangulation` with tetrahedra, oriented face gluings and signs,
optional `isosig`, selected joining disk, source correspondences, topology
checks and transported branching orders.

# Guarantees
1. The returned gluing record passes the closed normal-o-graph dual and
   capped sphere-link checks in the linked specification.
   [proved/agent: by construction]
2. Its chosen disk and source correspondences are retained, and the HC is
   not modified. [proved/agent: by construction]

# Not guaranteed
The computation does not independently prove the whole HC-to-manifold PL
identification. A branching order and an orientation sign are distinct data;
an unoriented isoSig does not identify them.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — forward dual tetrahedra and oriented-frame convention.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(-1,Dict(),Dict(),Dict(),Dict(),Dict(),nothing); hc = H.HeegaardCode(lc,Dict(),Dict(),Dict(),Dict(),nothing);

julia> try H.hc_triangulation(hc; signature=false) catch e; e isa H.KodaInputError end
true
```

Role: convert
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
@eval @doc """
    hc_triangulation(v::VirtualLink; signature=true, runner=nothing, timeout=120) -> HCClosedTriangulation

Return the oriented closed triangulation dual to a supported normal o-graph.

# Arguments
`v` must meet the closed normal-o-graph and dual-gluing conditions of
[the triangulation bridge](../specs/heegaard-snappy.md).
Dual construction and topology failures raise their validation exceptions
or `HCTriangulationError`. With `signature=true`, `runner` and positive
`timeout` have the domain of `hc_snappy`; backend errors propagate.

# Returns
`HCClosedTriangulation` with dual tetrahedra, oriented face gluings, source
correspondence and transported branching orders; `isosig` is absent when
`signature=false`.

# Guarantees
1. Successful return passes the dual face-pairing, orientation, Euler and
   capped sphere-link checks. [proved/agent: by construction]
2. The input is not modified. `signature=false` does not launch the backend.
   [proved/agent: by construction]

# Not guaranteed
The record is not an independently checked PL map to an externally named
manifold. Branching and orientation remain separate from unoriented isoSig.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — checked o-graph dual and positive vertex frame.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> v = Topology.lens_space_link(2,1); t = H.hc_triangulation(v; signature=false);

julia> (t.tetrahedra == length(Topology.real_gauss_code(v)[2]), length(t.gluings) == 2t.tetrahedra, t.isosig === nothing)
(true, true, true)
```

Role: convert
"""
function hc_triangulation(v::$(parentmodule(@__MODULE__).VirtualLink);
                               signature::Bool=true, runner=nothing, timeout::Real=120)
    result = _hc_closed_triangulation(v; signature=false, runner=runner, timeout=timeout)
    _hc_finish_frame(result, result.signs; signature=signature, runner=runner, timeout=timeout)
end

"""
    hc_snappy_runner() -> Cmd

Return the command selecting the pinned offline SnapPy Python environment.

# Arguments
`uv` must be on `PATH`, checked with `SnapPyUnavailable` if absent.
The caller supplies the already provisioned environment recorded in
`docs/ENVIRONMENT.md`; this function does not check the package cache.

# Returns
A `Cmd` for the existing Python 3.12 environment with SnapPy 3.3.2.

# Guarantees
1. The returned command requests the pinned backend in offline mode.
   [proved/agent: by construction]
2. No process is launched and no environment is modified.
   [proved/agent: by construction]

# Not guaranteed
A command value does not establish backend availability or compute a
mathematical invariant. Executing it can fail on an unprovisioned machine.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — external-runner boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.hc_snappy_runner() isa Cmd catch e; e isa H.SnapPyUnavailable end
true
```

Role: other
"""
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

"""
    hc_snappy(t::HCClosedTriangulation; kwargs...) -> Dict{String,Any}

Return the optional SnapPy bridge's analysis of a closed triangulation.

# Arguments
`t` must be a valid record in [the bridge domain](../specs/heegaard-snappy.md#入力と定義域);
manual construction does not certify it. Keywords include an existing Python
`runner::Cmd`, positive `timeout`, `signature_only` and `references`.
Invalid runner/timeout raises `ArgumentError`; availability and execution
failures raise `SnapPyUnavailable` or `SnapPyExecutionError` with diagnostics.

# Returns
A string-keyed backend reply with `ok` and isoSig information; further
analysis and comparison fields depend on the requested bridge mode.

# Guarantees
1. The reply is accepted only after the bridge reports success and the
   child process exits successfully. [proved/agent: by construction]
2. The triangulation is not modified. The call launches a child process and
   uses temporary files. [proved/agent: by construction]

# Not guaranteed
The evidence level is that of each backend field. Numerical isometry data
are not an independently checked PL certificate, and isoSig comparison does
not by itself preserve a chosen orientation or branching.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — backend analysis, evidence and orientation boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> t = H.HCClosedTriangulation(0, NamedTuple[], Int[], nothing, nothing, Dict{String,Any}(), (;));

julia> try H.hc_snappy(t; timeout=0) catch e; e isa ArgumentError end
true
```

Role: other
"""
hc_snappy(t::HCClosedTriangulation; kwargs...) = _hc_snappy_call(hc_triangulation_data(t); kwargs...)

"""
    hc_snappy(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, kwargs...) -> Dict{String,Any}

Return the optional SnapPy analysis of the triangulation obtained from an HC.

# Arguments
`hc` and disk-selection keywords satisfy `hc_triangulation` and
[the forward bridge domain](../specs/heegaard-snappy.md#入力と定義域).
Koda and triangulation errors propagate before the backend call. Remaining
keywords have the domain of `hc_snappy(::HCClosedTriangulation)`; invalid
runner/timeout and backend errors retain those exception types.

# Returns
The bridge dictionary for the resulting closed triangulation, including
isoSig and any requested backend analysis fields.

# Guarantees
1. The reply analyzes the triangulation returned by the same HC conversion
   with `signature=false`. [proved/agent: by construction]
2. The HC is not modified; successful conversion can be followed by one
   child-process invocation and temporary-file writes.
   [proved/agent: by construction]

# Not guaranteed
The reply does not upgrade the forward conversion's geometric identification
status. Backend invariants and unoriented identities are subject to their
individual evidence limits.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — HC forward conversion and backend evidence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(-1,Dict(),Dict(),Dict(),Dict(),Dict(),nothing); hc = H.HeegaardCode(lc,Dict(),Dict(),Dict(),Dict(),nothing);

julia> try H.hc_snappy(hc) catch e; e isa H.KodaInputError end
true
```

Role: other
"""
function hc_snappy(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, kwargs...)
    t = hc_triangulation(hc; joining_disk=joining_disk, disk_name=disk_name, disk_id=disk_id,
                         signature=false)
    hc_snappy(t; kwargs...)
end

"""
    hc_snappy_catalog(name::AbstractString; kwargs...) -> Dict{String,Any}

Request a catalog manifold's finite triangulation and dual spine from SnapPy.

# Arguments
`name` must be accepted by the catalog backend in
[the bridge domain](../specs/heegaard-snappy.md).
Keywords include an existing Python `runner::Cmd` and positive `timeout`;
invalid runner/timeout raises `ArgumentError`. Unknown catalog entries,
unavailable backends, malformed replies and execution errors propagate as
`SnapPyUnavailable` or `SnapPyExecutionError`. `extra_args` is reserved.

# Returns
A string-keyed bridge dictionary containing the exported finite
triangulation, topology, dual-spine information and backend identity data.

# Guarantees
1. A successful return is the backend's catalog-export reply for `name`.
   [proved/agent: by construction]
2. The call does not construct a branching from the manifold name.
   [proved/agent: by construction]
3. No input object is modified; a child process and temporary files are used.
   [proved/agent: by construction]

# Not guaranteed
Catalog naming and backend analysis alone do not identify a new HC with this
manifold. Branching, combing and a resolved HC remain separate constructions.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — catalog export and explicit branching boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.hc_snappy_catalog("S3"; timeout=0) catch e; e isa ArgumentError end
true
```

Role: other
"""
function hc_snappy_catalog(name::AbstractString; kwargs...)
    _hc_snappy_call(Dict("catalog_name"=>String(name)); extra_args=["--export-catalog"], kwargs...)
end

"""
    hc_snappy_compare(left::HCClosedTriangulation, right::HCClosedTriangulation; kwargs...) -> Dict{String,Any}

Ask the optional SnapPy bridge to compare the unoriented manifolds of two closed triangulations.

# Arguments
`left` and `right` are closed finite triangulations in the
[bridge domain](../specs/heegaard-snappy.md#入力と定義域); use `hc_triangulation`
to obtain checked records. Direct record construction does not certify them.
The optional `runner` is a `Cmd` ending in a Python executable with the bridge
dependencies available, and `timeout::Real=120` must be positive. Invalid
timeout or runner types raise `ArgumentError`. Missing environments or failure
before a JSON reply raise `SnapPyUnavailable`; malformed replies, timeouts
and reported backend failures raise `SnapPyExecutionError` with the reason.
Other bridge keywords are forwarded; `extra_args` is reserved by this wrapper.

# Returns
A string-keyed backend dictionary. In ordinary comparison mode it contains
`ok`, `same_manifold` (`true` or `nothing`), `evidence`, and `orientation`;
numerical comparison also supplies `is_isometric_to` with its status/value
or diagnostic. Use the evidence field to distinguish combinatorial identity
from a numerical isometry result.

# Guarantees
1. The ordinary comparison requests the unoriented comparison specified in
   [identification limits](../specs/heegaard-snappy.md#同定の限界).
   A failed numerical identification is returned as unresolved, not as a
   proof of distinct manifolds. [proved/agent: by construction]
2. The input records are not modified. The operation can launch an external
   process and create temporary files. [proved/agent: by construction]

# Not guaranteed
Numerical isometry is not a rigorous oriented homeomorphism certificate.
Successful execution and backend availability depend on the supplied environment.
Nonstandard forwarding options such as `signature_only=true` can select a
different backend reply schema.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> t = hc_triangulation(Topology.lens_space_link(2, 1); signature=false);

julia> try hc_snappy_compare(t, t; timeout=0); false catch e; e isa ArgumentError end
true
```

See also [`hc_triangulation`](@ref), [`hc_snappy`](@ref).

Role: condition
"""
function hc_snappy_compare(left::HCClosedTriangulation, right::HCClosedTriangulation; kwargs...)
    _hc_snappy_call(Dict("left"=>hc_triangulation_data(left), "right"=>hc_triangulation_data(right));
                    extra_args=["--compare"], kwargs...)
end

include("hc_triangulation_reverse.jl")
include("hc_triangulation_lc_cut.jl")
include("hc_triangulation_standard.jl")
export hc_find_branching, hc_triangulation_to_odata, hc_triangulation_to_heegaard_koda
