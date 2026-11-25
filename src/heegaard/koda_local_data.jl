# Source: docs/specs/heegaard-koda-conversion.md §7.0–7.2;
# data/heegaard/koda/*.json.
"""Load Koda local incidence tables; interpret by key lookup, not per-type methods."""

function _koda_data_dir()
    joinpath(dirname(@__DIR__), "..", "data", "heegaard", "koda") |> abspath
end

function load_koda_local_catalog()
    path = joinpath(_koda_data_dir(), "catalog.json")
    isfile(path) || throw(KodaImplementationError(:load_koda_catalog, "missing catalog.json", Dict("path" => path)))
    JSON3.read(read(path, String), Dict{String, Any})
end

function _koda_read_json(rel::AbstractString)
    path = joinpath(_koda_data_dir(), rel)
    isfile(path) || throw(KodaImplementationError(:load_koda_data, "missing data file", Dict("path" => path)))
    JSON3.read(read(path, String), Dict{String, Any})
end

"""
    load_koda_local_tables() -> NamedTuple

Load the installed local incidence and orientation tables for Koda conversion.

# Arguments
No arguments. The catalog and referenced JSON files must be installed and
readable. Missing files raise `KodaImplementationError` naming the path;
JSON, I/O and missing-catalog-key errors propagate.

# Returns
A named tuple with `catalog`, `vertex`, `edge`, `orientation`, `bw_boundary`
and `c_selection`. The last two fields are `nothing` if the catalog omits
their references. Field schemas are owned by
[the local-table specification, §7.0](../specs/heegaard-koda-conversion.md#70-専用データと共通の解釈手順).

# Guarantees
1. Each supplied field contains the data read from its catalog reference.
   [proved/agent: by construction]
2. Files are not modified and each call reads fresh table containers.
   [proved/agent: by construction]

# Not guaranteed
Loading does not validate every field, resolve pending entries, or certify
literature alignment. In particular, implementation readiness and source
alignment are distinct statuses in the catalog and §13.3.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§7.0, 13.3 — installed table schema and evidence boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> load_koda_local_tables().catalog["schema"]
"akaghef.heegaard.koda-local/2"
```

See also [`koda_table_decode`](@ref), [`koda_hc_from_boundary`](@ref).

Role: access
"""
function load_koda_local_tables()
    cat = load_koda_local_catalog()
    vertex = _koda_read_json(cat["vertex_models"])
    edge = _koda_read_json(cat["edge_gluings"])
    orient = _koda_read_json(cat["orientation"])
    bw = haskey(cat, "bw_boundary") ? _koda_read_json(cat["bw_boundary"]) : nothing
    csel = haskey(cat, "c_selection") ? _koda_read_json(cat["c_selection"]) : nothing
    (catalog=cat, vertex=vertex, edge=edge, orientation=orient, bw_boundary=bw, c_selection=csel)
end

"""
    koda_vertex_sign_from_frame(epsilon::Int, dalpha::Int, dbeta::Int) -> Int

Translate a local Heegaard frame and strand directions to the adopted PD vertex sign.

# Arguments
The three signs are the frame and strand-direction data of [the adopted local convention](../specs/heegaard-koda-conversion.md#73-向きの規約と局所符号). Each must be `1` or `-1`; otherwise `KodaInputError` with stage `:vertex_sign` is thrown. The caller establishes their geometric meaning.

# Returns
The PD vertex sign `1` or `-1` in the adopted convention.

# Guarantees
1. The returned sign is the local frame translation of the cited convention. [proved/agent: direct evaluation of the convention's sign relation]
2. The result depends only on the three supplied signs and has no side effects. [proved/agent: by construction]

# Not guaranteed
This operation does not recover a frame from a diagram, check global strand consistency, or establish a global Koda-conversion identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §7.3 — adopted local frame and PD sign translation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> (H.koda_vertex_sign_from_frame(1,1,1), H.koda_vertex_sign_from_frame(1,-1,1))
(-1, 1)
```

Role: construct
"""
function koda_vertex_sign_from_frame(epsilon::Int, dalpha::Int, dbeta::Int)
    epsilon in (-1, 1) && dalpha in (-1, 1) && dbeta in (-1, 1) ||
        throw(KodaInputError(:vertex_sign, "signs must be ±1", Dict()))
    -epsilon * dalpha * dbeta
end

function koda_passage_over_under_tables(vertex::AbstractDict)
    ports = vertex["pd_ports_ccw_from_incoming_under"]
    branching = vertex["branching"]
    Dict(
        "ports" => ports,
        "over_ports" => Int.(branching["over_passage_ports"]),
        "under_ports" => Int.(branching["under_passage_ports"]),
        "over_family" => "beta",
        "under_family" => "alpha",
    )
end

function koda_pending_table_fields(vertex::AbstractDict)
    pending = String[]
    nulls = get(vertex, "null_fields", nothing)
    nulls === nothing && return pending
    for (k, v) in nulls
        v === nothing && push!(pending, String(k))
    end
    pending
end
