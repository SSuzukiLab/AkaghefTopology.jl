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

function load_koda_local_tables()
    cat = load_koda_local_catalog()
    vertex = _koda_read_json(cat["vertex_models"])
    edge = _koda_read_json(cat["edge_gluings"])
    orient = _koda_read_json(cat["orientation"])
    bw = haskey(cat, "bw_boundary") ? _koda_read_json(cat["bw_boundary"]) : nothing
    csel = haskey(cat, "c_selection") ? _koda_read_json(cat["c_selection"]) : nothing
    (catalog=cat, vertex=vertex, edge=edge, orientation=orient, bw_boundary=bw, c_selection=csel)
end

"""PD vertex sign σ from candidate A (β over, α under) and local frame det."""
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
