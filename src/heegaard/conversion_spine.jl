"""K08 section 5: punctured Heegaard diagram -> singular graph -> normal o-graph.

Scope: cellular diagrams with an explicitly selected embedded joining disk.
The intermediate complex retains oriented sector attaching maps and gluing witnesses.
See docs/heegaard/joining-disk.md for conventions and the finite smoke test.
"""

_labelled_hash(H) = bytes2hex(SHA.sha256(_canonical_json(_json_dict(H))))

function _disk_json(D::JoiningDisk)
    Dict{String,Any}("labelled_diagram_hash"=>D.labelled_diagram_hash,
        "region"=>D.region,"boundary_darts"=>copy(D.boundary_darts),
        "boundary_orientation"=>"surface induced; disk on left")
end

"""
    joining_disk(H::CellularHeegaardDiagram, region::Integer) -> JoiningDisk

Select and validate a joining disk in a legacy cellular diagram.

# Arguments
`region` is a region index in the [legacy cellular representation](../specs/heegaard-terminology.md#legacy-cellular-records). Out-of-range and multiple-boundary regions raise `ArgumentError`; the selected disk must then pass `validate_joining_disk`, whose errors propagate. Hashing and incidence errors also propagate.

# Returns
A `JoiningDisk` bound to the labelled diagram, storing the selected region and a copied boundary-dart vector.

# Guarantees
1. The result passes `validate_joining_disk(H, result)` and has its implemented [joining-disk](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk) interpretation. [proved/agent: by construction]
2. The input diagram is not modified. [proved/agent: by construction]

# Not guaranteed
No full Heegaard cut-system certification, choice independence of the resulting branching, or general correctness of a subsequent Koda conversion is established by selecting this disk.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.2 — joining disk and its conversion domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); H.joining_disk(d,1).region
1
```

Role: convert
"""
function joining_disk(H::HeegaardDiagram, region::Integer)
    blocks = H.code.regions.boundary_cycle_blocks
    1 <= region <= length(blocks) || throw(ArgumentError("joining disk region is out of range"))
    length(blocks[region]) == 1 || throw(ArgumentError("joining region must be a disk, not a multi-boundary region"))
    D = JoiningDisk(_labelled_hash(H),Int(region),copy(boundary_cycles(H)[only(blocks[region])]))
    validate_joining_disk(H,D)
    D
end

"""
    validate_joining_disk(H::CellularHeegaardDiagram, D::JoiningDisk) -> Bool

Check a labelled legacy diagram's selected joining disk.

# Arguments
`H` and `D` must satisfy the [joining-disk conditions](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk)
in the [legacy cellular representation](../specs/heegaard-terminology.md#legacy-cellular-records).
Checked failures raise `ArgumentError` identifying an identity mismatch,
invalid diagram, out-of-range or nondisk region, boundary mismatch,
nonsimple boundary, unequal curve counts, or a failed curve-contact check.
Incidence and serialization errors propagate.

# Returns
`true` when the disk's implemented joining checks all succeed.

# Guarantees
1. The disk is bound to the supplied labelled diagram and its recorded
   region passes the current embedded-boundary and joining tests.
   [proved/agent: by construction]
2. Neither argument is modified. [proved/agent: by construction]

# Not guaranteed
This does not perform full Heegaard cut-system certification or prove the
subsequent Koda construction's general manifold identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.2 — joining disk and curve-contact domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); H.validate_joining_disk(d, H.joining_disk(d,1))
true
```

See also [`make_punctured`](@ref).

Role: condition
"""
function validate_joining_disk(H::HeegaardDiagram,D::JoiningDisk)
    D.labelled_diagram_hash == _labelled_hash(H) ||
        throw(ArgumentError("joining disk belongs to a different labelled diagram; reselect it explicitly"))
    validate(H).valid || throw(ArgumentError("joining disk requires a valid diagram"))
    blocks = H.code.regions.boundary_cycle_blocks
    1 <= D.region <= length(blocks) || throw(ArgumentError("joining disk region is out of range"))
    length(blocks[D.region]) == 1 || throw(ArgumentError("joining region is not a disk"))
    boundary = boundary_cycles(H)[only(blocks[D.region])]
    !isempty(boundary) && _cyclic_lexmin(D.boundary_darts) == _cyclic_lexmin(boundary) ||
        throw(ArgumentError("joining disk boundary does not match its region"))
    graph = ribbon_graph(H)
    points = [graph.darts[d].origin for d in boundary]
    length(unique(points)) == length(points) ||
        throw(ArgumentError("joining disk closure has a non-simple boundary (repeated intersection)"))
    edge_ids = [graph.darts[d].edge for d in boundary]
    length(unique(edge_ids)) == length(edge_ids) ||
        throw(ArgumentError("joining disk boundary repeats an edge"))
    arcs = [(graph.edges[e].family,graph.edges[e].curve_id) for e in edge_ids]
    expected = Set((c.family,c.id) for c in _curve_vector(H.code.code))
    length(H.code.code.alpha) == length(H.code.code.beta) ||
        throw(ArgumentError("joining disk requires equal numbers of alpha and beta curves"))
    Set(arcs) == expected || throw(ArgumentError("joining disk must meet every alpha and beta curve"))
    for key in expected
        count(i -> arcs[i] == key && arcs[mod1(i-1,length(arcs))] != key,eachindex(arcs)) == 1 ||
            throw(ArgumentError("joining disk must meet curve $key in one connected arc"))
    end
    true
end

const KODA_CONVENTION = "K08 left-hand ambient convention: beta over, alpha under; " *
    "singular edges follow the induced boundary of D (opposite to K08 slope orientations); " *
    "vertex sign from repository CCW PD frame"

"""
    make_punctured(H::CellularHeegaardDiagram; joining_disk=nothing, strategy=:koda, source_identity=:canonical) -> PuncturedHeegaardDiagram

Equip a legacy Heegaard diagram with a checked joining disk for Koda conversion.

# Arguments
`H` must be in [the punctured-diagram domain, §5](../specs/heegaard-koda-conversion.md).
`joining_disk` must be an explicit `JoiningDisk`; omission or another type
raises `ArgumentError`. `validate_joining_disk` checks its labelled-diagram
identity and joining conditions, with an explanatory `ArgumentError` on failure.
Only `strategy=:koda` and `source_identity` equal to `:canonical` or `:labelled`
are accepted; other values raise `ArgumentError`. Hashing and alpha-cut errors
propagate. The caller supplies geometric hypotheses beyond those finite checks.

# Returns
`PuncturedHeegaardDiagram` retaining `H`, a deep copy of the joining disk,
alpha-cut data and a `ConversionCertificate` with `evidence_status == :unverified`.
The certificate records the chosen source identity and joining data.

# Guarantees
1. The retained disk passes the joining-disk checks for the supplied diagram.
   [proved/agent: by construction]
2. The result records the punctured diagram of §5 for subsequent Koda
   conversion, without asserting that conversion has been verified.
   [proved/agent: by construction]
3. The input diagram and supplied disk are not modified. The diagram is
   shared with the result, while the disk is independently copied.
   [proved/agent: by construction]

# Not guaranteed
This does not choose a default disk, execute spine conversion or establish
a general manifold-identification theorem for the downstream finite encoding.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5–6 — explicit joining data and the punctured domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true));

julia> ph = H.make_punctured(d; joining_disk=H.joining_disk(d, 1)); Tuple(ph.joining_disk.boundary_darts)
(1, 8)
```

See also [`joining_disk`](@ref), [`heegaard_to_spine`](@ref), [`save_punctured`](@ref).

Role: convert
"""
function make_punctured(H::HeegaardDiagram;joining_disk=nothing,strategy=:koda,source_identity=:canonical)
    strategy == :koda || throw(ArgumentError("only :koda is supported"))
    joining_disk isa JoiningDisk || throw(ArgumentError("an explicit validated JoiningDisk is required; no default region is chosen"))
    validate_joining_disk(H,joining_disk)
    D = deepcopy(joining_disk)
    source_identity in (:canonical,:labelled) || throw(ArgumentError("unsupported source identity"))
    identity = source_identity==:canonical ? canonical_hash(H) : "labelled-sha256:"*_labelled_hash(H)
    cert = ConversionCertificate(:koda,identity,_disk_json(D),Any[],
        KODA_CONVENTION,Dict{String,Any}(),nothing,Dict{String,Any}(),:unverified)
    PuncturedHeegaardDiagram(H,D,cut_along_alpha(H;diagram_hash=identity),cert)
end

_base_diagram(h::HeegaardDiagram) = h
_base_diagram(h::PuncturedHeegaardDiagram) = h.diagram

"""
    save_punctured(path::AbstractString, PH::PuncturedHeegaardDiagram) -> AbstractString

Save a legacy punctured Heegaard diagram and its checked joining disk as JSON.

# Arguments
`path` is writable and its parent directory exists; file and serialization
errors propagate. Existing contents are overwritten. `PH` must pass
`validate_joining_disk`, whose `ArgumentError` identifies a disk/diagram
mismatch or a failed [joining condition](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).

# Returns
The supplied path after writing the `akaghef.heegaard.punctured/1` record.

# Guarantees
1. The saved record contains the labelled diagram and its selected joining
   disk, for reconstruction by `load_punctured`.
   [proved/agent: by construction]
2. `PH` is not modified; the path is written. [proved/agent: by construction]

# Not guaranteed
Writing is not atomic and does not create missing directories. Serialization
does not establish additional geometric source alignment.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); ph = H.make_punctured(d; joining_disk=H.joining_disk(d,1));

julia> mktempdir() do dir; p = joinpath(dir, "diagram.json"); H.save_punctured(p, ph) == p && isfile(p) end
true
```

See also [`load_punctured`](@ref).

Role: other
"""
function save_punctured(path::AbstractString,PH::PuncturedHeegaardDiagram)
    validate_joining_disk(PH.diagram,PH.joining_disk)
    value = _json_dict(PH)
    open(path,"w") do io
        println(io,_canonical_json(value))
    end
    path
end

_json_dict(PH::PuncturedHeegaardDiagram)=Dict("schema"=>"akaghef.heegaard.punctured/1",
    "diagram"=>_json_dict(PH.diagram),"joining_disk"=>_disk_json(PH.joining_disk))

"""
    load_punctured(path::AbstractString) -> PuncturedHeegaardDiagram

Load and revalidate a saved punctured legacy Heegaard diagram.

# Arguments
`path` names a readable JSON file from `save_punctured`. File and JSON errors
propagate. Unsupported schema or boundary orientation, or invalid joining
data, raise `ArgumentError`; failure of full diagram validation raises
`ValidationError`. Missing or ill-typed fields can raise native lookup or
conversion errors.

# Returns
A reconstructed `PuncturedHeegaardDiagram` with a certified legacy diagram,
validated joining disk, recomputed alpha-cut data and an unverified conversion
certificate using the default canonical source identity.

# Guarantees
1. Successful return requires both full legacy diagram certification and
   the saved joining disk's checks. [proved/agent: by construction]
2. The input file is not modified. [proved/agent: by construction]

# Not guaranteed
The saved format does not restore a prior conversion proof, an arbitrary
source-identity choice or cached downstream spine data. It records the
punctured input of [§5](../specs/heegaard-koda-conversion.md), not an independent
manifold-identification certificate.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); ph = H.make_punctured(d; joining_disk=H.joining_disk(d, 1));

julia> mktempdir() do dir; p = joinpath(dir, "example.json"); H.save_punctured(p, ph); Tuple(H.load_punctured(p).joining_disk.boundary_darts) end
(1, 8)
```

See also [`save_punctured`](@ref), [`make_punctured`](@ref).

Role: construct
"""
function load_punctured(path::AbstractString)
    value = JSON3.read(read(path,String))
    _punctured_from_json(value)
end

function _punctured_from_json(value)
    _spec_field(value,:schema) == "akaghef.heegaard.punctured/1" ||
        throw(ArgumentError("unsupported punctured Heegaard schema"))
    H = certify_heegaard(_heegaard_from_json(value.diagram))
    d = value.joining_disk
    _spec_field(d,:boundary_orientation) == "surface induced; disk on left" ||
        throw(ArgumentError("unsupported joining disk boundary orientation"))
    D = JoiningDisk(String(d.labelled_diagram_hash),Int(d.region),Int.(d.boundary_darts))
    make_punctured(H;joining_disk=D)
end

"""
    to_flow_spine(PH::PuncturedHeegaardDiagram) -> FlowSpine

Construct the legacy flow-spine record from a cellular diagram and selected joining disk.

# Arguments
The selected disk must pass `validate_joining_disk` and the source must be cellular. Failed disk checks, noncellular input, incoherent curve directions, unexpected singular degrees, reused source edges and vertex-free cases raise `ArgumentError`. The geometric conversion domain is [Koda construction](../specs/heegaard-koda-conversion.md#6-heegaard-diagramから-branched-spine-を作る).

# Returns
A `FlowSpine` retaining `PH`, with singular edges and source-dart paths, vertex signs, a branched-complex record and a conversion certificate whose evidence status is `:unverified`.

# Guarantees
1. The result records the directed singular data and cell correspondences of the implemented punctured-diagram construction in the adopted Koda convention. [proved/agent: by construction]
2. The source data are not modified and `spine.source === PH`. [proved/agent: by construction]

# Not guaranteed
The retained certificate does not prove the general manifold or flow identification. This legacy path rejects vertex-free singular components and does not guarantee all topology fields expected by `spine_topology`.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.1–6.4 — punctured-diagram construction and source-alignment limits.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> h = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); p = H.make_punctured(h; joining_disk=H.joining_disk(h,1));

julia> s = H.to_flow_spine(p); (length(s.vertex_signs), s.source === p)
(1, true)
```

Role: convert
"""
function to_flow_spine(PH::PuncturedHeegaardDiagram)
    H,D = PH.diagram,PH.joining_disk
    validate_joining_disk(H,D)
    is_cellular(H) || throw(ArgumentError("non-cellular source needs explicit standardization before conversion"))
    graph = ribbon_graph(H)
    removed_edges = Set(graph.darts[d].edge for d in D.boundary_darts)
    corners = Set(graph.darts[d].origin for d in D.boundary_darts)
    directions = Dict{Tuple{Symbol,Int},Int}()
    for d in D.boundary_darts
        dart = graph.darts[d]; edge = graph.edges[dart.edge]
        key = (edge.family,edge.curve_id); value = isodd(d) ? 1 : -1
        haskey(directions,key) && directions[key] != value &&
            throw(ArgumentError("incoherent joining boundary on a curve"))
        directions[key] = value
    end
    retained = [graph.darts[2e.id - (directions[(e.family,e.curve_id)] == 1 ? 1 : 0)]
                for e in graph.edges if !(e.id in removed_edges)]
    outgoing = Dict{Int,Vector{Int}}(); incoming = Dict{Int,Vector{Int}}()
    for d in retained
        push!(get!(outgoing,d.origin,Int[]),d.id)
        push!(get!(incoming,d.target,Int[]),d.id)
    end
    for p in graph.vertices
        degree = p in corners ? 1 : 2
        length(get(outgoing,p,Int[])) == degree == length(get(incoming,p,Int[])) ||
            throw(ArgumentError("unexpected singular degree at intersection $p"))
    end
    vertices = sort(setdiff(graph.vertices,collect(corners)))
    isempty(vertices) && throw(ArgumentError("vertex-free spine is not supported"))
    edges = SpineEdge[]; covered = Set{Int}()
    for p in vertices, start in outgoing[p]
        path = Int[]; d = start
        while true
            d in covered && throw(ArgumentError("smoothing reused a directed source edge"))
            push!(path,d); push!(covered,d)
            target = graph.darts[d].target
            target in corners || break
            d = only(outgoing[target])
        end
        first_dart,last_dart = graph.darts[first(path)],graph.darts[last(path)]
        push!(edges,SpineEdge(p,graph.edges[first_dart.edge].family,
            last_dart.target,graph.edges[last_dart.edge].family,path))
    end
    length(covered) == length(retained) || throw(ArgumentError("vertex-free singular component is not supported"))
    # beta=over, alpha=under. Repo sigma = det(over,under), hence the minus.
    signs = Dict(p.id => -p.sign * directions[(:alpha,p.alpha)] * directions[(:beta,p.beta)]
        for p in H.code.code.intersections if p.id in vertices)
    audit = Dict{String,Any}("removed_boundary_edges"=>sort(collect(removed_edges)),
        "smoothed_intersections"=>sort(collect(corners)),"remaining_intersections"=>vertices,
        "curve_directions"=>[Dict("family"=>String(k[1]),"curve_id"=>k[2],"direction"=>v)
            for (k,v) in sort(collect(directions);by=first)],
        "edge_paths"=>[Dict("source"=>[e.source_vertex,String(e.source_family)],
            "target"=>[e.target_vertex,String(e.target_family)],"source_darts"=>e.source_darts) for e in edges])
    complex = _koda_complex(H,D,edges,directions)
    audit["branched_complex"] = _complex_json(complex)
    # D's labelled digest was rechecked above, so its immutable-by-convention
    # source certificate remains applicable. Avoid repeating canonical search.
    cert = ConversionCertificate(:koda,PH.certificate.source_hash,_disk_json(D),Any[],
        KODA_CONVENTION,audit,nothing,Dict{String,Any}(),:unverified)
    FlowSpine(surface_genus(H),edges,signs,PH,cert,complex)
end

"""
    to_ograph(spine::FlowSpine) -> OData

Encode the directed singular strands and vertex signs of a flow spine as o-data.

# Arguments
`spine` supplies the consistent strand, sign and correspondence records of
[Koda encoding §7](../specs/heegaard-koda-conversion.md#7-branched-spine-と-o-data-の符号化復号).
Checked failures raise `ArgumentError` when singular edges do not define
oriented component cycles or a component does not close. Malformed metadata
and missing correspondence entries can raise lookup errors. The caller
establishes the geometric spine domain; use `spine_to_odata` when its sector
and vertex-free guards are required.

# Returns
An `OData` with signed Gauss words, orientation vector and conversion metadata,
including the source identity and strand/cell correspondences.

# Guarantees
1. The returned words encode the supplied closed strand cycles and vertex
   signs in the [fixed Koda convention](../specs/heegaard-koda-conversion.md#73-向きの規約と局所符号).
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This method does not separately check closed-normal C1–C3 or sector disk
conditions. Its metadata does not certify manifold identification; the
[Koda source-alignment limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合) remain.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§7.1–7.4 — strand encoding and orientation conventions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); ph = odata_to_heegaard_koda(g).value; s = heegaard_to_spine(ph).value;

julia> compare_odata_words(g, to_ograph(s))["equal_up_to_metadata"]
true
```

See also [`spine_to_odata`](@ref), [`odata_to_virtual_link`](@ref).

Role: convert
"""
function to_ograph(spine::FlowSpine)
    labels = sort(collect(keys(spine.vertex_signs)))
    label_map = Dict(p=>i for (i,p) in enumerate(labels))
    next = Dict((e.source_vertex,e.source_family)=>(e.target_vertex,e.target_family) for e in spine.edges)
    states = [(p,f) for p in labels for f in (:beta,:alpha)]
    Set(keys(next)) == Set(states) == Set(values(next)) ||
        throw(ArgumentError("singular edges do not define oriented component cycles"))
    seen = Set{Tuple{Int,Symbol}}(); gauss = Vector{Vector{Int}}()
    for start in states
        start in seen && continue
        word = Int[]; state = start
        while !(state in seen)
            push!(seen,state)
            p,f = state
            push!(word,(f == :beta ? 1 : -1)*label_map[p])
            state = next[state]
        end
        state == start || throw(ArgumentError("open o-graph component"))
        push!(gauss,word)
    end
    edge_of=Dict((e.source_vertex,e.source_family)=>i for (i,e) in enumerate(spine.edges))
    blocks=Dict{String,Any}[]; output_edges=Dict{String,Any}[]
    for (component,word) in enumerate(gauss)
        for (position,p) in enumerate(word)
            source=labels[abs(p)]; family=p>0 ? :beta : :alpha
            push!(output_edges,Dict("spine_edge"=>edge_of[(source,family)],
                "component"=>component,"position"=>position,"from"=>p,
                "to"=>word[mod1(position+1,length(word))]))
        end
        start=findfirst(i->sign(word[i])!=sign(word[mod1(i-1,length(word))]),eachindex(word))
        start===nothing && (start=1)
        ordered=[word[mod1(start+i,length(word))] for i in 0:length(word)-1]
        runs=Vector{Vector{Int}}()
        for p in ordered
            (isempty(runs) || sign(last(last(runs)))!=sign(p)) && push!(runs,Int[])
            push!(last(runs),p)
        end
        for run in runs
            push!(blocks,Dict("component"=>component,"family"=>first(run)>0 ? "beta" : "alpha",
                "word"=>run,"source_intersections"=>labels[abs.(run)]))
        end
    end
    OData(gauss,[spine.vertex_signs[p] for p in labels],Dict{String,Any}(
        "representation"=>"normal-o-graph","canonical"=>false,
        "source_hash"=>spine.certificate.source_hash,"convention"=>KODA_CONVENTION,
        "intersection_to_vertex"=>Dict(string(p)=>i for (p,i) in label_map),
        "edge_map"=>output_edges,"runs"=>blocks,
        "o_block_count"=>count(e->sign(e["from"])!=sign(e["to"]),output_edges)÷2))
end

"""
    to_heegaard(PH::PuncturedHeegaardDiagram; strategy=:retained_source) -> CellularHeegaardDiagram

Retrieve the cellular Heegaard diagram retained in a punctured diagram.

# Arguments
`strategy` must be `:retained_source` or its compatibility spelling `:block`, otherwise `ArgumentError` is thrown. No source validation is performed.

# Returns
The exact stored `PH.diagram` object.

# Guarantees
1. The result is the retained source diagram, not a newly reconstructed splitting. [proved/agent: by construction]
2. The input is not modified; the returned diagram is shared. [proved/agent: by construction]

# Not guaranteed
This does not invert o-data or construct a Heegaard diagram from an arbitrary spine. Use `odata_to_heegaard_koda` for that conversion domain.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.3 — retained source versus inverse construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> h = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); p = H.make_punctured(h; joining_disk=H.joining_disk(h,1));

julia> H.to_heegaard(p) === h
true
```

Role: access
"""
function to_heegaard(PH::PuncturedHeegaardDiagram;strategy=:retained_source)
    strategy in (:retained_source,:block) || throw(ArgumentError("this overload only retrieves the retained source"))
    PH.diagram
end
"""
    to_heegaard(spine::FlowSpine; strategy=:retained_source, display_cuts=false) -> CellularHeegaardDiagram

Retrieve the cellular Heegaard source retained by a flow spine.

# Arguments
`spine.source` must be a `PuncturedHeegaardDiagram`; absence raises `ArgumentError`. `display_cuts=true` raises `ArgumentError` naming the removed reverse-display option. `strategy` has the punctured overload's accepted values and errors.

# Returns
The source punctured diagram's stored cellular diagram.

# Guarantees
1. The result has the punctured overload's source-retrieval semantics. [proved/agent: by construction]
2. The input is not modified; source objects are shared. [proved/agent: by construction]

# Not guaranteed
The result is not reconstructed from the current spine cells, and editing spine records does not update the retained source.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.3 — source retrieval and inverse construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> h = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]], [[1,2,3]], [H.Intersection(i,1,1,s) for (i,s) in enumerate([1,-1,1])]); cellular=true)); p = H.make_punctured(h; joining_disk=H.joining_disk(h,1));

julia> H.to_heegaard(H.to_flow_spine(p)) === h
true
```

Role: access
"""
function to_heegaard(spine::FlowSpine;strategy=:retained_source,display_cuts::Bool=false)
    display_cuts && throw(ArgumentError("display_cuts belonged to the removed marked o-data reverse"))
    spine.source isa PuncturedHeegaardDiagram || throw(ArgumentError("no retained Heegaard source"))
    to_heegaard(spine.source;strategy=strategy)
end

_ograph_hash(G::OData) = bytes2hex(SHA.sha256(_canonical_json(Dict("gauss"=>G.gauss,"orientation"=>G.orientation))))

"""
    compare_representations(H::CellularHeegaardDiagram, G::OData) -> Dict{String,Any}

Compare retained source provenance between a cellular Heegaard diagram and o-data.

# Arguments
`H` is a legacy cellular diagram and `G` an `OData` record.
`G.metadata["source_hash"]`, when present, must be a string; otherwise string
operation errors can propagate. Canonicalization and serialization failures
of `H` propagate. Neither the o-data's validity nor a geometric conversion
between the inputs is checked.

# Returns
A dictionary containing `source_hash`, `ograph_hash`, `source_matches`,
`manifold_comparison == "not performed"` and `canonical_conversion == false`.
The source identity respects the labelled-hash form if recorded in `G`;
otherwise it uses the legacy canonical hash. The o-data fingerprint covers
its stored words and vertex signs.

# Guarantees
1. `source_matches` reports equality of the recorded provenance identity
   with the identity computed for `H`. [proved/agent: by construction]
2. No manifold comparison is performed, in the sense of
   [conversion evidence](../specs/heegaard-koda-conversion.md#103-多様体と付加構造).
   Neither input is modified. [proved/agent: by construction]

# Not guaranteed
Matching provenance can be supplied manually and is not an independent
mathematical certificate. This operation does not decide whether the inputs
represent the same manifold, nor whether their words are equivalent.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §10.3 — provenance and manifold-level correspondence are distinct evidence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> g = OData([[1,-1]], [1], Dict{String,Any}());

julia> compare_representations(d, g)["manifold_comparison"]
"not performed"
```

See also [`compare_odata_words`](@ref), [`Topology.HeegaardCodes.canonical_hash`](@ref).

Role: condition
"""
function compare_representations(H::HeegaardDiagram,G::OData)
    # Provenance is not an independently computed manifold invariant.
    identity=get(G.metadata,"source_hash","")
    source=startswith(identity,"labelled-sha256:") ? "labelled-sha256:"*_labelled_hash(H) : canonical_hash(H)
    _compare_representation_hash(source,G)
end
function _compare_representation_hash(source_hash,G)
    Dict{String,Any}("source_hash"=>source_hash,"ograph_hash"=>_ograph_hash(G),
        "source_matches"=>get(G.metadata,"source_hash",nothing) == source_hash,
        "manifold_comparison"=>"not performed","canonical_conversion"=>false)
end

"""
    to_odata(H::CellularHeegaardDiagram; joining_disk=nothing, strategy=:koda) -> ODataBundle

Convert a legacy cellular Heegaard diagram with a selected joining disk into an o-data bundle.

# Arguments
`H` and `joining_disk` must satisfy `make_punctured` and
[the forward Koda domain](../specs/heegaard-koda-conversion.md#5-heegaard-条件と-joining-disk-の入力検査).
The joining disk must be supplied explicitly; `strategy=:koda` is supported.
Puncture validation and spine-encoding exceptions propagate, including
`ArgumentError` for missing disk or unsupported strategy.

# Returns
An `ODataBundle` containing the punctured source, flow spine, o-data,
conversion certificate and representation-comparison record.

# Guarantees
1. On success, the bundle is the forward encoding of the checked punctured
   input in the fixed Koda convention. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The certificate has unverified evidence status and its comparison does not
perform a manifold comparison. General identification with the intended
thickened 3-manifold remains within [the specification's evidence limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合).

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5–7, 13.3 — punctured input, forward encoding and limits.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> ph = odata_to_heegaard_koda(OData([[1,-1]], [1], Dict{String,Any}())).value; b = to_odata(ph.diagram; joining_disk=ph.joining_disk); b.ograph isa OData
true
```

See also [`odata_to_virtual_link`](@ref), [`save_conversion`](@ref).

Role: convert
"""
function to_odata(H::HeegaardDiagram;joining_disk=nothing,strategy=:koda)
    to_odata(make_punctured(H;joining_disk=joining_disk,strategy=strategy))
end
"""
    to_odata(PH::PuncturedHeegaardDiagram) -> ODataBundle

Return the flow-spine encoding and correspondence bundle of a punctured legacy diagram.

# Arguments
`PH` lies in [the punctured Koda domain](../specs/heegaard-koda-conversion.md#5-heegaard-条件と-joining-disk-の入力検査).
Joining-disk, strand-cycle and conversion failures from `to_flow_spine` and
`to_ograph` propagate. The caller supplies compatible source and certificate
data; a manually constructed record is not certified by its type.

# Returns
`ODataBundle` with `punctured=PH`, `spine`, `ograph`, `certificate` and
`comparison`. The certificate records the source and generated identities,
joining disk, convention, label correspondence and `:unverified` status.

# Guarantees
1. The bundle's o-data encodes its flow spine, retaining the conversion
   correspondence for this punctured input. [proved/agent: by construction]
2. The input is not modified; the returned source is the supplied record.
   [proved/agent: by construction]

# Not guaranteed
Representation hashes and the comparison record do not prove manifold
equivalence. This is not a canonical conversion across diagram choices;
[the geometric identification limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合) still apply.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6–7, 10 — forward correspondence and round-trip scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> ph = odata_to_heegaard_koda(OData([[1,-1]], [1], Dict{String,Any}())).value; b = to_odata(ph); (b.ograph isa OData, b.certificate.evidence_status)
(true, :unverified)
```

See also [`save_conversion`](@ref), [`compare_representations`](@ref).

Role: convert
"""
function to_odata(PH::PuncturedHeegaardDiagram)
    spine = to_flow_spine(PH); G = to_ograph(spine)
    source_hash=spine.certificate.source_hash
    comparison = _compare_representation_hash(source_hash,G)
    cert = ConversionCertificate(:koda,source_hash,_disk_json(PH.joining_disk),
        Any[],KODA_CONVENTION,deepcopy(spine.certificate.label_map),_ograph_hash(G),comparison,:unverified)
    ODataBundle(PH,spine,G,cert,comparison)
end

function _conversion_json(c::ConversionCertificate)
    Dict{String,Any}("schema"=>"akaghef.heegaard.hgd2ograph/2","strategy"=>String(c.strategy),
        "source_diagram_hash"=>c.source_hash,"joining_disk"=>deepcopy(c.joining_disk),
        "isotopy_moves"=>deepcopy(c.isotopy_moves),"branch_convention"=>c.branch_convention,
        "label_map"=>deepcopy(c.label_map),"generated_hash"=>c.generated_hash,
        "comparison_invariants"=>deepcopy(c.comparison_invariants),"evidence_status"=>String(c.evidence_status))
end

function _conversion_bundle_json(bundle::ODataBundle)
    Dict("schema"=>"akaghef.heegaard.koda/1",
        "ograph"=>Dict("gauss"=>bundle.ograph.gauss,"orientation"=>bundle.ograph.orientation,
            "metadata"=>bundle.ograph.metadata),
        "complex"=>_complex_json(bundle.spine.complex),
        "certificate"=>_conversion_json(bundle.certificate))
end
"""
    save_conversion(path::AbstractString, value) -> AbstractString

Save the certificate record of a Heegaard-to-o-data conversion as JSON.

# Arguments
`value` is a `ConversionCertificate` or an `ODataBundle` containing one;
other values raise `ArgumentError` identifying the accepted types. `path`
must end in `.hgd2ograph.json`, checked with `ArgumentError`, and have a
writable existing parent. File and serialization errors propagate.

# Returns
The supplied path after writing the certificate JSON.

# Guarantees
1. The saved record preserves the certificate fields in the conversion
   serialization; a bundle contributes its certificate only.
   [proved/agent: by construction]
2. The value is not modified. Existing contents at the path are overwritten.
   [proved/agent: by construction]

# Not guaranteed
This does not save the bundle's complete spine or o-data and does not
validate the certificate's mathematical claims. Writing is not atomic.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.ConversionCertificate(:koda, "source", Dict(), Any[], "candidate_A", Dict(), nothing, Dict(), :unverified);

julia> mktempdir() do dir; p = joinpath(dir, "example.hgd2ograph.json"); H.save_conversion(p, c) == p && isfile(p) end
true
```

See also [`to_odata`](@ref).

Role: other
"""
function save_conversion(path::AbstractString,value)
    c = value isa ODataBundle ? value.certificate : value
    c isa ConversionCertificate || throw(ArgumentError("expected conversion certificate or bundle"))
    endswith(path,".hgd2ograph.json") || throw(ArgumentError("conversion path must end in .hgd2ograph.json"))
    open(path,"w") do io
        println(io,_canonical_json(_conversion_json(c)))
    end
    path
end
