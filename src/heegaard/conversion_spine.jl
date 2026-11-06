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

"""Select a disk face; check K08 section 3's embedded-boundary and joining conditions."""
function joining_disk(H::HeegaardDiagram, region::Integer)
    blocks = H.code.regions.boundary_cycle_blocks
    1 <= region <= length(blocks) || throw(ArgumentError("joining disk region is out of range"))
    length(blocks[region]) == 1 || throw(ArgumentError("joining region must be a disk, not a multi-boundary region"))
    D = JoiningDisk(_labelled_hash(H),Int(region),copy(boundary_cycles(H)[only(blocks[region])]))
    validate_joining_disk(H,D)
    D
end

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

"""Remove D-boundary edges and smooth the resulting bivalent corners.

Every singular edge retains its ordered source-dart path. Non-cellular and
vertex-free cases are rejected rather than silently standardized.
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

"""Retrieve retained source; this overload is not an inverse construction.
`:block` is retained as a compatibility spelling for `:retained_source`.
Use `to_heegaard(::OData)` for the independent dual-tetrahedron construction.
"""
function to_heegaard(PH::PuncturedHeegaardDiagram;strategy=:retained_source)
    strategy in (:retained_source,:block) || throw(ArgumentError("this overload only retrieves the retained source"))
    PH.diagram
end
function to_heegaard(spine::FlowSpine;strategy=:retained_source,display_cuts::Bool=false)
    # Source: OData reverse below. Display scaffolding belongs only to the new
    # independent splitting; never mutate a retrieved marked/source diagram.
    strategy==:dual_tetrahedra && return to_heegaard(to_ograph(spine);strategy=strategy,display_cuts=display_cuts)
    display_cuts && throw(ArgumentError("display_cuts requires the independent :dual_tetrahedra strategy"))
    spine.source isa PuncturedHeegaardDiagram || throw(ArgumentError("no retained Heegaard source"))
    to_heegaard(spine.source;strategy=strategy)
end

"""Independently construct a marked splitting from o-data, without source metadata.
The output represents the oriented dual-triangulated manifold; source marking,
joining D and branching/vector-field recovery are explicitly outside its claim.
`display_cuts=true` adds distinct alpha-left auxiliary meridians and an explicit
new infinity seed. Alpha and beta stay visible mathematical curves; their role
is never changed into a display cut. This option does not generate coordinates.
"""
function to_heegaard(G::OData;strategy=:dual_tetrahedra,display_cuts::Bool=false)
    strategy==:dual_tetrahedra || throw(ArgumentError("o-data reverse strategy must be :dual_tetrahedra"))
    _ograph_to_marked_heegaard(G;display_cuts=display_cuts)
end

"""Construct independently from the bundle's actual o-data.

Source: `to_heegaard(::OData)`; retained punctured diagrams, certificates and
comparison metadata are not consulted. This is a new splitting, not a marked
round-trip or recovery of the source joining disk or genus.
"""
function to_heegaard(bundle::ODataBundle;strategy=:dual_tetrahedra,display_cuts::Bool=false)
    to_heegaard(bundle.ograph;strategy=strategy,display_cuts=display_cuts)
end

_ograph_hash(G::OData) = bytes2hex(SHA.sha256(_canonical_json(Dict("gauss"=>G.gauss,"orientation"=>G.orientation))))

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

function to_odata(H::HeegaardDiagram;joining_disk=nothing,strategy=:koda)
    to_odata(make_punctured(H;joining_disk=joining_disk,strategy=strategy))
end
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
function save_conversion(path::AbstractString,value)
    c = value isa ODataBundle ? value.certificate : value
    c isa ConversionCertificate || throw(ArgumentError("expected conversion certificate or bundle"))
    endswith(path,".hgd2ograph.json") || throw(ArgumentError("conversion path must end in .hgd2ograph.json"))
    open(path,"w") do io
        println(io,_canonical_json(_conversion_json(c)))
    end
    path
end
