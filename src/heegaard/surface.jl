"""Surface reconstruction and compression-boundary summaries."""

function _surface_partition(object)
    object isa HeegaardDiagram && return object.code.regions
    object isa DecoratedHeegaardCode && return object.regions
    throw(ArgumentError("surface reconstruction requires explicit region decoration"))
end

"""
    complement_regions(object) -> Vector{ComplementRegion}

Return the planar complement-region records specified by a legacy ribbon decoration.

# Arguments
`object` is a `DecoratedHeegaardCode` or `CellularHeegaardDiagram`, with
[legacy region data](../specs/heegaard-terminology.md#legacy-cellular-records).
An object without an explicit decoration throws `ArgumentError` explaining
that requirement. Ribbon-construction errors propagate. The caller must
supply a valid boundary partition; partition validity is not checked here.

# Returns
A vector of `ComplementRegion` records, with region `id`, copied
`boundary_cycles` from the corresponding decoration block, and `genus == 0`.

# Guarantees
1. For a valid decorated surface with planar complementary regions, the
   records give its complement components and their boundary cycles.
   [proved/agent: by construction]
2. The input is unchanged and the returned boundary lists are independent
   copies. [proved/agent: by construction]

# Not guaranteed
This does not infer a missing decoration or discover positive-genus
complement regions. Invalid or overlapping blocks can be returned unchanged
in content; use `validate` to check the decorated input.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy cellular 記録 — boundary partitions and planar complement records.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> complement_regions(d)[1].boundary_cycles == [1]
true
```

See also [`boundary_cycles`](@ref), [`RegionPartition`](@ref), [`validate`](@ref).

Role: access
"""
function complement_regions(object)
    partition = _surface_partition(object)
    cycles = boundary_cycles(object)
    [ComplementRegion(index, copy(block), 0)
     for (index, block) in enumerate(partition.boundary_cycle_blocks)]
end

"""
    euler_characteristic(object) -> Int

Return the Euler characteristic of a decorated legacy ribbon surface with planar complements.

# Arguments
`object` is a `DecoratedHeegaardCode` or `CellularHeegaardDiagram` with valid
[legacy region data](../specs/heegaard-terminology.md#legacy-cellular-records).
Explicit decoration is required; otherwise `ArgumentError` explains that
requirement. A code with no intersections or any empty curve word throws
`ArgumentError` with `MARKED_INPUT_REQUIRED`, since the marked placement is
missing. Other ribbon errors propagate. The caller must supply a valid
partition and planar complement interpretation; these are not checked here.

# Returns
An `Int` equal to the Euler characteristic of that decorated surface on the
stated domain.

# Guarantees
1. The result is the Euler characteristic of the ribbon neighbourhood with
   the declared planar complementary regions attached.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not establish connectedness, surface validity, a complete cut
system or a Heegaard splitting. Invalid decorations can yield an integer
without representing the intended surface.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy cellular 記録 — ribbon surface and planar complementary regions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.euler_characteristic(d)
0
```

See also [`surface_genus`](@ref), [`complement_regions`](@ref), [`validate`](@ref).

Role: invariant
"""
function euler_characteristic(object)
    code = _undecorated(object)
    partition = _surface_partition(object)
    graph = ribbon_graph(code)
    cycles = boundary_cycles(graph)
    (isempty(code.intersections) || any(c->isempty(c.visits),_curve_vector(code))) &&
        throw(ArgumentError("MARKED_INPUT_REQUIRED: intersection-only data does not locate zero-intersection curves or specify their complementary surface"))
    region_euler = sum(2 - length(block) for block in partition.boundary_cycle_blocks)
    (length(graph.vertices) - length(graph.edges)) + region_euler
end

"""
    surface_genus(object) -> Int

Return the genus of the connected closed surface reconstructed from legacy region data.

# Arguments
`object` supplies a decorated legacy diagram, or a code whose surface
partition is otherwise available under [legacy records](../specs/heegaard-terminology.md#legacy-cellular-records).
Missing decoration, zero-intersection curves without marked surface data,
disconnected reconstruction, nonintegral genus or negative genus raise
`ArgumentError` identifying the obstruction. Incidence and ribbon errors
propagate; the caller supplies valid surface data.

# Returns
The nonnegative integer genus of the reconstructed surface.

# Guarantees
1. On valid connected orientable surface data, the result satisfies
   `2 - 2g == euler_characteristic(object)`.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is the genus of the represented surface, not the minimal Heegaard genus
of a 3-manifold. It does not locate unspecified empty curves.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy-cellular-records — complementary regions of the decorated ribbon surface.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]); H.surface_genus(d)
1
```

See also [`euler_characteristic`](@ref).

Role: invariant
"""
function surface_genus(object)
    code = _undecorated(object)
    χ = euler_characteristic(object)
    graph=ribbon_graph(object);blocks=_surface_partition(object).boundary_cycle_blocks
    owner=_ribbon_region_owners(graph,blocks);adj=Dict{Int,Set{Int}}()
    for e in graph.edges
        a,b=owner[2e.id-1],owner[2e.id]
        push!(get!(adj,a,Set{Int}()),b);push!(get!(adj,b,Set{Int}()),a)
    end
    length(_components(eachindex(blocks),adj))==1 ||
        throw(ArgumentError("the decorated ribbon surface is disconnected"))
    iseven(2 - χ) || throw(ArgumentError("reconstructed surface has non-integral genus"))
    genus = (2 - χ) ÷ 2
    genus >= 0 || throw(ArgumentError("reconstructed surface has negative genus"))
    genus
end

"""
    is_cellular(object) -> Bool

Test whether every region in a legacy planar decoration is a disk region.

# Arguments
`object` is a `DecoratedHeegaardCode` or `CellularHeegaardDiagram` in the
[legacy planar-region domain](../specs/heegaard-terminology.md#legacy-cellular-records).
An unsupported undecorated object raises `ArgumentError` requiring explicit
decoration. The caller guarantees a valid boundary partition; this function
does not check incidence, partition coverage or embeddedness.

# Returns
`Bool`, testing the disk-region condition for every stored block. An empty
list of blocks satisfies this universal test.

# Guarantees
1. On a valid planar complement decoration, the result states that all
   complementary regions are disks. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
A positive result does not certify the decoration or a closed Heegaard
splitting; arbitrary or incomplete partitions can satisfy the local test.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md, legacy cellular records — planar complement interpretation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.is_cellular(d)
true
```

See also [`complement_regions`](@ref), [`is_closed_heegaard`](@ref).

Role: condition
"""
function is_cellular(object)
    partition = _surface_partition(object)
    all(length(block) == 1 for block in partition.boundary_cycle_blocks)
end

function _ribbon_region_owners(graph,blocks)
    cycles=boundary_cycles(graph);owner=zeros(Int,length(graph.darts))
    for (r,block) in enumerate(blocks), c in block, d in cycles[c]
        owner[d]==0 || throw(ArgumentError("overlapping region boundary decoration"))
        owner[d]=r
    end
    all(>(0),owner) || throw(ArgumentError("incomplete region boundary decoration"))
    owner
end

"""Reconstruct the cut surface from genus-zero regions and transverse bands.

The RegionPartition contract attaches a genus-zero surface with b boundaries to
each block (Euler characteristic 2-b). Cutting alpha leaves beta edge strips:
glue their two oriented face intervals, subtracting one from Euler characteristic
per strip. The two oriented sides of each alpha circle give the remaining
boundaries. Exchange alpha/beta for the other compression. The independent
triangle-corner implementation in marked_topology.jl checks this calculation on
explicit embeddings; no result is deduced from the number of curves alone.
"""
function _ribbon_family_compression(object,role::Symbol)
    code=_undecorated(object);graph=ribbon_graph(code)
    blocks=_surface_partition(object).boundary_cycle_blocks
    owner=_ribbon_region_owners(graph,blocks);adj=Dict{Int,Set{Int}}()
    glued=[e for e in graph.edges if e.family!=role]
    for e in glued
        a,b=owner[2e.id-1],owner[2e.id]
        push!(get!(adj,a,Set{Int}()),b);push!(get!(adj,b,Set{Int}()),a)
    end
    groups=_components(eachindex(blocks),adj)
    membership=Dict(r=>i for (i,rs) in enumerate(groups) for r in rs)
    curves=role==:alpha ? code.alpha : code.beta
    sides=[Tuple{Int,Symbol}[] for _ in groups]
    for c in curves, (side,offset) in ((:left,-1),(:right,0))
        # Names follow the boundary traversal convention of ribbon.jl; only
        # the membership, not an arbitrary drawing-side convention, is used.
        roots=Set(membership[owner[2e.id+offset]] for e in graph.edges if e.family==role && e.curve_id==c.id)
        length(roots)==1 || throw(ArgumentError("curve side is not one connected boundary: $(c.label)"))
        push!(sides[only(roots)],(c.id,side))
    end
    map(enumerate(groups)) do (i,rs)
        strips=[e.id for e in glued if owner[2e.id-1] in rs]
        chi=sum(2-length(blocks[r]) for r in rs)-length(strips)
        b=length(sides[i]);twog=2-b-chi
        twog>=0 && iseven(twog) || throw(ArgumentError("invalid oriented cut-surface topology"))
        (;source_regions=sort(collect(rs)),glued_edges=strips,boundary_sides=sides[i],
          cut_euler_characteristic=chi,cut_boundary_components=b,genus=twog÷2,
          attached_caps=b,euler_characteristic=chi+b)
    end
end

"""
    compression_boundary_report(object) -> NamedTuple

Return the boundary surfaces obtained by compressing a decorated legacy surface along each curve family.

# Arguments
`object` supplies valid [legacy decorated surface data](../specs/heegaard-terminology.md#legacy-cellular-records).
Surface reconstruction, incidence and compression failures propagate as
`ArgumentError` or native lookup errors. Zero-intersection components without
marked embedding data are rejected. The caller supplies the valid orientable
surface and curve embedding required for the geometric interpretation.

# Returns
A named tuple with surface genus, alpha/beta component counts and genus
lists, separate completeness flags, a `closed` flag, detailed alpha/beta
component records and a method description.

# Guarantees
1. On the supported embedded-surface domain, the report describes the two
   compressed boundary surfaces of [the Heegaard construction](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
   [proved/agent: by construction]
2. The completeness flags distinguish the required minimal systems from
   merely spherical compression boundaries. [proved/agent: by construction]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The report does not locate empty curves or identify the resulting closed
3-manifold. A true `closed` flag alone does not imply both completeness flags.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 and lib/topology.jl/docs/specs/heegaard-terminology.md §legacy-cellular-records — compression boundaries and planar complementary regions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]); r = H.compression_boundary_report(d); (r.surface_genus, r.alpha_complete, r.beta_complete, r.closed)
(1, true, true, true)
```

See also [`certify_heegaard`](@ref).

Role: access
"""
function compression_boundary_report(object)
    code = _undecorated(object)
    genus = surface_genus(object)
    alpha=_ribbon_family_compression(object,:alpha)
    beta=_ribbon_family_compression(object,:beta)
    alpha_genera=Int[c.genus for c in alpha];beta_genera=Int[c.genus for c in beta]
    alpha_spherical=all(==(0),alpha_genera);beta_spherical=all(==(0),beta_genera)
    alpha_complete=alpha_spherical && length(alpha)==1 && length(code.alpha)==genus
    beta_complete=beta_spherical && length(beta)==1 && length(code.beta)==genus
    (; surface_genus=genus,
       alpha_components=length(alpha),beta_components=length(beta),
       alpha_boundary_genera=alpha_genera,beta_boundary_genera=beta_genera,
       alpha_complete, beta_complete,
       closed=alpha_spherical && beta_spherical,alpha,beta,
       method="genus-zero regions glued across transverse edge strips, then capped")
end

"Connected components of `nodes` under the adjacency `adj`, as sets."
function _components(nodes,adj)
    unseen=Set(nodes);out=Set{eltype(nodes)}[]
    while !isempty(unseen)
        seed=minimum(unseen);delete!(unseen,seed);comp=Set([seed]);stack=[seed]
        while !isempty(stack)
            v=pop!(stack)
            for w in get(adj,v,())
                if w in unseen
                    delete!(unseen,w);push!(comp,w);push!(stack,w)
                end
            end
        end
        push!(out,comp)
    end
    out
end
