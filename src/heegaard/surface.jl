"""Surface reconstruction and compression-boundary summaries."""

function _surface_partition(object)
    object isa HeegaardDiagram && return object.code.regions
    object isa DecoratedHeegaardCode && return object.regions
    throw(ArgumentError("surface reconstruction requires explicit region decoration"))
end

function complement_regions(object)
    partition = _surface_partition(object)
    cycles = boundary_cycles(object)
    [ComplementRegion(index, copy(block), 0)
     for (index, block) in enumerate(partition.boundary_cycle_blocks)]
end

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

function surface_genus(object)
    code = _undecorated(object)
    χ = euler_characteristic(object)
    graph=ribbon_graph(object);blocks=_surface_partition(object).boundary_cycle_blocks
    owner=_ribbon_region_owners(graph,blocks);adj=Dict{Int,Set{Int}}()
    for e in graph.edges
        a,b=owner[2e.id-1],owner[2e.id]
        push!(get!(adj,a,Set{Int}()),b);push!(get!(adj,b,Set{Int}()),a)
    end
    length(MarkedSurfaceDiagrams.components(eachindex(blocks),adj))==1 ||
        throw(ArgumentError("the decorated ribbon surface is disconnected"))
    iseven(2 - χ) || throw(ArgumentError("reconstructed surface has non-integral genus"))
    genus = (2 - χ) ÷ 2
    genus >= 0 || throw(ArgumentError("reconstructed surface has negative genus"))
    genus
end

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
    groups=MarkedSurfaceDiagrams.components(eachindex(blocks),adj)
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

"""Compute compression boundary components from the decorated ribbon surface.

Zero-intersection components require a full marked embedding. Redundant systems
are distinguished from complete minimal cut systems; either may compress to a
collection of spheres. Consumer algorithms requiring a basis must check the
separate alpha_complete/beta_complete predicates.
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
