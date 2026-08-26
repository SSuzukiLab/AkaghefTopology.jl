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
    if isempty(code.intersections)
        length(code.alpha) == length(code.beta) || throw(ArgumentError("alpha and beta curve counts differ"))
        isempty(partition.boundary_cycle_blocks) ||
            throw(ArgumentError("a zero-intersection code has no ribbon boundary cycles"))
        return 2 - 2length(code.alpha)
    end
    region_euler = sum(2 - length(block) for block in partition.boundary_cycle_blocks)
    (length(graph.vertices) - length(graph.edges)) + region_euler
end

function surface_genus(object)
    code = _undecorated(object)
    χ = euler_characteristic(object)
    iseven(2 - χ) || throw(ArgumentError("reconstructed surface has non-integral genus"))
    genus = (2 - χ) ÷ 2
    genus >= 0 || throw(ArgumentError("reconstructed surface has negative genus"))
    genus
end

function is_cellular(object)
    partition = _surface_partition(object)
    all(length(block) == 1 for block in partition.boundary_cycle_blocks)
end

"""Summarize the result of compressing along each meridian system.

The combinatorial input records disjoint systems.  A complete system on a
closed genus-`g` surface has a sphere as the compressed boundary; this report
keeps that check separate from the ribbon Euler calculation.
"""
function compression_boundary_report(object)
    code = _undecorated(object)
    genus = surface_genus(object)
    alpha_count = length(code.alpha)
    beta_count = length(code.beta)
    alpha_complete = alpha_count == genus
    beta_complete = beta_count == genus
    (; surface_genus=genus,
       alpha_components=alpha_complete ? 1 : missing,
       beta_components=beta_complete ? 1 : missing,
       alpha_boundary_genera=alpha_complete ? Int[0] : Int[],
       beta_boundary_genera=beta_complete ? Int[0] : Int[],
       alpha_complete, beta_complete,
       closed=alpha_complete && beta_complete && alpha_count == beta_count)
end
