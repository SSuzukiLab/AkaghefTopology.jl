"""Heegaard diagram subsystem for `Topology.jl`."""
module Heegaard

using SHA
using JSON3
using StructTypes
using LinearAlgebra
import Polyhedra

include("types.jl")
include("schema.jl")
include("ribbon.jl")
include("surface.jl")
include("validation.jl")
include("canonicalization.jl")
include("groups.jl")
include("moves.jl")
include("cut_surface.jl")
include("layout.jl")
include("mechanical_relaxation.jl")
include("feasible_router.jl")
include("planar_embedding.jl")
include("layout_transport.jl")
include("automatic_layout.jl")
include("plot_relaxation.jl")
include("surface_chart.jl")
include("joining_disk_layout.jl")
include("render.jl")
include("surface_chart_render.jl")
include("marked_display_plot.jl")
include("conversion_spine.jl")
include("koda_complex.jl")
include("../../vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl")
include("../../vendor/marked_surface_diagrams/src/builders.jl")
include("reference_surface.jl")
include("marked_display.jl")
include("marked_topology.jl")
include("marked_groups.jl")
include("marked_moves.jl")
include("marked_polygon_overlay.jl")
include("marked_subdivision.jl")
include("marked_annulus_product.jl")
include("marked_slides.jl")
include("marked_stabilization.jl")
include("marked_rechart.jl")
include("marked_twists.jl")
include("cellulated_conversion.jl")
include("ograph_reverse.jl")
include("marked_comparison.jl")

export MarkedSurfaceDiagrams, MarkedDiagramBuilders,
       lickorish_reference, prepare_marked_display, marked_topology_report,
       surface_group_presentation, marked_curve_words, marked_computation_report, cellulated_punctured,
       marked_display_source,
       marked_triangle_isotopy, marked_vertex_isotopy, marked_simplicial_action,
       replay_marked_move, inverse_marked_move,
       marked_subdivide, marked_collapse_subdivision,
       marked_dual_star, marked_collapse_dual_star,
       marked_parallel_copy, marked_remove_parallel,
       marked_handle_slide, marked_inverse_handle_slide,
       marked_stabilize, marked_destabilize,
       marked_dehn_twist, marked_inverse_dehn_twist,
       marked_comparison_report, odata_computation_report

export Curve, Intersection, HeegaardCode, RegionPartition, DecoratedHeegaardCode,
       HeegaardDiagram, JoiningDisk, PuncturedHeegaardDiagram, SpineEdge, SpineRegion,
       BranchedSpineComplex, FlowSpine, OData,
       ValidationIssue, ValidationReport, ValidationError, MoveCertificate,
       ConversionCertificate, HeegaardLayout, CutSurface, Puncture, Port,
       RouteCurve, PolylineCurve, CubicBezierSegment, BezierChain,
       RibbonGraph, GroupPresentation, HeegaardSurfaceChart,
       heegaard_code, decorate, certify_heegaard, load_heegaard, save_heegaard,
       canonicalize_representation, canonical_hash, source_manifest,
       validate, is_valid, repair_suggestions,
       ribbon_graph, boundary_cycles, complement_regions, surface_genus,
       euler_characteristic, is_cellular, compression_boundary_report,
       is_closed_heegaard, intersection_matrix, fundamental_group_presentation,
       relator, abelianization_matrix, smith_normal_form, first_homology,
       is_homology_sphere, group_probe, manifold_fingerprint,
       reverse_curve_orientation, add_canceling_pair, remove_canceling_pair,
       handle_slide, stabilize, destabilize, apply_move, inverse_move,
       enumerate_local_moves,
       cut_along_alpha, default_puncture_placement, place_ports,
       build_wiring_problem, gauss_layout, layout_heegaard, smooth_layout,
       mechanically_relax_layout,
       validate_layout, feasible_fixed_endpoint_routes,
       plot, surface_chart, render, render_tikz, inspect_html, save_layout, load_layout,
       joining_disk, validate_joining_disk, save_punctured, load_punctured,
       make_punctured, to_flow_spine, to_ograph, to_odata, to_heegaard,
       compare_representations, save_conversion

end
