"""Heegaard diagram subsystem for `Topology.jl`."""
module Heegaard

using SHA
using JSON3
using StructTypes
using LinearAlgebra

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
include("render.jl")
include("conversion_spine.jl")

export Curve, Intersection, HeegaardCode, RegionPartition, DecoratedHeegaardCode,
       HeegaardDiagram, PuncturedHeegaardDiagram, FlowSpine, OData,
       ValidationIssue, ValidationReport, ValidationError, MoveCertificate,
       ConversionCertificate, HeegaardLayout, CutSurface, Puncture, Port,
       RouteCurve, PolylineCurve, CubicBezierSegment, BezierChain,
       RibbonGraph, GroupPresentation,
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
       render, render_tikz, inspect_html, save_layout, load_layout,
       make_punctured, to_flow_spine, to_ograph, to_odata, to_heegaard,
       compare_representations, save_conversion

end
