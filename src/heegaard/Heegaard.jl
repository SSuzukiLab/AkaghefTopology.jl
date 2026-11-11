"""Heegaard diagrams in `Topology.jl`.

`Topology.HeegaardDiagram` is the resolved HC, the same struct as
`HeegaardLayer`. `Topology` exports that type beside `VirtualLink`. An
o-graph is one presentation held by `VirtualLink`. `to_odata` and
`hc_to_odata` map into that o-data; the image is not the diagram.

The cellulated certificate is `CellularHeegaardDiagram`. This module still
binds that struct to the name `HeegaardDiagram` for existing methods,
including `lescop_theta.jl`, and does not export that binding.
"""
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
include("conversion_spine.jl")
include("koda_complex.jl")
include("reference_surface.jl")
include("koda_conversion.jl")
include("hc_validate.jl")
# cf: docs/specs/heegaard-moves.md §3 — catalog rows applied by one move method
include("lc_moves.jl")
include("lc_draw.jl")
# cf: lib/topology.jl/src/heegaard/hc_compare.jl — α–β intersection, Smith factors, H₁, and compare_hc on HeegaardLayer
include("hc_compare.jl")

export lickorish_reference, LickorishConfiguration, HeegaardLayer, read_lc, read_hc,
       resolve_hc, is_standard,
       validate_diagram, movable, move, draw_heegaard,
       hc_intersection, hc_first_homology, compare_hc, hc_to_odata,
       odata_to_spine, spine_to_heegaard_koda, odata_to_hc,
       compare_odata_words, KodaInputError, KodaSpecificationPending

export Curve, Intersection, HeegaardCode, RegionPartition, DecoratedHeegaardCode,
       CellularHeegaardDiagram, JoiningDisk, PuncturedHeegaardDiagram, SpineEdge, SpineRegion,
       BranchedSpineComplex, FlowSpine, OData, ValidationIssue, ValidationReport,
       ValidationError, MoveCertificate, ConversionCertificate, CutSurface, Puncture,
       Port, RibbonGraph, GroupPresentation, heegaard_code, decorate, certify_heegaard,
       load_heegaard, save_heegaard, canonicalize_representation, canonical_hash,
       source_manifest, validate, is_valid, repair_suggestions, ribbon_graph,
       boundary_cycles, complement_regions, surface_genus, euler_characteristic,
       is_cellular, compression_boundary_report, is_closed_heegaard,
       intersection_matrix, fundamental_group_presentation, relator,
       abelianization_matrix, smith_normal_form, first_homology, is_homology_sphere,
       group_probe, manifold_fingerprint, reverse_curve_orientation,
       add_canceling_pair, remove_canceling_pair, handle_slide, stabilize, destabilize,
       apply_move, inverse_move, enumerate_local_moves, cut_along_alpha,
       default_puncture_placement, place_ports, build_wiring_problem, joining_disk,
       validate_joining_disk, save_punctured, load_punctured, make_punctured,
       to_flow_spine, to_ograph, to_odata, to_heegaard, compare_representations,
       save_conversion

include("lescop_theta.jl")
export lescop_intersection_terms, lescop_theta

end
