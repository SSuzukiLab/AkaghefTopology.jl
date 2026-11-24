"""HeegaardCodes: the HC layer of Heegaard diagrams in `Topology.jl` (docs/specs/module-structure.md).

`Topology.HeegaardDiagram` is the resolved HC, the same struct as
`HeegaardCode`. `Topology` exports that type beside `VirtualLink`. An
o-graph is one presentation held by `VirtualLink`. `to_odata` and
`hc_to_odata` map into that o-data; the image is not the diagram.

The cellulated certificate is `CellularHeegaardDiagram`. This module still
binds that struct to the name `HeegaardDiagram` for existing methods,
including `lescop_theta.jl`, and does not export that binding.
"""
module HeegaardCodes

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
include("bordered_lickorish_code.jl")
export BorderedLickorishCode, BorderedLCError, BorderedLCDiagnosis, BorderedCut,
       BLCFace, BLCPair, BLCPairHit,
       read_bordered_lc, bordered_lc_text, validate_bordered_lc, diagnose_bordered_lc,
       cut_bordered_lc, glue_bordered_lc, bordered_lc_faces,
       vertex_bordered_pair, digon_bordered_pair, swap_bordered_pair, apply_bordered_pair,
       digon_bordered_family, slide_bordered_family, apply_bordered_family,
       apply_slide_bordered, apply_slide_inverse, reverse_slide_band,
       apply_reversed_slide, apply_geometric_slide_inverse
include("hc_validate.jl")
# cf: docs/specs/heegaard-moves.md §3 — catalog rows applied by one move method
include("lc_moves.jl")
include("lc_isotopy_key.jl")
export lc_isotopy_key, LCIsotopyKey
# Joining-disk preparation calls movable/move. It stays out of lc_moves*.jl.
include("koda_joining_prepare.jl")
include("lc_draw.jl")
# cf: docs/specs/heegaard-plot.md — planar Julia entry on resolved HC
include("hc_plot.jl")
include("solid3d.jl")
export draw_heegaard_3d
export plot
# cf: lib/topology.jl/src/heegaard/hc_compare.jl — α–β intersection, Smith factors, H₁, and compare_hc on HeegaardCode
include("hc_compare.jl")
include("hc_invariants_linking.jl")
export HCLinkingForm, HCLinkingInputError, hc_linking_form, linking_pairing, Z_MOO
# cf: lib/topology.jl/docs/specs/module-structure.md §3.1 — HC search records on MoveGraphStoreDB
include("hc_search_db.jl")
export hc_search_state, put_hc_search_state!, load_hc_search_state,
       begin_hc_search_run!, finish_hc_search_run!, hc_search_run,
       hc_search_candidates, expand_hc_search_state!,
       hc_search_neighbors, query_hc_search_paths,
       migrate_hc_search_store!, hc_search_explore,
       hc_search_same_homology, hc_search_same_linking
include("hc_triangulation.jl")
include("hc_connected_sum.jl")
include("hc_destabilization.jl")
export hc_connected_sum, hc_destabilize
include("hc_invariants_tv.jl")
export Z_TuraevViro, hc_regina_runner, ReginaTVError
include("hc_invariants_tv_skeleton.jl")
export Z_TuraevViro_skeleton
# Lane-owned entry files (2026-11-14); each lane adds its own exports inside its file.
include("hc_standard_check.jl")
include("hc_invariants_extra.jl")
export HCClosedTriangulation, HCTriangulationError, SnapPyUnavailable, SnapPyExecutionError,
       hc_triangulation, hc_triangulation_data, hc_snappy, hc_snappy_runner,
       hc_snappy_compare, hc_snappy_catalog

# Old type names, kept as aliases for existing scripts (docs/specs/module-structure.md §2).
const LickorishConfiguration = LickorishCode
const HeegaardLayer = HeegaardCode
export LickorishConfiguration, HeegaardLayer

export lickorish_reference, LickorishCode, HeegaardCode, read_lc, read_hc,
       resolve_hc, is_standard,
       validate_diagram, movable, move, draw_heegaard,
       hc_intersection, hc_first_homology, compare_hc, hc_to_odata,
       odata_to_spine, spine_to_heegaard_koda, odata_to_hc,
       compare_odata_words, KodaInputError, KodaSpecificationPending

export Curve, Intersection, CombinatorialHeegaardCode, RegionPartition, DecoratedHeegaardCode,
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
export lescop_intersection_terms, lescop_theta, lescop_p1_difference,
       lescop_epsilon, lescop_spinc_label

include("hc_lescop.jl")
export hc_to_cellular, hc_lescop_theta, hc_lescop_p1_difference,
       hc_lescop_epsilon, hc_lescop_spinc_label

end
