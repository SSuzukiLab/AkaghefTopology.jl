# Source: docs/specs/heegaard-koda-conversion.md §2.2 (public boundary).
"""HC ↔ Koda branched spine ↔ o-data conversion (new LC/HC layer)."""

include("koda_conversion_errors.jl")
include("lickorish_configuration.jl")
include("lc_k_constants.jl")
include("lc_combinatorial.jl")
include("lc_realize.jl")
include("koda_local_data.jl")
include("koda_forget.jl")
include("koda_table_core.jl")
include("koda_interpreter.jl")
include("koda_spine_codec.jl")
include("koda_forward.jl")
include("koda_reverse.jl")

export LickorishConfiguration, HeegaardLayer, ForgottenHeegaard, CombinatorialPuncturedInput,
       read_lc, read_hc, validate_lc_input, validate_hc_input,
       realize_lc, realize_hc, forget_reference,
       enumerate_joining_disk_candidates, as_punctured_heegaard, heegaard_to_spine,
       odata_to_spine, spine_to_odata, spine_to_heegaard_koda,
       odata_to_heegaard_koda, mark_heegaard, odata_to_hc, hc_to_odata,
       compare_odata_words, koda_table_decode, koda_hc_from_boundary, koda_forward_gauss_from_hc,
       spine_topology, spine_pi1_presentation, joining_disk_failure,
       KodaConversionResult, KodaInputError, KodaSpecificationPending, KodaImplementationError, KodaDomainExcluded, KodaNotComputed,
       load_koda_local_tables, koda_vertex_sign_from_frame
