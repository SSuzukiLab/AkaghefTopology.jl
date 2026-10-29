using Topology, JSON3, Dates
const H = Topology.Heegaard
include(joinpath(@__DIR__,"../../examples/HeegaardAbalone.jl"))
include(joinpath(@__DIR__,"../../examples/VLExample.jl"))
out = joinpath(@__DIR__,"../../artifacts/heegaard/full_curves/abalone")
mkpath(out)
example = HeegaardAbaloneExample.fixture()
PH = H.make_punctured(example.diagram;joining_disk=example.disk)
H.save_punctured(joinpath(out,"S3.punctured.hgd.json"),PH)
# Exercise the persisted entry point, not only an in-memory construction.
bundle = H.to_odata(H.load_punctured(joinpath(out,"S3.punctured.hgd.json")))
converted = VirtualLink()
set_data!(converted;rgauss=bundle.ograph.gauss,orientation=bundle.ograph.orientation)
reference = only(filter(x->x.name=="abalone",run_vl_plot_calls())).link
actual = real_gauss_code(converted); expected = real_gauss_code(reference)
closed,criteria = is_closed(converted)
@assert actual == expected
@assert closed
@assert H.surface_genus(example.diagram)==1
@assert abs(only(H.intersection_matrix(example.diagram)))==1
H.save_conversion(joinpath(out,"S3.hgd2ograph.json"),bundle)
write(joinpath(out,"heegaard.svg"),HeegaardAbaloneExample.render_svg(example))
plot_svg(converted,joinpath(out,"ograph.svg");edge_labels=:none)
comparison = Dict("scope"=>"one explicit genus-one S3 diagram and selected disk; exact GC comparison",
    "date"=>string(now()),"julia"=>string(VERSION),"agent"=>"Codex",
    "gauss"=>actual[1],"vertex_signs"=>actual[2],"reference"=>"examples/VLExample.jl: abalone",
    "reference_gauss"=>expected[1],"reference_vertex_signs"=>expected[2],
    "exact_match"=>actual==expected,"closed"=>closed,"closedness_C1_C2_C3"=>collect(criteria),
    "genus"=>H.surface_genus(example.diagram),"algebraic_intersection"=>only(H.intersection_matrix(example.diagram)),
    "smoothed_intersections"=>bundle.certificate.label_map["smoothed_intersections"],
    "remaining_intersections"=>bundle.certificate.label_map["remaining_intersections"],
    "convention"=>H.KODA_CONVENTION,"general_conversion_verification"=>"not established",
    "input_hash"=>H._labelled_hash(example.diagram),"output_hash"=>bundle.certificate.generated_hash)
write(joinpath(out,"comparison.json"),H._canonical_json(comparison)*"\n")
println(JSON3.write(comparison))
