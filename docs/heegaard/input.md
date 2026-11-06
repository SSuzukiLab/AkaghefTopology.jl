# Mathematical input and lossless save

Use `Heegaard.load_heegaard(path)` and `Heegaard.save_heegaard(path, object)`
for JSON mathematical data. Legacy `akaghef.heegaard/1` returns its existing code
types. A native full-cellulation `version=1` document returns a defensive
`Dict{String,Any}`. Its genus, all family sizes, zero-intersection and disconnected
curves, ordered incidence, explicit marking, directed words, region marks and
metadata survive save/reload. Neither function infers a marking or requires the
curve counts to equal the genus. They validate the embedded source and D;
closed-Heegaard validity is a separate `validate`/compression question.

The native schema admits only its registered root fields: `version`, `surface`,
`patches`, `curves`, `crossings`, and optional `region_marks` and `metadata`.
Store source locators and other nonmathematical extras inside the metadata
object. Unknown mathematical fields and unsupported versions are rejected.
Array order, exact rational side-parameter strings and edge-word spelling are
retained. JSON whitespace and object-key order need not be retained.

## Construct and correct with existing tools

`examples/HeegaardInput.jl::heegaard_input_example(path)` builds two marked torus
curves with the existing `MarkedDiagramBuilders` and demonstrates edit/save/reload.
`grid_document` supports disk, annulus and torus; `lickorish_reference(g)` supplies
the existing closed reference surface with explicitly auxiliary reference curves
for the supported genera. This reference is a chosen marking; it does not recover
the marking of a supplied diagram. General surfaces can use the native patch
dictionary directly. Do not introduce alpha curves by renaming display cuts.

For a mathematical correction, copy the loaded dictionary, edit its curve word,
crossing declaration, gluing or D seed together with affected incidence, and call
`save_heegaard`. The full cellulation is checked again; invalid edits are rejected
before an existing file is opened. Kernel errors identify the failed constraint.
This editor operation is not an isotopy or a certified mathematical move.

The existing app's **データを開く** dialog imports the resulting JSON through its
existing `/api/render` endpoint. Rendering can additionally need an explicit
feasible layout and display cuts. A valid mathematical document can therefore
remain unsupported by an automatic layout. For a marked display, the app's
**図式＋配置を保存** link downloads the whole display envelope with its retained
source; its other download saves SVG. A `marked-cut-view/1` document is a display
envelope: call `marked_display_source(view)` to check consistency and extract its
retained source, then use `save_heegaard` for a separate mathematical source file.
Legacy views additionally expose separate mathematical and layout downloads.
Legacy joining-disk objects retain `save_punctured`/`load_punctured`.

## Transcribe literature

Consult the research workspace's `docs/LITERATURE.md` and the managed Work before
transcribing. Open the existing linked page/crop, record the Work, exact page and
figure locator and source asset in `metadata`, and distinguish visual inspection
from a checked mathematical transcription. Preserve the source original.

Specify the oriented surface patches and reversed side gluings first. Record
each labelled curve's family, orientation and full directed edge word, including
seam jumps and curves with no crossings. Declare ordered crossing pairs with
their signs; identify D by a region seed and its actual claimed condition. Mark
every route, gluing or orientation ambiguity explicitly instead of filling it in
from the drawing convention. Record any added auxiliary display cuts separately.

Validate, save, reload and compare the whole labelled oriented source with the
original crop: incidence, complementary regions, seam gluing, route class and D.
Then inspect its rendered realization. A successful render, representation digest
or transcription hash does not establish ambient isotopy with the source.
Unresolved figure details remain unresolved in metadata and the owning record.
