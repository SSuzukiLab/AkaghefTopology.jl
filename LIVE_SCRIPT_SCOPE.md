# Repository-wide Live Script scope

Audit date: 2026-08-13

This inventory covers every `.mlx` below the Execution repository.  Classification
is based on executable source extracted with MATLAB R2026a, not on directory or
filename alone.

- `T`: direct topology workflow. It creates or mutates `VirtualLink`/
  `URDiagram`, consumes Gauss/o-data, converts diagrams, computes a topology
  invariant, or renders an o-graph.
- `D`: required algebra/tensor dependency workflow reached by a `T` workflow.
- `R`: inspected research/example Live Script with no executable call into the
  topology dependency closure. It remains in the 60-file audit, but is not a
  topology-library acceptance test.
- `E`: empty or generated experiment scaffold.

All 60 Live Scripts are migration targets. The 17 `T` + `D` scripts are the
first implementation wave because they exercise the topology dependency
closure, but `R` is an implementation backlog, not an exclusion from the
user's scope. `E` scripts require a reproducible Julia placeholder or explicit
generated-artifact reproduction before they can be closed. `C250405Sweedler.mlx`
is included even though it depends on a pre-existing workspace variable `vl`;
excluding it due to the missing constructor would repeat the earlier scope
error.

| Class | Live Script | Evidence |
|---|---|---|
| R | `Docs/experiment/C260209SartoriRmatrix.mlx` | standalone symbolic R-matrix calculation |
| T | `Docs/experiment/C260225M_Oinv.mlx` | `VirtualLink`, o-data, weight, invariant, plot |
| E | `Results/Experiment1_Result1_20260125T131846/Snapshot/exFun.mlx` | generated experiment function |
| R | `algebra/Core/tensor/experiment_expression.mlx` | standalone tensor-expression experiment |
| R | `algebra/Examples/PolAlg/WeylAlgTest.mlx` | Weyl algebra example |
| R | `algebra/Examples/VectAlgebra/C260201KPalg.mlx` | vector algebra experiment; no topology API call |
| R | `algebra/Execution/C2410/A241030.mlx` | algebra scratch calculation |
| R | `algebra/Execution/C2410/Calc240908.mlx` | algebra scratch calculation |
| R | `algebra/Execution/C2410/Calc241006.mlx` | algebra scratch calculation |
| R | `algebra/Execution/C2410/Calc241107qModular.mlx` | modular algebra calculation |
| R | `algebra/Execution/C2410/Calc241128C2modalg.mlx` | modular algebra calculation |
| R | `algebra/Execution/C2410/Calc241217C2modalg.mlx` | modular algebra calculation |
| R | `algebra/Execution/C2410/Calc241228Cnmodalg.mlx` | modular algebra calculation |
| R | `algebra/Execution/C2410/Calc250127ModalgProof.mlx` | modular algebra proof experiment |
| R | `algebra/Execution/C2410/SL3.mlx` | Lie algebra experiment |
| R | `algebra/Execution/C2410/strEndVTest.mlx` | string endomorphism experiment |
| R | `algebra/Execution/C2508/C250820FpX.mlx` | polynomial algebra experiment |
| D | `algebra/Execution/C2508/C250825Dsl2.mlx` | Drinfeld/Heisenberg double and `Uqsl2BorelSmall` used by invariants |
| E | `algebra/Execution/C2508/C250909DH.mlx` | empty Live Script |
| R | `algebra/Execution/C2508/C250926QuasiHopfHn.mlx` | quasi-Hopf research; no topology workflow call |
| R | `algebra/Execution/C2511/Calc251112.mlx` | algebra calculation |
| R | `algebra/Execution/C2511/Calc251117Cnmodalg.mlx` | modular algebra calculation |
| R | `algebra/Execution/Docs/W241223QiitaMATLABAdvent.mlx` | algebra article source |
| R | `algebra/Execution/Docs/W241224QiitaMATLABAdvent.mlx` | algebra article source |
| R | `algebra/Execution/Docs/rewriteDocs1.mlx` | documentation rewrite helper |
| T | `algebra/Execution/MSTinv/C250405Sweedler.mlx` | mutates workspace `vl`, sets Gauss/o-data, calculates weights/invariants |
| R | `algebra/Execution/MSTinv/C250406SwAlg.mlx` | Sweedler/Heisenberg algebra experiment; no topology object |
| D | `algebra/Execution/MSTinv/C250406TensorNetwork.mlx` | o-graph tensor-expression and Hopf-algebra dependency |
| R | `algebra/Execution/MSTinv/C250410TraceOperator.mlx` | trace-operator research; no topology API call |
| R | `algebra/Execution/MSTinv/C250601FpX.mlx` | polynomial experiment |
| R | `algebra/Execution/MSTinv/QuantumMST/C250314.mlx` | numerical invariant-series research plots, not o-graph rendering |
| R | `algebra/Execution/MSTinv/QuantumMST/C250316.mlx` | numerical invariant-series research plots, not o-graph rendering |
| R | `algebra/Execution/MSTinv/QuantumMST/C250324cointegral_uqsl2.mlx` | cointegral research; no topology API call |
| R | `algebra/Execution/MSTinv/tensor_time.mlx` | tensor performance experiment |
| R | `algebra/Execution/SLnDual/sl2_dualbasis.mlx` | dual-basis derivation |
| R | `algebra/Execution/SLnDual/sl2borel_dualbasis.mlx` | dual-basis derivation |
| R | `algebra/Execution/SLnDual/sl3borel_dualbasis_.mlx` | dual-basis derivation |
| R | `algebra/Execution/qAnalog/qBinomTaylor.mlx` | q-binomial series research |
| R | `algebra/Execution/qAnalog/qLogarithm.mlx` | q-logarithm research |
| E | `algebra/untitled.mlx` | empty Live Script |
| E | `exFun.mlx` | generated experiment function |
| R | `experiment/ODE/Experiment1Function1.mlx` | ODE experiment with ordinary numeric plot |
| E | `experiment/ODE/Experiment1Function2.mlx` | generated experiment function |
| E | `experiment/ODE/Experiment1Initialization1.mlx` | generated experiment initialization |
| R | `experiment/PythonIntegrationWithSage.mlx` | Sage integration probe |
| R | `experiment/experiment_expression.mlx` | standalone tensor-expression experiment |
| T | `projects/invariant/C250811MSTinvUR.mlx` | `VirtualLink`, `URDiagram`, `VL2URD`, invariants, 2 o-graph plots |
| T | `projects/invariant/C250818MSTinv_uqsl2.mlx` | Gauss/o-data, weights, invariants, 6 o-graph plots |
| T | `projects/invariant/C250823MSTinv_extalg.mlx` | Gauss/o-data, weights, invariants, 2 o-graph plots |
| T | `projects/invariant/C251020KnotDsl2bfromUR.mlx` | `VirtualLink`, `URDiagram`, `VL2URD`, 2 o-graph plots |
| T | `projects/invariant/C251030WeylAlgTraceFormula.mlx` | `URDiagram` trace/simplification workflow |
| R | `projects/invariant/C251120KnotDUsl2bInvAndTraceFml2.mlx` | knot-series research; no topology-library object or conversion call |
| T | `projects/invariant/C260125MSTinv_uqsl2.mlx` | Gauss/o-data, weights, invariants, 10 o-graph plot calls |
| T | `projects/invariant/C260319MSTinv_unknot.mlx` | `VirtualLink`, Gauss code, weight and invariant |
| T | `projects/spine/example/C260220LensSpaces.mlx` | lens-space Gauss codes, disk table, weight, o-graph plot |
| T | `topology/Manifold/VLExample.mlx` | primary `VirtualLink` example; 19 o-graph plot calls |
| T | `topology/URDiagram/C240610URD.mlx` | `URDiagram` construction and mutation |
| T | `topology/URDiagram/C240623URD.mlx` | o-graph `URDiagram` reductions |
| T | `topology/URDiagram/URDforKnotCmpl.mlx` | knot-complement `URDiagram` conversion workflow |
| T | `topology/URDiagram/dev_URDiagram.mlx` | `URDiagram` development workflow and plot |

## Plot acceptance set

| Live Script | o-graph plot calls | MATLAB figures extracted |
|---|---:|---:|
| `Docs/experiment/C260225M_Oinv.mlx` | 1 | 1 |
| `projects/invariant/C250811MSTinvUR.mlx` | 2 | 2 |
| `projects/invariant/C250818MSTinv_uqsl2.mlx` | 6 | 6 |
| `projects/invariant/C250823MSTinv_extalg.mlx` | 2 | 2 |
| `projects/invariant/C251020KnotDsl2bfromUR.mlx` | 2 | 2 |
| `projects/invariant/C260125MSTinv_uqsl2.mlx` | 10 | 13 |
| `projects/spine/example/C260220LensSpaces.mlx` | 1 | 1 |
| `topology/Manifold/VLExample.mlx` | 19 | 10 |
| `topology/URDiagram/dev_URDiagram.mlx` | 1 | 0 |
| **Total** | **44** | **37** |

The mismatch between calls and embedded figures is expected: Live Scripts may
overwrite/reuse figures, save intermediate output, or contain cells whose output
was not stored. Acceptance therefore requires replaying all 44 calls; the 37
embedded figures alone are an incomplete oracle.
