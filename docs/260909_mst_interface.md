# MST 不変量インターフェース

本書は公開契約である。共通の設計方針は [AkaghefAlgebra の設計方針](../../algebra.jl/docs/260909_interface_design.md) を正本とする。

公開入口を `Z_MST(v, H; strategy=...)` とする。v は必要な整数重みなどを備えた入力、H は対象の MST 定義に必要な代数データを保持する。入力の重みを評価中に選び直さない。対象文献の仮定・規約・正規化定数を維持する。

出力型は H に依存する。Taft では根の指定を保持した厳密値を返す。簡約化戦略は評価方法を選び、完成した不変量の値を変更しない。途中の式を操作する利用者には、公開入口が使う構築・簡約化・評価・正規化の部品を提供する。

実装入口は `Z_MST`、段階別 API は `invariant_expression`、`simplify_invariant`、`evaluate_invariant`、`normalize_root` である。`hopf_invariant`、`taft_reduce`、`taft_evaluate` は互換入口として接続する。Taft の計算用型は内部 `TaftExact.Cyclotomic{N}`、完成値は公開 `Cyclotomic{N,M}`。使用方法は [README](../README.md#algebra-units) を参照する。

Taft 固有の型・移行・受入条件は [要件書](../../../threads/PJ260830/docs/reports/260909_mst_interface_requirements.md) に置く。

U(R) は `Z_MST(v, URAlgebra())` で厳密有理数を返す。
構築・簡約化・評価は同じ共通APIに接続し、URDiagramには依存しない。
対応範囲と規約は [U(R) backend](260909_ur_mst.md) を参照する。

値がないときの規則は [設計方針](../../algebra.jl/docs/260909_interface_design.md#値がないとき) に従う。不変量ごとの未定義の条件は、各 backend の文書に書く。
