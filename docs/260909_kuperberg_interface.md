# Kuperberg 不変量インターフェース（HC）

公開入口は `Z_Kuperberg(hc, H; strategy=...)` とする。`hc` は `Heegaard.read_hc` で読んだ `HeegaardLayer`、または `.hc` のパスである。代数 `H` は有限次元の**対合** Hopf 代数のみを定義域とする。積分規約は右積分・右余積分で `λ(Λ)=1`（[[PJ261007-LEDGER]]、`threads/PJ261007/docs/notes/2026-10-07-144500-mst-kuperberg-comparison.md`）。

段階別 API は MST と同型で `invariant_expression`、`simplify_invariant`、`evaluate_invariant` である。入力は joining disk / spine / MST の `_hopf_invariant_factors` を使わない。HC の α・β 語から α∩β だけを remapping した組合表現を作り、**lickorish-code.md §3 の passage sign** で crossing 符号（Kup91 §5 の (lower,upper) 組基に対応）を置く。交点は `positive_vertex` / `negative_vertex`、弧は無重みの port 同定、最後に Kup91 の `(dim H)^(g(S)−n_α−n_β)` を掛ける（pp. 12–13）。原典: G. Kuperberg, *Involutory Hopf algebras and 3-manifold invariants* (1991), PDF pool `900519_...pdf`（セミナー `.lnk` 経由）。

## 定義域

| 条件 | 失敗時 |
|---|---|
| `is_involutory(H)` | `DomainError`（非対合は framing 付きの別理論；未実装） |
| 右積分で `λ(Λ)=1` | `DomainError` |
| 種数 g に対し alpha・beta が各 g 本 | `DomainError` |
| `realize_hc` 後の圧縮が両族で complete cut system | `DomainError` |
| HC / LC の入力検査 | `KodaInputError` / `ArgumentError` |

値がないときの分類は [設計方針](../../algebra.jl/docs/260909_interface_design.md#値がないとき) に従う。未定義は `DomainError` で返す。

非対合の Hopf 代数（例：small Borel、Sweedler）は黙って受け取らない。
