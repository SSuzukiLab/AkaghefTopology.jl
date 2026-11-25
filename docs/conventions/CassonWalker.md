<!-- XQQQ1500000 -->
# Casson–Walker の規約

既存の幾何学的比較の状態は表のとおり draft である。下の Alexander 保存形式の辞書は実装とテストによって固定する。対象は `Invariants` の [casson_walker.jl](../../src/casson_walker.jl) で、手術係数と Seifert データから [[Les07#^anonnd70-blk-0138]] の正規化の λ を返す。HC は読まない。

| 規約 | 内容 | 状態 |
|---|---|---|
| 正規化 | λ = W₁∘Z₁、W₁(tata) = 2。λ_W の半分。残差の 24λ はこの λ | draft |
| レンズの名前 | p > 0 で L(p, q) = S³_U(p/q)。子午線の自己絡み数は −q/p。λ(L(p, q)) = −s(q, p)/2 | draft |
| S³ の結び目の (p, q) | SnapPy の `dehn_fill((p, q))` を、上の手術係数 p/q として読む | draft |
| 向きの反転 | `reverse_orientation()` のあと同じ (p, q) は、元の (−p, q) 充填の向きの反転。鏡像は反転のあと (−p, q)。λ は向きの反転で符号が変わる | 周辺系は台帳の規約。λ の符号は draft |
| Seifert の閉じた式 | 円周束の L(p, 1) と、Σ(2, 3, 5) と Σ(2, 3, 7) の [[Sav12#^sav12-sec-0544]] の値に合う式。本は Work に無い | draft |
| 符号を採らないもの | Regina の `isoSig`、Example3.lens の名前、HC の絡み形式の分数 | draft |

Dedekind 和 s(b, a) は [[Sav12#^sav12-sec-0540]]、印刷頁 200 の定義である。a > 0 で、k = 1 から a − 1 までの和である。a = 1 の和は空で、値は 0 になる。

手術の式は [[Les07#^anonnd70-blk-0314]] の補題である。S³ の結び目では、±1/24 が消えて、p > 0 のとき λ = (q/(2p))Δ''(1) − s(q, p)/2 になる。同じ補題の証明が引用する q > 0 の展開は、別の関数にした。正規化の文は [[Les07#^anonnd70-blk-0138]] で、PDF 頁 10 の図と頁 11 の図の間にある。補題の表示は PDF 頁 36 の図と頁 39 の図の間にある。絡み数の変化は [[Les07#^anonnd70-blk-0191]] で、PDF 頁 16 の図と頁 18 の図の間にある。
## Alexander normalization dictionary

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.

対象は `S³` の結び目の Alexander 多項式である。保存代表 `P` は
[Alexander coordinates](../specs/knot-classical-invariants.md#alexander-coordinates)
の最小指数 0、最高次係数正の規約を使う。手術の `delta2` は
[[Les07#^anonnd70-blk-0313]] の対称代表 `Δ(t)=Δ(t^-1)`、`Δ(1)=1` の
二階微分を使う。[入力条件 CW5](../specs/casson-walker.md#numbered-conditions)
はこの代表に関する条件である。

| 他の規約・位置 | こちらへの写像 | 確認 |
|---|---|---|
| 保存代表 `P`：[knot_classical_invariants.jl](../../src/knot_classical_invariants.jl) の `alexander_polynomial`。手術代表は [[Les07#^anonnd70-blk-0313]] | 有効な結び目多項式で `deg(P)=2d` のとき `Δ(t)=t^(-d)P(t)/P(1)`。従って `delta2=(P''(1)-2d P'(1)+d(d+1)P(1))/P(1)` | 下の五つの PD 入力で、実際の出力を中心化し、対称性、値 1、直接の二階微分とこの式の一致を確認 |

規則 a を決めるテストは
[knot_classical_invariants.jl](../../test/knot_classical_invariants.jl)
の `classical PD bridge and independent local invariants` と
[casson_walker.jl](../../test/casson_walker.jl) の `two surgery writings agree`。
2026-11-25 の direct driver は前者の実際の PD fixture から多項式を計算し、
次の値を得た（`verified-computationally`、この五入力に限る）。

| 結び目 | `deg(P)` | `P(1)` | `delta2` | `casson_walker_knot_surgery(1,1,delta2)` |
|---|---:|---:|---:|---:|
| `0_1` | 0 | 1 | 0 | 0 |
| `3_1` | 2 | 1 | 2 | 1 |
| `4_1` | 2 | −1 | −2 | −1 |
| `5_1` | 4 | 1 | 6 | 3 |
| `5_2` | 2 | 1 | 4 | 2 |

写像自体の微分公式は積の微分から従う（`proved/agent`）。全ての
入力で保存多項式が出典の結び目多項式と同じ `±t^k` 同値類を持つという
一般的同定は `unverified`。五つの独立の期待多項式と整合するため同定を
予想するが、Fox 行の符号・変数を出典と比較した証明が必要である。
非古典的 Gauss 入力や出典と一致しない多項式は、この同定を支持しない。
この変換を自動実行する公開関数は現時点ではない。

<!-- XQQQ1000000 -->
