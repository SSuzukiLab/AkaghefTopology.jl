<!-- XQQQ1500000 -->
# HeegaardCodes の規約

対象はモジュール `HeegaardCodes`（型 `HeegaardCode`、公開別名 `HeegaardDiagram`）である。このモジュールは `LickorishCodes` の上のレイヤーであり（[モジュール構成 §1](../specs/module-structure.md#1-構成)）、曲面・K・平面図の規約は [LickorishCodes の規約](LickorishCodes.md) に従う。

このページは、このモジュールの規約の一覧である。各項目の定義は、リンク先の仕様の節にある。

| 規約 | 定義の場所 | 状態 |
|---|---|---|
| 名称（LC、HC、Heegaard diagram）とレイヤーの区別 | [用語と記号](../specs/heegaard-terminology.md)、[HC 仕様 §1](../specs/heegaard-code.md#1-結論とレイヤー構造) | あかげふの決定 2026-11-09 |
| 役割を行で明示し、K 以外の成分すべてに役割を指定する | [HC 仕様 §2.1](../specs/heegaard-code.md#21-役割を行で明示する) | あかげふの確定 2026-11-09 |
| `point` と `weight` を HC に置く | [HC 仕様 §3](../specs/heegaard-code.md#3-hc-に置く追加情報) | あかげふの確定 2026-11-09 |
| 曲面の向きは ∂Hα。Hα の外向き法線を先に置くと右手系になる | [Koda 変換仕様 §6.1](../specs/heegaard-koda-conversion.md#61-b-の構成と向き) | 仕様 |
| K08 Fig.2 の局所規約。β が over、α が under（候補 A） | [Koda 変換仕様 §7.3](../specs/heegaard-koda-conversion.md#73-向きの規約と局所符号) | 採用規約は確定 |
| 細分領域番号は LC §3.1 の番号をそのまま使う | [Koda 変換仕様 §3.1](../specs/heegaard-koda-conversion.md#31-有限な回転系) | あかげふの決定 2026-11-16 |
| 比較用の辞書、語の開始位置、曲線名の比較順序 | [HC 比較用記録](../specs/heegaard-code.md#comparison-records) | fixed 2026-11-25 Codex lane lib-api-heegaard (rule a) |
| HC の行列と順序付き交差行列の符号翻訳 | [交差符号の対応](../specs/heegaard-code.md#intersection-sign-translation) | fixed 2026-11-25 Codex lane lib-api-heegaard (rule a) |
| legacy cellular 記録、切断境界・port の番号、表示の語 | [legacy 記録](../specs/heegaard-terminology.md#legacy-cellular-records) | fixed 2026-11-25 Codex lane lib-api-heegaard (rule a) |
| cellular basis の向き、行列の座標と基底変更 | [Koda 変換 §6.5](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates) | fixed 2026-11-25 Codex lane lib-api-heegaard (rule a) |
| Reidemeister torsion の指標座標と正規化 | [HC 指標座標](../specs/heegaard-code.md#reidemeister-character-coordinates) | fixed 2026-11-25 Codex lane lib-api-heegaard (rule a) |

<a id="computed-coordinate-conventions"></a>
## 計算・記録の座標規約

The following entries fix existing implementation choices by rule (a). Each linked specification names the tests that embody the choice and separates the finite comparisons from an unverified geometric or literature identification. Their status is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**.

| Convention | Definition and comparison scope |
|---|---|
| Alexander and Laurent free coordinates and normalization | [HC coordinates](../specs/heegaard-code.md#alexander-laurent-coordinates); the invariants-extra rank-one, rank-two and rank-three fixtures |
| Linking generators and boundary orientation | [HC linking coordinates](../specs/heegaard-code.md#linking-generator-conventions); named lens generator, mixed homology and two-primary tests |
| HC mirror | [HC passage mirror](../specs/heegaard-code.md#hc-mirror-convention); comparison, linking and torsion mirror tests |
| MOO and Turaev–Viro normalization | [HC state sums](../specs/heegaard-code.md#moo-tv-normalization); surgery sums and exact Regina/skeleton comparisons |
| Search identity, payload and snapshot schema | [HC search records](../specs/heegaard-code.md#hc-search-records); canonical-key and database-round-trip tests |
| Lescop matching, edge-side and affine Spinᶜ coordinates | [Lescop coordinates](../specs/heegaard-terminology.md#lescop-coordinates); edge-variation and matching/basepoint tests |
| CMR3 Euler obstruction and second-field records | [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing); frozen Euler classes and printed-weight equations |
| LC isotopy-key normalization | [LC key](../specs/lickorish-code.md#lc-isotopy-key-normalization); local relation and negative-control fixtures |
| Historical reference port records | [Historical translation](../specs/lickorish-code.md#historical-reference-records); modern display tests at genera 1–3, with the dedicated port-record regression still missing |
| Move route points | [Route coordinates](../specs/heegaard-moves.md#route-point-coordinates); embedded-path, ribbon-bank and inverse-slide tests |

<a id="comparison-conventions"></a>
## 比較の規約辞書

| 対応する記録 | 翻訳 | 照合 |
|---|---|---|
| [LC の負側・正側の交点順序](../specs/lickorish-code.md) と [HC 比較行列](../specs/heegaard-code.md#intersection-sign-translation) | 同じ向きと曲線順序で `A_HC = −I_(α,β)` | [hc_compare テスト](../../test/heegaard/hc_compare.jl) の `sign convention and presentation orientation`、両符号の direct driver |
| [legacy `Intersection.sign`](../specs/heegaard-terminology.md#legacy-cellular-records) と HC | 同じ配置の対応を指定して符号反転し、曲線 ID と名前で行・列を並べ替える | [core テスト](../../test/heegaard/core.jl) と上記の両符号照合 |
| HC 交差行列と群表示のアーベル化 | 生成元と関係の座標順が逆なので転置する | [hc_compare テスト](../../test/heegaard/hc_compare.jl) と direct driver |
| 入力の役割行の順序と比較出力の順序 | 名前の辞書式順への行・列の置換。正準形の改名は witness の辞書を使う | [hc_compare テスト](../../test/heegaard/hc_compare.jl) の改名・巡回移動と行順の照合 |
| cellular basis の順序・向きの変更 | [Koda 変換 §6.5](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates) の符号付き基底変更 | [koda_cmr3 テスト](../../test/heegaard/koda_cmr3.jl) の負係数、重複度、chain map と引き戻し |
| torsion の元の α 基底と Smith 基底、安定化前後 | [指標座標と正規化](../specs/heegaard-code.md#reidemeister-character-coordinates) の `x → x V` と指標の輸送。L(5,2) fixture では `[k] → [0,k]`、`τ → −τ` | [invariants extra テスト](../../test/heegaard/hc_invariants_extra.jl) のレンズ空間と安定化、異なる torsion 類の拒否 |
| [Laurent free-basis rows](../specs/heegaard-code.md#alexander-laurent-coordinates) and the exponent-column action | `e → M e` corresponds to `B → B Mᵀ`; compare values modulo a signed monomial. The univariate coefficient vector discards the lowest supported power and fixes the highest coefficient positive | [invariants-extra tests](../../test/heegaard/hc_invariants_extra.jl), variable swap, five surgeries and ranks 2–3 |
| [Named reference curves and Smith linking generators](../specs/heegaard-code.md#linking-generator-conventions) | A reference column `z` maps to `(U z)[J]`; generator change `X → X Q` gives `L → Qᵀ L Q mod 1`. For a cyclic group this is multiplication by a unit square | [linking tests](../../test/heegaard/hc_invariants_linking.jl), the named `m1` class on L(5,2) and two-primary forms |
| [Boundary side and mirror](../specs/heegaard-code.md#hc-mirror-convention) versus the source's unoriented lens identification | Changing boundary side or surface orientation gives `L → −L`; MOO is conjugated. An unoriented source comparison permits one common sign across all tested levels, not an independently selected sign per level | [linking tests](../../test/heegaard/hc_invariants_linking.jl), surgery normalization and orientation reversal |
| [MOO surgery sum and HC pairing sum](../specs/heegaard-code.md#moo-tv-normalization) | The same orientation gives the same value; reversing the identification conjugates it. No extra phase is selected per input. Exact source registration remains unresolved | [linking tests](../../test/heegaard/hc_invariants_linking.jl), levels 2–12 |
| [Regina TV coefficients and the r=3 skeleton sum](../specs/heegaard-code.md#moo-tv-normalization) | Embed the rational skeleton result in the returned cyclotomic field; no rescaling. Both give `TV(S³)=1/2` at r=3. A convention normalized to `TV(S³)=1` would divide by the displayed S³ value, but that rescaling is not performed here | [TV tests](../../test/heegaard/hc_invariants_tv.jl), both parities; the unimplemented rescaling is an algebraic comparison only |
| [Lescop test edge-side labels](../specs/heegaard-terminology.md#lescop-coordinates) and [the HC adapter's helper names](../../src/heegaard/hc_lescop.jl) | Test `left=owner[odd]` is the helper's `right`; test `right=owner[even]` is the helper's `left`. Translate through the dart-to-region map rather than matching these words | [Lescop tests](../../test/heegaard/lescop_theta.jl), `Lescop p1 difference across one beta edge` |
| [[OS04a]] matching/basepoint classes and [stored Spinᶜ coordinates](../specs/heegaard-terminology.md#lescop-coordinates) | The matching difference is represented by negative β-arrival coordinates modulo columns of `I`; the affine origin is the first matching and region 1. Source-to-diagram geometric identification is unverified | [Lescop tests](../../test/heegaard/lescop_theta.jl), `Spin^c labels from matchings and exterior regions` |
| [CMR3 Euler column and second-field weights](../specs/heegaard-koda-conversion.md#cmr3-euler-framing) | The right-hand side is `−c_P`; a convention writing `delta*x=c` uses `c=−c_P`. Mod-two solvability is not an integral framing | [CMR3 tests](../../test/heegaard/koda_cmr3.jl), printed weights and obstructed mod-two-soluble examples |
| [Historical port crossing sign and modern LC passage signs](../specs/lickorish-code.md#historical-reference-records) | Ordered crossing `(a,b,e)` becomes passages `a:−e`, `b:e`; named crossings map by the stated integer-label dictionary and curves are reordered into modern K order | [LC display tests](../../test/heegaard/hc_plot.jl) embody the modern side; the direct historical translation check and missing dedicated regression are stated in the specification |
| [Bare route corner, dictionary dart and face-boundary parameter](../specs/heegaard-moves.md#route-point-coordinates) | Bare corner `c` means dictionary dart `τ(c)` at `u=1/2`. Switching dictionary darts preserves oriented-edge `u`; boundary ranking uses `u` or `1−u` according to dart direction | [LC move tests](../../test/heegaard/lc_moves.jl), embedded paths and ribbon banks |

状態は上表の各定義と同じ **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**。これは既存のコードと試験が用いる選択を固定したものであり、文献にない HC の保存形式を文献の規約と同一視するものではない。legacy puncture の作業座標と LC の平面図、legacy 境界サイクル番号と LC 面番号は、そのまま同一視できない。必要な対応が記録されていない入力間の変換は未確立である。
<!-- XQQQ1000000 -->
