<!-- XQQQ1500000 -->
# HeegaardCode（HC）: Heegaard のレイヤーのデータ仕様

Status: LC・HC・Heegaard diagram の名称と二レイヤーの区別、§2.1 の行文法（role の明示と全指定）、point・weight を HC に置くことはあかげふが確定（2026-11-09）<br>
Scope: Lickorish Code に Heegaard 固有の情報を加える形式と入力検査

## 1. 結論とレイヤー構造

HeegaardCode（HC）は、Lickorish Code（LC）に α・β の役割、joining disk などを加えたデータである。数学の対象は Heegaard diagram と呼ぶ。語と `place` は [LC 仕様](lickorish-code.md) が定める。用語は [用語正本](heegaard-terminology.md)、spine との相互変換は [Koda 変換仕様](heegaard-koda-conversion.md) が定める。

公開型 `HeegaardDiagram` は `VirtualLink` と同格である。HC は、その図式が持つデータレイヤーである。o-data は `VirtualLink` の表示であり、Koda 変換は `HeegaardDiagram` からその o-data への射である。境界は [HeegaardDiagram の公開境界](heegaard-diagram.md) が定める。

HC は LC を一つ持つ。LC に含まれる K 以外の成分へ、HC が α・β の役割を与える。LC の成分名は文字列として扱う。たとえば LC の `alpha1` という名前も、LC 単独では役割を指定しない。

本書の形式は、HC が使う LC の非 K 成分すべてに α または β の役割を明示する形式である。追加の無役割の曲線を表示したい場合は LC の表示として扱う。HC の数学的対象に使う曲線集合が、形式から一意に決まる。

全指定により、K を忘れると残る曲線が α∪β に一致する。役割の全指定は、あかげふが確定した受入条件である。

## 2. HC の形式

### 2.1 役割を行で明示する

メモリ上の HC は、LC の値、役割の対応、disk、point、weight を持つ。テキスト形式は `.hc` とし、同じ場所を基準とする相対パスで `.lc` を参照する。読込時に参照先を一度読み、HC はその LC の値に結び付く。保存する際も LC と HC の組を保存する。

```text
lc : figure.lc
role c1 : alpha
role c2 : beta
disk D : F3 F4
point z : F3
weight F3 : -1
```

`lc` はちょうど一行とする。`role` は `role <成分名> : alpha` または `role <成分名> : beta` とする。一つの非 K 成分に役割をちょうど一度指定する。各族の曲線の順序は、参照した LC 内の成分の順序とする。K の成分は基準曲線である。

文字コード、空行、コメントの扱いは LC のテキスト形式に従う。`lc` 行のコロンの後は一つの相対パスとし、HC ファイルの親ディレクトリから解決する。役割の値は `alpha` と `beta` の二つである。各 disk 名と各 point 名は、それぞれの名前空間で一意とする。同じ領域への weight 行は一つとする。

`disk`、`point`、`weight` の参照は、その LC の細分領域と曲線の弧に結び付く。LC の値を変える操作では、HC の参照を同じ操作の対応写像で運ぶ。

<a id="22-比較案-成分名の規約から読む"></a>
### 2.2 既存入力の役割を明示する移行

HC を作る操作に限って `alpha1`、`alpha2`、`beta1` などから明示された命名規約に従って役割を作る。この操作は、LC と明示された命名規約を入力し、役割を明示した HC を返す。LC の読込と描画は同じままである。

成分名から読む方式は、既存例を手短に入力できる。明示する方式は、名前が `c1` や `gamma` でも同じ対象を扱え、曲線の改名と役割の変更を分けられる。HC の保存形式は明示方式である。命名規約から作る便利機能は、独立した HC 作成操作として定義する。

移行は LC の語の形式、語の符号、配置アルゴリズムを保持する。

## 3. HC に置く追加情報

### 3.1 disk

`disk` は HC に属する。既存の領域を用いて、曲線で分割された一枚の disk を指す。選択した disk を Koda 変換へ渡すときは、K を忘れた Heegaard diagram で joining 条件を入力検査する。

```text
disk D : F3
disk Dp : F3 F4
disk D : c1+
```

名前は任意とし、異なる名前の disk は別の指定である。複数の disk があるときは、変換要求で使用する名前を指定する。HC は disk を省略した Heegaard diagram も表す。

上の三行は、それぞれ独立した指定例である。同じ名前 D を一つの HC 内で重複定義すると入力エラーになる。

### 3.2 point の配置

`point` は HC に置く。従来の §1.1 の上のレイヤーという配置を保つ。点の支持は、LC の細分領域内、曲線の弧上、空語の曲線上で指定する。

```text
point z : F3
point w : alpha1/+8
point p : alpha1
```

`alpha1/+8` は、その曲線で +8 の次の crossing へ至る弧を指す。三行目は空語の曲線上の指定である。

数学的には、点を印した曲面は α・β を選ぶ前にも定義できる。点を HC に置く根拠は、今回の機能で基点付き Heegaard diagram と計算上の選択に用いることである。K08 §3 の based Heegaard diagram では各 β_i 上の点を選ぶ。HC の一般の `point` が、この特定の基点の条件を満たすかは、使用する計算の入力条件として確認する。

LC の語と配置は点を省略しても定まる。点を HC に置くことで、LC の土台を曲線配置として保てる。一般の点付き曲線配置の保存が必要になった場合は、点を LC の装飾へ置く案を改めて比較する。状態は、数学的な支持の説明が agent の論証、配置があかげふが確定である。

### 3.3 weight の配置

`weight` も HC に置く。指定した領域の整数係数として保存する。

```text
weight F3 : -1
weight c1+ : 0
```

書いたときだけ係数があり、省略は未指定とする。符号や値、支持の指定方法を保持する。全係数が定まった範囲では、向き付き領域 R に対する整数の和 `Σ n_R R` は 2 次元鎖であり、その境界は領域の貼合せから計算できる。

この鎖は一般の曲線配置にも定義できる。一方、Heegaard diagram の領域の係数として使う場合、その境界を α・β の曲線に沿って読む。今回 HC に置く根拠は、この用途と §1.1 のレイヤー構造である。K による細分上の係数は保持し、K を忘れた一つの補領域の係数として使う場合は、その内部で係数が一致することを入力条件とする。

LC の描画は係数から経路を選ばない。将来、一般の曲線配置の整数 2 次元鎖を土台の機能にする場合は、LC の装飾へ一般化する案を比較する。状態は agent の論証、配置があかげふが確定である。

**Evidence.**

- literature: [K08, ^k08-sec-0073]（§3 の based Heegaard diagram）

## 4. HC の入力検査

まず参照する LC の入力検査を一度行う。HC の検査は、役割と上のレイヤーの条件を続けて調べる。

| 段階 | 調べる条件 | 適用範囲 |
|---|---|---|
| HC の形式 | LC の参照が一つ、role の成分が実在し、非 K 成分へ役割が一意に指定される | HC の読込 |
| 追加行 | disk・point・weight の参照先が実在し、弧の符号付き crossing がその成分にある | 記載した行 |
| Heegaard diagram | α・β の各族が埋め込まれた円周の和で、同族内で交わらず、異族間は横断する | HC から数学的対象を作る入口 |
| 切断系 | 曲線系で切った曲面と圧縮後の境界型 | 完全切断系・閉多様体を要する計算 |
| joining disk | 元の曲面内の閉包の境界が単純で、各 α_i・β_j との交わりが一本の非空の連結な弧 | 選択した disk を用いる Koda 変換 |
| standard・flow 条件 | 補領域、spine の sector と正則近傍の条件 | 要求する変換の定義域 |

切断系と joining disk の詳細は [Koda 変換仕様 §5](heegaard-koda-conversion.md#5-heegaard-条件と-joining-disk-の入力検査) に従う。これらは LC 単独の受理条件に加えない。HC の基本的な受理も、すべての算法に必要な追加条件を一律には要求しない。公開入口は、どの条件を使用したかを返す。

既定経路では外部入力の検査だけを行う。構成した対象の独立な照合は DEBUG_MODE とテストで行う。

## 5. 保存と移行

LC は `.lc`、HC の参照と追加行は `.hc` とする。既存入力の語と `place` は `.lc` へ移す。役割を表していた成分名、`disk`、`point`、`weight` を持つ入力は、LC 本体と HC の指定に分ける。LC 内に残した成分名自体から役割を推測しない。

これは役割と追加行を実物のデータから分ける移行である。語と `place` の文字列は保持する。移行後も、HC が参照した LC の値と各追加行の支持を追跡できるよう記録する。

## 6. 確認事項

名称とレイヤーの区別、明示する role 行と `.hc` からの LC 参照、非 K 成分への役割の全指定、`point` と `weight` を HC に置くことは、あかげふが確定した。

Koda の局所表、C と基準曲面への同定の自動選択、spine の一般の対応の証明は、変換仕様の未解決欄に引き継ぐ。
<!-- XQQQ1000000 -->

<a id="comparison-records"></a>
## 比較用の記録と座標

[`HCCanonicalWitness` と `HCIdentityCorrespondence`](../../src/heegaard/hc_compare_canonical.jl) は向き付きの保存語の対応を記録する。`label_map` の定義域と値域は正の交点ラベルで、負の出現には同じ対応を符号付きで使う。`name_map` は入力曲線名から出力曲線名への辞書である。正準形の `k_words` 欄は K の名を保つ。witness の `name_map` は α と β の名をそれぞれ `a1,a2,…`、`b1,b2,…` にし、役割を持つ K 曲線にもこの改名を適用する。役割を持たない K とその他の曲線の名は保つ。したがって、この辞書を直接適用した HC が固定 K の名前を保持するとは限らない。二つの保存 HC の間の対応には `hc_identity_correspondence` を使う。`rotation_starts` は入力曲線名をキーに、対応を適用する前の保存語の巡回開始位置（1 始まり、空語は 1）を記録する。正準形 witness は入力から正準形へ、identity correspondence は source から target へ向く。対応を得られなかったことは非同一の証明ではなく、`HCNotComputed` になる場合がある。

比較の座標順序は §2.1 の入力行の宣言順と区別する。`hc_intersection` と `hc_presentation` は各役割の曲線名を Julia の文字列の辞書式順序で整列する（例えば `a10` は `a2` より先）。交差行列の行は α、列は β であり、提示の生成元は α、関係語は β の順である。正準形そのものの比較は [`HCCanonical` のフィールド順](../../src/heegaard/hc_compare.jl) の辞書式比較で、保持する witness の比較はソート済みの `rotation_starts`、`label_map`、`name_map` の順である。これらは格納形式の比較であって、多様体の標準形を定めない。探索上限に達した場合を異なる対象と判定してはならない。

<a id="intersection-sign-translation"></a>
### 交差行列の符号翻訳

同じ向き付き曲面・同じ向き付き α、β に対して、HC の行列を `A_HC`、順序付き接ベクトル `(Tα,Tβ)` による交差行列を `I_(α,β)` とすると、同じ行列座標では **`A_HC = −I_(α,β)`** である。HC は α 語で正、β 語で負の出現に `+1` を加える。一方、[LC の交点規約](lickorish-code.md) は負側、正側の順を曲面の正の向きとするためである。legacy `Intersection.sign` は後者の順序付き交差符号を用いる。同じ配置を両形式で表す場合は符号を反転し、曲線名と legacy 曲線 ID の対応で行・列を並べ替える。互いに独立に指定した二つの図式の行列へ、この変換を対応なしに適用してはならない。`hc_presentation` のアーベル化行列は行が関係、列が生成元なので `A_HC` の転置になる。

規約の状態: **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**。[HC 比較テスト](../../test/heegaard/hc_compare.jl) の `sign convention and presentation orientation` は両符号とアーベル化の転置を、正準形の改名・巡回移動の試験は辞書の輸送を確認する。[legacy core テスト](../../test/heegaard/core.jl) は正の `Intersection` の行列を確認する。今回の direct driver では両符号、転置、legacy 回転系の dart 被覆を照合し、`l1=[−1]` を β、`m1=[1]` を α とした標準 HC の witness が `l1→b1, m1→a1` を返すことも確認した。規約辞書の対応先は [HeegaardCodes の規約](../conventions/HeegaardCodes.md#comparison-conventions) にまとめる。一般の HC と legacy 図式を幾何的に同定する自動変換の正しさは、この局所符号照合だけでは確立しない。

<a id="reidemeister-character-coordinates"></a>
## Reidemeister torsion の指標座標と正規化

この節は [`hc_reidemeister_torsion` と `hc_torsion_character`](../../src/heegaard/hc_invariants_torsion.jl) が現在扱う、有限 H₁ の指標で評価した Fox 複体の値を定める。閉入力として α、β の本数がともに種数に等しいことと、関係行列の Smith 対角成分がすべて非零であることを検査する。これだけで任意の手製 HC の幾何的妥当性を証明するものではない。

### 指標のキーと基底変更

関係行列 `A` は `hc_presentation` の β 関係語を行、α 生成元を列とする。曲線順序は [比較用の記録と座標](#comparison-records) に従う。Smith 変換を `D = U A V` とする。`invariant_factors` は `D` の正の対角成分 `dᵢ` であり、1 も省かない。`order = lcm(dᵢ)` は H₁ の指数であって、一般には H₁ の位数 `prod(dᵢ)` と異なる。空の対角の場合は 1 とする。

`keys` は各座標 `0 ≤ kᵢ < dᵢ` の直積から零ベクトルを除いた列で、後ろの座標ほど速く変わる。`ζ` を `order` 次の標準の根とすると、キー `k` の指標は、元の α 生成元 `aⱼ` に `φₖ(aⱼ) = ∏ᵢ ζ^((order/dᵢ) kᵢ Vⱼᵢ)` を割り当てる。したがって元の指数和の行ベクトル `x` の Smith 座標は `x V` である。`A V = U⁻¹ D` より、この指標は関係を消す。[proved/agent: by construction、上の行列等式]

生成元を整数 unimodular 行列 `P` で `a′ᵢ = Σⱼ Pᵢⱼ aⱼ` と変更すると、関係行列は `A′ = A P⁻¹`、指標は `φ(a′ᵢ) = ∏ⱼ φ(aⱼ)^Pᵢⱼ` に輸送する。関係の基底変更は左からの unimodular 乗算である。別の入力の Smith キーを比較するには、この生成元の対応と各 Smith 変換を通して指標を同定する必要がある。返却記録は `V` を保存しないので、同じ整数キーだけでは図式間の同じ指標を指定できない。[proved/agent: 行座標の基底変更]

ベクトルによる accessor は格納キーと完全一致するものを選び、剰余を取らない。整数による accessor は非自明な巡回因子が一つだけのとき、その因子に整数を置き、他の座標を零にしたキーを選ぶ。零指標は producer が格納しない。[proved/agent: by construction]

### 値と正規化の曖昧さ

Fox 行列の行は β、列は α とする。現実装の双対類と 2-cell の符号選択に対して、値は
`τφ (φ(αᵢ)−1) εⱼ (φ(βⱼ♯)−1) = (−1)^(i+j+g+1) det φ(minorⱼᵢ)`
を満たす。零でない分母を持つ全 minor の値が一致すること、境界の合成が零であることを検査する。twisted 複体が非 acyclic の指標には `nothing` を格納する。これは計算された複体についての確認である。[proved/agent: by construction]

正規化は `τ₀ = +1`。代表は現実装の順序で最初に境界条件を満たす選択であり、全符号 `+1` の双対類、同じ符号の逆の類、それ以外の符号パターンの順に選ぶ。複数の通る選択の値は、一つの共通の符号 `s ∈ {±1}` と一つの `h ∈ H₁` による全指標に共通の変換 `τ′φ = s φ(h) τφ` で比較する。指標ごとに独立の単元を許す比較ではない。異なる類が残れば `HCTorsionError(:ambiguous_orientation)` となる。格納するのはこの代表であり、図式変更に不変な絶対値ではない。[proved/agent: by construction]

規約の状態: **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**。比較対象は [invariants extra テスト](../../test/heegaard/hc_invariants_extra.jl) の `Reidemeister torsion of lens spaces p ≤ 7`、`saved genus-one diagrams`、`stabilization of L(5, 2)`、`two 2-cell orientations outside one torsion class stop`。安定化の fixture は因子 `[1,5]` と、種数 1 の値に対する符号反転を要求する。種数 1 のキーから安定化後への対応は `[k] → [0,k]`、値の対応は `τ → −τ` である。これらはテストが指定する有限の照合であり、一般の安定化公式の証明ではない。このターンではテスト本文を照合し、実行は lead の fresh-process 検査に委ねる。現実装の関数を用いた direct driver は、関係行列 `[2]`、`[2 2;0 4]`、`[2 2 0;0 4 4;0 0 6]` について `D=UAV` と全 2、8、48 指標による関係の消滅を確認し、因子 `[1,5]` のキー順序と種数 2 の符号順序を確認した。[verified-computationally: 上記の有限な行列と全指標]

unverified: 一般の HC でこの Fox 複体を与えられた向き付き 3-manifold の cellular 複体と同定し、正規化を文献のものと同定する主張。レンズ空間と安定化 fixture の期待値は一致する同定を支持するが、任意の図式の cell の向きと双対類を運ぶ証明が必要である。反例となる対応付き図式、またはその証明で状態を更新する。実装先頭が挙げる Massuyeau, Lemma 5.3 は登録済み cite key による原典照合がなく、この節では `proved/literature` としない。

<a id="alexander-laurent-coordinates"></a>
## Alexander and Laurent coordinates

The coordinate convention is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**. It records the choices in [the Alexander producer](../../src/heegaard/hc_invariants_alexander.jl) and [the positive-Betti-number producers](../../src/heegaard/hc_invariants_positive.jl). Use the relation matrix and Smith convention of [the character coordinates](#reidemeister-character-coordinates): `D = U A V`, with relations in rows and the lexicographically ordered α generators in columns. Let `B` consist of the columns of `V` at zero Smith diagonal entries, changing each column's sign so that its first nonzero entry is positive. An exponent-sum row `x` maps to `x B`. These are coordinates on the free quotient of the presented first homology; torsion is sent to zero. The returned `free_basis` is this matrix of images of the original generators, not a matrix of lifts from the free quotient. [proved/agent: `A B = 0` and the unimodular Smith change of basis]

The standalone univariate Alexander producer uses `t`; the positive-torsion and multivariable producers use `t1,…,tb`, including `t1` when `b=1`. In the formulas below, `t` denotes that sole variable when the free rank is one. The Fox torsion uses the same minor sign and `τ₀ = +1` convention as the character producer. Defined minors must agree up to a signed Laurent monomial. For `b = 1`, the Alexander result is the torsion multiplied by `(t−1)^2`; its `coefficients` run from the lowest supported power to the highest, discard the common lowest power, and have positive last coefficient. Zero has the empty coefficient vector. This normalization does not impose value 1 at `t = 1` or a symmetric choice of powers. For positive torsion, the numerator and denominator are reduced by their gcd, zero is `0/1`, and the denominator sign is selected by the first coefficient in the polynomial library's iteration order. A Laurent monomial ambiguity remains. For `b ≥ 2`, the multivariable Alexander producer requires the reduced denominator to be a Laurent unit (or the numerator to vanish); the polynomial retains the returned free basis. [proved/agent: by construction]

`hc_laurent_gl_matrix` acts on exponent **columns** by `e → M e`. In the row convention above, the corresponding change is `B → B Mᵀ`. Its bounded search compares up to a common sign and monomial; failure at the bound is not inequivalence under all of `GL(b,Z)`. `hc_turaev_matches` compares in one polynomial parent: `numerator*(t−1)^2` with `Alexander*denominator` for `b = 1`, and `numerator` with `Alexander*denominator` for `b ≥ 2`, modulo those units. [proved/agent: by construction]

The [invariants-extra tests](../../test/heegaard/hc_invariants_extra.jl), `positive torsion of S1xS2`, `positive torsion of five 0-surgeries`, `multivariable Alexander of a 2-component 0-surgery and of T3`, and `nonconstant multivariable Alexander for b1 >= 2`, exercise ranks 1, 2 and 3. They require the two producers to return the same free basis, and compare with independent recorded Fox coefficients. `GL(2,Z) search swaps the variables` tests the exponent-coordinate translation. These are the finite comparison obligations embodied by the tests; their full package execution belongs to the lead's fresh run.

unverified: identification of every accepted HC Fox complex with the cellular complex of its claimed manifold, and exact alignment with the Nicolaescu and Massuyeau source normalizations cited in the source files. Agreement is expected on valid diagrams with the stated cell orientations, because the minor identities and the independent surgery fixtures agree with that interpretation; a corresponding diagram with a different torsion class would refute it. A registered-source comparison and a general cellular identification would establish the missing claim. This entry fixes the stored coordinates without promoting that claim to proved/literature.

<a id="linking-generator-conventions"></a>
## Linking generators and orientation

The convention is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**, following [the linking producer](../../src/heegaard/hc_invariants_linking.jl). The reference surface basis is `(l1,…,lg,m1,…,mg)`. Its matrix `Ω` uses the ordered surface intersection pairing from [the sign translation](#intersection-sign-translation). A curve `c` has coordinate column `Ω⁻¹ (I(k,c))ₖ`. The columns of `A` and `B` are the α and β meridians in their respective lexicographic name orders. The producer checks the alternating unimodular reference pairing, primitive isotropic meridian lattices and consistency of all stored curve intersections with these coordinates. [proved/agent: by construction]

For `C = [A B]` and `S = U C V`, let `r` be its rank and `J` the indices of diagonal entries `dᵢ > 1`. The torsion generators are the columns `X = U⁻¹[:,J]`; the free lifts are `U⁻¹[:,r+1:2g]`. A reference-coordinate column `z` has torsion coordinates `(U z)[J]`, reduced modulo the corresponding `dᵢ` when interpreted as a group element. Set `P = V[1:g,J]`, `R = V[g+1:2g,J]`, `E = A P`, `F = B R`. Then `E+F = X diag(d)`. The stored linking matrix is `Lᵢⱼ = (Eᵀ Ω X)ᵢⱼ/dᵢ mod 1` for `boundary_side=:alpha`, and its negative modulo 1 for `:beta`; representatives lie in `[0,1)`. `linking_pairing` evaluates `xᵀ L y mod 1`. Symmetry modulo integers, descent to the stated torsion orders and nonsingularity are checked. [proved/agent: matrix identities and by construction]

Changing generator lifts to `X Q` changes the pairing matrix to `Qᵀ L Q mod 1`, provided `Q` induces an automorphism of the indicated finite group. Thus a cyclic lens comparison with `q/p` allows multiplication by a unit square; comparison without a chosen manifold orientation also allows an overall minus sign. The [linking tests](../../test/heegaard/hc_invariants_linking.jl), `frozen identified L(5,2), available without local corpus`, evaluate the named `m1` class as `2/5` using `U`, rather than guessing the Smith generator. `mixed free and torsion, and two-primary forms` tests a rank-four surface lattice and distinct pairings on `(Z/2)^2`. `registered identified HC inputs` tests the Smith/decomposition identities and boundary-side sign. The test `orientation reversal` compares surface reversal with `L → −L` in the retained lifts. These fixture assertions specify the finite comparison; the optional corpus cases require their input files.

unverified: a general proof that these algebraic coordinates identify the linking pairing of every accepted HC with the geometrically claimed oriented 3-manifold. The primitive Lagrangian construction and independent lens comparisons support that identification; a diagram with an explicit surface/handlebody correspondence and a mismatched pairing would overturn it. The `Invariants.linking_form` HC dispatch uses this producer and inherits this boundary-side convention; it introduces no second generator convention.

<a id="moo-tv-normalization"></a>
## MOO and Turaev–Viro normalization

The stored normalization is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**. For the valid finite pairing `(T,L)` and free rank `b` of the preceding section, [the even-level MOO producer](../../src/heegaard/hc_invariants_linking.jl) evaluates
`Z_N = N^(b/2) / sqrt(|T|) * Σ_(x∈T) exp(2πi (N/2) L(x,x))`,
with positive real square roots and positive even integer `N`. It computes the same value in an exact cyclotomic field and may change the output conductor to a comparison field; the root metadata is not the input level. In particular, torsion-free input gives `N^(b/2)`, the empty genus-zero diagram gives 1, and negating `L` conjugates the sum. [proved/agent: finite-sum formula and by construction]

The [linking tests](../../test/heegaard/hc_invariants_linking.jl), `original MOO surgery normalization and sign`, compare with a separately evaluated surgery sum at levels 2,4,6,8,10,12, allowing one consistent orientation reversal across the levels for an unoriented source identification. The one-component surgery matrix `[3]` gives `−i` at level 2, `[3 1;1 2]` gives `−1`, and the specified three-component chain gives `i`. `free homology and Weeks` and `mixed free and torsion, and two-primary forms` test the free factor and noncyclic torsion. These do not identify orientations of every source table entry.

[The Regina adapter](../../src/heegaard/hc_invariants_tv.jl) retains the backend's exact Turaev–Viro coefficients without multiplying by a normalization factor. It uses `q = exp(2πi/c)` with `c = r` for odd `r` and `parity=false`, and `c = 2r` otherwise; `r ≥ 3`. The pinned backend is Regina 7.4.1. The [TV tests](../../test/heegaard/hc_invariants_tv.jl), `Regina exact Turaev–Viro adapter`, compare both parities at levels 3–6 with the independently recorded formulas `TV(S³)=−(q−q⁻¹)^2/(2r)` and `TV(S²×S¹)=1`, and with the specified lens formulas. Thus `TV(S³)=1/2` at `r=3`; this convention is not the rescaling that makes every `TV(S³)` equal to 1. [proved/agent: adapter preserves the returned coefficients; geometric comparisons are fixture obligations]

[The direct skeleton backend](../../src/heegaard/hc_invariants_tv_skeleton.jl) supports only `r=3`. Its faces are the LC surface faces and the α/β meridian disks; binary face colours satisfy even incidence at every edge, counting repeated incidences modulo two. For each permitted colouring, the selected face/edge/vertex counts give `χ`; its weight is `(-1)^χ` when `parity=true` and 1 otherwise. The value is one quarter of the sum. The empty genus-zero sphere has two states of Euler characteristics 0 and 2 and value `1/2`. The TV tests compare this value with the Regina backend for both parities. [proved/agent: finite combinatorial definition]

unverified: the complete literature identification of these formulas for arbitrary accepted HC inputs. The source files cite MOO/Deloup and Turaev–Viro (1992), but their registered-source alignment is unresolved. The independent surgery/closed-manifold fixtures support the stated intended normalizations; a discrepancy on a corresponding valid input would refute the identification. Backend agreement alone does not prove the HC-to-manifold conversion.

<a id="hc-mirror-convention"></a>
## HC mirror convention

The convention is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**. [The HC mirror](../../src/heegaard/hc_compare.jl) negates every signed passage and retains the traversal order of every curve. By the LC ordered-passage convention this reverses the encoded surface orientation. It preserves names, roles and weights, clears source paths, and transports direct face references by the unoriented cyclic sequence of their crossing labels. Only a unique such sequence supplies a face correspondence. Disk references without a unique target raise an error; a place or point token not recognized as a face key is retained unchanged, so arbitrary arc annotations have no general geometric transport guarantee. [proved/agent: by construction]

This is the operation used by the HC dispatch of `Invariants.mirror`; it differs from reversing selected curve traversals. In common reference coordinates, surface reversal sends `Ω → −Ω` while the reference-coordinate columns of curves remain the same. Consequently the linking construction sends `L → −L mod 1`, subject to identifying any differently chosen Smith generators. The passage negation sends `A_HC → −A_HC`. [proved/agent: ordered-pairing sign and the linking formula]

The [HC comparison tests](../../test/heegaard/hc_compare.jl), `orientation reversal and mirror are not identity`, keep these operations distinct. The [linking tests](../../test/heegaard/hc_invariants_linking.jl), `orientation reversal`, check the negative pairing and conjugate MOO value. The [invariants-extra tests](../../test/heegaard/hc_invariants_extra.jl), `mirror sends lens torsion to the conjugate at the inverse character` and `Alexander polynomial of a mirror`, compare the corresponding invariant transformations on their fixtures. A general manifold identification still requires the diagram and handlebody hypotheses of the relevant conversion; these tests do not certify arbitrary annotations or amphichirality.

<a id="hc-search-records"></a>
## HC search records

The format is **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**. Its exact database columns and constraints are owned by [schema version 1](../../db/hc_search_schema.sql); the serialization and producers are in [hc_search_db.jl](../../src/heegaard/hc_search_db.jl). A stored search state has representation `hc-search`, area 3, schema version 1 and complexity equal to genus. Its identity body is either `{mode: literal, writing: …}` or `{mode: canonical, form: …}`. Direct state construction in canonical mode does not silently fall back if canonicalization fails. The canonical form uses [the comparison convention](#comparison-records); literal identity uses the serialized writing, not original whitespace or file paths. [proved/agent: by construction]

The writing schema is `hc-search-writing/1`. It contains `genus`, an ordered `words` list of `{name,word}`, an ordered `roles` list of `{name,role}`, and the name/value pair lists `places`, `lc_disks`, `lc_points`, `lc_weights`, `disks`, `points`, `weights`. Component order is `hc_component_names`; pair-list keys are sorted as strings. Paths are omitted and load restores them as `nothing`. Stable JSON sorts dictionary keys by their string form, keeps sequence order, encodes symbols as strings, and requires integers to fit Int64. The payload preserves one writing; inserting another writing with the same identity does not replace the first `raw_json`. A canonical identity is an identity of the specified finite HC record, not a manifold classification. [proved/agent: by construction]

Runs store the script path and digest, label, caller scope, counts, UTC timestamps and status `running`, `completed`, `budget-exhausted` or `failed`. Operations store source/optional target states, run, depth, catalog row, engine version `hc-search-engine/1`, candidate JSON, outcome and inverse status. Outcomes are `success`, `rejected`, `error`; inverse statuses are independently `ok`, `mismatch`, `error`, `unchecked`. A successful move remains a successful edge when its inverse check fails. Operation identity is determined by run, source, engine version and candidate JSON; repeated insertion retains the existing row. The candidate records names, locations, vertices, before/after fragments, genus change and selected move data; composite sequences keep step summaries. This is a bounded-search receipt, not necessarily a complete replayable `LCMoveWitness`. [proved/agent: by construction]

An invariant snapshot has separate homology and linking dictionaries. Failure is `{ok:false,error:…}`. Successful homology has `ok`, `free_rank`, `torsion`, `invariant_factors`; the comparison ignores factors equal to 1 by using the normalized `torsion` field. Successful linking has `ok`, `free_rank`, `torsion` and `matrix`, whose exact rational strings are flattened **column first**. Matrix size is recovered as the square root of the entry count. Empty torsion needs no matrix comparison. One cyclic factor of order at most 10000 is compared up to a unit square; larger cyclic orders and multiple factors use literal matrix equality. No orientation-reversing sign is added to this comparison. Failure or a false comparison is not a certificate of different manifolds. [proved/agent: by construction]

The [HC-search tests](../../test/heegaard/hc_search.jl), `canonical key` and `search db round trip and one inverse`, embody canonical relabeling/cyclic-shift invariance, fixed reference names, payload round trip, deduplication by state identity and inverse status. The genus-one search fixture expects nine candidates. General search completeness, graph exhaustion beyond the recorded budgets, and identification of all geometrically equivalent writings are unverified; the catalog and finite budgets make a restricted search the expected interpretation, and a completeness proof would be needed to strengthen it.

<a id="simultaneous-curve-reversal"></a>
## Simultaneous curve-reversal question

question: On the standard one-crossing HC `l1=[−1]` (β), `m1=[1]` (α), selecting both curves in `hc_reverse_orientation(hc,["l1","m1"])` returns `l1=[1]`, `m1=[−1]`. [The operation](../../src/heegaard/hc_compare.jl) negates each passage once, even though both incident curves are selected, and reverses each selected traversal. A geometric reversal of both tangent vectors in the fixed oriented surface preserves their ordered intersection sign, hence would retain the passage signs under [the fixed sign translation](#intersection-sign-translation). The API's intended simultaneous-operation interpretation is not settled here. [verified-computationally: direct driver on this one-crossing HC; proved/agent: the two-vector sign calculation]

open: If the intended operation is simultaneous reversal of the selected curves on a fixed oriented surface, the expected correction is parity of the number of selected incident curves, rather than one negation whenever either is selected. The discrepancy is decisive on the input above. A documented interpretation that also changes the surface orientation or another datum would overturn that assessment and require an explicit translation rule. Until the intended operation is resolved, the current transformation is a recorded code operation, not a general certificate of geometric simultaneous reversal.

## LC の語と描画へのリンク

<a id="heegaardcodehd-一般化-heegaard-図式の正本データ仕様"></a>
[Lickorish Code（LC）: 曲面上の曲線配置の正本データ仕様](lickorish-code.md)

<a id="1-目的"></a>
[1. 目的](lickorish-code.md#1-目的)

<a id="11-レイヤー構造あかげふの決定2026-11-07"></a>
[1.1 レイヤー構造（あかげふの決定、2026-11-09）](lickorish-code.md#11-レイヤー構造あかげふの決定2026-11-09)

<a id="12-検証の方針あかげふの決定2026-11-07"></a>
[1.2 検証の方針（あかげふの決定、2026-11-07）](lickorish-code.md#12-検証の方針あかげふの決定2026-11-07)

<a id="2-lickorish-基準曲線網-k"></a>
[2. Lickorish 基準曲線網 K](lickorish-code.md#2-lickorish-基準曲線網-k)

<a id="3-データ形式gauss-code-形式"></a>
[3. LC のデータ形式（Gauss code 形式）](lickorish-code.md#3-lc-のデータ形式gauss-code-形式)

<a id="31-交点を持たない曲線joining-disk基点領域係数"></a>
[3.1 交点を持たない曲線と place](lickorish-code.md#31-交点を持たない曲線と-place)
[disk・point・weight の書式](heegaard-code.md#3-hc-に置く追加情報)

<a id="4-整合検査描画前すべて決定的"></a>
[4. LC の整合検査（描画前、すべて決定的）](lickorish-code.md#4-lc-の整合検査描画前すべて決定的)

<a id="5-意味"></a>
[5. 意味](lickorish-code.md#5-意味)

<a id="51-操作は前後の語"></a>
[5.1 操作は前後の語](lickorish-code.md#51-操作は前後の語)

<a id="6-描画アルゴリズム規範"></a>
[6. 描画アルゴリズム（規範）](lickorish-code.md#6-描画アルゴリズム規範)

<a id="61-原則"></a>
[6.1 原則](lickorish-code.md#61-原則)

<a id="62-表示の幾何種数-g--2"></a>
[6.2 表示の幾何（種数 g ≥ 2）](lickorish-code.md#62-表示の幾何種数-g--2)

<a id="63-工程"></a>
[6.3 工程](lickorish-code.md#63-工程)

<a id="64-実装への制約"></a>
[6.4 実装への制約](lickorish-code.md#64-実装への制約)

<a id="7-既存データの扱い"></a>
[7. 既存データの扱い](lickorish-code.md#7-既存データの扱い)

<a id="8-例"></a>
[8. 例](lickorish-code.md#8-例)

<a id="81-p28種数-2現行出力から読み取ったもの"></a>
[8.1 P28（種数 2、現行出力から読み取ったもの）](lickorish-code.md#81-p28種数-2現行出力から読み取ったもの)

<a id="82-描画規約の説明図種数-3"></a>
[8.2 描画規約の説明図（種数 3）](lickorish-code.md#82-描画規約の説明図種数-3)

<a id="9-作業分担と自走の条件"></a>
[9. 作業分担と自走の条件](lickorish-code.md#9-作業分担と自走の条件)

<a id="10-未確定事項"></a>
[10. 未確定事項](lickorish-code.md#10-未確定事項)
