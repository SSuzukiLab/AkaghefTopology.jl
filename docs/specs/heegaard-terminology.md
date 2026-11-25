<!-- XQQQ1500000 -->
# Heegaard の用語と記号

Status: 名称と表記の正本。2026-11-09 のあかげふの決定に基づく。LC と HC の区別は確定。HC の行文法はあかげふが確定。point・weight の HC 所属はあかげふが確定<br>
Scope: Lickorish Code、HeegaardCode、曲面上の曲線配置、Heegaard diagram、切断・描画、Koda 変換<br>
Audience: 仕様策定、実装、文書、試験、図の表示<br>
Owner context: 東京科学大学 情報理工学院 数理・計算科学系 鈴木咲衣研究室

## 1. 目的と正本

土台のデータを **Lickorish Code（LC）**、Heegaard のレイヤーのデータを **HeegaardCode（HC）** と呼ぶ。数学の対象は **Heegaard diagram** と呼び、略さない。この三つの名称と区別はあかげふが確定である。

ライブラリの公開型 `HeegaardDiagram` は `VirtualLink` と同格である。o-data は `VirtualLink` の表示である。Koda 変換は、`HeegaardDiagram` からその o-data への射である。公開境界は [HeegaardDiagram の公開境界](heegaard-diagram.md) が定める。

語と配置は [LC 仕様](lickorish-code.md)、役割と追加行は [HC 仕様](heegaard-code.md)、相互変換は [Koda 変換仕様](heegaard-koda-conversion.md) が定める。o-graph 側の用語は、研究リポジトリ基準の `lib/topology.jl/docs/specs/ograph-terminology.md` に従う。

## 2. LC と HC のレイヤー構造

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| Lickorish Code | LC | 曲面、Lickorish 基準曲線網 K、役割を持たない曲線の語、place を持つ土台のデータ |
| HeegaardCode | HC、`HeegaardCode` | LC に α・β の役割、joining disk など Heegaard 固有の情報を加えたデータ |
| Heegaard diagram | (S;α,β) | 閉連結向き付け可能曲面と、互いに横断する二つのコンパクトな 1 次元部分多様体。各族内の曲線は互いに交わらない。曲線本数は任意 |
| HeegaardDiagram | `HeegaardDiagram` | 認証した Heegaard diagram の公開型。`VirtualLink` と同格 |
| joining disk 付き Heegaard diagram | (S;α,β;D) | Heegaard diagram と選択した joining disk |
| punctured Heegaard diagram | (S;α,β;J) | Koda の joining union J を備える Heegaard diagram |
| 曲面 | S、Σ_g | 曲線を載せる向き付き閉連結曲面。g は種数 |
| 曲線配置 | C | 曲線の向きと交点の巡回順を含む、S 上の実際の配置 |

LC は α・β の役割を持たない。単なるトーラス上の曲線も LC である。語と `place` は LC、役割と `disk` は HC に属する。`point` と `weight` は、現在の上のレイヤーという配置を保つ HC 所属を [HC 仕様 §3](heegaard-code.md#3-hc-に置く追加情報) に記す。この所属はあかげふが確定した。

HC の確定した形式は、LC を参照する一行と、役割を明示する `role` 行、追加行からなる。LC の成分名は文字列であり、LC の読込と描画は名前から α・β を推測しない。

g=0 の S は球面、g=1 はトーラス、g≥2 は種数 g の閉曲面である。ベクトル場、twist front、framing は、これらに付随する別の構造に属する。

## 3. 基準曲線網と曲線

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| Lickorish 基準曲線網 | K | LC 仕様 §2 で固定する `l_1,…,l_g, m_1,…,m_g, n_1,…,n_(g−1)` の和集合。種数 0 では空 |
| 基準曲線 | l_i、m_i、n_i | K を構成する向き付き曲線 |
| 曲線 | curve | 配置の一つの向き付き閉曲線 |
| 成分 | component | 語の形式で曲線一本に対応する項目 |
| 役割 | role | 曲線を α、β などとして用いる上のレイヤーの指定 |
| α 曲線、β 曲線 | α_i、β_j | Heegaard diagramの各族に属する曲線 |
| 空語の曲線 | empty word | 交点を持たず、`place` が位置を指定する曲線 |
| 弧 | arc | 曲線を交点や切断境界で分けた連結な部分 |

Lickorish Code の Lickorish は、K が Lickorish の基準曲線網であることを指す。描く曲線は任意であり、その曲線が生成元であるという意味を持たない。`lickorish_reference(g)` は K を返す関数名として維持する。図の切断を K の l_i に沿って行う場合と、数学的操作として α に沿って切る場合は、切る曲線を明記する。

## 4. 交点、語、符号

文章での名称は **交点** とする。語の要素とその記録をプログラムやデータ形式で指すときは **crossing** を使う。

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| 交点 | x、crossing | 二つの曲線が横断する一つの点 |
| 交点ラベル | k | 一つの交点を指す正整数 |
| 巡回語 | cyclic word | 一つの曲線を向きに沿って辿った符号付き交点ラベルの巡回列 |
| occurrence | occurrence | 巡回語の中の一つの出現位置。交点一つは二つの occurrence を持つ |
| passage sign | +k、−k | LC では曲面上の枠の向きを記録する符号。−k の曲線を b、+k の曲線を a とすると (b,a) が正の枠 |
| 交点の局所符号 | ε(x) | 指定された順序の二曲線の接ベクトルが S の向きと一致すると +1、逆なら −1 |
| `orientation` | `+1` のベクトル | LC の語の保存では全交点で +1 に固定する欄 |

LC の passage sign は曲面上の向きを定める。o-data の passage sign は over/under を定める。それぞれの型と変換境界で意味を明示する。o-data の `orientation` は vertex sign の保存欄である。

曲線上の弧を指定するときは、成分名と符号付き crossing を使う。`alpha1/+8` は、alpha1 の +8 の occurrence から、その次の occurrence に至る向き付き弧を指す。

## 5. 面、細分した領域、補領域

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| 面 | K で切った面 | S\K の連結成分。種数 g≥1 の基準網では円板 |
| K と入力曲線による細分領域 | `F1`、`F2` など | LC 仕様 §3.1 が全曲線の回転系から番号を付ける領域 |
| 補領域 | R | Heegaard diagramでは S\(α∪β) の連結成分 |
| 曲線配置の補領域 | S\C の連結成分 | 役割を指定しない場合にも、どの曲線を除くかを明示して使う語 |
| 正規化した領域 | normalization | 境界の貼付けによる自己同定を解き、各境界の occurrence を分けた曲面 |
| 領域の閉包 | cl(R) | 元の S の中で閉包を取った部分集合 |

面は K で切った単位である。α・β を重ねると、一つの面の中がさらに分かれる。LC の `F` 番号はその細分した領域を指す。`F` 番号と K の面番号は、それぞれの番号体系として扱う。

K を忘れると、細分領域は貼り戻される。補領域は円板、円環、多重境界や正の種数を持つ曲面になり得る。正規化した領域が円板であることと、その閉包が埋め込まれた閉円板であることを区別する。

## 6. joining disk、印、完全切断系

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| joining union | J | #α=#β の図式の補領域の和。∂cl(J) が互いに素な単純閉曲線の和で、各 α_i・β_j との交わりがそれぞれ一本の非空の連結な弧になるもの |
| joining disk | D | joining union が一つの開円板成分である場合 |
| standard な punctured Heegaard diagram | standard | S\(α∪β∪J) の各成分が開円板である図式 |
| disk の指定 | `disk` | 既存領域の和で一枚の disk を指定する HC の行。Koda の joining 条件は使用時に判定する |
| 基点 | z、w、`point`（HC 所属） | 領域内または指定された曲線の弧上・空語曲線上の印 |
| 領域係数 | `weight` | 指定した領域に付く整数。HC 所属。省略時の扱いは HC 仕様 §3.3 に従う |
| 圧縮円板 | Δα_i、Δβ_j | Heegaard diagramに対応する 2-handle の core を、境界が S 上にあるよう延ばした円板 |
| 完全切断系 | complete cut system | 最小の場合、種数 g の曲面を g 本の曲線で切り、連結な種数 0・境界 2g 個の曲面を得る系 |

Heegaard diagramの定義、joining 条件、閉多様体を得る条件、完全切断系の条件は、入力条件として別々に記録する。joining disk と Heegaard は英語表記を使う。

## 7. 切断と描画

| 名称 | 表記 | 定義 |
|---|---|---|
| 切断境界 | cut boundary | 曲面を指定した曲線で切ると現れる境界成分 |
| 対になる切断境界 | l_i⁺、l_i⁻ など | 元の同じ曲線に戻すための貼付け対応を持つ二つの境界成分 |
| 穴 | hole | 図の上で切断境界を表す三角形や円の内部 |
| 切断境界上の通過点 | boundary point | 曲線の弧が切断境界に達する点。貼付け後の対応を持つ |
| 配置 | placement | LC の面ごとの語から決定的に計算した曲線の位置 |
| 描画座標 | coordinates | 配置を図として表示するための導出データ |
| 仕上げ | finishing | 交点構造と貼付け対応を保つ、角の丸め・穴の円化・力学的平滑化 |
| 破綻 | 破綻 | 図から読み取れる対象が入力と異なる、数学的に間違った図。あかげふの定義 |

見づらさは描画の品質の問題として扱う。破綻は、図が表す数学的対象と入力の相違で判定する。

図の説明には穴を使う。曲面の境界数や貼合せを論じるときは切断境界を使う。閉曲面そのものは境界を持たず、切断後の曲面が切断境界を持つ。

曲線上の crossing、切断境界上の通過点、o-graph の local port は、それぞれの対象を明示する。描画用の経路や制御点は、LC の位相的な正本から計算する。

## 8. spine と o-data の対応

| 名称 | 記号・表記 | 定義 |
|---|---|---|
| branched spine | P | 分岐を持つ二次元複体で、対象の三次元多様体がその上へ collapse するもの |
| 特異集合 | S(P) | 局所的に円板でない点の集合 |
| 真頂点 | V(P) | 特異集合の交点。局所 link が K4 の点 |
| 特異辺 | singular edge | 真頂点の間の特異集合の成分 |
| sector | sector | P\S(P) の連結成分 |
| sheet incidence | sheet incidence | 特異辺や真頂点に接する各局所曲面の接続。重複する接続も個別に保持する |
| branch direction | u | K08 の局所二枚側から一枚側へ向く、P の接方向内の法線。J の外向き。区間方向の横断向きとは別の量 |
| strand の向き | tΓ | 局所の二枚側の二つの sector が共通に誘導する境界向き。二枚側が strand の左、一枚側が右 |
| dα、dβ | ±1 | LC の曲線の向きに対する strand の向きの倍率。S の ∂J と LC の向きとの比を c とすると d=c。B の向きは S と同じにする |
| 分岐 | branching | sector の向きと各特異辺での一枚側・二枚側の整合 |
| flow-spine | flow-spine | 閉多様体の branched spine で、正則近傍の垂直境界が一つの円環になるもの |
| o-data | `OData` | `VirtualLink` の表示。normal o-graph の成分語と vertex sign を保存する |
| spine の対応 | correspondence | 真頂点・特異辺・sector・分岐・横断向きと、各貼付けを保つ明示的な写像 |

Fig.2 の局所規約では、strand は二枚側の誘導境界向きであり、branch direction は一枚側を向く。J に沿う branch direction は外向きである。採用する候補 A は B=S、d=c である。S の正の局所枠では β が over、α が under になる。根拠は [Koda 変換仕様 §7.3](heegaard-koda-conversion.md#73-向きの規約と局所符号) の Fig.19 の独立転写と全対応の比較に置く。局所翻訳の論証は agent の証明、14交点・28辺の照合は有限計算である。K08 p.210 の「∂J の逆」という文は、この採用規約と符号が反対になる。

全 strand の同時反転では、over/under を固定すれば vertex sign は同じで、各成分語の巡回方向が逆になる。同じ向きの局所図で over/under を交換する操作では、vertex sign が反転する。

sector は spine の領域である。Heegaard 側の面・補領域との対応は変換結果が返す。保存した入力の読戻し、spine の対応、表す多様体の同値は、それぞれの結果として記録する。

## 9. 操作、比較、検証

| 名称 | 定義 |
|---|---|
| LC の操作 | 曲線配置の before と after、および局所データと対応を返す操作 |
| HC の操作 | LC の操作と、役割・disk などの上のレイヤーの輸送を返す操作 |
| 曲線の同位 | 指定された固定曲面と参照枠の下での同位。比較が固定するものを明記する |
| 曲面同相による対応 | 向きを保つ曲面同相が曲線・交点・印を保存する対応 |
| 多様体の同値 | 構成した三次元多様体の向きを保つ同相による同値 |
| 構成時の対応 | 対象を作る工程と同時に計算したセル・曲線・印の写像 |
| 入力検査 | 公開入口で外部入力が定義域の条件を満たすかを調べる計算 |
| 実装の検査 | 作った結果の再実現や独立な対応照合。`DEBUG_MODE` とテストで実行する |

既定経路の検査は入力検査とする。数学的な正しさの根拠には、適用条件を照合した原典と構成の証明を用いる。図の品質は通常の入口から返す図で判断する。

### 9.1 群と不変量で使う用語

| 名称 | 記号・定義 |
|---|---|
| 交点行列 | A_ij は α_i と β_j の交点の局所符号の和。行を α、列を β とする |
| 群表示 | 生成元と関係語で記述した基本群の表示。切断系など、用いた構成の適用条件を伴う |
| 第一ホモロジー群 | H₁。群表示の可換化など、適用可能な構成から計算する群 |
| matching | 各 α 曲線と各 β 曲線をそれぞれ一回ずつ含む交点の選択 |
| 外部領域 | Lescop の計算などで、指定した点 w を含む補領域 |
| Lescop の Θ | 図式、matching、外部領域など、式が要求する選択に対して計算する量 |

外部領域の選択と、K の面のうち ∞ を含む面の表示は、各工程が選ぶ対象として記録する。入力した多様体の名前、ホモロジーの一致、spine の対応は、それぞれの主張の範囲を伴う。

### 9.2 move とカタログ

| 名称 | 定義 |
|---|---|
| move のデータ | source・target の有限な局所配置、境界の対応、適用条件、輸送、逆方向を持つカタログの行 |
| 適用箇所 | カタログの source を入力の曲面と語へ埋め込んだ具体的な位置 |
| witness | 適用箇所、向き、ラベル、支持、印の像を固定する記録 |
| 共通の書き換え | 固定した witness に従い、外部の語片と target を境界 port で再接続する操作 |
| 対応の輸送 | 曲線、交点、局所曲面と HC の印の支持を、構成時に計算した写像で新しいデータへ運ぶこと |

正本チェックアウトの研究リポジトリ基準 `lib/topology.jl/docs/specs/heegaard-moves.md` が、LC・HC の move の草稿を所有する。型名 MovePattern は有限な語の断片の共通表現を指す。Koda の模型のデータは、三 sheet と B-W の貼付けも持つ別の表である。

<a id="legacy-cellular-records"></a>
## legacy cellular 記録

[`CombinatorialHeegaardCode`](../../src/heegaard/types.jl) は、LC の K や `place` を持たない旧 cellular 経路の曲線・交点記録である。`Curve.family` は α または β、`id` はその族の正整数 ID、`visits` は向きに沿う交点 ID の巡回語である。ここでは交点 ID に交差符号を付けない。符号は `Intersection.sign` に独立して保存し、順序付き `(Tα,Tβ)` が曲面の向きに一致するとき正とする。符号付き LC 語への翻訳は [HC 比較の符号対応](heegaard-code.md#intersection-sign-translation) に従う。コンストラクタが族、正の ID、符号を検査することと、全訪問の incidence が整合することは別であり、後者には検証が必要である。

`RegionPartition.boundary_cycle_blocks` は `boundary_cycles(ribbon_graph(code))` の返却順に付けた 1 始まり番号のリストをまとめ、同じ補領域に属する境界を指定する。リスト内の番号は交点 ID でも LC の細分面番号でもない。回転系の辺は保存順の α 曲線、β 曲線と各語の位置から定まり、辺 `e` の dart `2e−1` は曲線の向き、`2e` は逆向きに対応する。境界サイクルは dart 列の最初の未訪問要素から順に列挙する。`DecoratedHeegaardCode` がこの境界の分割を加え、`CellularHeegaardDiagram` が検証報告を伴う。型の直接構築は幾何的認証を代替しない。公開の `HeegaardDiagram` 別名は現在の HC を指し、この legacy 型とは区別する。

`complement_regions` は各 block を genus 0 の補領域として読む。`ComplementRegion` の境界数はその block の長さであり、この経路の Euler 標数はこの平面補領域の解釈による。任意の境界分割に対して正の genus を持つ補領域まで指定できる記録ではない。

<a id="legacy-cut-indices"></a>
### 切断境界と port の番号

legacy `cut_along_alpha` の境界は保存された α 曲線の順に並ぶ。第 `k` 曲線の plus 側が `2k−1`、minus 側が `2k` で、`pair=k` はこの対の番号、`curve_id` は元の曲線 ID である。plus 側の交点順は保存された訪問順、minus 側はその逆順である。境界上の port はこの境界順と訪問順に正整数 ID を付ける。訪問数 `n` の境界の第 `r` port のパラメータは `(r−1/2)/n` である。各交点の β の出側は交差符号が正なら plus、負なら minus となり、入側は反対側である。β に沿う wire はある交点の出側から次の交点の入側へ向く。

既定の legacy puncture 配置は第 `k` 対の plus を `(k,1)`、minus を `(k,−1)` に置く。この配置は legacy の作業座標であり、[LC 平面図の規約](../conventions/LickorishCodes.md) の半平面や境界名との同定を与えない。両図式を対応させるには、境界の対応と曲面の向きを別に指定する必要がある。

<a id="presentation-word-encoding"></a>
### 群表示の語

`GroupPresentation` の `generators` は順序付き生成元名、各 `relators` の語は `(generator_index, exponent)` の有限列である。添字は生成元列への 1 始まり位置、指数は整数であり、空語は単位元を表す。`abelianization_matrix` の行は関係語、列は生成元で、成分は指数和である。型の直接構築では添字の整合性を認証しない。

legacy `relator` の返す対の先頭は α 曲線の ID であり、`fundamental_group_presentation` がこれを生成元の位置へ変換する。legacy の `beta_ids` は元の β 曲線 ID である。一方 `hc_presentation` の `beta_ids` は比較用 β 順序の位置を表す。群が表す空間は [Koda 変換仕様 §5.1](heegaard-koda-conversion.md#51-曲線系と多様体) の完全切断系の仮定のもとで定まり、任意の手製の語記録それ自体は多様体の認証を含まない。

これらの番号・語の規約の状態: **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**。[core テスト](../../test/heegaard/core.jl) の `Ribbon, surface, and exact algebra` と `Alpha-cut surface and ports` が境界被覆、表示のアーベル化、2 境界と wire の対応を照合する。今回の direct driver でも一交点の両符号の回転系と、正符号の境界 `(1,plus),(2,minus)`、port パラメータ `1/2`、wire `1→2` を照合した。LC の細分面番号との一般の対応は、この legacy 記録には含まれていない。

<a id="lescop-coordinates"></a>
## Lescop combing, edge-side and Spinᶜ coordinates

The stored conventions are **fixed 2026-11-25 Codex lane lib-api-heegaard (rule a)**, following [the cellular evaluator](../../src/heegaard/lescop_theta.jl) and [the HC adapter](../../src/heegaard/hc_lescop.jl). A matching chooses one crossing on every α and β curve. The Lescop diagram quantity additionally chooses the complement region of the exterior point. These data are intended to specify the infinity-combing used in Lescop's construction; they do not contain an explicit vector field or an absolute value of its `p1`. The evaluators for Θ and its edge variation require complete cut systems and an invertible ordered intersection matrix over the rationals. Their returned quantities are diagram calculations on this domain. [proved/agent: checked domain and data dependence]

The matrix `I` has α rows and β columns in the stored legacy curve order, and `J = I⁻¹` has β rows and α columns. Crossing signs are the ordered `(Tα,Tβ)` signs; the HC adapter applies [the sign translation](heegaard-code.md#intersection-sign-translation). Its cellular crossing IDs enumerate the α–β LC labels in ascending order. The adapter forgets roleless curves, merges adjacent LC faces, and numbers the resulting complement regions by the minimum numerical face ID in each merged block. An LC face number must therefore be translated through `face_region`, and an LC crossing label through `label_of_id`.

In the implemented arc sums, endpoints have half weight, interior crossings have weight one, and an arc from a matching crossing to itself is empty. For an oriented β edge from one crossing to the next in its stored traversal order, `lescop_p1_difference` returns `(lk, delta_p1)` with `delta_p1 = 8lk`. The test's side labels are operational: `left = owner[2e−1]`, `right = owner[2e]`, where `owner` maps a ribbon dart to its complement region. The tested identity is `delta_p1 = 4(Θ(right)−Θ(left))`. These names must not be substituted directly for the local variables of `_lescop_arc_regions`: the HC adapter assigns the odd dart to that helper's `right` region and the even dart to its `left` region. The correspondence, including its reversed naming, is retained in the [convention dictionary](../conventions/HeegaardCodes.md#comparison-conventions). [proved/agent: data-coordinate translation]

For the finite homology labels, the quotient is `Z^g / I Z^g`, with coordinates indexed by the stored α curves and relations the columns of `I`. It is not a Smith coordinate vector. If `h = |det I|`, the canonical representative is the lexicographically least equivalent vector in `{0,…,h−1}^g`; the implemented enumeration rejects a box larger than 20000 for genus greater than one. Equality of these coordinates across different diagrams requires a separately supplied basis correspondence. [proved/agent: integer-lattice definition]

For two matchings, traverse each β curve from the first chosen crossing to the second, exclude the departure crossing, include the arrival crossing, and sum the signed α-coordinate arrivals. `lescop_epsilon` is the class of the **negative** of this vector. The reference matching for the absolute stored label is the first lexicographic matching, and region 1 has shift zero. Across a dart-side edge from the operational left region to the operational right region, an α edge contributes `−eᵢ`; the β-edge classes are the unique assignments, if any, making these region shifts a coboundary in the finite quotient. A nonunique or absent assignment raises an error. `lescop_spinc_label` is the canonical class of the matching difference plus the selected region shift. Thus its zero is an affine record choice; it does not designate a canonical Spinᶜ structure of an unmarked manifold. [proved/agent: by construction]

The [Lescop tests](../../test/heegaard/lescop_theta.jl), `Lescop p1 difference across one beta edge`, assert both `8lk` and the Θ side-difference identity, while `Spin^c labels from matchings and exterior regions` checks matching differences, reversed differences, region shifts and HC/legacy coordinate translation, including a genus-two census diagram. `Lescop variation, rotation sign, and singular or multi-boundary regions` exercises the region interpretation beyond a one-boundary disk. These are finite comparison obligations; no absolute lens-space `p1` is assigned by them.

verified-computationally: a direct driver on 2026-11-25 checked matching differences, their reversals, affine labels and the β-edge Θ identity on the tests' genus-one `L(5,2)` input (5 matchings, 5 regions, 5 affine classes, 5 β edges) and genus-two census `L(3,1)` input (5 matchings, 4 regions, 3 affine classes, 6 β edges). This checks the stored coordinates on those two inputs, including the HC-to-cellular translation; it does not identify an absolute `p1`.

unverified: a general geometric identification of all these affine labels and side conventions with the source's infinity-combing and Spinᶜ labels. The matching-difference and basepoint tests support the intended interpretation, but a source-to-cellular correspondence including the side naming is required to establish it; a corresponding diagram with different source differences would refute it. The intended source is [[Les14]] §§1.3–1.5 for Θ, and Theorem 3.4 and Proposition 5.1 for combing variation, as cited by the HC adapter. The intended source for matching/basepoint variation is [[OS04a]] subsection `subsec:SpinCStructures`, lemma `VarySpinC`. These locators do not upgrade the coordinate translation to proved/literature. In particular, equality of Θ or of a `p1` difference is not used to define a Spinᶜ label.

## 10. 推奨する表記

| 使用する表記 | 使用場面 |
|---|---|
| Lickorish Code、LC | 土台のデータ |
| HeegaardCode、HC | Heegaard のレイヤーのデータ |
| Lickorish 基準曲線網 K | 基準曲線の網 |
| Heegaard diagram | 数学的な図式 |
| `HeegaardDiagram` | `VirtualLink` と同格の公開型 |
| 交点 | 本文と図の説明 |
| crossing | 語の要素とそのプログラム上の記録 |
| 面 | K で切った単位 |
| 補領域 | α∪β など、指定した曲線を除いた連結成分 |
| 穴 | 図の上の表示 |
| 切断境界 | 切断した曲面の位相 |
| joining disk | D の名称 |

プログラムの名称変更は、この用語とデータ仕様に合わせる。`lickorish_reference` は維持する。LC と HC の名称とレイヤーの区別は確定している。HC の形式はあかげふが確定した。point・weight の HC 所属もあかげふが確定した。
<!-- XQQQ1000000 -->
