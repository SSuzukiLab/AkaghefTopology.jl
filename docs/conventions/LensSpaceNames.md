<!-- XQQQ1500000 -->
# レンズ空間の名前と向きの規約

状態は fixed 2026-11-23 Claude である。自由度のある箇所は、リポジトリの既存の選択に合わせて機械的に一つに固定した。固定の規則（a：リポジトリの既存の選択、b：主な出典、c：引用文献で多数の選択、d：出典が先に挙げる選択）を各行に書く。他の規約との違いは末尾の辞書にまとめる。

対象は次の二つである。
- `VirtualLinks` の `lens_space_link(p, q)`。
- 登録例 `vl_example` のうち、レンズ空間の名前を持つもの（`lens_2_1`、`lens_3_1_mst2_example23`、`lens_3_2`）。

`VirtualLinks` の他の規約は、このページでは扱わない。向き付きの名前 L(p, q) の基準は、[Casson–Walker の規約](CassonWalker.md) の行「レンズの名前」にある。このページはその行を決め直さない。

## 用語

- **表示。** 閉じた normal o-graph を、実 Gauss 符号と頂点符号の組で与えたもの（`VirtualLink`）。表示 Γ が表す向き付き閉 3 次元多様体を M(Γ) と書く。M(Γ) は [[MST25]] の組合せ (M(Γ), v(Γ)) の多様体の部分である。
- **M^op。** M の向きを逆にしたもの（[[MST25#^mst25-blk-0187]]）。
- **向き付きの名前。** 向きを保つ同相の類を定める名前。**向きを持たない名前**は、同相の類だけを定める。
- **L(p, q)（向き付きの名前として）。** Casson–Walker の規約の行「レンズの名前」に従い、p > 0 で自明な結び目 U の p/q 手術 S³_U(p/q) を指す。[[Sav12]] は逆の規約をとる。L(p, q) を −p/q 手術とし（[[Sav12#^sav12-sec-0291]]）、二つの規約がどちらも使われていると注意している（[[Sav12#^sav12-sec-0273]]）。二つの対応は辞書にある。

## 規則

| 規約 | 内容 | 確かめた出典と検査 | 状態 |
|---|---|---|---|
| 出典の名前 | 出典から来たレンズ空間の名前は、向きを持たない名前として扱う。対象は三つある。(1) [[MST22]] が三葉結び目の 2 重巡回分岐被覆に付けた「M₂(K) = L(3,1)」（`lens_3_1_mst2_example23`）。(2) `lens_space_link(s, t)` の族の名前「L(s, t)」。(3) [[MST22]] と [[MST25]] が L(p, 1) の図に付けた名前 | (1) の出典 [[MST22#^mst22-blk-0154]] と図 [[MST22#^mst22-fig-0155]] は、三葉結び目の左右も L(p, q) の向きの規約も書かない。(2) はライブラリの docstring が MATLAB の `VLExample.mlx` と `C260220LensSpaces.mlx` を出典とし、向きを書かない。MATLAB のファイルは読んでいない。(3) [[MST22#^mst22-blk-0149]] と [[MST25#^mst25-fig-0060]] も向きを書かない | fixed 2026-11-23 Claude（規則 b：出典が向きを書かない） |
| 向きの比較 | 二つの表示の向きの関係は、次の二つだけで定める。一つは、出典が M^op を表すと定めた操作である（[[MST25#^mst25-blk-0188]] の命題 A.1 では、Γ×（`mst_sign_change`）と →Γ*（`mst_reflect_reverse`）がともに M^op を表す）。もう一つは、両方の表示に同じ手順で計算した、向きに敏感な値である。値の手順が使う全体の向きの符号は、関係に入らない | `test/lens_space_orientation.jl`。値はねじれ付き ℤ_N Dijkgraaf–Witten の状態和で、四面体の向きは over の対で与える。全体の符号を逆にすると、すべての値が複素共役になる | fixed 2026-11-23 Claude（規則 b） |
| o-data の向き | 表示 Γ の o-data が表す多様体の向きは、ライブラリの `hc_triangulation` が返す三角形分割の正の頂点枠で読む。Koda 変換仕様 §6.1 はこの枠を右手系とし（SnapPy と Regina の四面体の読み方と同じ）、分岐した四面体の符号 S と o-graph の頂点符号の関係を「頂点符号 = −S」とする | [HeegaardCodes の規約](HeegaardCodes.md) の行「曲面の向きは ∂Hα。Hα の外向き法線を先に置くと右手系になる」と [Koda 変換仕様 §6.1](../specs/heegaard-koda-conversion.md#61-b-の構成と向き) に従う。このページでは決め直さない。値の手順の over の対による枠は、この枠の逆であることを全例で確かめた（値がちょうど複素共役になる） | fixed 2026-11-23 Claude（規則 a：Koda 変換仕様 §6.1） |
| 向き付きの名前の付け方 | 族を基準にする。M(`lens_space_link(p, q)`) を L(p, q) と呼ぶ。出典の向きを持たない名前と違う場合は、辞書の行で読み替える | リポジトリには既にこの基準がある。`framed_lens_space_input` の docstring の「local all-positive `lens_space_link(p,q)` orientation convention」（PJ260829）である。下の O1 の判定により、D1 と D2 のもとで M(`lens_space_link(p, q)`) = S³_U(p/q) となり、Casson–Walker の規約の L(p, q) とも一致する | fixed 2026-11-23 Claude（規則 a：`framed_lens_space_input` と Casson–Walker の規約） |

採らなかった読みは二つある。一つは、O1 が決まるまで向き付きの名前を付けない読みである。もう一つは、[[MST22]] の例を基準にして M(`lens_3_1_mst2_example23`) を L(3, 1) と呼ぶ読みである。後者は族の `lens_space_link(p, q)` を L(p, −q) と呼ぶことになり、リポジトリの既存の基準と Casson–Walker の規約の両方と逆になる。

## 帰結

すべて `test/lens_space_orientation.jl` で確かめる（103 項目、Julia 1.12.6、`using Topology` のあと 34 秒）。スレッドの記録は `threads/PJ261120/code/checks/lens_orientation_dw.jl` と `lens_orientation_absolute.jl` にある。1 から 4 は向きの関係で、枠の選び方によらない。5 は行「o-data の向き」に依存する。

1. **例 23 は L(3, 2)。** M(`lens_3_1_mst2_example23`) ≅ M(`lens_space_link(3, 2)`) ≅ M(`lens_space_link(3, 1)`)^op で、同相は向きを保つ。よって例 23 の向き付きの名前は L(3, 2) である。[[MST22]] の向きを持たない名前 L(3, 1) とは矛盾しない。
2. **登録例 `lens_3_2`。** MATLAB の `VLExample.mlx` の引数をそのまま使ったこの登録例は、M(`lens_3_2`) ≅ M(`lens_space_link(3, 2)`) を満たす。
3. **族はレンズ空間の分類に従う。** qq' ≡ 1 (mod p) なら M(`lens_space_link(p, q)`) ≅ M(`lens_space_link(p, q')`)、また M(`lens_space_link(p, p − q)`) = M(`lens_space_link(p, q)`)^op である。p = 3, 6（ℤ₃）と p = 7（ℤ₇）で確かめた。これらの p では値が実数でなく、M と M^op を区別する。ℤ₃ の値は L(9, q) では 1 になり向きを区別しないので、p = 9 は検査に入れていない。
4. **MST25 の二つの操作。** Γ× と →Γ* は、登録閉例 9 件と上の族のすべてで値を複素共役にする。値が [[MST25]] の意味での向きの反転を検出することの確認である。
5. **標準の商との比較（行「o-data の向き」に依存）。** 標準の商を S³/⟨g_{p,q}⟩ と書く。g_{p,q}(z₁, z₂) = (ζz₁, ζ^q z₂)、ζ = e^{2πi/p} で、向きは ℂ² の単位球面を単位球の境界として外向き法線を先に置いて定める。この商の値を、二つの p 角形の結の三角形分割を作用で割ったものから同じ手順で計算した。行「o-data の向き」の枠のもとで、M(`lens_space_link(p, q)`) は S³/⟨g_{p,−q}⟩ と同じ値をとり、S³/⟨g_{p,q}⟩ とは異なる値をとる（p = 3, 6 で ℤ₃、p = 7 で ℤ₇、すべての q）。M(`lens_3_1_mst2_example23`) は S³/⟨g_{3,1}⟩ と同じ値をとる。商の値は、p = N で閉じた式 (1/N) Σ_b ζ_N^{q b²} に一致する。

## O1 の判定

O1 は「M(`lens_space_link(p, q)`) が S³_U(p/q) と S³_U(−p/q) のどちらか」という問いである。帰結 5 により、o-data から次のように決まる。依存は二つある。

- **D1：行「o-data の向き」。** `hc_triangulation` の正の頂点枠が右手系であること。リポジトリの仕様が定める規約で、このページはそれに従う（規則 a）。出典との比較は辞書の最後の三行にある（`threads/PJ261120/docs/notes/2026-11-23-frame-comparison.md`）。[[MST25]] の符号とは一致する。[[BP97]] 図 2.7 は二通りに読め、一方では一致し、他方では逆になる。どちらの読みかは、所蔵の図の解像度では矢印の向きが読み取れず、決まっていない。
- **D2：商と手術の対応。** S³/⟨g_{p,q}⟩ = S³_U(−p/q)（向きを保つ同相）。[[Sav12]] は L(p, q) を −p/q 手術とし（[[Sav12#^sav12-sec-0291]]）、Raymond と Hirzebruch–Neumann–Koch の向きの規約に従うと書く（[[Sav12#^sav12-sec-0273]]）。ただし、その L(p, q) が商 S³/⟨g_{p,q}⟩ であるという文は、所蔵の範囲に無い。見込み：成り立つ見込みが高い。商が −p/q 手術であることは標準的な事実で（記憶。Gompf–Stipsicz は未所蔵）、[[Sav12]] の規約とも一致する。覆るのは、所蔵した原典が商を +p/q 手術と書いている場合である。

D1 と D2 のもとでの結論：
- M(`lens_space_link(p, q)`) = S³/⟨g_{p,−q}⟩ = S³_U(p/q)。族の名前は Casson–Walker の規約の L(p, q) と一致する。
- M(`lens_3_1_mst2_example23`) = S³/⟨g_{3,1}⟩ = S³_U(−3)。これはこのページと Casson–Walker の規約では L(3, 2)、[[Sav12]] の規約では L(3, 1) である。

## 辞書（他の規約との対応）

| 他の規約 | 所在 | このページの規約への写像 | 検査 | 状態 |
|---|---|---|---|---|
| [[MST22]] の例 23 の「L(3, 1)」（向きを持たない名前） | [[MST22#^mst22-blk-0154]] | 向きを保って L(3, 2) | 帰結 1（`test/lens_space_orientation.jl`） | 確認済み |
| MATLAB の族の名前「L(s, t)」と、[[MST22]]・[[MST25]] の L(p, 1) の図の名前（向きを持たない名前） | `lens_space_link` の docstring、[[MST22#^mst22-blk-0149]]、[[MST25#^mst25-fig-0060]] | 同じ (s, t) の L(s, t)。向きは出典に無いので、写像は同相の類の上だけで定まる | 出典が向きを書かないことを確認 | 確認済み |
| 標準の商 S³/⟨g_{p,q}⟩（ℂ² の向き） | 帰結 5 | L(p, −q)（q ↦ −q） | 帰結 5（`threads/PJ261120/code/checks/lens_orientation_absolute.jl`、p = 3, 6, 7） | D1 のもとで確認済み |
| [[Sav12]] の L(p, q)（−p/q 手術） | [[Sav12#^sav12-sec-0291]] | L(p, −q)（q ↦ −q） | D2 による | 見込み高い（D2） |
| Casson–Walker の規約の L(p, q)（p/q 手術） | [Casson–Walker の規約](CassonWalker.md) の行「レンズの名前」 | 恒等 | O1 の判定による | D1 と D2 のもとで成り立つ |
| 値の手順の、over の対による四面体の枠 | `test/lens_space_orientation.jl` | ライブラリの枠の逆。値は複素共役になる | 登録閉例と族の全例 | 確認済み |
| [[BP97#^bp97-sec-0550]] 図 2.7 を「左手系のねじを双対四面体の向きとする」と読んだ向き | [[BP97#^bp97-sec-0550]] | 逆。枠つきの順 (0, 1, 2, 3) の符号が −1 になる | `threads/PJ261120/code/frames/compare_frames.jl`（読み bp-screw-is-hand）。登録閉例の 26 頂点と一辺ずらし 352 頂点で一様に逆 | 図の読みが未決（下の行と二者択一） |
| [[BP97#^bp97-sec-0550]] 図 2.7 を「水平面の通常の向きを持つ ℝ³ の中のねじ」と読んだ向き | [[BP97#^bp97-sec-0550]] | 恒等 | 同じスクリプト（読み bp-screw-in-R3）。同じ頂点で一様に一致 | 図の読みが未決（上の行と二者択一） |
| [[MST25]] の脚注「正の交差は負の向きの分岐四面体に対応する」 | [[MST25#^mst25-blk-0014]]、型の文 [[MST25#^mst25-blk-0050]]、図 [[MST25#^mst25-fig-0051]] | 恒等（GKS の頂点順での符号 = −交差符号） | 同じスクリプト（読み mst-footnote）。同じ頂点で一様に一致 | 確認済み |
<!-- XQQQ1000000 -->
