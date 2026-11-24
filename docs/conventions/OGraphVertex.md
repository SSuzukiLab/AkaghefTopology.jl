# o-graph の頂点・branching・流れの規約

状態は draft 2026-11-24 Claude である（行ごとの状態は表に書く）。書き方は [[Kod07]] §2.1 と [[K08]] §3・§5 に倣い、規約を図ではなく枠と向きで定める。符号の細部は [[MST25]]（MST3）に合わせる。固定の規則（a：リポジトリの既存の選択、b：主な出典、c：引用文献で多数の選択、d：出典が先に挙げる選択）を各行に書く。他の規約との違いは末尾の辞書にまとめる。

対象は `VirtualLinks` の閉じた normal o-graph の頂点符号、branching、流れの向きである。o-data が表す多様体の向きは [レンズ空間の名前と向きの規約](LensSpaceNames.md) の行「o-data の向き」に従い、このページでは決め直さない。Heegaard 図から branched spine を作る向きは [Koda 変換仕様 §6.1・§7](../specs/heegaard-koda-conversion.md#61-b-の構成と向き) に従う。

## 定義

以下、$M$ は向き付き閉 3 次元多様体で、$T_xM$ の基底 $(e_1,e_2,e_3)$ が $M$ の向きと両立するとき右手系と呼ぶ。$P\subset M$ は branched special spine である。

1. **branch locus と branch direction。** $S(P)$ を、どの近傍も円板でない点の集合とする。$S(P)\setminus V(P)$ の各成分に、$P$ の中の法線ベクトル $b$ を、局所的に一枚の側を向くように付ける。これを branch direction と呼ぶ（[[Kod07]] §2.1、[[K08]] §5 と同じ）。
2. **流れ。** $P$ は横断的に向き付けられている。すなわち、$P$ の正則近傍の区間束の繊維に大域的な向き $X$ が与えられている。$X$ を流れと呼ぶ。flow-spine では $X$ はその流れであり、Koda の構成では $X$ は Heegaard 曲面の正の方向（$\alpha$ 側から $\beta$ 側への向きを含め、[Koda 変換仕様 §6.1](../specs/heegaard-koda-conversion.md#61-b-の構成と向き) の外向き法線）で与えられる。
3. **sheet の上下。** 真の頂点 $v$ の近傍では、$P$ の sheet は $X$ に沿って高さの順に並ぶ。$X$ の下流側にある sheet を上の sheet と呼ぶ。
4. **strand の向き。** $S(P)\setminus V(P)$ の各成分には、隣接する 2-cell の向きが一致して誘導する向きを付ける。この向きの接ベクトルを $t$ と書く。
5. **over と under。** 真の頂点 $v$ を通る二本の特異線のうち、上の sheet が付く線を over、もう一本を under と呼ぶ。それぞれの $v$ での $t$ を $t_{\mathrm o}$、$t_{\mathrm u}$ と書く。
6. **頂点符号。** $\varepsilon(v)=\operatorname{sign}\det(t_{\mathrm o},t_{\mathrm u},X)$ と定める。行列式は $T_vM$ の右手系で測る。$\varepsilon=+1$ の頂点を型 $+$、$-1$ の頂点を型 $-$ と呼ぶ。この符号の向きは MST3 の型 $\pm$ に合わせて選んだ（下の行「頂点符号の正規化」）。

定義 6 は図の向きを使わない。上の sheet を $X$ の下流側で測るので、同じ spine の流れを反転すると over と under が入れ替わり、$\det(t_{\mathrm u},t_{\mathrm o},-X)=\det(t_{\mathrm o},t_{\mathrm u},X)$ によって頂点符号は変わらない。頂点符号が変わるのは、同じ見た目の図を、流れの約束が逆の二つの文献で読んだときである。その差は辞書の列「流れ」で扱う。

## 規則

| 規約 | 内容 | 確かめた出典と検査 | 状態 |
|---|---|---|---|
| branch direction | 定義 1。二枚側から一枚側へ | [[Kod07]] §2.1、[[K08]] §5（p. 209） | fixed 2026-11-24 Claude（規則 b） |
| MST3 の branching | $-b$（一枚側から二枚側）。MST3 の B3「strand を正方向に見て branching は右から左」は、定義 4 の $t$ と $-b$ の関係として読む | `docs/specs/ograph-source-conventions.md` §2 の B3 | fixed 2026-11-24 Claude（規則 a） |
| 流れ | 定義 2。flow-spine の流れ、Koda の構成では Heegaard 曲面の正の方向 | [Koda 変換仕様 §6.1・§7](../specs/heegaard-koda-conversion.md#61-b-の構成と向き)（β 側が手前、β が over）。MST3 の図はベクトル場を可視側から図の上へ出す向きに描く（[[MST25#^mst25-blk-0045]]、[[MST25#^mst25-fig-0046]]） | fixed 2026-11-24 Claude（規則 a） |
| 上の sheet | 定義 3。流れの下流側 | MST3 の B1「標準局所模型の upper local sheet ↔ over」の「上」は、流れが可視側へ出る図での上であり、流れの下流側に当たる（研究レーンの図の照合、`.research-scratch/PJ260816/grok/qhs2/flow_direction.md`） | fixed 2026-11-24 Claude（規則 a） |
| strand の向き | 定義 4 | MST3 の B3、[[BP97]] Figs. 3.4–3.5 | fixed 2026-11-24 Claude（規則 a） |
| 頂点符号の正規化 | $\det(t_{\mathrm o},t_{\mathrm u},X)>0$ を型 $+$ とする | MST3 Fig. 2.6（[[MST25#^mst25-fig-0051]]）の型 $+$ は、流れを読み手の側に向けて over が左、under が下で、$\det((-1,0,0),(0,-1,0),(0,0,1))=+1$。型 $-$ は under が上で $-1$。[[MST22]] Fig. 4.11（[[MST22#^mst22-fig-0131]]）も同じ。lib の自己検査（over $=(-1,0)$、under $=(0,-1)$ を型 $+$）と一致 | fixed 2026-11-24 Claude（規則 a） |
| 全正の abalone の $\theta_G$ | 語 $[+1,-1]$ の combing の [[Les12]] の $\theta_G$ は $+2$。流れの箱の各円板の横断の向きは joining disk が決める branching から取る（$D_\beta$ は $\beta$ に $-t_\beta$ を誘導する） | 独立の判定役 4 つ（2026-11-24、`.research-scratch/PJ260816/adjudicate2/`）。[[IIKN19a#^anonndu-blk-0390]] 補題 81 と系 82 | fixed 2026-11-24 Claude（規則 a） |

## 帰結

1. 頂点符号は流れの反転で変わらない（定義 6 の後の注）。
2. 全正の abalone の $\theta_G=+2$ と $C_H=q^2$ から、状態和の補正の法則は $C_H=q^{1+p_1/2}\,|H_1|\,\tau_{\mathrm{pert}}(M;q^2)$ の形になる（PJ260816 の要約 §4）。
3. 全正の abalone を $n_+$ 個、全負を $n_-$ 個、`connected_sum!` でつないだ $S^3$ の図の combing は、Gompf の接続和の公式 $\theta(\xi_1\#\xi_2)=\theta(\xi_1)+\theta(\xi_2)+2$ から $\theta_G=4n_+-2$ を持つと予言され、$C_H=q^{2n_+}$ がこれに合う（$n_++n_-\le2$ の 14 図、和を取る位置 3 通り、$\hbar^6$ まで。PJ260816 の要約 §4）。

## 辞書（他の規約との対応）

列「流れ」は、その文献の図で流れが読み手の手前を向くか奥を向くかである。空欄は照合待ちである。

| 他の規約 | 所在 | branch / branching | strand | 枠の手 | 流れ | このページへの写像 | 状態 |
|---|---|---|---|---|---|---|---|
| MST3（[[MST25]]）の型 $\pm$ | §2.1–2.3、Fig. 2.6 | branching $=-b$ | 同じ | 右手のねじ | 読み手の側（図の上） | 恒等 | 確認済み（図の照合） |
| [[Kod07]] の $l/r$ と branch direction | §2.1、Figs. 3–4 | $b$ と同じ | 逆 | | | strand を反転 | 照合待ち |
| [[K08]] の branch direction と Heegaard 図 | §3、§5、Fig. 12 | $b$ と同じ | 同じ | 左手系 $(\mathbf u,\mathbf v,\mathbf w)$ | Fig. 12 に横断の矢印は無い。どちらの向きで読んでも型は同じ（上の図が型 $-$、下の図が型 $+$） | 枠の手を反転（Koda 変換仕様 §6.1 の「頂点符号 $=-S$」） | 図は確認済み |
| [[IIKN23]] の $\ell/r$ | Fig. 1（[[IIKN23#^iikn23-fig-0048]]）、Fig. 3（[[IIKN23#^iikn23-fig-0071]]） | | | 右手系（[[IIKN23#^iikn23-blk-0047]]） | $\partial/\partial z$、可視側・図の上 | $r$ は型 $+$、$\ell$ は型 $-$ | owner 2026-11-23。図の照合で確認（研究レーン、見積もり 85%） |
| INS24 の $\ell/r$ と passage sign | `docs/specs/ograph-source-conventions.md` の INS24 行 | | | | Fig. 6 は前進側を図の下に置く | INS24 の $\ell\mapsto-1$ は chord の符号で、頂点符号ではない。passage の $p_{\mathrm{MST}}=-p_{\mathrm{INS24}}$（owner 2026-09-08）とは別の軸で、二つは両立する | 図は確認済み |
| [[Ish06]] の $\ell/r$ | Fig. 2、§1 | | 同じ | | $X=\partial_z$ | 照合待ち | 照合待ち |

## 検査

検査は lib のテストに置く。

1. [`test/ograph_vertex_convention.jl`](../../test/ograph_vertex_convention.jl)（integration 階）は、lib の 18 図（63 頂点）で、Koda の構成の座標から定義 6 の $\det(t_{\mathrm o},t_{\mathrm u},X)$ を計算し、o-data の頂点符号と一致することを確かめる。一致は二つの部分からなる。入力の符号が頂点の局所模型を選ぶので、$d_\alpha=d_\beta=+1$ での一致は表の書き換えであり、テスト名に `tautological` と書いた。独立な部分は、joining disk の境界から取った $d_\alpha d_\beta$ がすべての頂点で $+1$ であることで、このとき $\det=(\text{入力の符号})\,d_\alpha d_\beta$ が入力の符号に等しい。
2. [`test/ograph_abalone_theta.jl`](../../test/ograph_abalone_theta.jl)（unit 階）は、全正の abalone の流れの箱が箱であること（軌道がスパインに当たり、joining disk の branching に対する横断の符号が一定）と、二つの正則値の前像の絡み $k=+1$ から $\theta_G=-2+4k=+2$ を確かめる。較正は $x\mapsto ix$ の $k=0$（$\theta_G=-2$）、鏡像は $-2$。場は判定役が確かめた座標の模型である。その分岐曲面は lib の Koda の spine と同じ組合せ（辺 2 本、面 2 枚、境界語 $[-1,2,1,1,2]$ と $[-2]$、入射符号 $[-1,+1,+1]$）を持つことを、研究レーンが面ごとに照合した（2026-11-24、`.research-scratch/PJ260816/grok/qhs2/abalone_spine_match.md`）。この照合は lib のテストにはまだ入れていない。
