# Heegaard diagram から閉三角形分割と SnapPy への対応

公開入口は `Topology.Heegaard` に置く。パスの基準は研究リポジトリ根。
数学的入力は既存の [Koda 変換仕様](heegaard-koda-conversion.md) に従う。

## 入力と定義域

`Topology.HeegaardDiagram`（実体は `HeegaardLayer`、解決済みの HC）を
`hc_triangulation` に渡す。`hc_to_odata` の既存の定義域検査と joining disk
選択をそのまま使い、その出力 `VirtualLink` を `disk_tetrahedra` に渡す。
後者は非空の閉 normal o-graph を要求する。`VirtualLink` からの直接入口もある。
旧セル型を新しい公開型の土台にはしない。

joining disk は `joining_disk`、`disk_name`、`disk_id` で既存の選択器へ渡す。
省略時の選択・曖昧性の扱いも既存の選択器のままとする。単一 disk の経路であり、
joining union の一般の円環成分、真頂点なし、非円板 sector を補完しない。
入力を書き換えて joining disk を生成する機能ではない。

閉性は `disk_tetrahedra` の C1–C3 実装 witness に依存する。この witness は一般の
定理を機械証明したものではない。返却時には面の使用回数、商の頂点・辺・面、χを
計算する。ここで扱う closed normal o-graph の双対は、球面境界を塞いだ一頂点の
有限三角形分割である。頂点数1、χ=0を要求する。

HC の語の `hc_first_homology` は α、β の本数が種数に一致する条件を使うが、
完全切断系は単独では認証しない。カタログでの比較はまずその表示の可換化と
有限群への準同型数の比較であり、これだけを任意の HC の多様体認証と呼ばない。

## 面写像と向き

`src/ograph_tetrahedra.jl` が既に使う `sigma`、`source_rho`、`target_rho` から、
同じ検査ループ内で `face_gluings` を返す。既存の返却欄を保持し、面写像の欄を追加する。
橋はこの検査済み写像を読む。古い本体を読み込んだ場合は `HCTriangulationError` の
`face_gluings_unavailable` で止める。

四面体は o-graph の真頂点に一対一対応する。局所頂点は 0,1,2,3、辺は
(01,02,03,12,13,23)、向きは小さい頂点から大きい頂点。四面体の向きの符号は
既存規約の `-vertex_sign`。`permutation[i+1]` は頂点 i の像であり、省かれる
頂点も含めた4頂点の置換を保持する。面の残る3頂点の順序保存は既存の検査から
来る性質であって、橋で新しく仮定して作る写像ではない。

JSON の四面体番号・局所頂点番号は0始まり。元の o-graph の辺・頂点・disk 番号は
元ライブラリの1始まりを保持する。各貼り合わせを一方向だけ保存し、両側の面が
それぞれ一度だけ使われる。自己隣接と重複辺を潰さない。

## Julia の返却値と入口

`HCClosedTriangulation` の欄は次の通り。

| 欄 | 内容 |
|---|---|
| `tetrahedra` | 四面体数 |
| `gluings` | `tet, omitted, to_tet, to_omitted, permutation` と元の辺・strand・crossing の対応 |
| `signs` | 各四面体の向き |
| `isosig` | 有限三角形分割の Regina 互換 isoSig。明示的な `signature=false` のときだけ `nothing` |
| `joining_disk` | 実際に使った選択。o-data 直接入力では `nothing` |
| `correspondence` | Koda の対応・選択、四面体→真頂点、局所辺→disk、元の o-data |
| `topology` | 商の頂点・辺・面数と Euler 標数 |

`hc_triangulation(hc; joining_disk=..., signature=true)` が標準の入口。
`signature=false` では Julia だけの組合せ部分を返す。isoSig の新実装は作らず、
標準時は t3mlite の `Mcomplex.isosig()` を利用する。

`hc_snappy(hc; joining_disk=..., references=["m003(-3,1)"])` は Julia で
三角形分割を作り、Pythonを一度呼び、isoSig と解析結果の辞書を返す。
`hc_snappy(triangulation)`、`hc_snappy_compare(left,right)` も用意する。
`hc_triangulation_data` は JSON に渡す明示的なデータを返す。

## Python と任意依存

`lib/topology.jl/scripts/heegaard/hc_snappy.py --input <JSON>` が CLI。
`--reference` は複数指定可能、`--isosig-only` は構造だけ、`--export-catalog` は
`{"catalog_name":"m003(-3,1)"}` を読む。`--compare` は `left, right` の
二つの三角形分割を比較する。

環境は Python 3.12、SnapPy 3.3.2。通常の起動は
`uv run --no-project --python 3.12 --with snappy==3.3.2 python <script>`。
Juliaの既定 runner はさらに `--offline` を付ける。呼出し時にネットワークから
依存を取得しない。別の既存環境を使う場合は `runner::Cmd` を明示する。
`Project.toml` は変更しない。Regina 本体のインストールは不要。

有限三角形分割を t3mlite に読み、球面の頂点 link とχ=0を確認する。
`snappy_triangulation(remove_finite_vertices=true)` が持ち込む Dehn filling を
そのまま保持する。機械的に全カスプへ (1,0) を上書きしない。未充填カスプが
残れば対象が変わるので拒否する。

結果には `solution_type`、`volume`、H₁の elementary divisors・自由階数・torsion、
基本群の生成元・関係語・表示、Z₂・Z₃・Z₅・S₃への準同型数、`identify()`、参照との
`is_isometric_to` を含む。各計算は `status=ok/error/not_computed` と実際の値または
例外を持つ。準同型は共役類で割らずに数え、探索上限に達した値は返さない。
stdout は単一JSON、診断は stderr、全体の失敗は終了コード2。

Julia は欠落環境を `SnapPyUnavailable`、JSON破損・時間切れ・backend の失敗を
`SnapPyExecutionError` で返す。既存の Koda の例外は保持する。外部呼出しには
既定120秒の期限を置き、自分の呼出したプロセスだけを停止する。

## 同定の限界

数値的な solution type と体積は双曲性の厳密な証明ではない。負向きの四面体を含む
数値解も一律に非双曲とは判定しない。`identify()` の空配列、例外、等長探索の失敗を
「異なる多様体」へ変えない。有限群への準同型数と H₁ の一致も同相の証明ではない。
非双曲多様体は体積と等長判定では分類できず、未確定になり得る。

isoSig の一致は三角形分割の組合せ同型を表す。isoSig の不一致は多様体の相違を
意味しない。この版の isoSig と等長比較は向きを忘れた比較である。入力の四面体の
符号を保持することと、向きを保つ同相を認証することを区別する。

実装依存の確認先は [SnapPy Triangulation](https://snappy.computop.org/triangulation.html)
と [Manifold](https://snappy.computop.org/manifold.html)、および3.3.2の
`snappy.snap.t3mlite.Mcomplex.snappy_triangulation` と `Mcomplex.isosig`。

## 全入力カタログ

既存の `threads/PJ260922/data/261107-*/**/*.hc` と
`threads/PJ260922/code/lc-fixtures/**/*.hc` の全件、および transport の
`reason=digon_isotopy` に当たる新規入力を区別して記録する。後者は LC と HC の
本文を保存し、HC の参照ファイル名だけを作業用 `input.lc` に付け替える。
role、disk、point の不正な指定も修正しない。入力は本文の SHA256 で束縛する。

各入力を読み、既存の領域計算が列挙した全 joining disk について一行出す。
読込失敗・領域計算失敗・候補なしも一行残す。各群の計算と比較には欠損値を許し、
未実行を一致数に含めない。同じ入力の disk の全組について別表を作る。
不一致は `mismatches.json` に保存する。現象を消すための語・向き・期待値の調整はしない。

`hc_triangulation_catalog.jl` は独立実行の検証ドライバであり、通常のモジュールへは
includeしない。既存の `rows.jsonl` を上書きしない。直接ドライバ・提案差分読込は
通常のパッケージ実行と別の `mode` で記録する。
各実行の入力件数・本文hash・計算結果は出力先の `inputs.json`、`summary.json`、
`rows.jsonl` を実行終了時に読む。最新の実行先はカタログ根の `README.md` から辿る。

## 逆方向の実装範囲と未解決

`hc_snappy_catalog` は閉じたカタログ入力の `filled_triangulation()` と t3mlite の
既存 `save(format="spine")` を利用し、有限三角形分割と分岐を持たない双対spineの
語を返す。Python出口でのbranchingとo-graphは未計算とし、次のJulia入口へ渡す。
頂点が k 個なら、双対spineの対象は元の閉多様体から k 個の開球を除いたもの。
`removed_vertex_balls` に k を残し、複数の球面境界を一つの flow 条件へ読み替えない。

`hc_find_branching(data; orders=nothing, max_nodes=100_000)` は各四面体の24頂点順序を
候補とし、貼合せ面の残る3頂点の順序が一致する二項制約を解く。自己貼合せも同じ
変数への制約として扱う。候補削減後に最小候補数の変数から探索し、得た順序を再検査する。
`orders` を渡すと指定順序だけを検査する。探索尽くしの `not_branchable_on_input` と
上限到達の `branching_search_budget` を区別する。前者は同じ三角形分割上の結果である。
面写像による向きの整合・連結性を検査し、符号が入力にないときは一つの整合する向きを
選ぶ。選択した符号と由来を返す。この関数単独では頂点linkを認証しない。

`hc_triangulation_to_odata(data; orders=..., max_nodes=...)` は一頂点・χ=0を要求する。
選んだbranchingで頂点を順位に付け替え、四面体の符号を順列の符号で輸送する。
`_CP_GKS_FACE_OMITTED_VERTEX` を読み、incoming/outgoingと正負のpassageを復元する。
この既存表を複写せず、passageの後続置換の軌道からGauss語を作る。
既存の `disk_tetrahedra` がC1–C3を検査し、全4n面について元の順位付き貼合せと
再復号した貼合せが完全一致すること、四面体の向きと辺/sector数が一致することを確認する。
これが向きとbranchingを保存する明示的な組合せ同型の検査である。
返却値は `value::OData`、`virtual_link`、branchingの選択、四面体・面port・
局所順位辺からsheet germ/diskへの対応を持つ。大域辺の同値類と巡回linkは全ての
面写像の一致から同定される。

探索は一つのfull branchingを選ぶ。その選択がnormal条件を満たさないときは
`not_normal_ograph` を返す。他のbranchingでも不可能だとは主張しない。
多頂点入力は `multiple_finite_vertices`、貼合せの往復不整合は
`correspondence_failure` で止める。

`hc_each_branching(data; max_nodes=100_000, max_solutions=1024)` は同じCSPを再利用して
全解を列挙する。`solutions`、`complete`、`stop`、使用した上限を返す。上限到達時も
検査済みの解を保持する。`complete=false` は非分岐性の証明ではない。
`hc_find_branching` の `accept` callback は受理した解で止まり、全解を棄却した場合を
`branching_filter_exhausted`、解数上限を `branching_solution_budget` で区別する。

`hc_flow_certificate(data; orders=..., max_nodes=...)` は既存の `odata_to_spine` と
局所表を再利用する。球面境界が一つ、simpleな縫合円が一つ、B/Wの連結成分が各一つで
各χ=1であることを確認する。B/Wの共通辺を既存 `boundary_lift_map` で端点へ戻し、
縫合円の辺と完全一致すること、各頂点の次数が2であることを調べる。
この検査から各色の境界成分数を導く。decoderの `boundary_count` は別欄に保存し、
元の値を上書きしない。spineのχ=1、Koda曲線数=種数、joining diskの既存検査も要求する。

Kodaの構成は `koda_hc_from_boundary`、`koda_heegaard_from_hc`、
`_koda_select_joining_region`、`make_punctured` を使用する。最後の関数には既存の
`source_identity=:labelled` を明示し、標準化ハッシュの全列挙を行わない。
その図式を公開の `hc_to_odata` に渡し、向き付き語とvertex signの往復一致を検査する。
回転始点・成分順・vertex labelの付替えを許し、語の向きは保持する。
境界・縫合・B/W・往復の不一致は名前付きの `HCTriangulationError` となる。
この証明書は既存の局所表と境界interpreterの意味付けに依存する。

`hc_triangulation_to_heegaard_koda(data; ...)` はこの証明書と同じ構成を用い、
従来どおり `value` にpunctured cellular diagram、`conversion` にKoda結果型、
`encoding` に復元対応、`flow_witness` に境界検査を返す。

branchingは辺の大域的な向きと各四面体の頂点順序を要求する。pre-branchingと
weak branchingで代用しない。定義とbranchable subdivisionの存在は
[Benedetti–Petronio, AGT 14 (2014), pp.1006,1008–1009](https://msp.org/agt/2014/14-2/agt-v14-n2-p14-p.pdf)、
flow-spineの単一円環条件とB-W分解の縫合円は
[Koda, Toulouse 21 (2012), §1.3, §1.6](https://www.numdam.org/item/10.5802/afst.1349.pdf)を参照する。
細分は頂点を増やし得る。重心細分と一頂点化は、この版の実装範囲に含めない。

### 有界の移動探索

`hc_pachner_neighbors(data; runner, timeout=120)` はSnapPy 3.3.2の既存
`Mcomplex.two_to_three(unsafe_mode=false)` だけを使う。move後に `rebuild()` を行い、
全貼合せ・頂点link・χ・頂点数を検査する。合法な `three_to_two` を実行して元の
isoSigへ戻る候補だけを返す。move位置、前後のisoSig、残る四面体の番号対応、
新しい四面体、逆moveの辺を保存する。比較と重複除去の範囲は向きを忘れた組合せ同型である。
sourceからのPL写像全体を返すAPIではない。Python CLIでは `--pachner-neighbors` を使う。

`hc_prepare_inverse(data; max_depth=2, max_states=64, max_extra_tetrahedra=2,
max_nodes=100_000, max_solutions=1024, runner, timeout=120)` は、入力番号順の面を使う
幅優先探索を行う。各状態でbranchingを列挙し、flow条件を通る解を探す。
同じ入力と同じ版のbackendに対して順序を固定する。既定値の変更は新しい実行として
記録する。上限は状態数・深さ・四面体数・branchingのノード数と解数で固定し、
backendの120秒期限は実行上の停止として区別する。

返り値は `status`、`stop`、`limits`、`attempts`、`moves`、到達した三角形分割、
branching、o-data、対応、flow証明書、Koda境界コード、曲面切開を持つ辞書である。
現行版で `hc` は未構成を示す `nothing`。探索上限を多様体の不可能性に読み替えない。
多頂点では `multiple_finite_vertices_unreduced`、状態上限では `pachner_state_budget`、
branching上限では `branching_budget`、指定範囲の探索尽くしでは
`no_flow_branching_within_bounds` を返す。内部の対応不整合は
`inverse_certificate_failure` と詳細を残す。任意の辺の潰し、乱択、別入力への切替えは行わない。

### 曲面の切開と残る接続

`hc_surface_cut_system(cellular_diagram)` は既存のsigned ribbon graphから、
番号順にprimal spanning treeと、その補集合上のdual spanning treeを取る。
dual treeの辺を除いてfaceを一枚のdiskへまとめ、残余辺の数が2g、辺が正負一回ずつ、
primal treeを縮約した多角形の頂点が一つであることを検査する。元のcurve/arcへの
対応と、縮約前のdisk boundary dartsを返す。

この切開記録はPL同定の入力である。`hc_prepare_inverse` はここで `status=koda_only`、
`stop=polygon_normalization_and_reference_overlay_required` を返す。
punctured cellular diagramを解決済みHCとして保存しない。連結和のHCも、
二つの図式を曲面上で貼り、genus 2の基準網を運んで語を再抽出する工程を必要とする。
語の単純連結や、既知の不変量に合う語の作成はこの接続の代わりにならない。

有効な `OData` を得た後は既存の `odata_to_heegaard_koda` に接続できる。
それが返す punctured cellular diagram と、固定された Lickorish 基準網を持つ
公開の HC は別段階である。一般形式の種数1と標準多角形の種数2は次の入口で接続する。
αをKへ束縛する経路は、後節の `hc_triangulation_to_standard` を使う。
既存 `odata_to_hc` の `mark_heegaard` の実装範囲は変更しない。

### 種数1の基準網輸送と語の抽出

`hc_cut_to_resolved(diagram, cut=hc_surface_cut_system(diagram); joining_region)` は
Kodaの中間図式と、その図式から計算した切開記録を受け取り、`lc`、解決済みの `hc`、
`joining_disk_name="joining"`、`joining_faces`、`joining_region_id`、`certificate` を返す。
`hc_triangulation_to_resolved(input; orders, max_nodes=100_000)` は
`hc_flow_certificate` の後に同じ処理を行う。指定したbranchingと既存selectorの
joining regionを使い、失敗時に別のdiskを探索しない。

定義域は種数1の連結なcellular ribbon graphである。与えたcutがその図式からの
再計算と一致し、縮約語が `[a,b,-a,-b]`、指定した領域がjoining diskであることを要求する。
標準形の種数2は次節の範囲でも受け取る。種数3以上では
`general_polygon_slides_pending` で停止し、交換子ブロックでない語は
`polygon_not_commutator` で停止する。edge slide、各歩でのarcの輸送、連結和の
曲面構成はこの入口の実装範囲に含めない。

縮約前の境界dartを保持し、各retained edgeの二つの側と頂点sectorを展開する。
到着portから次のretained portまで回転のpredecessorを辿り、cotreeのhalf-dartを
順序付きで挿入する。展開した境界の点を有理凸多角形 `(i,i²)` に配置し、cotreeの
両端を直線で結ぶ。cotreeの端点が交互配置にならないことを確認する。
`l1` はa側から−a側へ、`m1` はb側から−b側へ向かう横断弧として置く。
これは境界語の巡回変更 `[b,-a,-b,a]=[l1,m1,-l1,-m1]` に対応し、
元の曲線と共有する辺を横断する位置へ基準網を分離する。

辺の端点パラメータは `u_l,u_m ∈ {1/(M+1),…,M/(M+1)}`、
`M=8(n_cotree+2)²+17` の有限格子で選ぶ。分子の和、次にl側の分子の昇順で試す。
判定は有理数で行い、K同士の一回の横断、portからの分離、交点パラメータの重複と
三重点の除外を要求する。上限到達は `overlay_grid_exhausted` とする。
入力名、ホモロジー、群、不変量、SnapPyの結果は候補選択に使わない。

有向接ベクトルC,Dの交点ではCのpassageを `−sign(det(C,D))`、Dをその逆とする。
旧交点ではαを `−epsilon`、βを `+epsilon` とし、旧頂点の回転順序がLCの回転と
巡回一致することを全頂点で確認する。各旧arcを交点順に細分し、Kとの交点を除けば
元の曲線語がそのまま戻る対応を保持する。Kの各segmentの両側のLC faceを合併し、
旧crossingの全cornerを用いて、合併した領域と元のribbon faceとの全単射を検査する。
LCのout cornerは元のoutgoing dartの反対向き、in cornerは元のincoming edgeの
順向きdartに対応する。指定領域に対応するLC faceの集合をHCの `disk joining` に書く。

証明書は境界lane、旧頂点からLC label、旧arcの細分、頂点回転検査、旧曲線語、
旧領域からLC face集合、corner対応、指定joining region、有限格子の上限・試行数・
選択パラメータを含む。各段階の不整合は名前付き `HCTriangulationError` で停止する。
保存は既存 `write_lc` と `write_hc(...; lc_filename=...)` を順に使い、
再読込後の順方向変換は `hc_triangulation(hc; joining_disk=result.joining_region_id)` を使う。
この番号は、輸送したface集合と完全一致する `forget_reference` の領域から取得する。
既存のnamed-disk解決に依存せず、同じ領域を番号で明示する。
同定の主張には、出どころとこの順方向出力を独立にSnapPyへ渡した比較を別途要求する。

### 標準形の種数2の基準網

同じ `hc_cut_to_resolved` は、縮約語が `[a1,b1,-a1,-b1,a2,b2,-a2,-b2]` の
種数2も受け取る。8辺の絶対値は4種類で、正負が一度ずつ現れることをcutの検査で
確認する。多角形の標準化は行わず、与えられた辺と元のjoining regionを保持する。

種数1と同じ縮約前の境界lane・有理凸多角形・cotree chordを使う。端点の順位から
別の座標を置き直さず、元の境界点とretained side上のパラメータをそのまま使う。
辺の対合は `(s,u) ↔ (-s,1-u)`。基準網の各有向arcを次のとおりに置く。

| 曲線 | 始点 | 終点 |
|---|---|---|
| l1 | (a1,u_l1) | (-a1,1-u_l1) |
| l2 | (-a2,1-u_l2) | (a2,u_l2) |
| m1 | (b1,u_m1) | (-b1,1-u_m1) |
| m2 | (-b2,1-u_m2) | (b2,u_m2) |
| n1の第1弧 | (b1,u_n1) | (-b2,1-u_n2) |
| n1の第2弧 | (b2,u_n2) | (-b1,1-u_n1) |

n1の2本の弧は端点対合で1本の円周となる。第2ハンドルのl2とm2をともに逆向きに
置くことで、Kだけに制限した語が `l1=[-1,-3]`、`l2=[-2,4]`、`m1=[1]`、
`m2=[2]`、`n1=[3,-4]` となる。Kの閉曲線の基点は巡回回転できる。
同じb側を使うmとnは、`u_ni<u_mi` の端点順で分離する。

`M=8(n_cotree+2)²+17` とし、`r=1,…,M` の昇順に `x=1/(r+1)` を試す。
各候補のパラメータを次で固定する。

```text
u_l1=(8+3x)/20    u_l2=(9+4x)/20
u_m1=(14+5x)/20   u_m2=(15+3x)/20
u_n1=(2+x)/20     u_n2=(3+2x)/20
```

Kの自己交差と余分なK同士の交差を除外し、4個のK交点と符号、端点対合、元のarcとの
横断性、共有支持、交点順の重複・三重点を有理数で検査する。共有支持を交点0として
処理せず、その候補を棄却する。上限到達は `overlay_grid_exhausted`。
ホモロジーやSnapPyの結果によるパラメータの調整は行わない。

元の交点のlabelはK交点4個の後に置き、元の各有向arc上の交点とKの各弧上の交点を
それぞれパラメータ順に読む。passage signと全cornerによる領域対応は種数1と同じ
処理を使う。証明書は共通の欄に加え、`reference_arcs` に各弧の符号付き辺と有理端点、
`even_handle_orientation_reversed`、`overlay_grid_rule` を持つ。

返り値の `joining_region_id` は、元のdiskを細分したface集合と一致する領域の番号である。
後続の既存 `forget_reference`・Koda入口の検査が失敗した場合は、その例外をそのまま
停止として扱う。別の領域へ切り替えない。HCの生成、HCの公開順方向入口の通過、
独立な三角形分割の同定は、それぞれの実行結果を区別して記録する。

SnapPyの `FundamentalGroup` は表示を返すが、3.3.2に群位数の公開APIはない。
H₁や有限商への準同型数から元の群の位数を代入しない。位数は未確認と記録し、
位数を決められないことだけを理由に逆方向の処理を止めない。

例の枠・選択条件は親の `.research-scratch/PJ260922/astra-lanes/example-table.md` を読む。
各実行後の表と停止は `threads/PJ260922/data/261112-example-table/table.csv` と
`stops.csv`、助言の採否・担当外の修正案は実行ターンのレーン報告 `REPORT-<turn>.md` を読む。

## αを基準曲線に束縛する逆方向

`hc_triangulation_to_standard(input; alpha_basis, orders=nothing, max_nodes=100_000)` は、
既存の `hc_flow_certificate` に続けて、同じ中間曲面上へKを置き、LCと解決済みHCを返す。
中間図式を直接渡す入口は
`hc_koda_to_standard(diagram; joining_region, alpha_basis)` である。
`alpha_basis` は必須で、`:l` と `:m` を受け付ける。あかげふの既定値を代入する処理は置かない。

定義域は、正の種数の連結なcellular曲面、完全な最小α系とβ系、指定済みのjoining disk。
完全性は `compression_boundary_report` による切開面の連結性・種数・境界本数から検査する。
三角形分割からの入口には、既存のbranching、有限頂点、closed、flow、Koda往復の条件も適用する。
この段でPachner探索の上限やjoining regionを変更する処理は行わない。

構成では元のα・βの符号と巡回語を保持する。αで切った2g穴の球面上で、順に双対曲線を
置く。各双対は指定したα辺の両側を、β辺だけを横切る単純な面双対経路で結ぶ。
すでに置いた基準曲線を避け、面の横断後は既存の `_lc_tau` に従う。
`:m` 指定では、この段階の二族を交換し、新しいlの向きと、その交点の両passageを反転する。
元のα・βの向きと順序は保たれる。

各 `n_i` は、m_iとm_(i+1)の指定したbankを結ぶ帯の正則近傍境界として作る。
既存の経路検査、barbellの境界軌道 `[1,3,-2,-3]`、境界語コンパイラ、共通の編集kernelを使う。
隣り合う帯は共有m上で反対側のbankを使う。面双対経路は各面を一度訪れる幅優先探索で求める。
交点の符号・辺上の挿入順は既存のdartとbankの規則から得る。
この内部の構成用編集記録は、登録済みのLC moveの実行証明書とは別のものである。

Kだけに制限した語を `_lc_k_constants` の定義と照合し、元の全曲線の巡回語を復元して検査する。
役割を持たないKを忘れた全領域が円板であること、元の全crossing cornerが対応すること、
元の全領域との全単射を検査する。指定joining diskの像は、その領域の全LC子faceをまとめた
`disk joining` とする。LCの子faceを一つ選び直す処理は行わない。

返り値は `lc`、`hc`、`joining_region_id`、`joining_faces`、`certificate`。
証明書はαの切開面、各双対の面経路、帯とbank、旧交点label対応、元の曲線の保存検査、
全cornerと全領域の対応、指定diskの輸送を保持する。
保存には `write_lc` と `write_hc(...; lc_filename=...)` を使う。
再読込したHCの順変換は `hc_triangulation(hc; disk_name="joining", signature=false)` で要求できる。
同定には元の三角形分割と、この公開順変換の像を独立にSnapPyへ渡す。
語の不変量は、既存の `hc_presentation`、`hc_first_homology`、`hc_hom_count` で別途要求する。

不正なbasis、非cellular、種数0、αまたはβの不完全性は、それぞれ
`standard_alpha_basis`、`standard_source_not_cellular`、`standard_positive_genus`、
`standard_alpha_not_complete`、`standard_beta_not_complete` で停止する。
経路・基準語の構成の不整合は `standard_reference_construction`、
元の曲線・領域・diskの対応の不整合は `standard_source_transport` を伴う
`HCTriangulationError`、既存のLC・Koda検査の失敗はその名前付き例外で返す。
仮定を満たす入力で経路が見つからない場合は、面またはbankの実装を調べる対象とする。

構成の接続性とbank順の監査は
`.research-scratch/PJ260922/astra-lanes/snappy-bridge/answer-9-2.md`、
有限入力での公開入口・独立同定・語の計算の実行状況は同フォルダの `REPORT-9.md` と
`threads/PJ260922/data/261112-example-table/` の当該実行記録を読む。
# HC の曲面操作

`hc_connected_sum(left, right; left_gap=nothing, right_gap=nothing)` は種数1のHC二つの曲面を連結和し、元のα・βを各成分の埋め込みに沿って運ぶ。各入力のα・βはそれぞれ一本で、全曲線は非空とする。点・係数・LC装飾と空語の曲線は現在の定義域に含まない。役割は入力から継承し、ℓ・mの既定値を選ばない。

左右のm1の平行な曲線を既存の帯境界コンパイラで作り、指定したgapで切った二本を向きを反転する首で接続してn1とする。gapの省略時はm1の語の最後の後ろを使う。固定Kの射影、曲面の種数、全入力曲線の符号付き語の埋め込みを検査する。返り値は `(lc, hc, certificate)`。証明書は入力の曲線・交点の写像、穿孔したfaceとgap、首の端点の対合を含む。入力のdisk印は消費して記録し、出力のHCにはdiskを付けない。α・βの補領域に生じる円環は円環として保持する。`certificate["status"] == "joining-disk-pending"` は、順方向への出口が完成したという意味ではない。

`hc_destabilize(hc; alpha, beta, alpha_basis, joining_disk=nothing, disk_name=nothing, disk_id=nothing)` は、指定したα・βが幾何的に一回交わるとき、そのハンドルを除く。α・βの番号はKを忘れた図式の順序である。入力と出力は完全な切断系を持つcellularな図式とし、出力の種数は正とする。点・係数・追加の印の輸送は現在の定義域に含まない。`alpha_basis` は `:l` または `:m` を必須とする。

交点pで二曲線を切り、残る各曲線を支持の帯に沿ってslideしてから、一穴トーラスを切除して円板で塞ぐ。pを除いたα・βの巡回列をA・Bとし、q∈Bとr∈Aが作る新交点の符号は `-epsilon(q)*epsilon(r)/epsilon(p)`。残るαでqを置換するAの向きは `-epsilon(q)/epsilon(p)`、残るβでrを置換するBの向きは `-epsilon(r)/epsilon(p)` とする。種数と両切断系を再検査し、固定Kは既存のstandard入口で作り直す。

指定したjoining diskの全ての生き残ったcornerが一つの新faceへ写ることを要求し、そのfaceでjoining diskの条件を検査する。cornerが残らない場合、写像が複数になる場合、輸送先がjoining diskでない場合は停止する。別のdiskや別のキャンセル対への選び替えは行わない。返り値はstandard入口と同じ形式で、`certificate["destabilization"]` に交点の置換表・曲線の出所・種数・diskのcorner対応を加える。

これらの操作の定義域外と支持の不成立は、理由を持つ既存の `LCMoveError` で返す。destabilizationは指定された支持の一回の操作であり、最小種数の探索は行わない。三角形分割と独立な出どころとの同定は、従来どおりSnapPyまたはisoSigの一致を用いる。H₁・有限群への準同型数の一致だけで同定済みとは扱わない。
