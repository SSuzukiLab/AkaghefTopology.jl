# Topology のモジュール構成

Status: あかげふが 2026-11-14 に推奨案を承認した構成。移行中。<br>
Scope: `Topology` パッケージのモジュール分割、型名、総称関数の持ち主、依存の向き<br>
Readers: 人間は §1 の図と表だけを読めば足りる。AI agent は最後まで読む。

## 1. 構成

```
Topology                 パッケージ。再エクスポートだけを持つ
│
├─ VirtualLinks          型 VirtualLink（o-graph とリンク）
├─ KnotGauss             型 KnotGaussDiagram（古典結び目の Gauss 語）
├─ URDiagrams            型 URDiagram
├─ HeegaardCodes         型 HeegaardCode（HC）。LickorishCodes の上のレイヤー
│     └ 取り込む ─ LickorishCodes   型 LickorishCode（LC）。曲面・K・曲線の語
│
├─ MoveAPI               総称名と契約だけ：movable, move, 探索 DB の操作
├─ Invariants            総称名と契約だけ：表示によらない不変量の名前
├─ MoveGraphStore        移動グラフ DB の器（既存）
└─ Algebra               表示によらない代数：有限代数、Hopf 代数、テンソルネットワーク、Laurent 多項式
```

依存は下から上へ向かう。`Algebra`・`MoveGraphStore`・`MoveAPI`・`Invariants` は表示のモジュールを知らない。表示のモジュールは、それらを取り込み、自分の型のメソッドを足す。

| モジュール | 主な型 | 担うもの |
|---|---|---|
| `Algebra` | `FiniteAlgebra`、`FiniteHopfAlgebra`、`SparseFactor` | 表示によらない代数とテンソル縮約 |
| `MoveAPI` | なし | `movable`、`move`、探索 DB の総称名と契約 |
| `Invariants` | なし | 不変量の総称名と契約 |
| `MoveGraphStore` | `MoveGraphStoreDB` | 状態・遷移・探索ランの永続化 |
| `VirtualLinks` | `VirtualLink` | o-graph とリンク、移動、重み、MST、o-graph 側の Kuperberg |
| `KnotGauss` | `KnotGaussDiagram` | 古典結び目、Reidemeister 移動、分類 |
| `URDiagrams` | `URDiagram` | U(R) 図式と MST の trace 規約 |
| `LickorishCodes` | `LickorishCode` | 曲面と K と曲線の語、実現、描画、LC の移動 |
| `HeegaardCodes` | `HeegaardCode`（公開別名 `HeegaardDiagram`） | 役割・円板・点、検査、変換、不変量、HC の移動 |

### 1.1 Root exports

bug: `Topology` exports `neighbors` at [src/Topology.jl:104](../../src/Topology.jl#L104),
but the root module cannot resolve the binding: its imported module `Graphs`
and its submodule `MoveGraphStore` export different functions named `neighbors`
(`using Graphs` at [line 6](../../src/Topology.jl#L6), `using .MoveGraphStore` at
[line 165](../../src/Topology.jl#L165); the store export is at
[src/move_graph_store.jl:17](../../src/move_graph_store.jl#L17)). After
`using Topology`, the observed contract-example error is
``UndefVarError: `neighbors` not defined in `Main.APIExample` ``, with Julia's hint
that two or more modules export different bindings, resulting in ambiguity,
and that a binding also exists in `Graphs`. The generated index records
`Topology.neighbors: exported binding is undefined.` The root export needs an
explicit choice of binding; the store query is accessible as
`Topology.MoveGraphStore.neighbors`.

## 2. 名前の規則

1. モジュール名は、それが持つ主な型の複数形にする。Julia では、同じ親からモジュールと同名の型を両方見せると衝突するためである（`DataFrames` が `DataFrame` を持つのと同じ形）。
2. 型名は仕様の正式名に揃える。LC は `LickorishCode`、HC は `HeegaardCode` である（[用語集](heegaard-terminology.md) §2）。`HeegaardDiagram` は `HeegaardCode` の公開別名として残す。
3. 旧名は移行期間だけ `const` の別名として残す。旧名は `LickorishConfiguration`、`HeegaardLayer`、モジュール `Heegaard` である。新しいコードは旧名を使わない。
4. HC 導入前の位置を持たない型（`types.jl` の旧 `HeegaardCode`）は名前を空け、`CombinatorialHeegaardCode` と呼ぶ。2026-11-09 のあかげふの決定により、`HeegaardCode` の名前は HC が上書きする。

## 3. 総称関数の持ち主

総称関数は、一つのモジュールだけが宣言する。他のモジュールは `import ..MoveAPI: movable` のように取り込み、自分の型のメソッドを足す。こうすると全表示で一つの関数になり、引数の型で振り分けられる。Julia の振り分けは型で行われ、モジュールでは行われない。別々のモジュールが同名の関数を `function` で定義すると、それは別の関数になる。

### 3.1 MoveAPI

| 名前 | 契約 |
|---|---|
| `movable(x)` | 現在の一引数メソッドは LC/HC の既定の有限候補集合を返す。逆安定化はこの集合に含めない。[契約](../api/heegaard-code/Topology.MoveAPI.movable.md)を参照 |
| `movable(x, kind; 選択...)` | 候補の範囲と省略時の意味は表示別の[契約](../api/heegaard-code/Topology.MoveAPI.movable.md)による。KnotGauss の空の選択はその種類の全局所候補。VirtualLink は局所頂点の選択を必要とし、省略は全頂点部分集合の列挙ではない。HC の行指定はその行の候補に限る |
| `move(x, witness)` | witness を適用した新しい値を返す。`x` は変えない |
| `move!(x, ...)` | 破壊的な版。可変な表示（`VirtualLink`）だけが持つ |
| `search_state(x)` | 探索 DB の状態の鍵。同じ表示の値は同じ鍵になる |
| `put_search_state!`、`load_search_state` | 状態の保存と復元 |
| `begin_search_run!`、`finish_search_run!` | 探索ランの記録 |
| `expand_search_state!(store, x; kinds)` | `movable` の全候補を適用し、遷移として保存する |
| `search_neighbors`、`query_search_paths` | 保存した遷移の読み出し |
| `verify_search_store` | 現在の o-graph 版はスキーマ版・外部キーを検査し、保存行を数える。遷移の再生は行わない。[契約](../api/move-graph/verify_search_store.md)を参照 |

既存の `ograph_search_*` は `VirtualLink` のメソッドの実装として残す。

探索 DB の関数のうち、`load_search_state`、`begin_search_run!`、`finish_search_run!`、`search_neighbors`、`query_search_paths` は、o-graph の版と `hc_search_*` の版で位置引数が同じ（`store` と状態の鍵など）である。Julia はキーワード引数で振り分けないので、このままでは一つの総称関数のメソッドにできない。表示の型を `store` の次の位置引数として渡す形にする。例えば `load_search_state(store, HeegaardCode, state_id)`、`search_neighbors(store, VirtualLink, state_id)` とする。これは `parse(T, s)` と同じ書き方である。o-graph の版の今の呼び方は残す。（2026-11-16、統括のセッションが決めた。あかげふが変えてよい。）

### 3.2 Invariants

表示を引数に取る名前は、このモジュールにメソッドを持たない。各名前の docstring が次の契約を書く。Casson–Walker の式は手術係数と Seifert データから計算する関数で、`src/casson_walker.jl` にある。HC は読まない。規約は [Casson–Walker の規約](../conventions/CassonWalker.md) に draft で置いた。

- 入力は閉 3 次元多様体を表す値である。
- 出力の型。
- `MoveAPI` の、多様体を保つと宣言された移動の前後で値が等しい。

総称名にするのは、二つ以上の表示で同じ多様体の量として意味を持つものだけである。現在の候補は `first_homology`、`linking_form`、`Z_Kuperberg` である。`Z_MST` は Heegaard 側に無いので、`VirtualLinks` だけの関数として残す。二つ目の表示に実装ができた時点で総称名に昇格させる。

上の不変性の契約は、移動の実装を確かめる探索の判定にそのまま使う。

## 4. 移行の順序

1. 型の改名と旧名の別名。`LickorishCode`、`HeegaardCode`、`CombinatorialHeegaardCode`。
2. `Heegaard` を `LickorishCodes` と `HeegaardCodes` に分ける。`Heegaard` は旧名の別名として残す。
3. 最上位のファイルを `Algebra`、`VirtualLinks`、`KnotGauss`、`URDiagrams` に分ける。
4. `MoveAPI` と `Invariants` を作り、既存のメソッドをそこへ付け替える。

2026-11-14 の状態：第 1 段は済んだ。第 2 段のうち、モジュール `Heegaard` を `HeegaardCodes` に改名する部分は済んだ（`Topology.Heegaard` は別名）。`LickorishCodes` の切り出しは保留している。`lc_moves*.jl` が `Union{LickorishCode,HeegaardCode}` を受け、HC のレイヤーの処理を含むため、LC のレイヤーから HC のレイヤーへの依存を先に解く必要がある。

2026-11-16 の計画：切り出しは四段で行う。参照の行と段ごとの受入は、研究リポジトリの根から `threads/PJ260922/data/261114-module-rename/lc-heegaard-split.md` にある。段 1 は `lickorish_configuration.jl` を分け、`LickorishCode` の定義と読み書きを下に、`HeegaardCode` と `resolve_hc`、`read_hc`、`validate_hc_input` を `HeegaardCodes` に残す。段 2 は `hc_text`、`write_hc`、`realize_hc`、`draw_heegaard(::HeegaardCode)` と、役割を読む移動を `HeegaardCodes` 側のメソッドに移す。段 3 は `Union{LickorishCode,HeegaardCode}` を分け、LC の関数は `LickorishCode` だけを受ける。段 4 で HC を名指ししなくなったファイルをモジュール `LickorishCodes` に移し、`HeegaardCodes` がそれを取り込む。各段の受入は移行前と同じ試験結果で、公開名は残す。`lc_moves.jl` と `lc_moves_support.jl` は移動のレーンが、`koda_conversion.jl` の include は変換のレーンが持っているので、そのファイルを触る段はレーンが空いてから行う。

2026-11-16 の状態：第 3 段は、最上位の実装ファイルを `Algebra`、`URDiagrams`、`KnotGauss`、`VirtualLinks` に入れ終えた。`orthogonal_layout.jl` と `layout_invariants.jl` は `VirtualLink` だけを取るので `VirtualLinks` にある。`svg_plot.jl` は `URDiagram` と `VirtualLink` の両方を描くので `Topology` の直下に残してある。2026-11-16 に中身を見て、前の二ファイルは `VirtualLink` 専用、`svg_plot.jl` は共有と決めた。共有ファイルは `Topology` 直下に残し、描画モジュールへの分離は §5 のとおり保留する。`kuperberg_invariant.jl` は `HeegaardCode` から `VirtualLink` を組み立てるので、`HeegaardCodes` の後に `Topology` が取り込んでいる。`move_graph_store.jl` は §1 の `MoveGraphStore` のままである。

第 4 段は、`MoveAPI` と `Invariants` の総称名と契約がある。総称名に付いているのは、`VirtualLink` と `KnotGaussDiagram` の `movable` と `move`、`VirtualLink` の `move!`、`VirtualLink` の探索 DB（型引数のない呼び方は既存の `ograph_search_*` を呼ぶ）、`HeegaardCode` の `movable`（引数なし、行、行 id）と `move(::HeegaardCode, ::LCMoveWitness)`、`Z_Kuperberg`、`mirror`（`HeegaardCode` は `hc_mirror`、`VirtualLink` は `mst_sign_change`、`KnotGaussDiagram` は図式の鏡像）、`connected_sum`（二つの `HeegaardCode` は `hc_connected_sum`、`VirtualLink` は成分の選択を取る `connected_sum!`）である。`load_search_state`、`begin_search_run!`、`finish_search_run!`、`search_neighbors`、`query_search_paths` は、`store` の次に `HeegaardCode` または `VirtualLink` を取る。HC の版は `load_hc_search_state`、`begin_hc_search_run!`、`finish_hc_search_run!`、`hc_search_neighbors`、`query_hc_search_paths` へ、o-graph の版は既存の `ograph_search_*` へ、それぞれ一行で転送する。`search_state(::HeegaardCode)`、`put_search_state!(store, ::HeegaardCode)`、`expand_search_state!(store, run_id, ::HeegaardCode)` は `hc_search_state`、`put_hc_search_state!`、`expand_hc_search_state!` へ転送する。`first_homology(::HeegaardCode)` は `hc_first_homology`、`linking_form(::HeegaardCode)` は `hc_linking_form` へ転送する。`ograph_h1` は `OGraphPi1Presentation` を取るので `VirtualLink` のメソッドにはしていない。`VirtualLink` の絡み形式の関数はリポジトリに無い。別の関数のまま残っているのは、`LickorishCode` の行指定の `movable` と `move`、`knot_canonical_state`、`expand_reidemeister_neighbors!`、`HeegaardCodes.first_homology`（群の表示を取る）である。`verify_search_store` は o-graph のメソッド一つである。`Z_MST` は `VirtualLinks` の関数のままである。

2026-11-16 に、統括がサブモジュール `lib/topology.jl` の HEAD（`af90b62`）の `export` を `names(Topology)` と内側のモジュールの名前に照らした。消えた公開名は無かった。このターンは公開名を足さず、上のメソッドだけを足した。

各段の受入は、移行前と同じ試験結果である。移行前に失敗していた試験は、同じ失敗のまま残ってよい。

## 5. 未決

- `orthogonal_layout.jl` と `layout_invariants.jl` は `VirtualLink` 専用である。前者は `orthogonal_layout(::VirtualLink)`（`src/orthogonal_layout.jl` 294 行）と、曲げ数の補助 `minimal_bending_numbers(::AbstractMatrix)`（30 行）を持つ。後者は `check_layout(::LinkLayout)`（`src/layout_invariants.jl` 194 行）と `check_layout(::VirtualLink)`（206 行）を持つ。`LinkLayout` は `VirtualLink` の配置である。両方とも `VirtualLinks.jl` が取り込んでおり、そこに置く。
- `svg_plot.jl` は共有である。`plot_svg(::URDiagram)`（193 行）、`plot_urd_structure_svg(::URDiagram)`（209 行）、`plot_svg(::VirtualLink)`（511 行）を持つ。`Topology.jl` 193 行が直下で取り込んでいる。描画の別モジュールへ出す作業は、`Topology.jl` をこのレーンだけで編集できる段まで保留する。
- `heegaard/moves.jl` と `heegaard/types.jl` の旧 API を残すか消すか。使用箇所は `threads/PJ260922/data/261114-module-rename/lc-heegaard-split.md` の後半にある。消す場合はあかげふの確認を経る。この一覧を書いた時点では消していない。
- `MoveAPI` と `Invariants` を一つのモジュールにまとめるかどうか。
- `HeegaardCodes` の中の名前 `HeegaardDiagram` は `CellularHeegaardDiagram` を指し、`Topology.HeegaardDiagram`（`HeegaardCode`）とは別の型である。`Z_Kuperberg` のメソッドはこの二つの別々の型に付いている。内側の名前を改めるかどうか。改名はあかげふの確認を経る。
