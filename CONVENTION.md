<!-- XQQQ1500000 -->
# Topology.jl の規約

規約はモジュールごとのページに置く。このファイルはその目次である。モジュールの一覧は[モジュール構成 §1](docs/specs/module-structure.md#1-構成)に従う。

| モジュール | 規約のページ |
|---|---|
| `LickorishCodes` | [LickorishCodes の規約](docs/conventions/LickorishCodes.md)（平面図の規約は全文、ほかは一覧） |
| `HeegaardCodes` | [HeegaardCodes の規約](docs/conventions/HeegaardCodes.md)（一覧） |
| `VirtualLinks` の Z₃ 辺色 | [OGraphColours の規約](docs/conventions/OGraphColours.md)（全行 confirmed 2026-11-22） |
| `VirtualLinks` のレンズ空間の名前と向き | [レンズ空間の名前と向きの規約](docs/conventions/LensSpaceNames.md)（fixed 2026-11-23） |
| `VirtualLinks` の頂点符号・branching・流れの向き | [o-graph の頂点・branching・流れの規約](docs/conventions/OGraphVertex.md)（draft。符号の正規化、流れ、全正の abalone の $\theta_G$ は fixed 2026-11-24） |
| `VirtualLinks` の残り、`KnotGauss`、`URDiagrams`、`Algebra`、`MoveAPI`、`MoveGraphStore` | 未作成。規約の棚卸しはまだしていない |
| `Invariants` の Casson–Walker | [Casson–Walker の規約](docs/conventions/CassonWalker.md)（draft） |
| 手術後の図式でのスライドの逆 | 下の「手術後の図式での逆」（draft） |

## 手術後の図式での逆

状態は draft である。あかげふの受理はまだ無い。長い文と比較は `docs/specs/bordered-lickorish-code.md` の「案：手術後の図式での逆」にある。関数は `apply_geometric_slide_inverse` で、`apply_reversed_slide` は手術前の図式で帯を認定する別の関数のままである。

| 規約 | 内容 | 状態 |
|---|---|---|
| 往路の端点 | 始点は \(a_{\mathrm{new}}\) の補弧の先頭の in 半辺。終点は前向きの終点と同じラベルと同じ符号で、in と out を入れ替えた半辺 | draft |
| 弧の向き | 往路と復路の入れ替えは各弧の向きの反転を含む。新しい往路は \(a_{\mathrm{new}}\) から \(b\) へ向かう。平行弧の反転は §7.7.3 にある | draft |
| 矩形 | 手術後の閉じた図式の一つの面の上の道として認定する。手術が面を分けたときは手術後の面で読む | draft |
| 付着辺 | 補弧の長さを \(n\) とすると、長くなった語の番号は \(n+1\)。回し戻す量は前向きの始点の辺番号で、witness に残す | draft |
| 対応の記録 | 逆の band の宿主半辺、付着辺の番号、振ったラベル、置き換えた曲線と保つ曲線を、返す記録に残す | draft |

別の読みで結果が変わる点は、仕様の案の小節に書いた。前向きの始点を手術後も往路の始点にすると、別の面に分かれたときに矩形の認定で止まり、語は \(a_{\mathrm{new}}\) のまま残る。終点の in と out を保つと、二角形を含む例では同じ理由で止まる。端点だけを入れ替えて向きを保つと、往路が \(b\) から始まる。付着辺を 1 番と読むと、長さ 1 の \(a\) では元のラベルが残らない。対応を残さず番号を選び直すと、消すラベルが変わり得る。
<!-- XQQQ1000000 -->
