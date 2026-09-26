# 旧Heegaard coreと配置からの移行

## 失われた情報は推測で埋めない

現在把握している旧coreは、genus、交点のα/β所属・符号、各曲線上のcyclic order、joining disk境界などである。この情報は重要であり、捨てない。しかし、固定されたmarked surface上での経路を指定する情報が別に必要な場合がある。

旧文書から何か一つを自動描画できたことは、その出力が元の固定曲面上の意図した同位類である証拠ではない。旧coreだけを新形式へ移す際に、routerが選んだ経路を「復元した」と呼ばない。

## 移行状態

| 状態 | 意味 | 認証表示 |
|---|---|---|
| `legacy_unmarked` | cyclic orders等だけがある | marked同位類は未指定 |
| `adopted_realization` | 現在の保存配置を新たな明示経路として採用した | 採用した図のクラスを指定。元の意図との一致は別 |
| `marked_verified` | 基準曲面・埋込み・幾何を検証した | 本仕様の同位類保存を表示可能 |
| `transport_verified` | 旧atlasとの明示的な同じmarkingの輸送も検証した | 該当する旧viewとの対応を主張可能 |

これらはアプリ側の状態enumであり、新数学文書の未定義キーとして勝手に追加しない。非数学metadataまたはwrapperに置く。

## 実際のsaved layoutがある場合の採用手順

本配布物は任意の旧Bezier/SVGを自動取り込むimporterを実装していない。以下はimporterを追加する際の要求仕様であり、既存データをこの環境で移行済みと扱わない。

まず旧chartのorientation、boundary gluing、portのparameter、routeの点列、annulusでの継続を固定する。カメラ投影や画面上の交差を数学上の交点と混同しない。

元が浮動小数点なら、値を確定した小数として有理数化し、その**採用したPL経路**を明示する。Bezierをサンプルした折れ線を使う場合は、原curveとの同位に関する区間保証がなければ「原Bezierと認証済み同位」とは言わない。

曲線、true boundaries、基準seamsの完全arrangementを作る。意図する交点のみをラベル付き頂点にし、各辺を分割する。真の交点以外の交差が見つかれば、黙って新交点にせず止める。

各基準円板の全領域を細分・三角形分割する。無交点loopは削除せず、必要な補助頂点を追加する。非連結な成分がどの補領域に入るかは、採用した埋込みから確定する。

変換後`validate`を通し、`legacy_projection`と旧coreのcyclic order・sign・curve labelsを照合する。これは必要条件だが、marking保存の十分条件ではない。旧atlas→新基準曲面の同じmarkingの明示対応を別に保存する。

Dは新図式の領域として選び直し、旧D境界の各辺・cornerへの対応を記録する。領域が円板でない、同じ旧idが異なる新境界sectorへ重複する、方向規約が一致しない場合は止める。

## 証跡の内容

wrapperには、元ファイルのhash、採用した時点の明示同意／プロジェクト決定、元chartと新sideの対応、portのparameter変換、旧route→新dart列、旧crossing→新vertex、旧D→新region witness、符号変換、幾何輸送の検査結果を保持する。

ただし、これらの保存だけで証明が完成するわけではない。例えば「old_route_7→new_curve_7」という名前の対応だけでなく、実際の辺と向き・局所写像が検査可能である必要がある。

## 教科書をコードで書く際の運用

新規の図は最初から基準曲面と語をソースにする。見た目の変更はlayout編集とし、数学文書を変えない。数学的なmoveを行う場合には、前後の文書と対応・局所証跡を別のrewriteとして残す。

meridianやlongitudeの基準図を固定しても、それとの交点**数**を追加しただけで図式全体の配置が指定されたとは限らない。基準線との順序付き交点と、各disk内での配置を含めて本形式へコンパイルする。Lickorishの生成元はmapping classを操作する入力言語としては有用だが、現在の教科書図式を忠実に描くrendererの内部で任意のtwistを選ばせない。
