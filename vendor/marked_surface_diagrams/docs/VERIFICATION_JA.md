# 検証状態と受け入れ試験

## 実際に実行した範囲

Python 3.13.5で、標準ライブラリだけの参照実装を実行した。`tools/verify.py`は36のunit test、12個のJSON/TOML例の構造検査、厳密な有理座標での描画、SVG小数化後の厳密再検査を行う。実際のログは`artifacts/python-tests.txt`、機械可読結果は`artifacts/verification-status.json`。

| 検証群 | 内容 |
|---|---|
| 表現 | JSON/TOMLの同内容、語のparse/print、旧coreへの符号正規化射影 |
| 基準曲面 | sphere、closed genus 1–4、annulus、disk、Euler数・link・boundary |
| 経路 | 同じ円環端点でwind −1/0/+1、同じmarking識別子、異なるseam語 |
| 図式 | 正負交点、両曲線のcyclic order、bigon、無交点平行曲線、contractibleな入れ子 |
| 領域 | 開円板と閉包が埋込み円板であることの区別、annulus領域、family cut-system |
| Joining | 幾何的DとKoda条件の区別、明示curve方向audit、監査で元文書を変えない |
| 不正入力 | 不正向き、誤genus、二重貼合せ、欠落交点、誤符号、α同士交差、辺共有、実在しない経路 |
| 幾何 | 有理数以外を拒否、反転／重複を拒否、一内部頂点の線形同位、SVG再検査 |
| 状態管理 | 入力の防御的コピー、数学文書の変更で既存Stateを失効、未知semantic keyを拒否 |
| 補助細分 | 12回のstellar refinementで交点・領域を保ち、文書hashと同位を混同しない |

12例は、円環3、mer/lonトーラス1、無交点平行トーラス1、(+,-,+)三交点トーラス1、disk内の入れ子1、closed genus1–4の空図4、sphere空図1である。

`torus_three_crossings`は新しく定義した固定トーラス上の図式であり、既存アプリのAbalone配置から抽出したものではない。既存Poincaréのルートデータは含めていない。

## 概要HTMLの検査

Chromiumで自己完結HTMLを実際に表示し、7例の切替、mesh表示切替、Dの領域塗りつぶし、JavaScriptエラーがないことを検査した。`artifacts/html-test-status.json`と画面画像を保存している。JSON Schema 2020-12による全12例の外形検査も実行した。これらは既存Heegaardアプリの試験とは別である。

## 実行していない範囲

Julia runtimeが配布環境にないため、Juliaネイティブ実装は実行していない。ネイティブ版の意味が参照版と一致することは、提供した`test/runtests.jl`をユーザー環境で通す必要がある。静的に記述したコードを「実行検証済み」とは扱わない。

既存Windows作業ツリーの型定義・未コミット差分は取得できなかった。従来アプリへのdispatch、保存再読込、UI、sphere-plus-annuli atlasへの輸送、Koda converter接続も、この配布環境では試験していない。

Lean等での形式証明は行っていない。数学的根拠は`PROOFS_JA.md`の仮定付き証明であり、実行テストはその証明の代替ではない。

## 参照実装とJulia版の独立性

同じ設計を二言語に書き分けたものであり、完全に独立した二つの理論実装ではない。共通する設計上の誤りは、相互照合だけでは見つからない。そのため、Euler数、既知の円環wind、交点順序、region closureなどの外部的に予測できる例を含めている。

Julia側はPythonの厳密な有理座標を含む期待値と直接比較する。`expected.toml`を書き換えて失敗を消すことを避け、差があればまず意味規約と移植処理を調べる。

## 表示保証の限界

`certified_centerlines=true`は、出力小数を正確な有理数として解釈したSVG中心線の埋込み保証である。ブラウザの内部丸め、画素解像度、stroke幅、文字ラベル同士の重なりは保証外。見やすさの目視検査と数学的経路認証を別軸で扱う。

初期rendererは凸disk atlasであり、既存の穴を円形に配置したviewと同じ美観ではない。経路の意味を失わないことを先に保証し、美観の変換は明示的なPL輸送で追加する。
