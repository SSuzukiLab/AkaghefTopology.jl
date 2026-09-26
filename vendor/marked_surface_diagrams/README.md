# 固定された曲面上のHeegaard図式：経路保存・描画カーネル

**この束の対象は、曲線それぞれの自由ホモトピー類ではなく、固定された向き付き曲面上の、交点付き図式全体の埋込みである。**

保存するのは「固定した曲面の円板分解・辺の貼合せ」「各円板の向き付き細分」「その細分の辺をたどる曲線語」「交点ラベル・符号」「選択した領域」。描画器は既存の組合せ埋込みを実現するだけで、経路探索をやり直さない。

## 最初に読むもの

| ファイル | 内容 |
|---|---|
| `overview.html` | 図と実行結果付きの概要。ブラウザで単独表示可能 |
| `docs/SPECIFICATION_JA.md` | 意味論、データ型、記号列、API、描画、制約 |
| `docs/PROOFS_JA.md` | 表現の一意性・表現可能性・描画保存の証明と仮定 |
| `docs/CODEX_HANDOFF.md` | 既存Juliaアプリへの具体的な接続順序と受入条件 |
| `docs/MIGRATION_JA.md` | 旧core/layout・joining disk・marking移行 |
| `docs/VERIFICATION_JA.md` | 実行済み／未実行の区別、試験範囲と限界 |
| `src/MarkedSurfaceDiagrams.jl` | Juliaネイティブ実装。標準ライブラリのみ |
| `src/builders.jl` | 整数格子上の記号的経路の入力補助 |
| `tools/reference.py` | 実行済みPython参照実装。標準ライブラリのみ |
| `test/expected.toml` | 実行済み参照実装の厳密な有理座標・符号・位相の期待値 |

## 実行

Julia 1.10以降を想定している。Julia版はこの配布環境では**未実行**。

```sh
julia --startup-file=no test/runtests.jl
julia --startup-file=no tools/run.jl examples/torus_three_crossings.toml --svg artifacts/julia-three.svg --report artifacts/julia-three.json
julia --startup-file=no examples/author_annulus.jl
```

実行確認済みの参照経路はPython 3.11以降で動く。

```sh
python tools/verify.py
python tools/reference.py examples/annulus_winding_p1.json --svg artifacts/annulus.svg --report artifacts/annulus.json
```

Pythonは既存Juliaアプリの実行依存にしない。移植の照合オラクルと、Juliaの接続前に仕様を検証するためのものである。

## 保証の読み方

有効な一つの文書からは、一つの **marked surface 上の同時ambient isotopy class** が定まる。同じ文書を別の適合描画で実現してもこの類は変わらない。

ただし、同じ同位類が一つの文書に正規化されるわけではない。`representation_digest` は文書識別子であり、同位類の完全不変量ではない。独立した曲線の同位、bigon除去、handleslide、Dehn twist、stabilizationは、この保証と別の操作である。

## 現在の検証状態

Python参照実装の36試験、および12例の厳密描画・小数化後の再検査を実行済み。Juliaの同一期待値照合、およびユーザー環境の既存アプリへの据付けは未実行。`artifacts/verification-status.json` に機械可読な状態を保存している。

既存の未コミット作業ツリーを読めていないため、この束には既存リポジトリのファイルを書き換えるパッチを含めない。旧coreから失われたmarkingを復元した、あるいは既存Poincaré例と同じ固定曲面上の経路を取得した、とは主張しない。
