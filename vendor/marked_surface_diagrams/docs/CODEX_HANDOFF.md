# Codexへの据付け指示

## 作業対象と変更してはいけない意味

対象は既存の`topology-heegaard-full-curves`内のHeegaard可視化である。ただしこの束の作成環境では、現在の未コミット作業ツリー・型定義を取得できなかった。したがって以下の接続点は、会話に現れたファイルとAPI名をもとにした契約であり、現物の引数型を確認せず機械的に上書きしてはならない。

数学的設計は固定する。signed cyclic orderだけからfixed-surface routingを推定しない。保存された経路を描画器が変更しない。Dを円板性と向き条件の未検証な塗りつぶしとして扱わない。既存の未コミット変更を破棄・reset・上書きしない。

## 導入順序

### 環境と差分の確認

作業ツリー状態を記録し、既存`types.jl`, `cut_surface.jl`, `surface_chart.jl`, `surface_chart_render.jl`, `joining_disk_layout.jl`, `automatic_layout.jl`, `layout_transport.jl`, `plot_relaxation.jl`, `conversion_spine.jl`を読む。型の実在と向き規約を確認する。新しい数学を設計し直す作業ではなく、型とI/Oの適合箇所を特定する作業である。

本束の`src/MarkedSurfaceDiagrams.jl`を独立submoduleとして追加し、標準ライブラリTOMLとSHAだけで読み込む。既存`Heegaard.jl`の親moduleに合わせて`include`位置のみ調整する。namespace衝突を避け、既存`plot`の意味を変更しない。

### ネイティブ照合を先に通す

```sh
julia --startup-file=no test/runtests.jl
julia --startup-file=no examples/author_annulus.jl
```

テストはPythonで実行済みの`test/expected.toml`に対し、全12例の有理数座標・符号・cycle・seam語・領域位相・ハッシュを照合する。Julia版は配布環境では未実行なので、この段階を省略しない。失敗時にはJulia移植の不整合を修正し、期待値を合わせるために数学データを変更しない。

Pythonはアプリの依存にしない。照合が必要なら`python tools/verify.py`を使う。GitHub Actions用のworkflowは同梱するが、未実行であり、既存CI設定へ勝手に置き換えない。

### 保存と読み込み

新しい文書を別のversioned形式として読み込ませる。既存JSONライブラリでdictを得た後`validate`へ渡す。Juliaファイル単独ではTOMLで試せる。

既存数理coreとは別に、`marked_embedding`のpayloadとして本形式全文を持つ。source-of-truthはこのpayload。既存coreは`legacy_projection(st)`から生成する互換射影とし、保存時に両者を照合する。元のcoreへ後付けで無根拠なmarkingを割り当てない。

既存のlayoutファイルには描画座標と証跡だけを保存する。数理文書のrepresentation hashを必ず紐付ける。別hashのlayoutを読んだら拒否する。数理文書の更新時にはStateとlayoutを作り直す。

### 最初の画面接続

最初は本カーネルの円板atlas SVGを既存画面の新しいviewとして表示する。呼び出しは次で完結する。

```julia
st  = MarkedSurfaceDiagrams.validate(payload)
pos = MarkedSurfaceDiagrams.realize(st)
MarkedSurfaceDiagrams.write_svg(st, pos, output_path)
```

`automatic_layout`、force-directed router、端点だけを受け取る`_route_path`をこの経路に挟まない。word内の局所dart列がrouteである。

既存`surface.js`にはview切替、curve IDによるhighlight、mesh表示切替、数学文書保存、report保存を接続する。SVGの`data-curve`は元のcurve ID。crossing IDはreportのcyclesから取る。数学上の交差にover/underの空白を勝手に加えない。

### 従来の球面＋annulus view

従来viewの形へ戻すことは、初回導入の必須条件にしない。使用する場合には、本marked cellulationと従来atlasの間の明示的な同じmarkingの対応を作る。両方を共通細分し、各triangle上の向きを保つaffine mapとして輸送する。点の対応と辺列を明示し、図式とDも同じ写像で運ぶ。

端点の対応と符号が一致するだけでは不十分。holesの位置を変更するmap、seamを越えるmap、環を一周ねじるmapは、どのmapping classなのかを区別する必要がある。輸送証跡を実装するまでは、旧viewを「認証済み」と表示しない。

### Joining diskとKoda

Dの選択は`region_marks`で保存する。領域idの並び順に依存する数字を、文書版をまたいだIDとして用いない。

画面表示の「Dを左に見る」とKodaのoutward-positive規約を分離する。`check_joining_disk`がfalseならKoda変換を認証済みとして実行しない。仮のcurve反転auditがtrueでも、元の文書の認証がtrueになったわけではない。

`legacy_projection`が返すcrossing signは(α,β)順に直してある。既存coreの符号規約が逆なら、その変換をadapter一箇所に書き、逆変換を含めて試験する。曲線向きの反転によるcycleの逆順も同時に扱う。

### エラーの挙動

エラーはcodeとdetailをそのまま見せる。無効入力を適当な経路に置換して成功扱いにしない。上限超過、precision不足、marking欠落、同位輸送証跡欠落は区別する。写真・スクリーンショットの見た目だけで数学的認証フラグを立てない。

## 接続箇所の契約

| 既存側 | 接続する内容 | 禁止事項 |
|---|---|---|
| `types.jl` | marked payloadを保持するwrapper／別型 | signed ordersを完全なmarked encodingとみなす |
| `render_request.jl` | marked/legacyを明示dispatch | legacy不足時の自動推測 |
| `surface_payload.jl` | hash、curve/crossing/regionの対応を出力 | 配置座標だけを数学正本にする |
| `surface_chart_render.jl` | verified atlasを表示 | 頂点間の再routing |
| `joining_disk_layout.jl` | region trianglesを同じmapで輸送 | Dだけ別座標・別経路で描く |
| `plot_relaxation.jl` | 認証された全体PL変形に限定 | sampled nonintersectionだけで同位を主張 |
| `conversion_spine.jl` | 検証済みcore射影とKoda条件を受領 | 方向修正を複数箇所で重複適用 |
| `surface.js` | view・highlight・数理/幾何保存の分離 | 表示の都合で数理dataを更新 |

## 受入基準

初回据付けの完了条件は、Julia照合テスト全通過、既存アプリで12例のファイル読込・保存再読込、marking/hash mismatchの拒否、円環±1/0の表示区別、交点ゼロ二曲線の保持、Dの認証状態の表示分離である。

既存Abalone・Poincaréへの移行は別の受入段階である。実在の保存routeとatlasの対応を採用し、同位類の由来を残した上で移行する。添付の`torus_three_crossings`は同じ(+,-,+)cycleをもつ明示例だが、ユーザーの従来Abaloneのmarkingを復元したものではない。Poincaréの固定経路fixtureはこの束には含めていない。

## 完了報告に必ず含めるもの

変更ファイル、Julia/アプリの実行結果、数学hash、未対応の輸送・import、旧fixtureを採用した際の出典と規約変換を記載する。「同位類を完全分類した」「教科書画像を自動で完全復元した」「Kodaを形式証明した」とは書かない。
