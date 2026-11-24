<!-- XQQQ1500000 -->
# HeegaardCodes の規約

対象はモジュール `HeegaardCodes`（型 `HeegaardCode`、公開別名 `HeegaardDiagram`）である。このモジュールは `LickorishCodes` の上のレイヤーであり（[モジュール構成 §1](../specs/module-structure.md#1-構成)）、曲面・K・平面図の規約は [LickorishCodes の規約](LickorishCodes.md) に従う。

このページは、このモジュールの規約の一覧である。各項目の定義は、リンク先の仕様の節にある。

| 規約 | 定義の場所 | 状態 |
|---|---|---|
| 名称（LC、HC、Heegaard diagram）とレイヤーの区別 | [用語と記号](../specs/heegaard-terminology.md)、[HC 仕様 §1](../specs/heegaard-code.md#1-結論とレイヤー構造) | あかげふの決定 2026-11-09 |
| 役割を行で明示し、K 以外の成分すべてに役割を指定する | [HC 仕様 §2.1](../specs/heegaard-code.md#21-役割を行で明示する) | あかげふの確定 2026-11-09 |
| `point` と `weight` を HC に置く | [HC 仕様 §3](../specs/heegaard-code.md#3-hc-に置く追加情報) | あかげふの確定 2026-11-09 |
| 曲面の向きは ∂Hα。Hα の外向き法線を先に置くと右手系になる | [Koda 変換仕様 §6.1](../specs/heegaard-koda-conversion.md#61-b-の構成と向き) | 仕様 |
| K08 Fig.2 の局所規約。β が over、α が under（候補 A） | [Koda 変換仕様 §7.3](../specs/heegaard-koda-conversion.md#73-向きの規約と局所符号) | 採用規約は確定 |
| 細分領域番号は LC §3.1 の番号をそのまま使う | [Koda 変換仕様 §3.1](../specs/heegaard-koda-conversion.md#31-有限な回転系) | あかげふの決定 2026-11-16 |
<!-- XQQQ1000000 -->
