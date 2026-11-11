<!-- XQQQ1500000 -->
# HeegaardDiagram の公開境界

Status: ライブラリの型の位置。2026-11-11 の仕様追記<br>
Scope: `HeegaardDiagram` と `VirtualLink` の同格、および Koda 変換の射

`HeegaardDiagram` は、このライブラリが Heegaard diagram を表す公開型であり、`VirtualLink` と同格である。

その実体は解決済みの HC であり、構造体 `HeegaardLayer` と同じ型である。曲面と Lickorish の語は LC に、α・β の役割と joining disk はその HC の欄にある。`standard` は α の支持を既存の meridian `m` に束縛する入力略記であり、別の計算型ではない。Koda の standard punctured diagram でもない。

セル分割で認証した図式は公開名 `HeegaardDiagram` を使わない。その構造体名は `CellularHeegaardDiagram` である。語の配列を `VirtualLink` の `gauss` と `orientation` に置くことは、語の保存である。

o-data は `VirtualLink` の表示である。normal o-graph の成分語と vertex sign は、この表示の中にある。

Koda 変換は、`HeegaardDiagram` から `VirtualLink` の o-data への射である。定義域の検査、joining disk、spine の対応は [Koda 変換仕様](heegaard-koda-conversion.md) が定める。用語は [Heegaard の用語と記号](heegaard-terminology.md) に従う。
<!-- XQQQ1000000 -->
