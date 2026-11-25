# Casson–Walker

Status: draft。規約は [Casson–Walker の規約](../conventions/CassonWalker.md)。<br>
Scope: 手術係数と Seifert データからの λ。HC からの計算は範囲外。<br>
Readers: 関数名と、二つの手術の式が一致する範囲。

## 1. 関数

`Invariants` が次を返す。引数は整数と有理数で、HC は受け取らない。

| 演算名 | 対象となる値 |
|---|---|
| `dedekind_sum(b, a)` | s(b, a)。a > 0 |
| `casson_walker_lens(p, q)` | λ(L(p, q))。p > 0、gcd(p, q) = 1 |
| `casson_walker_knot_surgery(p, q, delta2)` | S³ の結び目の p/q 手術。`delta2` は Δ''(1) |
| `casson_walker_surgery_pup(p, q, delta2)` | q > 0 の展開 |
| `casson_walker_brieskorn(p, q, r)` | [[Sav12#^sav12-sec-0544]] の式 |
| `casson_walker_seifert(b, fibers)` | `(lambda, h1, E)`。E = 0 のときは `nothing` |

`casson_walker_seifert` の式は、E = −b − Σ β/α、S(β/α) = 12·sign(α)·s(β, |α|)、Q = 2 − n + Σ 1/α² として、λ = (3·sign(E) − E − Σ S − Q/E)/24 である。|H₁| = |E|·Π|α| である。

ここで結び目手術の微分は [[Les07#^anonnd70-blk-0313]] の Alexander 多項式のものを指す。[規約の辞書](../conventions/CassonWalker.md#alexander-normalization-dictionary)が保存形式からの変換を所有する。Seifert の一般式の幾何学的同定は `unverified`：Brieskorn と円周束の校正値は整合するが、全入力を覆う出典と向きの比較が必要である。

### Numbered conditions

| 条件 | 入力の条件 | 検査の境界 |
|---|---|---|
| CW1 | 有理数は有限で、途中の整数・有理数演算が表現可能である。 | 有限性と固定幅 overflow は一般には検査しない。大きい演算には `BigInt` を使う。 |
| CW2 | Dedekind 和の `a` は正の整数、`b` は整数。 | `a > 0` は検査する。互いに素でない入力にも同じ有限和を評価する。 |
| CW3 | レンズの座標は整数 `p > 0`、`gcd(p,abs(q)) = 1`。 | 両条件を検査する。名前と向きは規約ページのもの。 |
| CW4 | `S³` の結び目の有限手術係数は整数の組 `(p,q)` で `p != 0`、`q != 0`、`gcd(abs(p),abs(q)) = 1`。 | `q != 0` と互いに素は検査する。`p != 0` は有理ホモロジー球の意味のため呼ぶ側が保証する。負の `q` は許す。 |
| CW5 | `delta2` はその結び目の対称かつ値が 1 で正規化された Alexander 多項式の `Δ''(1)`。 | 整数または有理数を受ける。結び目との対応、対称性、正規化は呼ぶ側が保証する。 |
| CW6 | 正分母展開は CW4 に加えて `q > 0`。 | `q > 0` と互いに素を検査する。 |
| CW7 | Brieskorn の整数重複度 `p,q,r` は正。ホモロジー球の意味には [[Sav12#^sav12-sec-0544]] Theorem 19.4 の互いに素の仮定と向きも必要。 | 正値のみ検査する。定理の残りの仮定は呼ぶ側が保証する。 |
| CW8 | Seifert データは整数 `b` と `Vector{Tuple{Int,Int}}` の組 `(alpha,beta)` で `alpha != 0`、`gcd(abs(alpha),abs(beta)) = 1`。 | 各組の二条件を検査する。`E = 0` は `nothing` に対応する。幾何学的同定は一般には未検証。 |

条件は [casson_walker.jl](../../src/casson_walker.jl) の検査境界を記述する。演算の例外・返り値・保証は各 [契約](../api/casson-walker/INDEX.md)が所有する。

## 2. 一致

S³ の結び目では、補題の形と q > 0 の展開は同じ値を返す。レンズ空間は Δ''(1) = 0 の結び目手術なので、Dedekind 和の値とその展開は同じ値を返す。整数ホモロジー球 Σ(p, q, r) では、[[Sav12#^sav12-sec-0544]] の λ がこの正規化の λ である。

表示を引数に取る `casson_walker(x)` は、まだ名前にしていない。
