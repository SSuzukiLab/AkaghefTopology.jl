# Casson–Walker

Status: draft。規約は [Casson–Walker の規約](../conventions/CassonWalker.md)。<br>
Scope: 手術係数と Seifert データからの λ。HC からの計算は範囲外。<br>
Readers: 関数名と、二つの手術の式が一致する範囲。

## 1. 関数

`Invariants` が次を返す。引数は整数と有理数で、HC は受け取らない。

| 関数 | 値 |
|---|---|
| `dedekind_sum(b, a)` | s(b, a)。a > 0 |
| `casson_walker_lens(p, q)` | λ(L(p, q))。p > 0、gcd(p, q) = 1 |
| `casson_walker_knot_surgery(p, q, delta2)` | S³ の結び目の p/q 手術。`delta2` は Δ''(1) |
| `casson_walker_surgery_pup(p, q, delta2)` | q > 0 の展開 |
| `casson_walker_brieskorn(p, q, r)` | [[Sav12#^sav12-sec-0544]] の式 |
| `casson_walker_seifert(b, fibers)` | `(lambda, h1, E)`。E = 0 のときは `nothing` |

`casson_walker_seifert` の式は、E = −b − Σ β/α、S(β/α) = 12·sign(α)·s(β, |α|)、Q = 2 − n + Σ 1/α² として、λ = (3·sign(E) − E − Σ S − Q/E)/24 である。|H₁| = |E|·Π|α| である。

## 2. 一致

S³ の結び目では、補題の形と q > 0 の展開は同じ値を返す。レンズ空間は Δ''(1) = 0 の結び目手術なので、Dedekind 和の値とその展開は同じ値を返す。整数ホモロジー球 Σ(p, q, r) では、[[Sav12#^sav12-sec-0544]] の λ がこの正規化の λ である。

表示を引数に取る `casson_walker(x)` は、まだ名前にしていない。
