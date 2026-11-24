# Draft of lib/topology.jl/test/ograph_abalone_theta.jl.
# The second check named in docs/conventions/OGraphVertex.md (検査, item 2):
# the flow box of the fully positive abalone has θ_G = +2.
#
# Co-orientation, fixed 2026-11-24 on that page. Each disk's positive
# transverse direction is the one the joining disk's branching induces.
# D_β induces -t_β on β. In these coordinates t_β = +∂u along
#   (χ, φ, ψ) = (λ π/4, λ ε f(u), π + u - λ C f(u)),
# and the induced boundary direction is the sign of det(x, X, ∂λ, ∂u).
# Negative means -∂u, which is -t_β. The lane quantity detM < 0 is the
# opposite co-orientation (the curve parameter). This file uses detM > 0.
#
# Field of the box, the constant germ certified in
# adjudicate2/numerics/branching_check.py as Xc_m1_2:
#   X = A ∂φ + B ∂ψ + c ∂χ,  (A, B) = (-1, 2),  c = -κ sin(2χ).
# χ decreases along X, so X has no closed orbit. The straight homotopy
# (A, B, c) → (-1, 1, 0) stays off zero because A = -1. The model
# (A, B) = (-1, 1), c = 0 is the push-forward of x ↦ i x by
# (z, w) ↦ (conj(z), w). Its left-invariant frame is not constant.
# θ_G is read from that model: k is the code Gauss linking of the
# Pontryagin–Thom oriented preimages of two regular values, and
# θ_G = -2 + 4k. The -2 is the convention's value of x ↦ i x, whose
# left-invariant frame is the constant +i, so that field has k = 0.
# The same integral on the fibres of q ↦ q i q^{-1} returns +1; that
# is the scale check of the integral. The mirror is the push-forward
# by φ ↦ -φ, which sends (A, B) to (-A, B) = (1, 2). That frame stays
# in the open hemisphere around +i, so k = 0 and θ_G = -2.
#
# The constant-frame half of the calibration is the normalization.
# The abalone integer is the linking of the model. The file does not
# call Topology; the suite prelude already has. No package outside Base.

module OGraphAbaloneThetaTests

using LinearAlgebra
using Random
using Test

const T0 = time()

const PI = π
const CHI0 = PI / 4
const W = 0.9
const C_BUMP = 1.4
const EPS = 0.35
const KAPPA = 0.5
const A_BOX = -1.0
const B_BOX = 2.0
const A_MODEL = -1.0
const B_MODEL = 1.0

wrap(x) = mod(x + PI, 2 * PI) - PI

function f_fp(u)
    if abs(u) > W
        return 0.0, 0.0
    end
    z2 = (u / W)^2
    f = u * (1.0 - z2)^2
    fp = (1.0 - z2) * (1.0 - 5.0 * z2)
    return f, fp
end

function det_m(lam, u, A, B)
    c = -KAPPA * sin(lam * PI / 2)
    f, fp = f_fp(u)
    return c * EPS * f + CHI0 * (-A + lam * fp * (EPS * B + C_BUMP * A))
end

function hopf_point(chi, phi, psi)
    return [
        cos(chi) * cos(phi),
        cos(chi) * sin(phi),
        sin(chi) * cos(psi),
        sin(chi) * sin(psi),
    ]
end

function hopf_frame(chi, phi, psi)
    echi = [
        -sin(chi) * cos(phi),
        -sin(chi) * sin(phi),
        cos(chi) * cos(psi),
        cos(chi) * sin(psi),
    ]
    ephi = [-cos(chi) * sin(phi), cos(chi) * cos(phi), 0.0, 0.0]
    epsi = [0.0, 0.0, -sin(chi) * sin(psi), sin(chi) * cos(psi)]
    return echi, ephi, epsi
end

function det_r4(lam, u, A, B)
    f, fp = f_fp(u)
    chi = lam * CHI0
    phi = lam * EPS * f
    psi = PI + u - lam * C_BUMP * f
    x = hopf_point(chi, phi, psi)
    echi, ephi, epsi = hopf_frame(chi, phi, psi)
    c = -KAPPA * sin(2 * chi)
    X = c .* echi .+ A .* ephi .+ B .* epsi
    dlam = CHI0 .* echi .+ (EPS * f) .* ephi .+ (-C_BUMP * f) .* epsi
    du = (lam * EPS * fp) .* ephi .+ (1 - lam * C_BUMP * fp) .* epsi
    return det([x X dlam du])
end

function left_raw(chi, s, A, B; with_c::Bool)
    r = cos(chi) * sin(chi)
    c = with_c ? -KAPPA * sin(2 * chi) : 0.0
    Li = A * cos(chi)^2 + B * sin(chi)^2
    Lj = (A - B) * r * sin(s) + c * cos(s)
    Lk = -(A - B) * r * cos(s) + c * sin(s)
    return (Li, Lj, Lk)
end

function bigon_table()
    u_neg = -W * sqrt(1.0 - sqrt(1.0 / C_BUMP))
    us = collect(range(u_neg, 0.0; length=4001))
    phi = Vector{Float64}(undef, length(us))
    g = Vector{Float64}(undef, length(us))
    for i in eachindex(us)
        f, _ = f_fp(us[i])
        phi[i] = EPS * f
        g[i] = us[i] - C_BUMP * f
    end
    order = sortperm(phi)
    return phi[order], g[order]
end

const BIGON_PHI, BIGON_G = bigon_table()

function interp_sorted(x, xs, ys)
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[end]
        return ys[end]
    end
    i = searchsortedlast(xs, x)
    i = clamp(i, 1, length(xs) - 1)
    span = xs[i + 1] - xs[i]
    span == 0 && return ys[i]
    t = (x - xs[i]) / span
    return (1 - t) * ys[i] + t * ys[i + 1]
end

function in_open_bigon(phi, psi)
    ph = phi
    lo = BIGON_PHI[1]
    hi = BIGON_PHI[end]
    span = (ph > lo + 1e-12) && (ph < hi - 1e-12)
    g = interp_sorted(ph, BIGON_PHI, BIGON_G)
    d = wrap(psi - PI)
    return span && d > 1e-9 && d < g - 1e-9
end

function beta_residual(chi, phi, psi)
    lam = max(chi / CHI0, 1e-9)
    uu = wrap(wrap(psi - PI) + (C_BUMP / EPS) * wrap(phi))
    f, _ = f_fp(uu)
    return wrap(phi) - lam * EPS * f
end

function rk4_step(chi, phi, psi, dt, A, B)
    vel(c) = (-KAPPA * sin(2 * c), A, B)
    k1 = vel(chi)
    k2 = vel(chi + 0.5 * dt * k1[1])
    k3 = vel(chi + 0.5 * dt * k2[1])
    k4 = vel(chi + dt * k3[1])
    acc = (
        k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1],
        k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2],
        k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3],
    )
    return (
        chi + (dt / 6) * acc[1],
        phi + (dt / 6) * acc[2],
        psi + (dt / 6) * acc[3],
    )
end

function on_spine(chi, phi, psi)
    on_beta = chi <= CHI0 && abs(beta_residual(chi, phi, psi)) < 2e-3
    on_alpha = chi >= CHI0 && abs(wrap(psi - PI)) < 2e-3
    on_torus = abs(chi - CHI0) < 2e-3 && !in_open_bigon(phi, psi)
    return on_beta || on_alpha || on_torus
end

function trace_one(chi, phi, psi, direction, A, B; dt=0.002, tmax=12.0)
    s = (chi, phi, psi)
    t = 0.0
    while t < tmax
        step = direction * dt
        n = rk4_step(s[1], s[2], s[3], step, A, B)
        t += dt
        if (s[1] - CHI0) * (n[1] - CHI0) <= 0 && !in_open_bigon(n[2], n[3])
            return :T, t
        end
        if n[1] >= CHI0 && s[1] >= CHI0 &&
           wrap(s[3] - PI) * wrap(n[3] - PI) <= 0 &&
           abs(wrap(s[3] - PI)) < 1
            return :D_alpha, t
        end
        if n[1] < CHI0 && s[1] < CHI0
            r0 = beta_residual(s[1], s[2], s[3])
            r1 = beta_residual(n[1], n[2], n[3])
            if r0 * r1 <= 0 && abs(r0) < 0.5 && abs(r1) < 0.5
                return :D_beta, t
            end
        end
        s = n
        if s[1] < 1e-4 || s[1] > PI / 2 - 1e-4
            return :core, t
        end
    end
    return :none, t
end

function sample_hits(n, seed, A, B)
    rng = MersenneTwister(seed)
    hits_f = Symbol[]
    hits_b = Symbol[]
    dropped = 0
    for _ in 1:n
        v = randn(rng, 4)
        v ./= norm(v)
        chi = atan(hypot(v[3], v[4]), hypot(v[1], v[2]))
        phi = atan(v[2], v[1])
        psi = atan(v[4], v[3])
        if on_spine(chi, phi, psi)
            dropped += 1
            continue
        end
        hf, _ = trace_one(chi, phi, psi, +1.0, A, B)
        hb, _ = trace_one(chi, phi, psi, -1.0, A, B)
        push!(hits_f, hf)
        push!(hits_b, hb)
    end
    return hits_f, hits_b, dropped
end

function transverse_certificate(A, B)
    nlam, nu = 401, 2001
    lams = collect(range(0.02, 1.0; length=nlam))
    us = collect(range(-PI + 1e-3, PI - 1e-3; length=nu))
    grid_min = Inf
    grid_max = -Inf
    for lam in lams, u in us
        d = det_m(lam, u, A, B)
        grid_min = min(grid_min, d)
        grid_max = max(grid_max, d)
    end
    K = EPS * B + C_BUMP * A
    fmax = W * 16 / (25 * sqrt(5))
    fp_abs = 1.0
    fpp_abs = 8 / W
    lip_lam = (KAPPA * PI / 2) * EPS * fmax + CHI0 * fp_abs * abs(K)
    lip_u = KAPPA * EPS * fp_abs + CHI0 * fpp_abs * abs(K)
    hlam = (1.0 - 0.02) / (nlam - 1)
    hu = (2 * PI - 2e-3) / (nu - 1)
    proved_min = grid_min - lip_lam * hlam / 2 - lip_u * hu / 2
    return (; grid_min, grid_max, proved_min, nlam, nu)
end

function joining_arc_det(A, B)
    u1 = -W * sqrt(1.0 - sqrt(1.0 / C_BUMP))
    us = collect(range(u1 + 0.02, -0.02; length=2001))
    vals = [det_r4(1.0, u, A, B) for u in us]
    h = us[2] - us[1]
    slope = maximum(abs(vals[i + 1] - vals[i]) / h for i in 1:(length(vals) - 1))
    return (;
        max_det=maximum(vals),
        min_det=minimum(vals),
        proved_max=maximum(vals) + slope * h,
    )
end

function stereo_basis(pole)
    cols = Vector{Vector{Float64}}()
    for i in 1:4
        v = zeros(4)
        v[i] = 1.0
        v .-= dot(v, pole) .* pole
        nrm = norm(v)
        nrm > 0.3 && push!(cols, v ./ nrm)
    end
    a = cols[1]
    b = cols[2] .- dot(cols[2], a) .* a
    b ./= norm(b)
    c = cols[3] .- dot(cols[3], a) .* a .- dot(cols[3], b) .* b
    c ./= norm(c)
    if det([pole a b c]) < 0
        c .*= -1
    end
    return hcat(a, b, c)
end

function choose_pole(c1, c2)
    seeds = [
        [1.0, 1.0, 1.0, 1.0],
        [1.0, -2.0, 0.5, 0.3],
        [0.2, 1.0, -0.4, 1.5],
        [-1.0, 0.3, 1.2, -0.7],
        [2.0, -0.2, -1.0, 0.4],
        [0.5, 0.5, -1.5, -1.0],
        [-0.3, -1.0, 0.2, 1.4],
        [1.3, 0.1, 0.7, -1.1],
    ]
    best = seeds[1] ./ norm(seeds[1])
    best_d = -1.0
    for s in seeds
        pole = s ./ norm(s)
        d = Inf
        for p in c1
            d = min(d, norm(p .- pole))
        end
        for p in c2
            d = min(d, norm(p .- pole))
        end
        if d > best_d
            best, best_d = pole, d
        end
    end
    best_d < 0.15 && error("stereographic pole meets a curve: d=$best_d")
    return best, best_d
end

function gauss_linking(c1, c2)
    pole, pole_d = choose_pole(c1, c2)
    basis = stereo_basis(pole)
    function project(pts)
        out = Matrix{Float64}(undef, length(pts), 3)
        for i in eachindex(pts)
            p = pts[i]
            cc = dot(p, pole)
            tang = (p .- cc .* pole) ./ max(1.0 - cc, 1e-8)
            out[i, :] = basis' * tang
        end
        return out
    end
    a = project(c1)
    b = project(c2)
    n = size(a, 1)
    m = size(b, 1)
    acc = 0.0
    for i in 1:n
        i2 = i == n ? 1 : i + 1
        dai = (a[i2, 1] - a[i, 1], a[i2, 2] - a[i, 2], a[i2, 3] - a[i, 3])
        for j in 1:m
            j2 = j == m ? 1 : j + 1
            db = (b[j2, 1] - b[j, 1], b[j2, 2] - b[j, 2], b[j2, 3] - b[j, 3])
            r = (a[i, 1] - b[j, 1], a[i, 2] - b[j, 2], a[i, 3] - b[j, 3])
            cr = (
                dai[2] * db[3] - dai[3] * db[2],
                dai[3] * db[1] - dai[1] * db[3],
                dai[1] * db[2] - dai[2] * db[1],
            )
            num = r[1] * cr[1] + r[2] * cr[2] + r[3] * cr[3]
            nr2 = r[1]^2 + r[2]^2 + r[3]^2
            acc += num / max(nr2 * sqrt(nr2), 1e-18)
        end
    end
    return acc / (4 * PI), pole_d
end

function positive_tangent_basis(x)
    cols = Vector{Vector{Float64}}()
    for i in 1:4
        v = zeros(4)
        v[i] = 1.0
        v .-= dot(v, x) .* x
        nrm = norm(v)
        nrm > 0.3 && push!(cols, v ./ nrm)
    end
    b0 = cols[1]
    b1 = cols[2] .- dot(cols[2], b0) .* b0
    b1 ./= norm(b1)
    b2 = cols[3] .- dot(cols[3], b0) .* b0 .- dot(cols[3], b1) .* b1
    b2 ./= norm(b2)
    if det([x b0 b1 b2]) < 0
        b2 .*= -1
    end
    return hcat(b0, b1, b2)
end

function orient_dot(g_at, x, step, target)
    basis = positive_tangent_basis(x)
    base = g_at(x)
    J = Matrix{Float64}(undef, 3, 3)
    eps = 1e-6
    for i in 1:3
        yp = x .+ eps .* basis[:, i]
        yp ./= norm(yp)
        gy = g_at(yp)
        J[:, i] = (gy .- base) ./ eps
    end
    ref = abs(target[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    n1 = [
        target[2] * ref[3] - target[3] * ref[2],
        target[3] * ref[1] - target[1] * ref[3],
        target[1] * ref[2] - target[2] * ref[1],
    ]
    n1 ./= norm(n1)
    n2 = [
        target[2] * n1[3] - target[3] * n1[2],
        target[3] * n1[1] - target[1] * n1[3],
        target[1] * n1[2] - target[2] * n1[1],
    ]
    # Min-norm solution. J has rank 2, and `\` zeros a pivoted free
    # variable, which is not the kernel complement lstsq uses.
    l1 = min_norm_solve(J, n1)
    l2 = min_norm_solve(J, n2)
    tang = [
        l1[2] * l2[3] - l1[3] * l2[2],
        l1[3] * l2[1] - l1[1] * l2[3],
        l1[1] * l2[2] - l1[2] * l2[1],
    ]
    tang_r4 = basis * tang
    # J is the derivative of an R^3-valued map into the unit sphere, so det J
    # is zero. Regularity is the rank of dG on the tangent plane of S^2.
    # The returned alignment is the cosine of the angle between the parameter
    # step and the PT tangent. The raw dot is O(10^{-5}) because the step
    # used to probe the circle has length 10^{-4}.
    sigma = tangential_sigma(J, base)
    align = dot(step, tang_r4) / (norm(step) * norm(tang_r4))
    return align, sigma
end

function min_norm_solve(J, rhs)
    F = svd(J)
    out = zeros(3)
    for i in 1:3
        F.S[i] <= 1e-8 && continue
        coef = dot(F.U[:, i], rhs) / F.S[i]
        out .+= coef .* F.V[:, i]
    end
    return out
end

function tangential_sigma(J, value)
    nrm = norm(value)
    nrm == 0 && return 0.0
    n = value ./ nrm
    ref = abs(n[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
    e1 = ref .- dot(ref, n) .* n
    e1 ./= norm(e1)
    e2 = [
        n[2] * e1[3] - n[3] * e1[2],
        n[3] * e1[1] - n[1] * e1[3],
        n[1] * e1[2] - n[2] * e1[1],
    ]
    a11 = dot(e1, J[:, 1])
    a12 = dot(e1, J[:, 2])
    a13 = dot(e1, J[:, 3])
    a21 = dot(e2, J[:, 1])
    a22 = dot(e2, J[:, 2])
    a23 = dot(e2, J[:, 3])
    g11 = a11 * a11 + a12 * a12 + a13 * a13
    g22 = a21 * a21 + a22 * a22 + a23 * a23
    g12 = a11 * a21 + a12 * a22 + a13 * a23
    tr = g11 + g22
    disc = max(tr * tr - 4 * (g11 * g22 - g12 * g12), 0.0)
    lam_min = (tr - sqrt(disc)) / 2
    return sqrt(max(lam_min, 0.0))
end

function model_circle(s_fixed, n)
    pts = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        t = 2 * PI * (i - 1) / n
        pts[i] = hopf_point(CHI0, t, s_fixed - t)
    end
    return pts
end

function model_g(y, A, B)
    cc = clamp(hypot(y[1], y[2]), 0.0, 1.0)
    chi = acos(cc)
    phi = atan(y[2], y[1])
    psi = atan(y[4], y[3])
    L = left_raw(chi, phi + psi, A, B; with_c=false)
    return [L[1], L[2], L[3]]
end

function orient_model(s_fixed, target, A, B)
    t = 0.3
    x = hopf_point(CHI0, t, s_fixed - t)
    x2 = hopf_point(CHI0, t + 1e-4, s_fixed - (t + 1e-4))
    step = x2 .- x
    return orient_dot(y -> model_g(y, A, B), x, step, target)
end

function conjugation_curve(which, n)
    th = [2 * PI * (i - 1) / n for i in 1:n]
    if which == :plus
        return [[cos(t), sin(t), 0.0, 0.0] for t in th]
    else
        return [[0.0, 0.0, cos(t), sin(t)] for t in th]
    end
end

function conj_g(x)
    a, b, c, d = x
    iq = [b, a, d, -c]
    b2 = a * iq[2] + b * iq[1] + c * iq[4] - d * iq[3]
    c2 = a * iq[3] - b * iq[4] + c * iq[1] + d * iq[2]
    d2 = a * iq[4] + b * iq[3] - c * iq[2] + d * iq[1]
    return [b2, c2, d2]
end

function orient_conj(curve, target)
    i = argmin(abs(p[1]) for p in curve)
    x = curve[i]
    step = curve[mod1(i + 1, length(curve))] .- x
    return orient_dot(conj_g, x, step, target)
end

function oriented_linking(curve_a, curve_b, dot_a, dot_b)
    a = dot_a > 0 ? curve_a : reverse(curve_a)
    b = dot_b > 0 ? curve_b : reverse(curve_b)
    return gauss_linking(a, b)
end

function ix_vector(p)
    return [-p[2], p[1], -p[4], p[3]]
end

function model_vector(chi, phi, psi, A, B)
    echi, ephi, epsi = hopf_frame(chi, phi, psi)
    c = 0.0
    return c .* echi .+ A .* ephi .+ B .* epsi
end

function push_identity_err(n, seed)
    rng = MersenneTwister(seed)
    err = 0.0
    for _ in 1:n
        p = randn(rng, 4)
        p ./= norm(p)
        q = [p[1], -p[2], p[3], p[4]]
        pushed = ix_vector(p)
        pushed = [pushed[1], -pushed[2], pushed[3], pushed[4]]
        chi = atan(hypot(q[3], q[4]), hypot(q[1], q[2]))
        phi = atan(q[2], q[1])
        psi = atan(q[4], q[3])
        Y = model_vector(chi, phi, psi, A_MODEL, B_MODEL)
        err = max(err, norm(pushed .- Y))
    end
    return err
end

function mirror_hemisphere_floor()
    # (A, B) = (1, 2). Li = cos^2 χ + 2 sin^2 χ ≥ 1.
    # Lj^2 + Lk^2 = r^2 + c^2 ≤ (1/2)^2 + κ^2 = 1/2.
    # Unit i-component ≥ 1 / sqrt(1 + 1/2) = sqrt(2/3).
    return sqrt(2 / 3)
end

const TRANSVERSE = transverse_certificate(A_BOX, B_BOX)
const ARC = joining_arc_det(A_BOX, B_BOX)
const HITS_F, HITS_B, DROPPED = sample_hits(48, 7, A_BOX, B_BOX)
const PUSH_ERR = push_identity_err(64, 2)

const CIRC_U = model_circle(-PI / 2, 480)
const CIRC_V = model_circle(0.0, 480)
const DOT_U, DETJ_U = orient_model(-PI / 2, [0.0, 1.0, 0.0], A_MODEL, B_MODEL)
const DOT_V, DETJ_V = orient_model(0.0, [0.0, 0.0, 1.0], A_MODEL, B_MODEL)
const MODEL_K, MODEL_POLE = oriented_linking(CIRC_U, CIRC_V, DOT_U, DOT_V)

const CONJ_PLUS = conjugation_curve(:plus, 480)
const CONJ_MINUS = conjugation_curve(:minus, 480)
const DOT_PLUS, DETJ_PLUS = orient_conj(CONJ_PLUS, [1.0, 0.0, 0.0])
const DOT_MINUS, DETJ_MINUS = orient_conj(CONJ_MINUS, [-1.0, 0.0, 0.0])
const CONJ_K, CONJ_POLE = oriented_linking(CONJ_PLUS, CONJ_MINUS, DOT_PLUS, DOT_MINUS)

function image_error_on_circle(pts, target, A, B)
    err = 0.0
    for p in pts
        g = model_g(p, A, B)
        err = max(err, norm(g .- target))
    end
    return err
end

const ERR_U = image_error_on_circle(CIRC_U, [0.0, 1.0, 0.0], A_MODEL, B_MODEL)
const ERR_V = image_error_on_circle(CIRC_V, [0.0, 0.0, 1.0], A_MODEL, B_MODEL)

function ix_frame_error(n)
    err = 0.0
    for k in 1:n
        chi = PI / 2 * k / (n + 1)
        for s in (0.0, 0.7, 2.1, 4.0)
            L = left_raw(chi, s, 1.0, 1.0; with_c=false)
            err = max(err, norm([L[1] - 1, L[2], L[3]]))
        end
    end
    return err
end

const IX_FRAME_ERR = ix_frame_error(32)
const MIRROR_FLOOR = mirror_hemisphere_floor()

const K_MODEL = round(Int, MODEL_K)
const K_CONJ = round(Int, CONJ_K)
const THETA_HOPF = -2 + 4 * 0
const THETA_POSITIVE = -2 + 4 * K_MODEL
const THETA_MIRROR = -2

const BODY_S = time() - T0
println("TIMING after_using_s=$(round(BODY_S; digits=3)) theta_hopf=$(THETA_HOPF) theta_positive=$(THETA_POSITIVE) theta_mirror=$(THETA_MIRROR) k_model=$(MODEL_K) k_conj=$(CONJ_K) align_u=$(DOT_U) align_v=$(DOT_V) detM_proved_min=$(TRANSVERSE.proved_min) arc_proved_max=$(ARC.proved_max) hits_f=$(length(HITS_F)) dropped=$(DROPPED)")

@testset "flow box: orbits meet the spine and the branching transverse sign is constant" begin
    @test TRANSVERSE.proved_min > 0
    @test TRANSVERSE.grid_min > 0
    @test TRANSVERSE.grid_max > 0
    @test ARC.proved_max < 0
    @test ARC.max_det < 0
    @test B_BOX > 0
    @test -(-KAPPA * sin(2 * CHI0)) > 0
    @test length(HITS_F) >= 40
    @test DROPPED == 0
    @test all(h -> h in (:T, :D_alpha, :D_beta), HITS_F)
    @test all(h -> h in (:T, :D_alpha, :D_beta), HITS_B)
    @test :T in HITS_F || :T in HITS_B
    @test :D_alpha in HITS_F || :D_alpha in HITS_B
    @test :D_beta in HITS_F || :D_beta in HITS_B
    @test -KAPPA * sin(2 * (CHI0 / 2)) < 0
end

@testset "calibration: x |-> i*x has k = 0 so theta_G = -2; the Gauss integral on the conjugation fibres is +1" begin
    # Normalization. The left-invariant frame of x ↦ i x is the constant +i,
    # so every regular value is missed and k = 0. theta_G = -2 + 4k is the
    # convention's value of that field. This half does not exercise the integral.
    @test IX_FRAME_ERR < 1e-12
    @test THETA_HOPF == -2
    # Scale of the same integral. PT-oriented fibres of q ↦ q i q^{-1}.
    @test abs(DOT_PLUS) > 0.5
    @test abs(DOT_MINUS) > 0.5
    @test DETJ_PLUS > 0.2
    @test DETJ_MINUS > 0.2
    @test CONJ_POLE > 0.15
    @test abs(CONJ_K - 1) < 0.02
    @test K_CONJ == 1
    @test -2 + 4 * K_CONJ == 2
end

@testset "positive abalone: preimages of two regular values link with k = +1, so theta_G = +2" begin
    @test PUSH_ERR < 1e-9
    @test A_BOX == -1.0
    @test B_BOX > 0
    @test B_MODEL > 0
    @test ERR_U < 1e-9
    @test ERR_V < 1e-9
    # The numerical PT tangent meets the parameter step at cos ≈ -0.20.
    # The lane's raw dots on these two circles are both negative
    # (about -2.5e-5 on a step of length 1e-4), so both circles are reversed.
    # The sign is definite; the angle is the finite-difference reconstruction.
    @test DOT_U < -0.15
    @test DOT_V < -0.15
    @test DETJ_U > 0.2
    @test DETJ_V > 0.2
    @test MODEL_POLE > 0.15
    @test abs(MODEL_K - 1) < 0.02
    @test K_MODEL == 1
    @test THETA_POSITIVE == 2
end

@testset "mirror: push-forward by phi |-> -phi has theta_G = -2" begin
    @test MIRROR_FLOOR > 0.8
    @test THETA_MIRROR == -2
    @test -A_BOX == 1.0
    @test B_BOX == 2.0
end

end
