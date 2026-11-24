# Consequences of docs/conventions/LensSpaceNames.md (fixed 2026-11-23): relative orientations of lens-space presentations.
#
# The value used is a twisted Z_N Dijkgraaf–Witten state sum on the dual triangulation of a closed normal o-graph
# (N prime, cocycle ω(a,b,c) = ζ^{a·carry(b,c)}, ζ = exp(2πi/N), carry(b,c) = 1 if b + c ≥ N). The tetrahedron of a
# vertex is oriented by its over pair (slot order 1,2,3,4 positive when the over pair is (2,4)), and the vertex order is
# the claim 3.1 local rule of PJ261120. Every assertion compares two presentations (equal or complex conjugate), so the
# global orientation sign of this rule does not enter; reversing it conjugates every value.
isdefined(@__MODULE__, :vl_example) || include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))

const _LS_VL = Topology.VirtualLinks
const _LS_SIG_PLUS = ((:under, :head), (:over, :tail), (:under, :tail), (:over, :head))
const _LS_SIG_MINUS = ((:under, :head), (:over, :head), (:under, :tail), (:over, :tail))

function _ls_permsign(o)
    s = 1
    for i in 1:length(o), j in i+1:length(o)
        o[i] > o[j] && (s = -s)
    end
    s
end

function _ls_find!(parent, rel, i)
    parent[i] == i && return (i, 0)
    (r, p) = _ls_find!(parent, rel, parent[i])
    parent[i] = r; rel[i] = mod(rel[i] + p, 2)
    (r, rel[i])
end

"Nullspace basis of A x = 0 over Z/N (N prime)."
function _ls_nullspace_mod(A::Matrix{Int}, N::Int)
    A = mod.(A, N); m, n = size(A); pivots = Int[]; row = 1
    for col in 1:n
        row > m && break
        p = findfirst(r -> A[r, col] != 0, row:m)
        p === nothing && continue
        p += row - 1
        A[row, :], A[p, :] = A[p, :], A[row, :]
        A[row, :] = mod.(A[row, :] .* invmod(A[row, col], N), N)
        for r in 1:m
            r != row && A[r, col] != 0 && (A[r, :] = mod.(A[r, :] .- A[r, col] .* A[row, :], N))
        end
        push!(pivots, col); row += 1
    end
    free = setdiff(1:n, pivots)
    basis = Vector{Int}[]
    for f in free
        v = zeros(Int, n); v[f] = 1
        for (i, pc) in enumerate(pivots)
            v[pc] = mod(-A[i, f], N)
        end
        push!(basis, v)
    end
    basis
end

"Twisted Z_N Dijkgraaf–Witten value of the closed normal o-graph presented by `link`."
function _ls_dw(link; N::Int=3, k::Int=1)
    nc = _LS_VL.normal_coloured_ograph(ograph_state(link))
    g = nc.graph
    n = length(g.vertex_halfedges)
    own = Dict{Int,Tuple{Int,Int}}()
    for (v, cyc) in enumerate(g.vertex_halfedges), (slot, h) in enumerate(cyc)
        own[Int(h)] = (v, slot)
    end
    chart = _LS_VL._leg_chart(g, :b)
    colour = Dict(minmax(Int(r.halfedge), Int(g.mate[r.halfedge])) => mod(Int(r.colour), 3) for r in g.edge_colours)
    gluings = []
    for (h, m) in enumerate(g.mate)
        h < Int(m) || continue
        c = colour[(h, Int(m))]
        phi = Dict(own[chart[h][i]][2] => own[chart[Int(m)][_LS_VL.fig25_frame_reversal(_LS_VL.fig25_tau(c, i))]][2] for i in 1:3)
        push!(gluings, (v1=own[h][1], k1=own[h][2], v2=own[Int(m)][1], k2=own[Int(m)][2], phi=phi))
    end
    # orientation from the over pairs, checked against every face gluing
    s = [Tuple(sort(collect(g.over_positions[v]))) == (2, 4) ? 1 : -1 for v in 1:n]
    for x in gluings
        img = [x.phi[q] for q in 1:4 if q != x.k1]
        @assert s[x.v1] * (-1)^(x.k1 - 1) * _ls_permsign(img) == -s[x.v2] * (-1)^(x.k2 - 1) "gluing does not reverse orientation"
    end
    # vertex order: claim 3.1 local rule by slot-type signature
    kind = Dict{Int,Symbol}(); for e in nc.edges; kind[Int(e.tail)] = :tail; kind[Int(e.head)] = :head; end
    order = Vector{NTuple{4,Int}}(undef, n)
    for (v, cyc) in enumerate(g.vertex_halfedges)
        sig = Tuple((slot in g.over_positions[v] ? :over : :under, kind[Int(h)]) for (slot, h) in enumerate(cyc))
        order[v] = sig == _LS_SIG_PLUS ? (3, 1, 2, 4) : sig == _LS_SIG_MINUS ? (1, 3, 2, 4) : error("unexpected slot types $sig")
    end
    rank(v, q) = findfirst(==(q), order[v])
    for x in gluings      # branching: every face gluing preserves the vertex order
        face = sort([q for q in 1:4 if q != x.k1]; by=q -> rank(x.v1, q))
        @assert issorted([rank(x.v2, x.phi[q]) for q in face]) "the local rule is not a branching here"
    end
    # edge classes of the triangulation (signed union-find on oriented tetrahedron edges)
    keys_ = [(v, a, b) for v in 1:n for a in 1:4 for b in a+1:4]
    idx = Dict(key => i for (i, key) in enumerate(keys_))
    parent = collect(eachindex(keys_)); rel = zeros(Int, length(keys_))
    for x in gluings
        f = [q for q in 1:4 if q != x.k1]
        for i in 1:3, j in i+1:3
            a, b = f[i], f[j]; A, B = x.phi[a], x.phi[b]
            (ri, pi_) = _ls_find!(parent, rel, idx[(x.v1, a, b)]); (rj, pj) = _ls_find!(parent, rel, idx[(x.v2, min(A, B), max(A, B))])
            want = A < B ? 0 : 1
            if ri == rj
                @assert mod(pi_ + pj + want, 2) == 0 "inconsistent edge orientation"
            else
                parent[ri] = rj; rel[ri] = mod(pi_ + pj + want, 2)
            end
        end
    end
    roots = Dict{Int,Int}(); cls = Dict{Tuple{Int,Int,Int},Tuple{Int,Int}}()
    for (i, key) in enumerate(keys_)
        (r, p) = _ls_find!(parent, rel, i)
        haskey(roots, r) || (roots[r] = length(roots) + 1)
        cls[key] = (roots[r], p == 0 ? 1 : -1)
    end
    E = length(roots)
    # vertex classes of the triangulation
    vparent = collect(1:4n)
    vfind(i) = vparent[i] == i ? i : (vparent[i] = vfind(vparent[i]))
    for x in gluings, q in 1:4
        q == x.k1 && continue
        a = vfind(4(x.v1 - 1) + q); b = vfind(4(x.v2 - 1) + x.phi[q]); a != b && (vparent[a] = b)
    end
    V = length(unique(vfind(i) for i in 1:4n))
    # flat Z_N labellings: kernel of the coboundary on every triangle of every tetrahedron
    coeff(v, a, b) = a < b ? (cls[(v, a, b)][1], cls[(v, a, b)][2]) : (cls[(v, b, a)][1], -cls[(v, b, a)][2])
    rows = Vector{Int}[]
    for v in 1:n, a in 1:4, b in a+1:4, c in b+1:4
        r = zeros(Int, E)
        for (p, q, sgn) in ((a, b, 1), (b, c, 1), (a, c, -1))
            (cl, o) = coeff(v, p, q); r[cl] += sgn * o
        end
        push!(rows, r)
    end
    basis = _ls_nullspace_mod(permutedims(reduce(hcat, rows)), N)
    label(x, v, p, q) = (cl_o = coeff(v, p, q); mod(cl_o[2] * x[cl_o[1]], N))
    carry(b, c) = b + c >= N ? 1 : 0
    total = 0.0 + 0.0im
    for coefs in Iterators.product(ntuple(_ -> 0:N-1, length(basis))...)
        x = isempty(basis) ? zeros(Int, E) : mod.(sum(coefs[i] .* basis[i] for i in eachindex(basis)), N)
        e = 0
        for v in 1:n
            o = order[v]
            e += s[v] * _ls_permsign(o) * label(x, v, o[1], o[2]) * carry(label(x, v, o[2], o[3]), label(x, v, o[3], o[4]))
        end
        total += cis(2pi * mod(k * e, N) / N)
    end
    total / N^V
end

_ls_near(a, b) = abs(a - b) < 1e-9

"""The same state sum on a Delta-complex: tetrahedra 1..n with vertex labels 0..3; a pairing (s, a, t, b, perm) sends
vertex i of tetrahedron s (i ≠ a) to vertex perm[i+1] of tetrahedron t; orders[t] lists the vertices from first to last
and must form a branching; eps[t] is the sign of orders[t] in the oriented frame."""
function _ls_dw_delta(n, pairings, orders, eps; N::Int=3, k::Int=1)
    keys_ = [(t, a, b) for t in 1:n for a in 0:3 for b in a+1:3]
    idx = Dict(key => i for (i, key) in enumerate(keys_))
    parent = collect(eachindex(keys_)); rel = zeros(Int, length(keys_))
    for (s, a, t, b, perm) in pairings
        f = [i for i in 0:3 if i != a]
        for x in 1:3, y in x+1:3
            i, j = f[x], f[y]; I, J = perm[i+1], perm[j+1]
            (ri, pi_) = _ls_find!(parent, rel, idx[(s, i, j)]); (rj, pj) = _ls_find!(parent, rel, idx[(t, min(I, J), max(I, J))])
            want = I < J ? 0 : 1
            if ri == rj
                @assert mod(pi_ + pj + want, 2) == 0 "inconsistent edge orientation"
            else
                parent[ri] = rj; rel[ri] = mod(pi_ + pj + want, 2)
            end
        end
        rs(i) = findfirst(==(i), orders[s]); rt(i) = findfirst(==(i), orders[t])
        @assert issorted([rt(perm[i+1]) for i in sort(f; by=rs)]) "the vertex orders are not a branching"
    end
    roots = Dict{Int,Int}(); cls = Dict{Tuple{Int,Int,Int},Tuple{Int,Int}}()
    for (i, key) in enumerate(keys_)
        (r, p) = _ls_find!(parent, rel, i)
        haskey(roots, r) || (roots[r] = length(roots) + 1)
        cls[key] = (roots[r], p == 0 ? 1 : -1)
    end
    E = length(roots)
    vp = collect(1:4n); vf(i) = vp[i] == i ? i : (vp[i] = vf(vp[i]))
    for (s, a, t, b, perm) in pairings, i in 0:3
        i == a && continue
        x = vf(4(s - 1) + i + 1); y = vf(4(t - 1) + perm[i+1] + 1); x != y && (vp[x] = y)
    end
    V = length(unique(vf(i) for i in 1:4n))
    coeff(t, a, b) = a < b ? cls[(t, a, b)] : (cls[(t, b, a)][1], -cls[(t, b, a)][2])
    rows = Vector{Int}[]
    for t in 1:n, a in 0:3, b in a+1:3, c in b+1:3
        r = zeros(Int, E)
        for (p, q, sg) in ((a, b, 1), (b, c, 1), (a, c, -1))
            (cl, o) = coeff(t, p, q); r[cl] += sg * o
        end
        push!(rows, r)
    end
    basis = _ls_nullspace_mod(permutedims(reduce(hcat, rows)), N)
    lab(x, t, p, q) = (co = coeff(t, p, q); mod(co[2] * x[co[1]], N))
    carry(b, c) = b + c >= N ? 1 : 0
    total = 0.0 + 0.0im
    for co in Iterators.product(ntuple(_ -> 0:N-1, length(basis))...)
        x = isempty(basis) ? zeros(Int, E) : mod.(sum(co[i] .* basis[i] for i in eachindex(basis)), N)
        e = sum(eps[t] * lab(x, t, orders[t][1], orders[t][2]) *
                carry(lab(x, t, orders[t][2], orders[t][3]), lab(x, t, orders[t][3], orders[t][4])) for t in 1:n)
        total += cis(2pi * mod(k * e, N) / N)
    end
    total / N^V
end

"""Value on the library's oriented triangulation of a presentation (`hc_triangulation`, positive vertex frame, which the
Koda conversion spec §6.1 reads as right-handed), with its recorded branching."""
function _ls_dw_library_frame(link; N::Int=3, k::Int=1)
    tri = Topology.Heegaard.hc_triangulation(link; signature=false)
    @assert all(==(1), tri.signs) "expected the positive vertex frame"
    pairings = [(g.tet + 1, g.omitted, g.to_tet + 1, g.to_omitted, collect(g.permutation)) for g in tri.gluings]
    orders = [collect(o) for o in tri.branch_orders]
    _ls_dw_delta(tri.tetrahedra, pairings, orders, [_ls_permsign(o) for o in orders]; N=N, k=k)
end

"""Value of the lens space S³/⟨(z₁,z₂) ↦ (ζz₁, ζ^q z₂)⟩, ζ = exp(2πi/p), oriented as the boundary of the unit ball of ℂ²
(outward normal first, coordinates (x₁,y₁,x₂,y₂)). Triangulation: the join of two p-gons u_i = (ζ^i, 0), v_j = (0, ζ^j),
tetrahedra [u_i, u_{i+1}, v_j, v_{j+1}] modulo the action; the order u_i < u_{i+1} < v_j < v_{j+1} is a branching."""
function _ls_dw_quotient_lens(p::Int, q::Int; N::Int=3, k::Int=1)
    θ = 2pi / p
    u(i) = [cos(i * θ), sin(i * θ), 0.0, 0.0]; v(j) = [0.0, 0.0, cos(j * θ), sin(j * θ)]
    function sign_of(j)
        P = (u(0), u(1), v(j), v(j + 1)); c = (P[1] + P[2] + P[3] + P[4]) / 4
        M = hcat(c, P[2] - P[1], P[3] - P[1], P[4] - P[1])
        d = M[1,1] * _ls_det3(M[2:4, 2:4]) - M[1,2] * _ls_det3(M[2:4, [1,3,4]]) + M[1,3] * _ls_det3(M[2:4, [1,2,4]]) - M[1,4] * _ls_det3(M[2:4, 1:3])
        d > 0 ? 1 : -1
    end
    rep(i, j) = mod(j - q * i, p) + 1                       # T_{i,j} is the image of T_{0, j-qi}
    pairings = Tuple{Int,Int,Int,Int,Vector{Int}}[]
    for j in 0:p-1
        push!(pairings, (j + 1, 0, rep(1, j), 1, [1, 0, 2, 3]))      # [u_1, v_j, v_{j+1}] of T_{0,j} = face opposite u_2 of T_{1,j}
        push!(pairings, (j + 1, 2, rep(0, j + 1), 3, [0, 1, 3, 2]))  # [u_0, u_1, v_{j+1}] of T_{0,j} = face opposite v_{j+2} of T_{0,j+1}
    end
    _ls_dw_delta(p, pairings, [[0, 1, 2, 3] for _ in 1:p], [sign_of(j) for j in 0:p-1]; N=N, k=k)
end
_ls_det3(A) = A[1,1] * (A[2,2] * A[3,3] - A[2,3] * A[3,2]) - A[1,2] * (A[2,1] * A[3,3] - A[2,3] * A[3,1]) + A[1,3] * (A[2,1] * A[3,2] - A[2,2] * A[3,1])

@testset "lens space orientation (LensSpaceNames)" begin
    registered = ("s3_abalone", "lens_2_1", "lens_3_1_mst2_example23", "lens_3_2", "poincare", "quaternion",
                  "s2xs1", "s2xs1_2vertex", "odata_variant")
    family = [(p, q, p == 7 ? 7 : 3) for p in (3, 6, 7) for q in 1:p-1 if gcd(p, q) == 1]

    @testset "MST25 Proposition A.1: Γ× and →(Γ*) conjugate the value" begin
        for name in registered
            link = vl_example(name); z = _ls_dw(link)
            @test _ls_near(_ls_dw(mst_sign_change(link)), conj(z))
            @test _ls_near(_ls_dw(mst_reflect_reverse(link)), conj(z))
        end
        for (p, q, N) in family
            link = lens_space_link(p, q); z = _ls_dw(link; N=N)
            @test abs(imag(z)) > 0.1                                   # the value sees the orientation
            @test _ls_near(_ls_dw(mst_sign_change(link); N=N), conj(z))
            @test _ls_near(_ls_dw(mst_reflect_reverse(link); N=N), conj(z))
        end
    end

    @testset "lens family: L(p,q) ≅ L(p,q⁻¹) and L(p,-q) = L(p,q)^op" begin
        for (p, q, N) in family
            z = _ls_dw(lens_space_link(p, q); N=N)
            @test _ls_near(z, _ls_dw(lens_space_link(p, invmod(q, p)); N=N))
            @test _ls_near(z, conj(_ls_dw(lens_space_link(p, p - q); N=N)))
        end
    end

    @testset "registered examples against the family" begin
        z23 = _ls_dw(vl_example("lens_3_1_mst2_example23"))
        @test _ls_near(z23, conj(_ls_dw(lens_space_link(3, 1))))
        @test _ls_near(z23, _ls_dw(lens_space_link(3, 2)))
        @test _ls_near(_ls_dw(vl_example("lens_3_2")), _ls_dw(lens_space_link(3, 2)))
    end

    # Absolute orientation. Depends on the Koda conversion spec §6.1 (the positive vertex frame of hc_triangulation is
    # right-handed); that row is fixed by the HeegaardCodes conventions (rule a), so these consequences hold in the library frame.
    @testset "o-data against the standard quotient S³/⟨(ζz₁, ζ^q z₂)⟩ (Koda §6.1 frame)" begin
        for (p, q, N) in family
            link = lens_space_link(p, q)
            zl = _ls_dw_library_frame(link; N=N)
            @test _ls_near(zl, conj(_ls_dw(link; N=N)))                     # the over-pair rule is the opposite frame
            @test _ls_near(zl, _ls_dw_quotient_lens(p, p - q; N=N))         # lens_space_link(p,q) presents L(p,-q)
            @test !_ls_near(zl, _ls_dw_quotient_lens(p, q; N=N))
        end
        z23 = _ls_dw_library_frame(vl_example("lens_3_1_mst2_example23"))
        @test _ls_near(z23, _ls_dw_quotient_lens(3, 1))                   # MST22's example presents L(3,1)
        @test _ls_near(_ls_dw_quotient_lens(3, 1), (1 + 2cis(2pi / 3)) / 3)   # closed form (1/N) Σ_b ζ^{q b²} for p = N
    end
end
