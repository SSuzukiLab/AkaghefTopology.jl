# cf: lib/topology.jl/docs/specs/lickorish-code.md §3 — passage sign; −k on b and +k on a means the frame (b, a) is positive
# cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — intersection rows are α, columns are β; H₁ is the abelianization of the Heegaard presentation
# cf: threads/PJ260922/code/hc_matrix.py — (α, β) entry is +1 when + is on α and − is on β
# cf: lib/topology.jl/src/heegaard/groups.jl — smith_normal_form, abelianization_matrix, first_homology(::GroupPresentation); α generates, β is the relator
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3.1 — face numbers come from the rotation system, not from stored label integers

"""計算と比較。入口の引数は HeegaardLayer だけである。"""

struct HCClosedDiagramExcluded <: Exception
    condition::Symbol
    witness::Dict{String,Any}
end

function Base.showerror(io::IO, err::HCClosedDiagramExcluded)
    print(io, "HCClosedDiagramExcluded(", err.condition, ")")
end

struct HCNotComputed <: Exception
    reason::Symbol
    detail::Dict{String,Any}
end

function Base.showerror(io::IO, err::HCNotComputed)
    print(io, "HCNotComputed(", err.reason, ")")
end

struct HCIntersection
    rows::Vector{String}
    columns::Vector{String}
    matrix::Matrix{Int}
    row_order::Symbol
    column_order::Symbol
    sign_convention::String
end

struct HCInvariantFactors
    factors::Vector{Int}
    source::Symbol
    manifold_homology::Bool
end

struct HCPresentationResult
    presentation::GroupPresentation
    generator_family::Symbol
    relator_family::Symbol
    generator_names::Vector{String}
    relator_names::Vector{String}
end

struct HCFirstHomology
    free_rank::Int
    torsion::Vector{Int}
    invariant_factors::Vector{Int}
    generator_family::Symbol
    relator_family::Symbol
    condition::Symbol
end

struct HCCanonical
    genus::Int
    k_words::Vector{Pair{String,Vector{Int}}}
    alpha_words::Vector{Vector{Int}}
    beta_words::Vector{Vector{Int}}
    other_words::Vector{Pair{String,Vector{Int}}}
    places::Vector{Pair{String,String}}
    lc_disks::Vector{Pair{String,Vector{String}}}
    hc_disks::Vector{Pair{String,Vector{String}}}
    lc_points::Vector{Pair{String,String}}
    hc_points::Vector{Pair{String,String}}
    lc_weights::Vector{Pair{String,Int}}
    hc_weights::Vector{Pair{String,Int}}
end

function Base.:(==)(a::HCCanonical, b::HCCanonical)
    all(getfield(a, i) == getfield(b, i) for i in 1:fieldcount(HCCanonical))
end

function Base.hash(canonical::HCCanonical, h::UInt)
    for i in 1:fieldcount(HCCanonical)
        h = hash(getfield(canonical, i), h)
    end
    hash(HCCanonical, h)
end

struct HCMoveStep
    kind::Symbol
    detail::Dict{String,Any}
end

struct HCMoveSearchResult
    steps::Vector{HCMoveStep}
end

struct HCVerdict
    level::Symbol
    status::Symbol
    reason::Symbol
    evidence::Dict{String,Any}
end

struct HCComparison
    identity::HCVerdict
    isotopy::HCVerdict
    move_equivalence::HCVerdict
    manifold::HCVerdict
end

struct HCNamedRelations
    identical::Bool
    orientation_reversal::Symbol
    mirror::Symbol
end

const _HC_SIGN = "+1 when + is on alpha and - is on beta; -1 for the opposite; labels that are not alpha-beta do not contribute"
const _HC_MOVE_KINDS = (:isotopy, :handle_slide, :stabilization, :destabilization)
const _HC_TIE_LIMIT = 2_000_000

function _hc_k_set(genus::Integer)
    Set(_lc_k_names(genus))
end

function _hc_owners(words::AbstractDict)
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in words
        for x in word
            x == 0 && throw(KodaInputError(:hc_compare, "label 0 is not a passage", Dict("curve" => String(name))))
            lab = abs(x)
            dest = x > 0 ? plus : minus
            haskey(dest, lab) && throw(KodaInputError(
                :hc_compare, "passage sign is not unique", Dict("label" => lab, "curve" => String(name))))
            dest[lab] = String(name)
        end
    end
    owners = Dict{Int,Tuple{String,String}}()
    for lab in keys(plus)
        haskey(minus, lab) || throw(KodaInputError(
            :hc_compare, "label is missing a minus passage", Dict("label" => lab)))
        owners[lab] = (plus[lab], minus[lab])
    end
    for lab in keys(minus)
        haskey(plus, lab) || throw(KodaInputError(
            :hc_compare, "label is missing a plus passage", Dict("label" => lab)))
    end
    owners
end

function _hc_curve_lists(hc::HeegaardLayer)
    genus = hc.lc.genus
    kn = _hc_k_set(genus)
    alphas = String[]
    betas = String[]
    for (name, role) in hc.roles
        String(name) in kn && continue
        haskey(hc.lc.words, name) || throw(KodaInputError(
            :hc_compare, "role names a missing component", Dict("curve" => String(name))))
        if role == :alpha
            push!(alphas, String(name))
        elseif role == :beta
            push!(betas, String(name))
        else
            throw(KodaInputError(:hc_compare, "invalid role", Dict("curve" => String(name), "role" => String(role))))
        end
    end
    others = String[]
    for name in keys(hc.lc.words)
        sname = String(name)
        sname in kn && continue
        haskey(hc.roles, name) && continue
        push!(others, sname)
    end
    sort!(alphas)
    sort!(betas)
    sort!(others)
    alphas, betas, others
end

function _hc_index(names::Vector{String})
    Dict(name => i for (i, name) in enumerate(names))
end

"""
α 行、β 列の交点行列を返す。並びは成分名の辞書順で、結果の rows と columns がそのまま行と列である。

符号は passage sign から読む。+ が α にあり − が β にある交点は +1、逆は −1。
これは曲面上の枠 (β, α) が正のとき +1 になる hc_matrix.py の規約であり、接ベクトルの順 (Tα, Tβ) そのものではない。
K どうし、同族どうしの交点は足さない。種数が曲線本数と違っても行列は返す。これは多様体のホモロジーではない。
"""
function hc_intersection(hc::HeegaardLayer)
    alphas, betas, _ = _hc_curve_lists(hc)
    owners = _hc_owners(hc.lc.words)
    ai = _hc_index(alphas)
    bi = _hc_index(betas)
    matrix = zeros(Int, length(alphas), length(betas))
    for (_, (plus_c, minus_c)) in owners
        if haskey(ai, plus_c) && haskey(bi, minus_c)
            matrix[ai[plus_c], bi[minus_c]] += 1
        elseif haskey(ai, minus_c) && haskey(bi, plus_c)
            matrix[ai[minus_c], bi[plus_c]] -= 1
        end
    end
    HCIntersection(alphas, betas, matrix, :sorted_component_name, :sorted_component_name, _HC_SIGN)
end

"""
交点行列の Smith 不変因子。1 を含む正の対角成分で、0 は除く。多様体のホモロジーではない。

抽出は groups.jl の first_homology と同じ式で、対角化そのものは smith_normal_form を使う。
行列とその転置の不変因子は一致する。
"""
function hc_invariant_factors(hc::HeegaardLayer)
    matrix = hc_intersection(hc).matrix
    smith = smith_normal_form(matrix)
    diagonal = [abs(smith[i, i]) for i in 1:min(size(smith)...)]
    HCInvariantFactors(filter(!iszero, diagonal), :intersection_matrix, false)
end

function _hc_relator(word::Vector{Int}, owners, ai)
    letters = Tuple{Int,Int}[]
    for x in word
        plus_c, minus_c = owners[abs(x)]
        if haskey(ai, plus_c) && !haskey(ai, minus_c)
            push!(letters, (ai[plus_c], 1))
        elseif haskey(ai, minus_c) && !haskey(ai, plus_c)
            push!(letters, (ai[minus_c], -1))
        end
    end
    letters
end

"""
π₁ の表示を返す。生成元は α、関係式は β である。

β の語を保存された順に歩き、α との交点だけを文字にする。指数は交点行列と同じで、
+ が α なら +1、+ が β なら −1。可換化は交点行列の転置になる。
交点のない β は空語の関係式、交点のない α は関係式に現れない生成元、種数 0 で曲線が無いときは自明群である。
語の始点は保存された語の先頭であり、巡回で取り替えた表示とはタプルとして違う。群としては共役で同値である。
"""
function hc_presentation(hc::HeegaardLayer)
    alphas, betas, _ = _hc_curve_lists(hc)
    owners = isempty(hc.lc.words) ? Dict{Int,Tuple{String,String}}() : _hc_owners(hc.lc.words)
    ai = _hc_index(alphas)
    relators = [_hc_relator(hc.lc.words[name], owners, ai) for name in betas]
    meta = Dict{String,Any}(
        "generator_family" => "alpha",
        "relator_family" => "beta",
        "generator_names" => alphas,
        "relator_names" => betas,
        "generator_order" => "sorted_component_name",
        "relator_order" => "sorted_component_name",
        "sign_convention" => _HC_SIGN,
        "beta_ids" => "positions in relator_names, not HeegaardCode curve ids",
    )
    presentation = GroupPresentation(alphas, relators, collect(1:length(betas)), meta)
    HCPresentationResult(presentation, :alpha, :beta, alphas, betas)
end

"""
閉じた図式の H₁ を返す。α の本数と β の本数がどちらも種数と一致しないときは
HCClosedDiagramExcluded(:alpha_beta_count_ne_genus) を投げる。

返す群は、α を生成元、β を関係式とした表示の可換化である。完全切断系であることはここでは調べない。
本数が種数と一致することだけを、閉じた 3 次元多様体の図式として H₁ を返す条件にする。
"""
function hc_first_homology(hc::HeegaardLayer)
    alphas, betas, _ = _hc_curve_lists(hc)
    genus = hc.lc.genus
    if length(alphas) != genus || length(betas) != genus
        throw(HCClosedDiagramExcluded(
            :alpha_beta_count_ne_genus,
            Dict{String,Any}(
                "genus" => genus,
                "n_alpha" => length(alphas),
                "n_beta" => length(betas),
                "alpha_names" => alphas,
                "beta_names" => betas,
            ),
        ))
    end
    presented = hc_presentation(hc)
    homology = first_homology(presented.presentation)
    HCFirstHomology(
        homology.free_rank,
        homology.torsion,
        homology.invariant_factors,
        :alpha,
        :beta,
        :alpha_beta_count_eq_genus,
    )
end

function _hc_s3_tables()
    perms = [(1, 2, 3), (2, 1, 3), (1, 3, 2), (3, 2, 1), (2, 3, 1), (3, 1, 2)]
    mul = zeros(Int, 6, 6)
    inv = zeros(Int, 6)
    for i in 1:6, j in 1:6
        p, q = perms[i], perms[j]
        composed = ntuple(k -> q[p[k]], 3)
        mul[i, j] = findfirst(==(composed), perms)
    end
    for i in 1:6
        p = perms[i]
        inverse = ntuple(k -> findfirst(==(k), p), 3)
        inv[i] = findfirst(==(inverse), perms)
    end
    mul, inv
end

const _HC_S3_MUL, _HC_S3_INV = _hc_s3_tables()

function _hc_eval_s3(relator, assign)
    acc = 1
    for (gen, exp) in relator
        step = exp >= 0 ? assign[gen] : _HC_S3_INV[assign[gen]]
        for _ in 1:abs(exp)
            acc = _HC_S3_MUL[acc, step]
        end
    end
    acc
end

function _hc_count_s3(relators, n::Int)
    n == 0 && return BigInt(1)
    used = falses(n)
    for relator in relators
        for (gen, _) in relator
            1 <= gen <= n || throw(KodaInputError(
                :hc_hom_count, "relator names an unknown generator", Dict("generator" => gen)))
            used[gen] = true
        end
    end
    active = [i for i in 1:n if used[i]]
    remap = Dict(old => i for (i, old) in enumerate(active))
    mapped = Vector{Tuple{Int,Int}}[]
    for relator in relators
        isempty(relator) && continue
        push!(mapped, [(remap[gen], exp) for (gen, exp) in relator])
    end
    free = n - length(active)
    length(active) == 0 && return BigInt(6)^free
    assign = zeros(Int, length(active))
    nodes = Ref(0)
    count = Ref(0)
    function rec(i::Int)
        if i > length(active)
            count[] += 1
            return
        end
        for g in 1:6
            nodes[] += 1
            nodes[] > _HC_TIE_LIMIT && throw(HCNotComputed(
                :hom_search_budget,
                Dict{String,Any}("group" => "S3", "generators" => n, "budget" => _HC_TIE_LIMIT),
            ))
            assign[i] = g
            ok = true
            for relator in mapped
                any(letter -> letter[1] > i, relator) && continue
                if _hc_eval_s3(relator, assign) != 1
                    ok = false
                    break
                end
            end
            ok && rec(i + 1)
        end
    end
    rec(1)
    count[] * BigInt(6)^free
end

"""
|Hom(π₁, G)| を返す。G は :Z2, :Z3, :Z5, :S3。表示は hc_presentation と同じで、α が生成元、β が関係式。

巡回群への準同型は可換化を通るので、first_homology の自由階数とねじれから数える。S₃ は語をそのまま評価する。
これは表示された群の準同型の個数であり、本数が種数と違う図式でも返す。多様体の H₁ ではない。
"""
function hc_hom_count(presented::HCPresentationResult, group::Symbol)
    presentation = presented.presentation
    if group in (:Z2, :Z3, :Z5)
        modulus = group == :Z2 ? 2 : group == :Z3 ? 3 : 5
        homology = first_homology(presentation)
        exponent = homology.free_rank + count(d -> d % modulus == 0, homology.torsion)
        return BigInt(modulus)^exponent
    elseif group == :S3
        return _hc_count_s3(presentation.relators, length(presentation.generators))
    else
        throw(ArgumentError("unsupported hom target: $group"))
    end
end

hc_hom_count(hc::HeegaardLayer, group::Symbol) = hc_hom_count(hc_presentation(hc), group)

function _hc_permutations(items::Vector{String})
    n = length(items)
    n == 0 && return [String[]]
    n > 8 && throw(HCNotComputed(:canonical_tie_group, Dict{String,Any}("n" => n)))
    out = Vector{String}[]
    used = falses(n)
    current = Vector{String}(undef, n)
    function rec(depth::Int)
        if depth > n
            push!(out, copy(current))
            return
        end
        for i in 1:n
            used[i] && continue
            used[i] = true
            current[depth] = items[i]
            rec(depth + 1)
            used[i] = false
        end
    end
    rec(1)
    out
end

function _hc_emit_rotation(word::Vector{Int}, start::Int, label_map::Dict{Int,Int}, next_id::Int)
    n = length(word)
    n == 0 && return Int[], label_map, next_id
    local_map = copy(label_map)
    nid = next_id
    emitted = Vector{Int}(undef, n)
    for i in 1:n
        x = word[mod1(start + i - 1, n)]
        lab = abs(x)
        if !haskey(local_map, lab)
            local_map[lab] = nid
            nid += 1
        end
        emitted[i] = sign(x) * local_map[lab]
    end
    emitted, local_map, nid
end

include(joinpath(@__DIR__, "hc_compare_form.jl"))
include(joinpath(@__DIR__, "hc_compare_canonical.jl"))


hc_identical(a::HeegaardLayer, b::HeegaardLayer) = hc_canonical(a) == hc_canonical(b)

function _hc_replace_words(hc::HeegaardLayer, words::Dict{String,Vector{Int}})
    lc = hc.lc
    new_lc = LickorishConfiguration(lc.genus, words, copy(lc.places), copy(lc.disks), copy(lc.points), copy(lc.weights), nothing)
    HeegaardLayer(new_lc, copy(hc.roles), copy(hc.disks), copy(hc.points), copy(hc.weights), nothing)
end

"""
指定した曲線の向きを取り替える。語を逆順にし、その曲線上の交点ラベルの両側の符号を反転する。

両側を反転するのは、passage sign の + 一つと − 一つを保つためである。片方だけ反転すると語が LC の形を離れる。
向きを取り替えた曲線どうしの交点は、それぞれの曲線で一度だけ反転する。
"""
function hc_reverse_orientation(hc::HeegaardLayer, curves)
    names = Set(String.(curves))
    for name in names
        haskey(hc.lc.words, name) || throw(KodaInputError(
            :hc_reverse_orientation, "missing curve", Dict("curve" => name)))
    end
    labels = Set{Int}()
    for name in names
        for x in hc.lc.words[name]
            push!(labels, abs(x))
        end
    end
    words = Dict{String,Vector{Int}}()
    for (name, word) in hc.lc.words
        if String(name) in names
            words[String(name)] = Int[-x for x in reverse(word)]
        else
            words[String(name)] = Int[abs(x) in labels ? -x : x for x in word]
        end
    end
    _hc_replace_words(hc, words)
end

"""鏡像。全部の passage sign を反転し、語の順序は逆にしない。曲面の向きを反転した枠の取り替えである。"""
function hc_mirror(hc::HeegaardLayer)
    words = Dict{String,Vector{Int}}(String(name) => Int[-x for x in word] for (name, word) in hc.lc.words)
    _hc_replace_words(hc, words)
end

function _hc_orientation_related(a::HeegaardLayer, b::HeegaardLayer)
    names = sort!(collect(String.(keys(a.lc.words))))
    length(names) > 12 && return :undetermined
    n = length(names)
    for mask in 1:(1 << n) - 1
        chosen = String[names[i] for i in 1:n if (mask >> (i - 1)) & 1 == 1]
        hc_identical(hc_reverse_orientation(a, chosen), b) && return :yes
    end
    :no
end

"""同一性、向きの取り替え、鏡像を別々に判定する。向きの探索は曲線が 12 本を超えると :undetermined である。"""
function hc_named_relations(a::HeegaardLayer, b::HeegaardLayer)
    HCNamedRelations(hc_identical(a, b), _hc_orientation_related(a, b), hc_identical(hc_mirror(a), b) ? :yes : :no)
end

function _hc_move_levels(move_result::HCMoveSearchResult)
    isempty(move_result.steps) && return Symbol[]
    kinds = [step.kind for step in move_result.steps]
    all(kind -> kind in _HC_MOVE_KINDS, kinds) || return Symbol[]
    levels = Symbol[:move_equivalence, :manifold]
    all(kind -> kind == :isotopy, kinds) && push!(levels, :isotopy)
    levels
end

_hc_move_levels(::Nothing) = Symbol[]

function _hc_invariant_snapshot(hc::HeegaardLayer)
    presented = hc_presentation(hc)
    hom = Dict{Symbol,BigInt}()
    hom_error = Dict{Symbol,Symbol}()
    for group in (:Z2, :Z3, :Z5, :S3)
        try
            hom[group] = hc_hom_count(presented, group)
        catch err
            err isa HCNotComputed || rethrow()
            hom_error[group] = err.reason
        end
    end
    homology = nothing
    homology_error = nothing
    try
        homology = hc_first_homology(hc)
    catch err
        err isa HCClosedDiagramExcluded || rethrow()
        homology_error = err.condition
    end
    (; hom, hom_error, homology, homology_error)
end

function _hc_invariants_differ(left, right)
    for group in (:Z2, :Z3, :Z5, :S3)
        haskey(left.hom, group) && haskey(right.hom, group) && left.hom[group] != right.hom[group] &&
            return true, Symbol(:hom_, group)
    end
    if left.homology !== nothing && right.homology !== nothing
        if left.homology.free_rank != right.homology.free_rank || left.homology.torsion != right.homology.torsion
            return true, :h1
        end
    end
    false, :none
end

function _hc_agreeing(left, right)
    agreeing = String[]
    if left.homology !== nothing && right.homology !== nothing &&
            left.homology.free_rank == right.homology.free_rank &&
            left.homology.torsion == right.homology.torsion
        push!(agreeing, "h1")
    end
    for group in (:Z2, :Z3, :Z5, :S3)
        haskey(left.hom, group) && haskey(right.hom, group) && left.hom[group] == right.hom[group] &&
            push!(agreeing, string("hom_", group))
    end
    agreeing
end

function _hc_invariant_evidence(left, right)
    evidence = Dict{String,Any}(
        "agreeing_invariants" => _hc_agreeing(left, right),
        "hom" => Dict(string(group) => haskey(left.hom, group) ? left.hom[group] : "not_computed" for group in (:Z2, :Z3, :Z5, :S3)),
    )
    if left.homology !== nothing
        evidence["h1"] = Dict("free_rank" => left.homology.free_rank, "torsion" => left.homology.torsion)
    else
        evidence["h1_condition"] = left.homology_error
    end
    evidence
end

function _hc_equal(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :equal, reason, evidence)
end

function _hc_different(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :different, reason, evidence)
end

function _hc_undetermined(level::Symbol, reason::Symbol, evidence=Dict{String,Any}())
    HCVerdict(level, :undetermined, reason, evidence)
end

"""
4 つの水準を返す。各水準の status は :equal、:different、:undetermined のいずれか。

同一性は hc_canonical で決める。同位と移動同値は、同一なら :equal、不変量が違えば :different、
それ以外は :undetermined。move_result は探索結果を受け取る形だけで、この関数は move を探さない。
同一性が決まらないときは HCNotComputed を投げ、:different にはしない。
表す多様体は、不変量が一致しても :equal にしない。:equal は同一性か、与えられた move の列だけを根拠にする。
"""
function compare_hc(a::HeegaardLayer, b::HeegaardLayer; move_result::Union{Nothing,HCMoveSearchResult}=nothing)
    identical = hc_identical(a, b)
    if identical
        evidence = Dict{String,Any}("basis" => "canonical_form")
        return HCComparison(
            _hc_equal(:identity, :canonical_form, evidence),
            _hc_equal(:isotopy, :identical, evidence),
            _hc_equal(:move_equivalence, :identical, evidence),
            _hc_equal(:manifold, :identical, evidence),
        )
    end
    identity = _hc_different(:identity, :canonical_form)
    left = _hc_invariant_snapshot(a)
    right = _hc_invariant_snapshot(b)
    differ, differ_reason = _hc_invariants_differ(left, right)
    witnessed = _hc_move_levels(move_result)
    move_evidence = Dict{String,Any}(
        "steps" => move_result === nothing ? Any[] : [Dict("kind" => String(step.kind)) for step in move_result.steps],
    )
    function coarse(level::Symbol)
        differ && return _hc_different(level, differ_reason, _hc_invariant_evidence(left, right))
        level in witnessed && return _hc_equal(level, :move_sequence, move_evidence)
        _hc_undetermined(level, :no_move_search, _hc_invariant_evidence(left, right))
    end
    manifold = if differ
        _hc_different(:manifold, differ_reason, _hc_invariant_evidence(left, right))
    elseif :manifold in witnessed
        _hc_equal(:manifold, :move_sequence, move_evidence)
    else
        _hc_undetermined(:manifold, :invariants_agree_not_identity, _hc_invariant_evidence(left, right))
    end
    HCComparison(identity, coarse(:isotopy), coarse(:move_equivalence), manifold)
end

export HCClosedDiagramExcluded, HCNotComputed,
    HCIntersection, HCInvariantFactors, HCPresentationResult, HCFirstHomology, HCCanonical,
    HCCanonicalWitness, HCIdentityCorrespondence,
    HCMoveStep, HCMoveSearchResult, HCVerdict, HCComparison, HCNamedRelations,
    hc_intersection, hc_invariant_factors, hc_presentation, hc_first_homology, hc_hom_count,
    hc_canonical, hc_canonical_witness, hc_identity_correspondence, hc_apply_correspondence,
    hc_identical, hc_reverse_orientation, hc_mirror, hc_named_relations, compare_hc
