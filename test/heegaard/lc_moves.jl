# cf: lib/topology.jl/docs/specs/heegaard-moves.md §8 — inverse, checker, genus, cyclic wrap
using Test
using Topology
using SHA

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const CHECKER = joinpath(ROOT, "threads", "PJ260922", "code", "lc_check.py")

@assert isdefined(H, :load_heegaard_moves) "Topology must load the move API through its package entrypoint"

function parse_lc(text)
    H._lc_parse_text(text)
end

function checker(lc)
    dir = mktempdir()
    path = joinpath(dir, "probe.lc")
    H.write_lc(path, lc)
    out = read(Cmd(`uv run --no-sync python $CHECKER $path`; ignorestatus=true), String)
    occursin("\"status\": \"pass\"", out) || occursin("\"status\":\"pass\"", out), out
end

function cyclic_equal(left::Vector{Int}, right::Vector{Int})
    length(left) == length(right) || return false
    isempty(left) && return true
    any(i -> left == vcat(right[i:end], right[1:i - 1]), 1:length(right))
end

function relabel(word, mapping)
    [sign(entry) * mapping[abs(entry)] for entry in word]
end

function same_up_to(left, right)
    H._component_names(left) == H._component_names(right) || return false
    labels = sort(unique(abs.(reduce(vcat, values(left.words); init=Int[]))))
    right_labels = sort(unique(abs.(reduce(vcat, values(right.words); init=Int[]))))
    length(labels) == length(right_labels) || return false
    # K labels are pinned. Remaining labels are matched by trying the identity first,
    # then the bijection that sends left's non-K labels onto right's in sorted order
    # only when both sides use the same crossing pattern. The composed move maps are
    # checked by the caller; this predicate searches permutations of the non-K labels.
    limit = left.genus <= 0 ? 0 : 3 * left.genus - 2
    fixed = [label for label in labels if label <= limit]
    all(label in right_labels for label in fixed) || return false
    free_left = [label for label in labels if label > limit]
    free_right = [label for label in right_labels if label > limit]
    length(free_left) == length(free_right) || return false
    # Propagate labels along named cyclic components. No factorial label search
    # and no size cutoff: crossings shared with the next word constrain its turn.
    names = sort(collect(keys(left.words)); by=n -> length(left.words[n]))
    function visit(i, mapping, used)
        i > length(names) && return true
        a, b = left.words[names[i]], right.words[names[i]]
        length(a) == length(b) || return false
        isempty(a) && return visit(i+1,mapping,used)
        for shift in 0:length(b)-1
            next, occupied = copy(mapping), copy(used)
            ok = true
            for j in eachindex(a)
                x,y = a[j],b[mod1(j+shift,length(b))]
                if sign(x) != sign(y) || (haskey(next,abs(x)) ? next[abs(x)] != abs(y) : abs(y) in occupied)
                    ok=false; break
                end
                next[abs(x)]=abs(y); push!(occupied,abs(y))
            end
            ok && visit(i+1,next,occupied) && return true
        end
        false
    end
    left.genus == right.genus && visit(1,Dict(label=>label for label in fixed),Set(fixed))
end

function permutations(values)
    isempty(values) && return [Int[]]
    out = Vector{Vector{Int}}()
    for i in eachindex(values)
        rest = vcat(values[1:i - 1], values[i + 1:end])
        for tail in permutations(rest)
            push!(out, vcat(values[i], tail))
        end
    end
    out
end

function algebraic_intersections(lc, alphas, betas)
    words = lc.words
    owners = Dict{Int,Tuple{String,Int}}()
    for (name, word) in words
        for entry in word
            owners[sign(entry) * abs(entry) == entry ? abs(entry) + 0 : abs(entry)] = (name, sign(entry))
        end
    end
    # rebuild properly
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in words
        for entry in word
            entry > 0 ? (plus[abs(entry)] = name) : (minus[abs(entry)] = name)
        end
    end
    matrix = zeros(Int, length(alphas), length(betas))
    for (i, alpha) in enumerate(alphas), (j, beta) in enumerate(betas)
        for label in keys(plus)
            pair = (plus[label], minus[label])
            if pair == (alpha, beta)
                matrix[i, j] += 1
            elseif pair == (beta, alpha)
                matrix[i, j] -= 1
            end
        end
    end
    matrix
end

const DIGON = """
genus: 1
l1 : -1 2 -3
m1 : 1
c : -2 3
"""

const DIGON_WRAP = """
genus: 1
l1 : -1 2 -3
m1 : 1 5
c : 3 -5 -2
"""

const TRIANGLE = """
genus: 1
l1 : 3 -1
m1 : 2 1
c : -3 -2
"""

@testset "lc moves" begin
@testset "catalog and isotopy" begin
    ids = H.heegaard_move_ids()
    @test "lc.isotopy.digon.remove.00" in ids
    @test "lc.twist.annulus" in ids
    @test "hc.slide.pants" in ids
    @test "hc.stabilization.once" in ids

    digon = parse_lc(DIGON)
    hits = H.movable(digon, "lc.isotopy.digon.remove.00")
    @test length(hits) == 1
    forward = H.move(digon, only(hits))
    ok, report = checker(forward.output)
    @test ok
    @test forward.output.genus == digon.genus
    back = H.move(forward.output, forward.inverse)
    @test back.output.genus == digon.genus
    @test same_up_to(digon, back.output)

    wrapped = parse_lc(DIGON_WRAP)
    whits = H.movable(wrapped, "lc.isotopy.digon.remove.00")
    @test length(whits) == 1
    wforward = H.move(wrapped, only(whits))
    wok, wreport = checker(wforward.output)
    @test wok
    wback = H.move(wforward.output, wforward.inverse)
    @test same_up_to(wrapped, wback.output)

    tri = parse_lc(TRIANGLE)
    thits = H.movable(tri, "lc.isotopy.vertex.000")
    @test length(thits) == 1
    tforward = H.move(tri, only(thits))
    tok, treport = checker(tforward.output)
    @test tok
    tback = H.move(tforward.output, tforward.inverse)
    @test same_up_to(tri, tback.output)
    @test tforward.output.genus == 1
end

@testset "marks and unresolved" begin
    digon = parse_lc(DIGON)
    marked = Topology.Heegaard.LickorishCode(
        digon.genus, digon.words, digon.places, Dict("D" => ["F1"]), digon.points, digon.weights, nothing)
    # Owner decision 2026-11-12: free marks no longer gate geometric moves.
    mr=H.move(marked,first(H.movable(marked, "lc.isotopy.digon.remove.00")))
    @test checker(mr.output)[1]
    @test all(r in H._move_mark_cells(mr.output).faces for r in mr.output.disks["D"])
    @test_throws H.LCMoveError H.movable(digon, "lc.fixed.unresolved")
end

@testset "digon and parallel corpus properties" begin
    files = ["code/lc-fixtures/operations/g1-base.lc", "code/lc-fixtures/p28-after.lc",
        "data/261107-lickorish-data/group-18/case-01.lc",
        "data/261107-ha9-curve-figures/p0254-i174-pair.lc"]
    for rel in files
        lc=H.read_lc(joinpath(ROOT,"threads","PJ260922",rel))
        names=H._component_names(lc)
        original=algebraic_intersections(lc,names,names)
        hits=H.movable(lc,"lc.isotopy.digon.create.ab")
        @test !isempty(hits)
        for w in (length(hits)<=24 ? hits : hits[unique([1,cld(length(hits),2),length(hits)])])
            @test sum(sign.(w.after_frags[1])) == 0
            r=H.move(lc,w)
            @test checker(r.output)[1]
            @test r.output.genus == lc.genus
            @test algebraic_intersections(r.output,names,names) == original
            back=H.move(r.output,r.inverse)
            @test checker(back.output)[1]
            @test same_up_to(lc,back.output)
        end
        for w in H.movable(lc,"lc.isotopy.parallel")
            r=H.move(lc,w)
            @test checker(r.output)[1]
            @test length(r.output.words) == length(lc.words)+1
            @test algebraic_intersections(r.output,names,names) == original
            c=w.extra["curve"]; copy_name=only(w.extra["created_components"])
            @test algebraic_intersections(r.output,[copy_name],names) == algebraic_intersections(lc,[c],names)
            @test same_up_to(lc,H.move(r.output,r.inverse).output)
        end
    end
end

@testset "shared support and embedded path" begin
    lc=parse_lc("genus: 1\nl1 : -1 -2\nm1 : +1\nb : +2\n")
    hc=H.resolve_hc(lc,Dict("b"=>:beta);standard=true)
    shared=only(H.movable(hc,"lc.twist.annulus";q=1,support="b",acted_components=["m1"]))
    result=H.move(hc,shared)
    @test haskey(shared.extra,"sequence")
    @test checker(result.output.lc)[1]
    @test same_up_to(lc,H.move(result.output,result.inverse).output.lc)
    base=H.read_lc(joinpath(ROOT,"threads","PJ260922","code","lc-fixtures","operations","g1-base.lc"))
    active=H._lc_nonempty_words(base.words)
    faces,right=H._lc_faces(active); occ=H._lc_rotation_occ(active)
    left=Dict(H._lc_tau(active,occ,c)=>id for (c,id) in right)
    face=first(faces)
    boundary=H._move_left_boundary(active,occ,face)
    p=H.lc_move_path(base;start=boundary[1].corner,finish=boundary[2].corner)
    @test p["faces"] == [face.id]
    @test isempty(p["crossings"])
    c=boundary[1].corner; other=H._lc_tau(active,occ,c)
    destination=only(f for f in faces if f.id==left[other])
    finish=first(e.corner for e in H._move_left_boundary(active,occ,destination) if e.corner!=other)
    p=H.lc_move_path(base;start=boundary[2].corner,crossings=[c],finish=finish)
    @test length(p["faces"]) == 2
    @test_throws H.LCMoveError H.lc_move_path(base;start=boundary[2].corner,crossings=[c],finish=finish,avoid=[boundary[1].component])
    @test_throws H.LCMoveError H.lc_move_path(base;start=c,finish=c)
end

# Spec §6.2: F_q(s,t)=(s+q*t mod 1, t), support a along s. q=2 sends b to 2a+b,
# so <a,b> stays ±1. b meets m1, which crosses a, once per turn. The §8
# determinant-2 example needs support c=a+b and q=1, tested separately below.
@testset "partial twist determinant" begin
    torus = parse_lc("""
    genus: 1
    l1 : -2 -1
    m1 : 1 3
    a : -4 -3
    b : 2 4
    """)
    before = algebraic_intersections(torus, ["a"], ["b"])[1, 1]
    @test abs(before) == 1
    @test algebraic_intersections(torus, ["m1"], ["b"])[1, 1] == 0
    hits = H.movable(torus, "lc.twist.annulus"; q=2, support="a", acted_components=["b"])
    @test !isempty(hits)
    result = H.move(torus, only(hits))
    ok, report = checker(result.output)
    @test ok
    after = algebraic_intersections(result.output, ["a"], ["b"])[1, 1]
    @test abs(after) == 1
    transverse = algebraic_intersections(result.output, ["m1"], ["b"])[1, 1]
    @test abs(transverse) == 2
    @test result.output.genus == torus.genus
end

# Independent flat-torus fixture: a=(1,0) at y=1/5, b=(0,1) at x=2/7,
# c=(1,1) through (1/11,1/13); l1,m1 are the coordinate circles.
@testset "partial twist changes determinant for c=a+b" begin
    lc=parse_lc("genus: 1\nl1 : -1 -3 -2\nm1 : 1 4 5\na : -4 -7 -6\nb : 2 6 8\nc : 7 -8 -5 3\n")
    @test checker(lc)[1]
    @test abs(algebraic_intersections(lc,["a"],["b"])[1,1])==1
    w=only(H.movable(lc,"lc.twist.annulus";q=1,support="c",acted_components=["b"]))
    r=H.move(lc,w)
    @test checker(r.output)[1]
    @test abs(algebraic_intersections(r.output,["a"],["b"])[1,1])==2
    @test same_up_to(lc,H.move(r.output,r.inverse).output)
end

@testset "handle slide empty band" begin
    lc = H.read_lc(joinpath(ROOT, "threads", "PJ260922", "code", "lc-fixtures", "p28-after.lc"))
    hc = H.HeegaardCode(lc, Dict("alpha1" => :alpha, "alpha2" => :alpha, "beta1" => :beta, "beta2" => :beta),
        Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
    slides = H.movable(hc, "hc.slide.pants"; curves=["alpha1", "alpha2"], band=Any[])
    @test !isempty(slides)
    slid = H.move(hc, first(slides))
    sok, sreport = checker(slid.output.lc)
    @test sok
    @test slid.output.lc.genus == lc.genus
    back = H.move(slid.output, slid.inverse)
    @test back.output.lc.genus == lc.genus
    @test same_up_to(lc,back.output.lc)
end

@testset "stabilization genus" begin
    base = H.read_lc(joinpath(ROOT, "threads", "PJ260922", "code", "lc-fixtures", "operations", "g1-base.lc"))
    g1 = H.HeegaardCode(base, Dict("alpha1" => :alpha, "beta1" => :beta),
        Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing)
    stabs = H.movable(g1, "hc.stabilization.once"; disk="F1", k_identification="ladder_end")
    @test !isempty(stabs)
    stabilized = H.move(g1, only(stabs))
    @test stabilized.output.lc.genus == base.genus + 1
    tok, treport = checker(stabilized.output.lc)
    @test tok
    undone = H.move(stabilized.output, stabilized.inverse)
    @test undone.output.lc.genus == base.genus
    @test same_up_to(base, undone.output.lc)
end

@testset "group-15 annulus twist inverse restores words" begin
    hc = H.read_hc(joinpath(ROOT, "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc"))
    original = Dict(name => copy(word) for (name, word) in hc.lc.words)
    for (support, acted) in (("m1", ["alpha1"]), ("alpha1", ["beta1"]), ("beta1", ["alpha1"]))
        hit = only(H.movable(hc, "lc.twist.annulus"; q=1, support=support, acted_components=acted))
        forward = H.move(hc, hit)
        back = H.move(forward.output, forward.inverse)
        @test back.output.lc.genus == hc.lc.genus
        @test back.output.lc.words == original
    end
end

@testset "ribbon banks and corpus relations" begin
    for path in ("code/lc-fixtures/p28-after.hc","data/261107-lickorish-data/group-18/case-01.hc")
        hc=H.read_hc(joinpath(ROOT,"threads","PJ260922",path));lc=hc.lc
        names=H._component_names(lc);original=algebraic_intersections(lc,names,names)
        faces,right=H._lc_faces(lc.words);occ=H._lc_rotation_occ(lc.words)
        for role in (:alpha,:beta)
            fam=sort([n for(n,r) in hc.roles if r==role]);a,b=fam[1:2]
            bands=H._move_face_bands(lc,a,b)
            @test !isempty(bands)
            chosen=Any[first(bands)]
            # A specified route leaves a face and comes back across the same edge
            # at two distinct parameters; all face chords must remain disjoint.
            found=false
            for band in bands
                f=only(f for f in faces if f.id==right[H._move_corner(band["start"]["dart"])] )
                for c in f.corners
                    occ[c[1]][c[2]][1] in fam && continue
                    trial=deepcopy(band)
                    trial["crossings"]=[Dict("dart"=>collect(c),"u"=>"1/3"),
                        Dict("dart"=>collect(H._lc_tau(lc.words,occ,c)),"u"=>"2/3")]
                    try
                        H._move_certify_route(lc,trial;avoid=fam)
                        push!(chosen,trial);found=true;break
                    catch e
                        e isa H.LCMoveError || rethrow()
                    end
                end
                found && break
            end
            @test found
            for band in chosen
                w=only(H.movable(hc,"hc.slide.pants";curves=[a,b],band=band))
                r=H.move(hc,w);back=H.move(r.output,r.inverse)
                @test checker(r.output.lc)[1]
                @test same_up_to(lc,back.output.lc)
                coef=w.extra["support"]["coefficient"];ia=findfirst(==(a),names);ib=findfirst(==(b),names)
                expected=[original[i,j]+(i==ia)*coef*original[ib,j]+(j==ia)*coef*original[i,ib] for i in eachindex(names),j in eachindex(names)]
                @test algebraic_intersections(r.output.lc,names,names)==expected
            end
        end
    end
    for si in (-1,1),sj in (-1,1)
        @test H._move_barbell_boundary(si,sj)==[1,3,2si*sj,-3]
    end
end

@testset "common affine atlas and corpus transvections" begin
    for n in (1,2,5,17),q in (0,-1,1,-3,3)
        atlas=Dict{String,Any}("runs"=>[(component="c$(i%2)",core_label=i,
            p0=(BigInt(i)//BigInt(2n+1),BigInt(0)//BigInt(1)),
            p1=(BigInt(i)//BigInt(2n+1),BigInt(1)//BigInt(1))) for i in 1:n])
        acted=Set(["c0","c1"])
        @test H._move_annulus_map(H._move_annulus_map(atlas,acted,q),acted,-q)==atlas
        @test H._move_annulus_map(H._move_annulus_map(atlas,acted,q),acted,2)==H._move_annulus_map(atlas,acted,q+2)
    end
    for path in ("code/lc-fixtures/operations/g1-base.lc","code/lc-fixtures/p28-after.lc","data/261107-lickorish-data/group-18/case-01.lc")
        lc=H.read_lc(joinpath(ROOT,"threads","PJ260922",path));names=H._component_names(lc)
        original=algebraic_intersections(lc,names,names);acted=[first(n for n in names if !H._is_k(lc,n))]
        for core in (first(names),last(names)),q in (-3,1)
            w=only(H.movable(lc,"lc.twist.annulus";q=q,support=core,acted_components=acted))
            r=H.move(lc,w);back=H.move(r.output,r.inverse);c=findfirst(==(core),names)
            expected=[original[i,j]-q*((names[i] in acted)*original[c,i]*original[c,j]+(names[j] in acted)*original[c,j]*original[i,c]) for i in eachindex(names),j in eachindex(names)]
            @test checker(r.output)[1]
            @test algebraic_intersections(r.output,names,names)==expected
            @test same_up_to(lc,back.output)
            @test r.inverse.extra["annulus_recompute"]
        end
    end
end

@testset "surface inclusion on corpus representatives" begin
    for path in ("code/lc-fixtures/operations/g1-base.hc","code/lc-fixtures/p28-after.hc","data/261107-lickorish-data/group-18/case-01.hc","data/261107-lickorish-data/group-08/case-01.hc")
        hc=H.read_hc(joinpath(ROOT,"threads","PJ260922",path));lc=hc.lc
        faces,right=H._lc_faces(lc.words);disk=right[(lc.genus,"-","out")]
        w=only(H.movable(hc,"hc.stabilization.once";disk=disk,k_identification="ladder_end"))
        r=H.move(hc,w);back=H.move(r.output,r.inverse)
        @test checker(r.output.lc)[1]
        @test same_up_to(lc,back.output.lc)
        @test r.output.lc.genus==lc.genus+1
        names=H._component_names(lc)
        @test algebraic_intersections(r.output.lc,names,names)==algebraic_intersections(lc,names,names)
        newa=only(n for(n,role) in r.output.roles if role==:alpha && !haskey(hc.roles,n))
        newb=only(n for(n,role) in r.output.roles if role==:beta && !haskey(hc.roles,n))
        @test length(intersect(abs.(r.output.lc.words[newa]),abs.(r.output.lc.words[newb])))==1
        for n in keys(hc.roles),new in (newa,newb)
            @test isempty(intersect(abs.(r.output.lc.words[n]),abs.(r.output.lc.words[new])))
        end
        for (n,word) in H._lc_k_constants(lc.genus+1)
            @test cyclic_equal(filter(x->abs(x)<=3*(lc.genus+1)-2,r.output.lc.words[n]),word)
        end
        @test !isempty(w.extra["support"]["handle_triangles"])
        @test !isempty(w.extra["support"]["boundary_gluing"])
        @test !isempty(r.correspondence["k_identification"]["faces"])
        @test back.correspondence["support_direction"] == -1
    end
end


@testset "punctured torus template without old non-K crossings" begin
    lc=parse_lc("genus: 1\nl1 : -1\nm1 : 1\n")
    empty_layer=H.resolve_hc(lc,Dict("l1"=>:beta);standard=true)
    w=only(H.movable(empty_layer,"hc.stabilization.once";disk="F1",k_identification="ladder_end"))
    r=H.move(empty_layer,w)
    @test checker(r.output.lc)[1]
    @test length(unique(abs.(reduce(vcat,values(r.output.lc.words)))))==8
    @test same_up_to(lc,H.move(r.output,r.inverse).output.lc)
end

# cf: answer-3-full.md §§1-3 — independently verified oriented collars, normal annulus certificates and fixed marks.
function embedded_collar_fixture(lc,core)
    word=lc.words[core];n=length(word);boundaries=[Any[],Any[]];runs=Any[]
    for (i,x) in enumerate(word)
        partner=H._owner_of(lc.words,abs(x),x>0 ? "-" : "+");w=lc.words[partner];p=findfirst(==(-x),w);s=i//(2n+1)
        for level in 0:1
            side=level==0 ? -1 : 1;before=side*sign(x)==1
            edge=before ? mod1(p-1,length(w)) : p
            # Both boundaries follow the core orientation, hence passage sign sign(x).
            y=sign(x)==1 ? w[edge] : w[mod1(edge+1,length(w))]
            dart=[abs(y),y>0 ? "+" : "-",sign(x)==1 ? "out" : "in"]
            push!(boundaries[level+1],Dict("id"=>"p$(i)_$level","dart"=>dart,"u"=>(before ? "2/3" : "1/3"),"s"=>string(s)))
        end
        a=x<0 ? 0 : 1;b=1-a
        push!(runs,Dict("component"=>partner,"start_port"=>"p$(i)_$a","finish_port"=>"p$(i)_$b","points"=>[[string(s),a],[string(s),b]]))
    end
    push!(runs,Dict("component"=>core,"closed"=>true,"points"=>[[0,"1/2"],[1,"1/2"]]))
    Dict("kind"=>"embedded_annulus","boundary0"=>Dict("events"=>boundaries[1]),"boundary1"=>Dict("events"=>boundaries[2]),"runs"=>runs)
end

@testset "standard alpha support separation" begin
    for g in 1:3
        lc=parse_lc("genus: $g\n"*join([n*" : "*join(w," ") for(n,w) in H._lc_k_constants(g)],"\n"))
        hc=H.resolve_hc(lc,Dict("l$i"=>:beta for i in 1:g);standard=true)
        automatic=only(H.movable(hc,"lc.twist.annulus";support="l1",q=1,acted_components=["m1"]))
        ar=H.move(hc,automatic)
        @test same_up_to(lc,H.move(ar.output,ar.inverse).output.lc)
        ws=H.movable(hc,"hc.role-support.separate";curve="m1")
        @test length(ws)==2
        for w in ws
            r=H.move(hc,w);back=H.move(r.output,r.inverse).output;copyname=only(w.extra["created_components"])
            @test checker(r.output.lc)[1]
            @test back.roles==hc.roles
            @test same_up_to(lc,back.lc)
            @test !haskey(r.output.roles,"m1") && r.output.roles[copyname]==:alpha
            names=H._component_names(lc);renamed=[n=="m1" ? copyname : n for n in names]
            @test algebraic_intersections(lc,names,names)==algebraic_intersections(r.output.lc,renamed,renamed)
            twist=H.move(r.output,only(H.movable(r.output,"lc.twist.annulus";support="l1",q=1,acted_components=[copyname])))
            @test checker(twist.output.lc)[1]
            @test same_up_to(r.output.lc,H.move(twist.output,twist.inverse).output.lc)
        end
    end
end

@testset "normal annulus with full source overlay" begin
    for path in ("code/lc-fixtures/operations/g1-base.lc","code/lc-fixtures/p28-after.lc","data/261107-lickorish-data/group-18/case-01.lc")
        lc=H.read_lc(joinpath(ROOT,"threads/PJ260922",path));names=H._component_names(lc);m=algebraic_intersections(lc,names,names)
        core="m1";c=findfirst(==(core),names);acted=[first(n for n in names if !H._is_k(lc,n))]
        raw=embedded_collar_fixture(lc,core)
        for q in (-3,1,3)
            w=only(H.movable(lc,"lc.twist.annulus";support=raw,q=q,acted_components=acted));r=H.move(lc,w)
            @test checker(r.output)[1]
            @test same_up_to(lc,H.move(r.output,r.inverse).output)
            expected=[m[i,j]-q*((names[i] in acted)*m[c,i]*m[c,j]+(names[j] in acted)*m[c,j]*m[i,c]) for i in eachindex(names),j in eachindex(names)]
            @test algebraic_intersections(r.output,names,names)==expected
            # Apply a second twist to the wound source, not to a fresh thin collar.
            wound=w.extra["embedded_annulus_after"]
            w2=only(H.movable(r.output,"lc.twist.annulus";support=wound,q=-q,acted_components=acted))
            @test same_up_to(lc,H.move(r.output,w2).output)
        end
        missing=deepcopy(raw);pop!(missing["runs"])
        @test_throws H.LCMoveError H.movable(lc,"lc.twist.annulus";support=missing,q=1,acted_components=acted)
        bad=deepcopy(raw);bad["runs"][1]["points"][1][1]="17/19"
        @test_throws H.LCMoveError H.movable(lc,"lc.twist.annulus";support=bad,q=1,acted_components=acted)
    end
    # A returning cap and a closed circle intersect twice with opposite signs.
    cap=(component="cap",closed=false,points=[(1//5,0//1),(1//2,1//2),(4//5,0//1)])
    arc=(component="arc",closed=true,points=[(0//1,1//5),(1//1,1//5)])
    @test_throws H.LCMoveError H._move_pl_events([cap,merge(arc,(closed=false,))])
    ev=H._move_pl_events([cap,arc])
    @test length(ev.hits)==2
    @test sum(x[2] for x in ev.events[1])==0
    for q in (-3,1,3)
        mapped=[merge(r,(points=[(s+q*t,t) for(s,t) in r.points],)) for r in (cap,arc)]
        @test length(H._move_pl_events(mapped).hits)==2
    end
end

@testset "all empty-band occurrences" begin
    hc=H.read_hc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/p28-after.hc"));lc=hc.lc
    faces,right=H._lc_faces(lc.words);occ=H._lc_rotation_occ(lc.words)
    a,b="alpha1","alpha2";names=H._component_names(lc);m=algebraic_intersections(lc,names,names);ia=findfirst(==(a),names);ib=findfirst(==(b),names)
    expected=Set((x,y) for f in faces for x in f.corners for y in f.corners if occ[x[1]][x[2]][1]==a && occ[y[1]][y[2]][1]==b)
    ws=H.movable(hc,"hc.slide.pants";curves=[a,b],band=Any[])
    actual=Set((w.extra["support"]["route"].start.dart,w.extra["support"]["route"].finish.dart) for w in ws)
    @test actual==expected
    @test length(ws)==length(expected)>1
    for w in ws
        r=H.move(hc,w);coef=w.extra["support"]["coefficient"]
        @test checker(r.output.lc)[1]
        @test same_up_to(lc,H.move(r.output,r.inverse).output.lc)
        @test algebraic_intersections(r.output.lc,names,names)==[m[i,j]+(i==ia)*coef*m[ib,j]+(j==ia)*coef*m[i,ib] for i in eachindex(names),j in eachindex(names)]
    end
end

@testset "returning runs in a certified source annulus" begin
    lc=H.read_lc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/operations/g1-base.lc"))
    raw=embedded_collar_fixture(lc,"m1");e=last(raw["boundary0"]["events"]);d=H._move_corner(e["dart"]);occ=H._lc_rotation_occ(lc.words)
    push!(raw["boundary0"]["events"],Dict("id"=>"cap_in","dart"=>collect(H._lc_tau(lc.words,occ,d)),"u"=>"1/2","s"=>"3/5"))
    push!(raw["boundary0"]["events"],Dict("id"=>"cap_out","dart"=>collect(d),"u"=>"2/3","s"=>"4/5"))
    push!(raw["runs"],Dict("component"=>"alpha1","start_port"=>"cap_in","finish_port"=>"cap_out","points"=>[["3/5",0],["7/10","1/4"],["4/5",0]]))
    names=H._component_names(lc);m=algebraic_intersections(lc,names,names);c=findfirst(==("m1"),names);acted=["alpha1"]
    for q in (-3,-1,1,3)
        w=only(H.movable(lc,"lc.twist.annulus";support=raw,q=q,acted_components=acted));r=H.move(lc,w)
        @test checker(r.output)[1]
        @test same_up_to(lc,H.move(r.output,r.inverse).output)
        @test algebraic_intersections(r.output,names,names)==[m[i,j]-q*((names[i] in acted)*m[c,i]*m[c,j]+(names[j] in acted)*m[c,j]*m[i,c]) for i in eachindex(names),j in eachindex(names)]
        altered=deepcopy(w);altered.after_frags[1]=Int[]
        err=try H.move(lc,altered);nothing catch e;e end
        @test err isa H.LCMoveError && err.kind==:witness_plan_mismatch
    end
end

@testset "fixed marks on unchanged cells" begin
    hc=H.read_hc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/p28-after.hc"));lc=hc.lc
    ws=H.movable(hc,"lc.isotopy.digon.create.ab")
    w=first(w for w in ws if !isempty(first(H._move_untouched_cells(lc,w))));safe=first(first(H._move_untouched_cells(lc,w)))
    marked=H.HeegaardCode(lc,hc.roles,Dict("D"=>[safe]),Dict("z"=>safe),Dict(safe=>0),nothing)
    selected=H.movable(marked,w.row_id;mark_transport="fixed_cells")
    w2=only(x for x in selected if x.locations==w.locations && x.after_frags==w.after_frags)
    r=H.move(marked,w2);back=H.move(r.output,r.inverse).output
    @test back.disks==marked.disks && back.points==marked.points && back.weights==marked.weights
    @test length(r.output.weights)==1 && only(values(r.output.weights))==0
    f=only(f for f in first(H._lc_faces(lc.words)) if f.id==safe)
    g=only(f for f in first(H._lc_faces(r.output.lc.words)) if f.id==r.output.points["z"])
    @test Set(g.corners)==Set((r.correspondence["label_map"][c[1]],c[2],c[3]) for c in f.corners)
    changed=H.HeegaardCode(lc,hc.roles,marked.disks,Dict("z"=>safe),Dict(safe=>1),nothing)
    changed_result=H.move(changed,w2)
    @test sum(values(changed_result.output.weights))==1
    inside=w.extra["face"]
    ambiguous=H.HeegaardCode(lc,hc.roles,Dict{String,Vector{String}}(),Dict("z"=>inside),Dict{String,Int}(),nothing)
    aw=only(x for x in H.movable(ambiguous,w.row_id;mark_transport="fixed_cells") if x.locations==w.locations && x.after_frags==w.after_frags)
    moved=H.move(ambiguous,aw)
    @test moved.output.points["z"] in H._move_mark_cells(moved.output.lc).faces
    @test any(r.kind==:point && r.mark=="z" for r in moved.correspondence["mark_cells"].refs)
end

@testset "free mark transport and geometry inverse" begin
    for path in ("code/lc-fixtures/operations/g1-base.hc","code/lc-fixtures/p28-after.hc","data/261107-lickorish-data/group-18/case-01.hc")
        hc=H.read_hc(joinpath(ROOT,"threads/PJ260922",path));w=first(H.movable(hc,"lc.isotopy.digon.create.ab"));face=w.extra["face"]
        marked=H.HeegaardCode(hc.lc,hc.roles,Dict("D"=>[face]),Dict("z"=>face),Dict(face=>0),nothing)
        plain=H.move(hc,w);r=H.move(marked,w);cells=H._move_mark_cells(r.output.lc)
        @test r.output.lc.words==plain.output.lc.words
        @test checker(r.output.lc)[1]
        @test all(x in cells.faces for x in r.output.disks["D"])
        @test r.output.points["z"] in vcat(cells.faces,cells.edges,cells.empties)
        @test !isempty(r.output.weights) && all(==(0),values(r.output.weights))
        # A geometric inverse remains usable after the user changes free marks.
        edited=H.HeegaardCode(r.output.lc,r.output.roles,r.output.disks,Dict("z"=>last(cells.faces)),r.output.weights,nothing)
        back=H.move(edited,r.inverse).output
        @test same_up_to(hc.lc,back.lc)
        @test back.roles==hc.roles
    end
    weights=Dict{String,Int}();H._move_add_weight!(weights,"F1",2);H._move_add_weight!(weights,"F1",-2)
    @test haskey(weights,"F1") && weights["F1"]==0 && !haskey(weights,"F2")
    source=H.read_hc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/operations/g1-base.hc"))
    source=H.HeegaardCode(source.lc,source.roles,Dict("D"=>["F1"]),Dict("z"=>"F1"),Dict("F1"=>0),nothing)
    empty_target=H.LickorishCode(0,Dict{String,Vector{Int}}(),Dict{String,String}(),Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
    deferred,receipt=H._move_transport_marks(source,empty_target,nothing,Dict{Int,Int}())
    @test isempty(deferred.words) && isempty(deferred.disks) && isempty(deferred.points) && isempty(deferred.weights)
    @test length(receipt.refs)==3 && all(r.mode==:deferred for r in receipt.refs)
end

@testset "data-defined mark edits" begin
    for path in ("data/261107-lickorish-data/group-15/case-01.hc","data/261107-lickorish-data/group-18/case-01.hc")
        hc=H.read_hc(joinpath(ROOT,"threads/PJ260922",path));fs=H._move_mark_cells(hc.lc).faces
        h=H.HeegaardCode(hc.lc,hc.roles,copy(hc.disks),Dict("z"=>first(fs)),Dict(first(fs)=>2,last(fs)=>-2),nothing)
        js=[r.id for r in H.forget_reference(H.realize_hc(h)).regions if r.joining_disk]
        jobs=[("hc.mark.point.move",(;name="z",target=last(fs))),
            ("hc.mark.disk.move",(;name="D",target=Dict("joining_region"=>last(js)))),
            ("hc.mark.weight.move",(;name=first(fs),target=last(fs)))]
        for (row,kw) in jobs
            w=only(H.movable(h,row;kw...));r=H.move(h,w);back=H.move(r.output,r.inverse).output
            @test r.output.lc===h.lc && r.output.roles==h.roles
            @test H._marks(back)==H._marks(h)
            @test checker(r.output.lc)[1]
            if row=="hc.mark.weight.move";@test haskey(r.output.weights,last(fs)) && r.output.weights[last(fs)]==0;end
        end
        @test_throws H.LCMoveError H.movable(h,"hc.mark.point.move";name="z",target="F999999")
        @test_throws H.LCMoveError H.movable(h,"hc.mark.disk.move";name="D",target=Dict("joining_region"=>999999))
    end
end

@testset "automatic shared-support slides" begin
    for g in 2:3
        lc=parse_lc("genus: $g\n"*join([n*" : "*join(w," ") for(n,w) in H._lc_k_constants(g)],"\n"))
        h=H.resolve_hc(lc,Dict("l$i"=>:beta for i in 1:g);standard=true)
        alpha=["m$i" for i in 1:g];beta=["l$i" for i in 1:g];m=algebraic_intersections(lc,alpha,beta)
        ws=H.movable(h,"hc.slide.pants";curves=["m1","m2"],band=Any[])
        @test !isempty(ws)
        for w in ws
            r=H.move(h,w);back=H.move(r.output,r.inverse).output;rename=w.extra["preparation"]["role_component"]
            a=[get(rename,n,n) for n in alpha];coef=last(w.extra["sequence"]).extra["support"]["coefficient"]
            expected=copy(m);expected[1,:]+=coef*m[2,:]
            @test checker(r.output.lc)[1]
            @test same_up_to(lc,back.lc) && back.roles==h.roles
            @test algebraic_intersections(r.output.lc,a,beta)==expected
        end
    end
end

@testset "free attachment disk across K faces" begin
    h=H.read_hc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/p28-after.hc"));lc=h.lc
    faces=first(H._lc_faces(lc.words));reassigned=0
    for f in faces
        marked=H.HeegaardCode(lc,h.roles,Dict("D"=>[f.id]),Dict("z"=>f.id),Dict(f.id=>0),nothing)
        w=only(H.movable(marked,"hc.stabilization.once";disk="D",k_identification="ladder_end"))
        r=H.move(marked,w);back=H.move(r.output,r.inverse).output
        reassigned+=w.extra["preparation"]["reassigned"]
        @test checker(r.output.lc)[1]
        @test r.output.lc.genus==lc.genus+1
        @test same_up_to(lc,back.lc) && back.roles==h.roles
    end
    @test reassigned>0
end

end

# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-4-3-harvest.md —
# the creation band misses an old corner sector; actual nested circles follow their root.
@testset "digon creation with actual nested placed circles" begin
    base=H.read_lc(joinpath(ROOT,"threads/PJ260922/code/lc-fixtures/operations/g1-base.lc"))
    words=deepcopy(base.words)
    for n in ("c","r","d","e");words[n]=Int[];end
    places=Dict("c"=>"F1","r"=>"F1","d"=>"c+","e"=>"F2")
    lc=H.LickorishCode(base.genus,words,places,base.disks,base.points,base.weights,nothing)
    ws=[w for w in H.movable(lc,"lc.isotopy.digon.create.ab") if w.extra["face"]=="F1"]
    @test !isempty(ws)
    names=H._component_names(lc);m=algebraic_intersections(lc,names,names)
    for w in ws
        r=H.move(lc,w);after=r.output;back=H.move(after,r.inverse).output
        @test checker(after)[1]
        @test same_up_to(lc,back) && back.places==places
        @test after.places["c"]==after.places["r"] && after.places["d"]=="c+"
        target=only(f for f in first(H._lc_faces(after.words)) if f.id==after.places["c"])
        @test any(c[1] in values(r.correspondence["label_map"]) for c in target.corners)
        @test algebraic_intersections(after,names,names)==m
    end
    probe=H.move(lc,first(ws));created=Set(probe.correspondence["created_labels"])
    bigons=[f.id for f in first(H._lc_faces(probe.output.words)) if length(f.corners)==2 && all(c[1] in created for c in f.corners)]
    @test length(bigons)==1
    blocked_words=deepcopy(probe.output.words);blocked_words["block"]=Int[]
    blocked_places=copy(probe.output.places);blocked_places["block"]=only(bigons)
    blocked=H.LickorishCode(lc.genus,blocked_words,blocked_places,lc.disks,lc.points,lc.weights,nothing)
    @test checker(blocked)[1]
    @test all(Set(w.vertices)!=Set(probe.inverse.vertices) for w in H.movable(blocked,probe.inverse.row_id))
end

include("lc_moves_inverse_walk.jl")

println("LC_MOVES_TESTS_DONE")
