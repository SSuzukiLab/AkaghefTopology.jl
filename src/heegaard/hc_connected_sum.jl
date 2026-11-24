# cf: docs/specs/heegaard-snappy.md — genus-one surface connected sum and explicit pending joining disk.
# cf: .research-scratch/PJ260922/astra-lanes/snappy-bridge/answer-10-3-full-main.md, Implementation and Splice-order judgment — oriented collar compiler and neck gluing; reconstructed.

"""Connect two genus-one HC surfaces, retaining the input alpha/beta roles.
Returns (lc,hc,certificate). The marked punctures and reversed neck are recorded;
input disk marks are consumed and the output has no joining disk.
Empty curves and point/weight/LC decorations are outside this operation scope.
"""
function hc_connected_sum(left::HeegaardCode, right::HeegaardCode;
        left_gap::Union{Nothing,Int}=nothing, right_gap::Union{Nothing,Int}=nothing)
    row=Dict{String,Any}("id"=>"construction","inverse_id"=>nothing,"digest"=>"",
        "relation"=>"insertion","layer"=>"LC")
    target="__cs_parallel"
    labels(ws)=Set(abs(x) for w in values(ws) for x in w)
    rw(w,M)=[sign(x)*M[abs(x)] for x in w]
    fail(k,m,d=Dict{String,Any}())=throw(LCMoveError(k,m,Dict{String,Any}(d)))
    function kcheck(lc,ks,limit)
        p=Dict(k=>[x for x in lc.words[k] if abs(x)<=limit] for k in ks)
        q=LickorishCode(count(k->startswith(k,"m"),ks),p,
            Dict{String,String}(),Dict{String,Vector{String}}(),
            Dict{String,String}(),Dict{String,Int}(),nothing)
        validate_lc_input(q); _lc_validate_k_constants!(q); p
    end
    function compile_side(hc,tag,gap)
        realize_hc(hc); lc=hc.lc; side_name=string(tag)
        lc.genus==1 || fail(:connected_sum_input,"each input must have genus one",
            Dict("side"=>side_name,"genus"=>lc.genus))
        all(r->count(==(r),values(hc.roles))==1,(:alpha,:beta)) ||
            fail(:connected_sum_input,"each input needs exactly one alpha and one beta",
                Dict("side"=>side_name,"roles"=>copy(hc.roles)))
        isempty(lc.places) && all(w->!isempty(w),values(lc.words)) ||
            fail(:connected_sum_input,"placed or empty circles are unsupported",
                Dict("side"=>side_name))
        isempty(lc.disks) && isempty(lc.points) && isempty(lc.weights) &&
            isempty(hc.points) && isempty(hc.weights) ||
            fail(:connected_sum_input,"point/weight/LC decorations are unsupported",
                Dict("side"=>side_name))
        haskey(lc.words,target) && fail(:connected_sum_input,
            "reserved component name is already used",Dict("side"=>side_name,"name"=>target))
        kcheck(lc,["l1","m1"],1)
        m=lc.words["m1"]; n=length(m)
        n>0 && length(unique(abs.(m)))==n ||
            fail(:connected_sum_seed_support,"m1 is not a nonempty simple core",
                Dict("side"=>side_name))
        e=gap===nothing ? n : gap
        1<=e<=n || fail(:connected_sum_seed_support,"seed gap is out of range",
            Dict("side"=>side_name,"gap_after"=>e,"length"=>n))
        side,dir=tag===:left ? (-1,1) : (1,-1)
        j=side==-1 ? e : mod1(e+1,n); x=m[j]
        corner=(abs(x),x>0 ? "+" : "-",side==-1 ? "out" : "in")
        _,right_face=_lc_faces(lc.words); face=get(right_face,corner,nothing)
        face===nothing && fail(:connected_sum_seed_support,
            "chosen collar side has no incident input face",
            Dict("side"=>side_name,"gap_after"=>e,"corner"=>collect(corner)))
        run=(kind=:curve,component="m1",edge=e,side=side,direction=dir,survives=false)
        w=_move_compile_boundary(lc,row,target,[run],
            Dict("kind"=>"connected_sum_collar","side"=>side_name,
                 "face"=>face,"gap_after"=>e))
        r=_kernel_apply(w.names,lc.words,w.before_frags,w.after_frags,
            w.locations,w.vertices,w.capacity;label_map=w.extra["label_map"])[1]
        length(r)==length(w.names) ||
            fail(:connected_sum_parallel_support,"compiler returned a partial word bank")
        ws=Dict(w.names[i]=>r[i] for i in eachindex(w.names))
        tmp=LickorishCode(1,ws,Dict{String,String}(),
            Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
        validate_lc_input(tmp); _lc_validate_surface!(tmp)
        p=ws[target]; p==w.after_frags[1] || fail(:connected_sum_parallel_support,
            "compiler moved the recorded seed cut",Dict("side"=>side_name))
        inds=dir==1 ? [mod1(e+k,n) for k in 1:n] :
                      [mod1(e-k,n) for k in 0:n-1]
        length(p)==n || fail(:connected_sum_parallel_support,
            "parallel copy omitted source crossings",Dict("side"=>side_name))
        old=labels(lc.words)
        all(x->!(abs(x) in old),p) && length(unique(abs.(p)))==n ||
            fail(:connected_sum_parallel_support,
                "parallel copy did not give one fresh crossing per m1 visit",
                Dict("side"=>side_name))
        hits=[i for i in eachindex(inds) if abs(m[inds[i]])==1]
        length(hits)==1 || fail(:connected_sum_parallel_support,
            "copied K crossing is not unique",Dict("side"=>side_name))
        corr=[(source_position=inds[i],source_entry=m[inds[i]],copy_entry=p[i])
              for i in eachindex(inds)]
        (;hc,lc,ws,p,old,face,corner,gap=e,side,dir,
          kcopy=abs(p[only(hits)]),corr)
    end
    L=compile_side(left,:left,left_gap); R=compile_side(right,:right,right_gap)
    LM=Dict(1=>1,L.kcopy=>3); RM=Dict(1=>2,R.kcopy=>4); fresh=5
    for x in sort!(collect(labels(L.ws))); haskey(LM,x)||(LM[x]=fresh;fresh+=1); end
    for x in sort!(collect(labels(R.ws))); haskey(RM,x)||(RM[x]=fresh;fresh+=1); end
    LC=Dict(n=>n=="l1" ? "l1" : n=="m1" ? "m1" : "L__"*n
            for n in keys(L.lc.words))
    RC=Dict(n=>n=="l1" ? "l2" : n=="m1" ? "m2" : "R__"*n
            for n in keys(R.lc.words))
    words=Dict{String,Vector{Int}}("l1"=>rw(L.ws["l1"],LM),
        "m1"=>rw(L.ws["m1"],LM),"l2"=>rw(R.ws["l1"],RM),
        "m2"=>rw(R.ws["m1"],RM),"n1"=>vcat(rw(L.p,LM),rw(R.p,RM)))
    for n in keys(L.lc.words); n in ("l1","m1")||(words[LC[n]]=rw(L.ws[n],LM)); end
    for n in keys(R.lc.words); n in ("l1","m1")||(words[RC[n]]=rw(R.ws[n],RM)); end
    roles=Dict{String,Symbol}()
    for (n,r) in L.hc.roles; roles[LC[n]]=r; end
    for (n,r) in R.hc.roles; roles[RC[n]]=r; end
    lc=LickorishCode(2,words,Dict{String,String}(),
        Dict{String,Vector{String}}(),Dict{String,String}(),
        Dict{String,Int}(),nothing)
    validate_lc_input(lc)
    try
        _lc_validate_surface!(lc)
    catch e
        fail(:connected_sum_neck_order,
            "left-then-right splice does not realize a genus-two surface",
            Dict("cause"=>sprint(showerror,e)))
    end
    kp=kcheck(lc,["l1","l2","m1","m2","n1"],4)
    function source_check(S,C,M)
        keep=Set(M[x] for x in S.old)
        for (n,w) in S.lc.words
            got=[x for x in words[C[n]] if abs(x) in keep]
            got==rw(w,M) || fail(:connected_sum_embedding,
                "deleting new labels did not recover a source word",
                Dict("component"=>C[n],"expected"=>rw(w,M),"got"=>got))
        end
    end
    source_check(L,LC,LM); source_check(R,RC,RM)
    hc=resolve_hc(lc,roles;disks=Dict{String,Vector{String}}()); realize_hc(hc)
    pc(S,M)=[Dict("source_position"=>c.source_position,
        "source_entry"=>sign(c.source_entry)*M[abs(c.source_entry)],
        "n1_entry"=>sign(c.copy_entry)*M[abs(c.copy_entry)]) for c in S.corr]
    certificate=Dict{String,Any}(
        "status"=>"joining-disk-pending","joining_disk"=>nothing,
        "forward_koda_snappy"=>"not run or claimed",
        "punctures"=>Dict("left"=>Dict("face"=>L.face,"gap_after"=>L.gap,
            "side"=>L.side,"direction"=>L.dir,"corner"=>collect(L.corner)),
            "right"=>Dict("face"=>R.face,"gap_after"=>R.gap,
            "side"=>R.side,"direction"=>R.dir,"corner"=>collect(R.corner))),
        "neck"=>Dict("orientation"=>"reversing","parameter_map"=>"s_right=1-s_left",
            "n1_order"=>["left_parallel","right_parallel"],
            "endpoint_pairs"=>[["left_end","right_start"],["right_end","left_start"]]),
        "source_embeddings"=>Dict(
            "left"=>Dict("source_path"=>L.lc.source_path,"components"=>LC,"labels"=>Dict(x=>LM[x] for x in L.old)),
            "right"=>Dict("source_path"=>R.lc.source_path,"components"=>RC,"labels"=>Dict(x=>RM[x] for x in R.old))),
        "parallel_crossings"=>Dict("left"=>pc(L,LM),"right"=>pc(R,RM)),
        "discarded_input_disks"=>Dict("left"=>deepcopy(left.disks),
                                      "right"=>deepcopy(right.disks)),
        "checks"=>Dict("genus"=>lc.genus==2,"K_projection"=>kp,
            "source_words"=>true,"alpha_count"=>count(==(:alpha),values(roles)),
            "beta_count"=>count(==(:beta),values(roles))))
    (;lc,hc,certificate)
end
