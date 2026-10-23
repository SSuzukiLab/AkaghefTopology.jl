"""
    lickorish_reference(g)

Named oriented disk decomposition for Fig. 12.3 (Ahara--Sakasai, p. 244).
The curves l/m/n are REFERENCE curves, never the input alpha/beta families.
The surface is the oriented ribbon neighbourhood capped by its named face disks.
This constructs a chosen marked model, not an automatic marking of legacy data.

Conventions: l_i follows the arrow in the figure; m_i crosses l_i positively.
n_i is oriented with signs +(l_i,n_i), -(l_{i+1},n_i). The n_i orientations
are our choice (the figure does not assign arrows to them).
"""
function lickorish_reference(g::Integer)
    M=MarkedSurfaceDiagrams
    M.need(!(g isa Bool) && 1<=g<=16,"REFERENCE_GENUS","1 <= genus <= 16 (kernel patch limit)")
    cycles=Dict{String,Vector{String}}(); pairs=Dict(); signs=Dict()
    for i in 1:g
        v="lm_$i"; cycles["m$i"]=[v]; pairs[v]=("l$i","m$i"); signs[v]=1
        cycles["l$i"]=[v]
        i<g && push!(cycles["l$i"],"ln_$(i)_$i")
        i>1 && push!(cycles["l$i"],"ln_$(i)_$(i-1)")
    end
    for i in 1:g-1
        a,b="ln_$(i)_$i","ln_$(i+1)_$i"
        cycles["n$i"]=[a,b]; pairs[a]=("l$i","n$i"); pairs[b]=("l$(i+1)","n$i")
        signs[a]=1; signs[b]=-1
    end
    rotations=Dict(); sigma=Dict(); twin=Dict()
    for v in sort(collect(keys(pairs)))
        a,b=pairs[v]; e=signs[v]
        r=[(a,v,1),(b,v,e),(a,v,-1),(b,v,-e)]; rotations[v]=r
        for i in 1:4; sigma[r[i]]=r[mod1(i+1,4)]; end
    end
    for c in sort(collect(keys(cycles)))
        vs=cycles[c]
        for (x,y) in zip(vs,circshift(vs,-1)); twin[(c,x,1)]=(c,y,-1); twin[(c,y,-1)]=(c,x,1); end
    end
    inverse=Dict(v=>k for (k,v) in sigma); unseen=Set(keys(twin)); faces=[]
    while !isempty(unseen)
        start=minimum(unseen); z=start; f=[]
        while true
            M.need(z in unseen,"REFERENCE_FACE",z)
            delete!(unseen,z);push!(f,z);z=inverse[twin[z]];z==start && break
        end
        push!(faces,f)
    end
    # Vertex disks plus face disks. Ports on a vertex disk are glued directly
    # along each ribbon edge. Side OCCURRENCES remain distinct even on loops.
    patches=[];gluings=[];ports=Dict();corners=Dict()
    for v in sort(collect(keys(rotations)))
        sides=Dict();order=[];bd=[]
        for (j,d) in enumerate(rotations[v])
            i=j-1;a="a$i";s="g$i"
            sides[a]=["b$(2i)","p$i","b$(2i+1)"]
            sides[s]=["b$(2i+1)","b$(mod(2i+2,8))"]
            append!(order,[a,s]);append!(bd,sides[a]);ports[d]=["V_$v",a,"p$i"];corners[d]=["V_$v",s]
        end
        push!(patches,Dict("face"=>"V_$v","side_order"=>order,"sides"=>sides,
            "triangles"=>[["z",bd[i],bd[mod1(i+1,length(bd))]] for i in eachindex(bd)]))
    end
    for d in sort(collect(keys(twin)))
        e=twin[d];isless(d,e) && push!(gluings,[ports[d][1:2],ports[e][1:2]])
    end
    for (j,face) in enumerate(faces)
        f="F$j";n=length(face);order=["e$k" for k in 1:n]
        M.need(n>=3,"REFERENCE_FACE_SIZE",(g,j,n))
        push!(patches,Dict("face"=>f,"side_order"=>order,
            "sides"=>Dict("e$k"=>["t$k","t$(mod1(k+1,n))"] for k in 1:n),
            "triangles"=>[["z","t$k","t$(mod1(k+1,n))"] for k in 1:n]))
        for (k,d) in enumerate(face);push!(gluings,[[f,"e$k"],corners[d]]);end
    end
    curves=[Dict("id"=>c,"role"=>"aux","closed"=>true,
        "word"=>join(["V_$v{"*ports[(c,v,-1)][3]*" z "*ports[(c,v,1)][3]*"}" for v in cycles[c]]," | "))
        for c in sort(collect(keys(cycles)))]
    doc=Dict("version"=>1,"surface"=>Dict("marking_id"=>"ahara-sakasai-fig12.3-g$g-v1",
        "genus"=>g,"gluings"=>gluings,"boundary"=>[]),"patches"=>patches,"curves"=>curves,
        "crossings"=>[Dict("id"=>v,"at"=>["V_$v","z"],"curves"=>collect(pairs[v]),"sign"=>signs[v]) for v in sort(collect(keys(pairs)))],
        "metadata"=>Dict("source"=>"AnonND60, Fig. 12.3, printed p.244 / PDF p.251",
            "source_work_id"=>"w-a3c2a004-c587-4cfd-ae05-b73a924d15d7",
            "source_status"=>"figure visually inspected; oriented combinatorial transcription, not owner-checked",
            "reference_only"=>true,"display_cuts"=>["l$i" for i in 1:g]))
    st=M.validate(doc);rs=M.regions(st)
    M.need(length(rs)==g && all(r["normalized_disk"] for r in rs),"REFERENCE_NOT_CELLULAR",g)
    cuts=Set("l$i" for i in 1:g)
    cut=M.regions(st;barriers=Set(e for (e,cs) in st.used if !isdisjoint(Set(cs),cuts)))
    M.need(length(cut)==1 && cut[1]["genus"]==0 && cut[1]["boundary_components"]==2g,"REFERENCE_DISPLAY_CUT",g)
    doc
end
