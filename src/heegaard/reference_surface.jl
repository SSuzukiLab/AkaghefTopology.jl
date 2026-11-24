"""
    lickorish_reference(g)

Named oriented disk decomposition.
cf:[[AS13#^anonnd60-sec-0740]] Fig. 12.3, printed p.244 / PDF p.251 — linked OCR
crop references/derived/n8n-jobs/01d1e1652c07435ea1f51e46a2933bfc/assets/p0251-0001.jpeg.
The crop was visually inspected; the mathematical transcription is agent-owned.
The curves l/m/n are REFERENCE curves, never the input alpha/beta families.
The surface is the oriented ribbon neighbourhood capped by its named face disks.
Returns the reference crossings and K curves only. The surface check is
V - E + F = 2 - 2g on the traced ribbon faces (E = 2V for 4-valent crossings).

Conventions: l_i follows the arrow in the figure; m_i crosses l_i positively.
n_i is oriented with signs +(l_i,n_i), -(l_{i+1},n_i). The n_i orientations
are our choice (the figure does not assign arrows to them).
"""
function lickorish_reference(g::Integer)
    !(g isa Bool) && 1<=g<=16 || throw(ArgumentError("REFERENCE_GENUS: 1 <= genus <= 16"))
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
            z in unseen || throw(ArgumentError("REFERENCE_FACE: $z"))
            delete!(unseen,z);push!(f,z);z=inverse[twin[z]];z==start && break
        end
        push!(faces,f)
    end
    nv=length(rotations)
    nv-2nv+length(faces)==2-2g || throw(ArgumentError("REFERENCE_NOT_GENUS_$g"))
    all(length(f)>=3 for f in faces) || throw(ArgumentError("REFERENCE_FACE_SIZE: genus $g"))
    ports=Dict(d=>"p$(j-1)" for v in keys(rotations) for (j,d) in enumerate(rotations[v]))
    curves=[Dict("id"=>c,"role"=>"aux","closed"=>true,
        "word"=>join(["V_$v{"*ports[(c,v,-1)]*" z "*ports[(c,v,1)]*"}" for v in cycles[c]]," | "))
        for c in sort(collect(keys(cycles)))]
    Dict("genus"=>g,"curves"=>curves,
        "crossings"=>[Dict("id"=>v,"at"=>["V_$v","z"],"curves"=>collect(pairs[v]),"sign"=>signs[v]) for v in sort(collect(keys(pairs)))],
        "metadata"=>Dict("source"=>"AS13, Fig. 12.3, printed p.244 / PDF p.251",
            "source_work_id"=>"w-a3c2a004-c587-4cfd-ae05-b73a924d15d7",
            "source_figure"=>Dict("cite_key"=>"AS13","anchor"=>"anonnd60-sec-0740",
                "title"=>"パズルゲームで楽しむ写像類群入門",
                "figure"=>"12.3","printed_page"=>244,"pdf_page"=>251,
                "asset"=>"references/derived/n8n-jobs/01d1e1652c07435ea1f51e46a2933bfc/assets/p0251-0001.jpeg",
                "asset_sha256"=>"c385ea4b13a46c869d40c25e96742f761d9b7121f1413931dac28e96e031dc73"),
            "source_status"=>"figure visually inspected; oriented combinatorial transcription, not owner-checked",
            "reference_only"=>true,"display_cuts"=>["l$i" for i in 1:g]))
end
