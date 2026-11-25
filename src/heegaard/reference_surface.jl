"""
    lickorish_reference(g::Integer) -> Dict{String,Any}

Construct the named oriented reference-curve record for genus g.

# Arguments
`g` must be a non-Boolean integer from 1 through 16; otherwise `ArgumentError` with `REFERENCE_GENUS` is thrown. Internal face-orbit, genus or face-size inconsistencies also raise `ArgumentError`.

# Returns
A dictionary with `genus`, `curves`, `crossings` and `metadata`. The curve records have auxiliary roles; this is reference data, not an assigned alpha/beta diagram or a `LickorishCode`.

# Guarantees
1. The returned curve arrangement carries a connected oriented ribbon surface whose capped face decomposition has genus `g`. [proved/agent: connected reference graph and the checked Euler characteristic]
2. The record contains only the reference-curve arrangement and source metadata, without assigning a Heegaard splitting. [proved/agent: by construction]

# Conventions
Historical reference-record translation: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/lickorish-code.md#historical-reference-records).

# Not guaranteed
The historical record's port words are not asserted to be the modern LC interchange format. The source figure's geometric identification and the translation of its auxiliary orientation choices to every modern LC consumer are not independently certified by this constructor.

# Sources
cf: lib/topology.jl/docs/specs/lickorish-code.md — reference curves and the distinction from input curves.
cf: [[AS13#^anonnd60-sec-0740]] Fig. 12.3 — source of the reference arrangement.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> r = H.lickorish_reference(1); (length(r["curves"]), length(r["crossings"]))
(2, 1)
```

Role: construct
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
