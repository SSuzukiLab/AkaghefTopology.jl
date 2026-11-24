# cf: Turaev–Viro (1992), DOI 10.1016/0040-9383(92)90015-A, §9.1 and §9.2.A — r=3 surface-colouring sum; Work registration unresolved, unverified.
# cf: lib/topology.jl/src/heegaard/lc_combinatorial.jl::_lc_faces / _lc_tau — existing oriented surface CW cells.

"""Exact r=3 TV sum on Σ ∪ meridian disks, without triangulation or Regina.

The full LC cell structure subdivides the skeleton. Reference-only edges impose
equal face colours. Meridian edges impose the admissible binary triple. For
`parity=true`, each selected closed surface contributes (-1)^χ; for `false`,
each contributes 1. Both use the original normalisation 1/4 (χ(skeleton)=2).
Returns Rational{BigInt}; r>=4 is outside this backend. Completeness of the
meridian systems is checked by hc_linking_form; its value is unused.
"""
function Z_TuraevViro_skeleton(hc::HeegaardCode,r::Integer=3;
        parity::Bool=true,certificate::Bool=false,max_states::Integer=1_000_000)
    r==3 || throw(ArgumentError("the direct skeleton backend currently supports r=3"))
    max_states>0 || throw(ArgumentError("max_states must be positive"))
    hc_linking_form(hc) # Reuse existing complete primitive meridian validation.
    if hc.lc.genus==0 && isempty(hc.lc.words)
        result=big(1)//2
        cert=(method="TV-1992-9.2.A",r=3,parity=parity,states=big(2),skeleton_chi=2,
              surface_euler_characteristics=[0,2],special_case="empty genus-zero sphere")
        return certificate ? (;value=result,certificate=cert) : result
    end
    cell=realize_hc(hc).value["realization"]
    faces=cell.faces
    faceindex=Dict(f.id=>i for (i,f) in enumerate(faces))
    meridians=vcat(hc_curve_names(hc,:alpha),hc_curve_names(hc,:beta))
    nf=length(faces);nv=length(cell.occurrence)
    disks=Dict(c=>nf+i for (i,c) in enumerate(meridians))
    edges=Vector{Int}[]
    for c in sort!(collect(keys(hc.lc.words)))
        word=hc.lc.words[c]
        for x in word
            corner=(abs(x),_lc_sign_char(x),"out")
            opposite=_lc_tau(hc.lc.words,cell.occurrence,corner)
            incident=[faceindex[cell.corner_face[corner]],faceindex[cell.corner_face[opposite]]]
            haskey(disks,c) && push!(incident,disks[c])
            push!(edges,incident)
        end
    end
    vertexcells=Vector{Int}[]
    for label in sort!(collect(keys(cell.occurrence)))
        incident=[faceindex[cell.corner_face[(label,s,w)]] for (s,w) in _LC_ROT_SEQ]
        for (c,_) in values(cell.occurrence[label])
            haskey(disks,c) && push!(incident,disks[c])
        end
        push!(vertexcells,unique(incident))
    end
    chi=nv-length(edges)+nf+length(meridians)
    chi==2 || throw(ArgumentError("closed Heegaard skeleton must have Euler characteristic 2"))
    field=AbstractAlgebra.GF(2)
    constraints=AbstractAlgebra.zero_matrix(field,length(edges),nf+length(meridians))
    for (i,incident) in enumerate(edges), j in incident
        constraints[i,j]+=one(field) # Repeated incidences cancel modulo two.
    end
    dim,kernel=AbstractAlgebra.nullspace(constraints)
    states=big(2)^dim
    states<=min(max_states,typemax(Int)) || throw(ArgumentError("direct state sum needs $states states, above max_states=$max_states or the Int range"))
    eulers=Int[];total=big(0)
    for mask in 0:Int(states)-1
        color=falses(nf+length(meridians))
        for j in 1:dim
            isodd(mask>>(j-1)) || continue
            for i in eachindex(color)
                color[i] = xor(color[i],!iszero(kernel[i,j]))
            end
        end
        surface_chi=count(color)-count(e->any(color[j] for j in e),edges)+
                    count(v->any(color[j] for j in v),vertexcells)
        push!(eulers,surface_chi)
        total += parity && isodd(surface_chi) ? -1 : 1
    end
    result=total//big(4)
    cert=(method="TV-1992-9.2.A",r=3,parity=parity,states=states,kernel_dimension=dim,
          surface_euler_characteristics=eulers,skeleton_chi=chi,
          cells=(vertices=nv,edges=length(edges),surface_faces=nf,meridian_disks=length(meridians)))
    certificate ? (;value=result,certificate=cert) : result
end
