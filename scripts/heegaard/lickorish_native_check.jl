# Source: docs/specs/heegaard-lickorish-data.md sections 3, 4, 6.4;
# src/virtual_link.jl::{set_data!,pd_code,regions} and
# src/heegaard/reference_surface.jl::lickorish_reference.
# JSON-lines worker: retain loaded Julia state across multiple checked inputs.
using Topology, JSON3

function native_check(d)
    g = Int(d.genus)
    names = String.(collect(keys(d.words)))
    words = Dict(String(k) => Int[x for x in v] for (k,v) in pairs(d.words))
    labels = sort!(unique(abs.(reduce(vcat, values(words); init=Int[]))))
    compact = Dict(v=>i for (i,v) in enumerate(labels))
    v = VirtualLink()
    # Section 3.1 / ld_check.py::nonempty_words: empty curves have placement
    # metadata and do not contribute to the active crossing rotation system.
    set_data!(v; gauss=Vector{Int}[Int[sign(x)*compact[abs(x)] for x in words[k]]
                                  for k in names if !isempty(words[k])],
              orientation=ones(Int,length(labels)))
    pd = pd_code(v)
    rs = regions(pd)
    # Connectivity is input validation, as in ld_check.py::surface. Regions
    # themselves are computed only by the existing VirtualLink implementation.
    parent = collect(1:length(labels))
    root(i) = parent[i]==i ? i : (parent[i]=root(parent[i]))
    for word in v.gauss
        isempty(word) && continue
        for label in word
            parent[root(abs(label))]=root(abs(first(word)))
        end
    end
    components = length(Set(root(i) for i in eachindex(parent)))
    if g == 0
        return Dict("vertices"=>length(labels), "faces"=>length(rs), "chi"=>length(rs)-length(labels),
                    "components"=>components, "constants"=>Dict(), "regions"=>rs)
    end
    ref = Topology.Heegaard.lickorish_reference(g)
    # The reference owns K cycles and crossing signs. Label numbering is the
    # spec order: lm_1..lm_g, then the two ln crossings for each n_i.
    ids = ["lm_$i" for i in 1:g]
    for i in 1:g-1
        append!(ids,["ln_$(i)_$i","ln_$(i+1)_$i"])
    end
    ix = Dict(id=>i for (i,id) in enumerate(ids))
    crosses = Dict(c["id"]=>c for c in ref["crossings"])
    constants = Dict{String,Vector{Int}}()
    for curve in ref["curves"]
        name = curve["id"]
        constants[name] = Int[]
        for visit in split(curve["word"]," | ")
            id = split(split(visit,"{")[1],"V_";limit=2)[2]
            cr = crosses[id]
            passage = (name==cr["curves"][1] ? -1 : 1)*cr["sign"]
            push!(constants[name],passage*ix[id])
        end
    end
    Dict("vertices"=>length(labels), "faces"=>length(rs),
         "chi"=>length(rs)-length(labels), "components"=>components, "constants"=>constants,
         "pd"=>[collect(row) for row in eachrow(pd)], "regions"=>rs)
end

println("READY"); flush(stdout)
for line in eachline(stdin)
    isempty(strip(line)) && continue
    try
        println(JSON3.write(native_check(JSON3.read(line))))
    catch error
        println(JSON3.write(Dict("error"=>sprint(showerror,error))))
    end
    flush(stdout)
end
