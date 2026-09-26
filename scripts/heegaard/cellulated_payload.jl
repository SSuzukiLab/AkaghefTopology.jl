# Shared by the app endpoint and fixture build. No automatic routing dependency.
if !isdefined(@__MODULE__, :MarkedSurfaceDiagrams)
    include("../../vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl")
end
const MSD = MarkedSurfaceDiagrams
const CELL_PROJECT = "akaghef.heegaard.cellulated-project/1"
const CELL_LAYOUT = "akaghef.heegaard.cellulated-layout/1"

function cellulated_payload(input::Dict; layout=nothing, id="imported", title="Marked cellulated diagram")
    if get(input,"schema",nothing)==CELL_PROJECT
        layout === nothing || error("A project already supplies its layout; do not override it implicitly")
        layout=input["layout"]; input=input["diagram"]
    end
    st=MSD.validate(input)
    digest=MSD.representation_digest(st.doc)
    if layout === nothing
        pos=MSD.realize(st)
    else
        get(layout,"schema",nothing)==CELL_LAYOUT || error("CELL_LAYOUT_SCHEMA")
        get(layout,"representation_sha256",nothing)==digest || error("CELL_LAYOUT_HASH_MISMATCH")
        pos=Dict(f=>Dict(v=>Tuple(MSD.qnum.(p)) for (v,p) in ps) for (f,ps) in layout["positions"])
        MSD.certify_layout(st,pos)
        # Fixed chart boundary; saved interior vertex adjustments are allowed.
        for (f,bd) in st.boundaries
            all(pos[f][v]==p for (v,p) in zip(bd,MSD.circle_boundary(length(bd)))) || error("CELL_LAYOUT_BOUNDARY_CHANGED")
        end
    end
    report=MSD.report(st,pos)
    report["joining_disks"]=[merge(Dict("id"=>m["id"],"claim"=>m["claim"]),MSD.check_joining_disk(st,m["seed"])) for m in get(st.doc,"region_marks",[])]
    svg,certificate=mktempdir() do dir
        path=joinpath(dir,"diagram.svg")
        cert=MSD.write_svg(st,pos,path)
        read(path,String),cert
    end
    report["svg_certificate"]=certificate
    saved=Dict("schema"=>CELL_LAYOUT,"representation_sha256"=>digest,"positions"=>report["rational_positions"])
    result=Dict("view"=>"cellulated","id"=>id,"title"=>title,"diagram"=>st.doc,
        "layout"=>saved,"report"=>report,"svg"=>svg,"genus"=>st.surface_report["genus"],
        "intersections"=>length(st.doc["crossings"]))
    # Optional locally audited source projection, bound to the WHOLE mathematical
    # representation. Never trust imported metadata as an SVG/route witness.
    witness=joinpath(@__DIR__,"../../artifacts/heegaard/cellulated/projections",digest*".json")
    if isfile(witness)
        projection=JSON3.read(read(witness,String),Dict{String,Any})
        projection["representation_sha256"]==digest || error("projection digest mismatch")
        result["source_projection"]=projection
    end
    result
end
