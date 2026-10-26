# Include after MarkedSurfaceDiagrams.jl. Integer grid coordinates below are
# AUTHORING indices, not saved drawing positions and not a shortest-path request.
module MarkedDiagramBuilders
using ..MarkedSurfaceDiagrams
export grid_document, waypoint_word, add_curve!, declare_crossing!
vertex(x,y)="v$(x)_$(y)"
function grid_document(n::Int=4;surface::Symbol=:torus,diagonal::Symbol=:ascending)
    n>=2 || throw(ArgumentError("n must be at least 2"))
    surface in (:torus,:annulus,:disk) || throw(ArgumentError("unknown surface"))
    diagonal in (:ascending,:descending) || throw(ArgumentError("unknown diagonal"))
    sides=Dict("south"=>[vertex(x,0) for x in 0:n],"east"=>[vertex(n,y) for y in 0:n],
      "north"=>[vertex(x,n) for x in n:-1:0],"west"=>[vertex(0,y) for y in n:-1:0])
    triangles=[]
    for y in 0:n-1,x in 0:n-1
        a,b,c,d=vertex(x,y),vertex(x+1,y),vertex(x+1,y+1),vertex(x,y+1)
        append!(triangles,diagonal==:ascending ? [[a,b,c],[a,c,d]] : [[a,b,d],[b,c,d]])
    end
    gl=[];bd=[]
    if surface in (:torus,:annulus);push!(gl,[["P","east"],["P","west"]]);end
    if surface==:torus
        push!(gl,[["P","north"],["P","south"]])
    elseif surface==:annulus
        append!(bd,[Dict("side"=>["P","south"],"label"=>"inner"),Dict("side"=>["P","north"],"label"=>"outer")])
    else
        append!(bd,[Dict("side"=>["P",s],"label"=>"outer") for s in ["south","east","north","west"]])
    end
    Dict{String,Any}("version"=>1,"surface"=>Dict("marking_id"=>"$(surface)-square-v1","genus"=>(surface==:torus ? 1 : 0),"gluings"=>gl,"boundary"=>bd),
      "patches"=>[Dict("face"=>"P","side_order"=>["south","east","north","west"],"sides"=>sides,"triangles"=>triangles)],"curves"=>[],"crossings"=>[])
end
function waypoint_word(blocks;face::String="P")
    # blocks=[[(x0,y0),...],[(...),...]]; each seam jump is an explicit block break.
    words=String[]
    for waypoints in blocks
        length(waypoints)>=2 || throw(ArgumentError("at least two waypoints per block"));vs=String[]
        for (a,b) in zip(waypoints[1:end-1],waypoints[2:end])
            dx,dy=b[1]-a[1],b[2]-a[2];n=max(abs(dx),abs(dy))
            n>0 && (dx==0 || dy==0 || abs(dx)==abs(dy)) || throw(ArgumentError("explicit segments must follow grid edges"))
            sx,sy=sign(dx),sign(dy)
            seg=[vertex(a[1]+k*sx,a[2]+k*sy) for k in 0:n]
            append!(vs,isempty(vs) ? seg : seg[2:end])
        end
        push!(words,face*"{"*join(vs," ")*"}")
    end
    join(words," | ")
end
function add_curve!(d,id,role,closed,blocks;face="P")
    push!(d["curves"],Dict("id"=>id,"role"=>role,"closed"=>closed,"word"=>waypoint_word(blocks;face=face)));d
end
function declare_crossing!(d,id,x,y,first,second,sign;face="P")
    push!(d["crossings"],Dict("id"=>id,"at"=>[face,vertex(x,y)],"curves"=>[first,second],"sign"=>sign));d
end
end
