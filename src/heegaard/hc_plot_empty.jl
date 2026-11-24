# cf: lib/topology.jl/scripts/heegaard/lc_empty.py add_empty — circles in named rotation faces
function _hp_empty!(runs,lc,roles,arcs,intervals,map_point,triangles,holes,W,H,square,names)
    words=lc.words;empty=[c for c in names if isempty(words[c])]
    isempty(empty) && return
    for a in arcs
        labels=vcat(a.first,a.mids)
        for (j,label) in enumerate(labels)
            signed=only(x for x in words[a.curve] if abs(x)==label)
            count=a.face=="Finf" ? 300 : 1;model=a.eventpaths[j]
            path=[map_point(a.face,_hp_lerp(model[i],model[i+1],k/count)) for i in 1:length(model)-1 for k in 0:count-1]
            push!(path,map_point(a.face,model[end]));push!(get!(intervals,(a.curve,signed),Vector{HCPlotPoint}[]),path)
        end
    end
    segments=Tuple{HCPlotPoint,HCPlotPoint}[(p[i],p[i+1]) for paths in values(intervals) for p in paths for i in 1:length(p)-1]
    for r in runs
        append!(segments,[(r.points[i],r.points[i+1]) for i in 1:length(r.points)-1])
        r.closed && push!(segments,(r.points[end],r.points[1]))
    end
    faces,_=_lc_faces(words);face_map=Dict(f.id=>f.corners for f in faces);occ=_lc_rotation_occ(words)
    circles=Dict{String,Vector{HCPlotPoint}}();disks=Dict{String,Tuple{HCPlotPoint,Float64}}();inside=Dict{String,Set{String}}()
    required=Set{String}()
    allowed(p)=all((_hp_dist(p,c)<r)==(name in required) for (name,(c,r)) in disks)
    function free_disk(paths)
        best=-Inf;bestp=(0.0,0.0)
        for path in paths,i in 1:length(path)-1
            a,b=path[i:i+1];len=_hp_dist(a,b);len<1e-7 && continue
            for t in (0.5,0.38196601125,0.61803398875)
                p=_hp_lerp(a,b,t);clearance=min(len*t,len*(1-t))
                for (c,d) in segments
                    dist=_hp_nearest(p,c,d)[1];dist>1e-7 && (clearance=min(clearance,dist))
                end
                offset=clearance/2;c=(p[1]-(b[2]-a[2])*offset/len,p[2]+(b[1]-a[1])*offset/len)
                allowed(c) || continue
                square && !(0<c[1]<W && 0<c[2]<H) && continue
                any(_hp_dist(c,h)<=r for (_, (h,r)) in holes) && continue
                in_triangle=false
                for tri in triangles
                    cross=[_hp_cross(_hp_sub(tri[mod1(k+1,3)],tri[k]),_hp_sub(c,tri[k])) for k in 1:3]
                    if minimum(cross)>=0 || maximum(cross)<=0;in_triangle=true;break;end
                end
                in_triangle && continue
                radius=offset*.8
                for (a,b) in segments;radius=min(radius,_hp_nearest(c,a,b)[1]*.8);end
                if radius>1e-7 && radius>best;best=radius;bestp=c;end
            end
        end
        best>0 || throw(ArgumentError("no nondegenerate local disk in named face"))
        bestp,best
    end
    function face_paths(fid)
        paths=Vector{HCPlotPoint}[]
        for (label,sign,way) in face_map[fid]
            c,i=occ[label][sign];w=words[c]
            if way=="out";append!(paths,get(intervals,(c,w[i+1]),Vector{HCPlotPoint}[]))
            else;append!(paths,reverse.(get(intervals,(c,w[mod1(i,length(w))]),Vector{HCPlotPoint}[])));end
        end
        paths
    end
    pending=copy(empty)
    while !isempty(pending)
        advanced=false
        for c in copy(pending)
            ref=lc.places[c];parent=(endswith(ref,"+") || endswith(ref,"-")) ? ref[1:end-1] : ""
            parent in empty && !haskey(circles,parent) && continue
            required=copy(get(inside,parent,Set{String}()))
            haskey(circles,parent) && endswith(ref,"+") && push!(required,parent)
            if ref=="outer"
                center,radius=isempty(circles) ? ((260.0,260.0),160.0) : free_disk(reverse.(collect(values(circles))))
            elseif haskey(face_map,ref)
                center,radius=free_disk(face_paths(ref))
            elseif haskey(circles,parent)
                path=circles[parent]
                if endswith(ref,"+")
                    siblings=[name for name in empty if lc.places[name]==ref];origin,extent=disks[parent];n=length(siblings)
                    angle=2pi*(findfirst(==(c),siblings)-1)/n
                    center=n==1 ? origin : _hp_add(origin,(extent*.5*cos(angle),extent*.5*sin(angle)))
                    radius=min(n==1 ? extent/3 : extent*sin(pi/n)/3,minimum(_hp_nearest(center,a,b)[1] for (a,b) in segments)*.8)
                    if radius<1e-7 || !allowed(center);center,radius=free_disk([path]);end
                else;center,radius=free_disk([reverse(path)]);end
            elseif haskey(words,parent) && !isempty(words[parent])
                paths=[p for ((name,_),parts) in intervals if name==parent for p in parts]
                center,radius=free_disk(endswith(ref,"+") ? paths : reverse.(paths))
            else;continue;end
            pts=[_hp_add(center,(radius*cos(2pi*k/64),radius*sin(2pi*k/64))) for k in 0:63]
            push!(runs,HCPlotRun(c,c,get(roles,c,:none),true,pts,Int[]))
            disks[c]=(center,radius);inside[c]=copy(required);circles[c]=vcat(pts,[pts[1]])
            append!(segments,[(circles[c][i],circles[c][i+1]) for i in 1:64])
            deleteat!(pending,findfirst(==(c),pending));advanced=true
        end
        advanced || throw(ArgumentError("unsupported placement references: $pending"))
    end
end
