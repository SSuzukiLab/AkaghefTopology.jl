# cf: docs/specs/solid3d.md — the square chart, role-free words and HTML output.
@testset "solid cut chart" begin
    lc = H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\na : -3 +4 -5 +6\nb : -4 +5 -6 -2\n")
    hc = H.resolve_hc(lc,Dict("a"=>:alpha,"b"=>:beta))
    planar_svg = H.plot(lc).svg
    mktempdir() do dir
        for input in (lc,hc)
            data = H._solid_data(input,true)
            @test data.square
            @test data.kinds["a"] == (input === lc ? :curve : :alpha)
            for (name,word) in lc.words
                @test sort([e.label for e in data.paths if e.name == name && !e.copy_on_cut]) == sort(word)
            end
            path = joinpath(dir,"diagram.html")
            @test H.draw_heegaard_3d(input,path;crossing_labels=true) == path
            html = read(path,String)
            @test occursin("data-labels=\"true\"",html)
            @test !occursin("__SOLID_CHART__",html)
        end
    end
    @test H.plot(lc).svg == planar_svg
end
# cf: docs/specs/solid3d.md — arclength smoothing pins ports and rounds a corner.
@testset "solid curve smoothing" begin
    pts = [(0.,0.), (60.,0.), (60.,60.)]
    r = H.HCPlotRun("c","c",:none,false,pts,Int[])
    q = H._solid_smooth_run(r,18.).points
    @test first(q) == first(pts) && last(q) == last(pts)
    @test minimum(H._hp_dist(v,pts[2]) for v in q) > 4
    turns = [begin
        a = H._hp_sub(q[i],q[i-1])
        b = H._hp_sub(q[i+1],q[i])
        acosd(clamp(H._hp_dot(a,b)/(H._hp_norm(a)*H._hp_norm(b)),-1,1))
    end for i in 2:length(q)-1]
    @test maximum(turns) < 12
    @test r.points == pts
end
