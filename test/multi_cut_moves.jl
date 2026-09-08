using Test
using Topology

@testset "closed BMP retains cyclic reconnection compatibility" begin
    cyclic_key(words)=sort([isempty(w) ? () : minimum(Tuple(circshift(w,k)) for k in 0:length(w)-1) for w in words])
    for words in ([[1,2],[-1],[-2]],[[1,2],[-1,-2]]), p in (0,1)
        source=VirtualLink(); set_data!(source;rgauss=words,orientation=[-1,1])
        forward=move_boundary_transport(source,"BMP";param=p,vertices=[1,2],strict=true)
        reverse=move_boundary_transport(forward.target,"BMP";param=p,vertices=[1,2,3],strict=true)
        @test cyclic_key(reverse.target.gauss)==cyclic_key(source.gauss)
        @test reverse.target.orientation==source.orientation
        @test isempty(forward.endpoints) && isempty(reverse.endpoints)
    end
end

# These are source-side fixtures.  Their expected targets below are written
# from the MATLAB fragment/capacity construction, rather than obtained by
# running the Julia move and copying its result.
function _bmp_different_opposite_fixture()
    v=VirtualLink()
    set_data!(v;rgauss=[[1,2,0],[-1,0],[-2,0]],orientation=[-1,1])
    v
end

function _bmp_same_opposite_fixture()
    v=VirtualLink()
    set_data!(v;rgauss=[[1,2,0],[-1,-2,0]],orientation=[1,-1])
    v
end

function _mp_three_cut_fixture(;weighted=false)
    v=VirtualLink()
    set_data!(v;rgauss=[[-1,2,0],[1,0],[-2,0]],orientation=[1,1])
    weighted && set_weight!(v,[[11,0,13],[21,22],[31,32]])
    v
end

_endpoint_key(endpoint)=(endpoint.strand,endpoint.side,endpoint.lane)

function _boundary_inventory(v)
    rows=circuit_b_boundary(v)
    sort!(vcat([_endpoint_key(row.start) for row in rows],
               [_endpoint_key(row.stop) for row in rows]))
end

function _transport_inventory(transport)
    sort!([(row.source.strand,row.source.side,row.target.strand,row.target.side)
           for row in transport.endpoints])
end

const _BMP_DIFF_ENDPOINTS = [
    (1,:in,1,:in), (1,:out,1,:out),
    (2,:in,2,:in), (2,:out,3,:out),
    (3,:in,3,:in), (3,:out,2,:out),
]

const _MP_ENDPOINTS = [
    (1,:in,1,:in), (1,:out,1,:out),
    (2,:in,2,:in), (2,:out,2,:out),
    (3,:in,3,:in), (3,:out,3,:out),
]

@testset "BMP multi-cut forward preserves every endpoint" begin
    source=_bmp_different_opposite_fixture()
    source_inventory=_boundary_inventory(source)
    @test length(source_inventory)==18

    expected_by_parameter=Dict(
        0 => ([[3,0],[2,0],[1,0],[-2,-1,-3]],[-1,1,1]),
        1 => ([[3,0],[2,0],[1,0],[-1,-2,-3]],[1,-1,-1]),
    )
    for parameter in (0,1)
        transport=move_boundary_transport(source,"BMP";
            param=parameter,vertices=[1,2],strict=true)
        target=transport.target
        expected_code,expected_orientation=expected_by_parameter[parameter]
        @test target.gauss==expected_code
        @test target.orientation==expected_orientation
        @test strand_cuts(target)==[true,true,true,false]
        @test _transport_inventory(transport)==sort(_BMP_DIFF_ENDPOINTS)
        @test _boundary_inventory(target)==source_inventory

        inverse=move_boundary_transport(target,"BMP";
            param=parameter,vertices=[1,2,3],strict=true)
        @test inverse.target.gauss==source.gauss
        @test inverse.target.orientation==source.orientation
        @test _boundary_inventory(inverse.target)==source_inventory
        inverse_map=Dict(row.source=>row.target for row in inverse.endpoints)
        @test all(inverse_map[row.target]==row.source for row in transport.endpoints)
    end
end

@testset "B02 preserves distinct empty open intervals" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[1,0],[2,0],[-1,-2]],orientation=[1,-1])
    transport=move_boundary_transport(source,"B02";vertices=[1,2],strict=true)
    @test transport.target.gauss==[[0],[0]]
    @test isempty(transport.target.orientation)
    @test length(real_edge_table(transport.target))==2
    @test length(transport.endpoints)==4
    @test length(unique(row.target for row in transport.endpoints))==4
    @test _boundary_inventory(transport.target)==_boundary_inventory(source)
end

@testset "BMP signed same-component opposite visits" begin
    source=_bmp_same_opposite_fixture()
    expected_by_parameter=Dict(
        0 => ([[-2,0],[-3,0],[-1],[2,1,3]],[1,-1,-1]),
        1 => ([[-2,0],[-3,0],[-1],[1,2,3]],[-1,1,1]),
    )
    for parameter in (0,1)
        transport=move_boundary_transport(source,"BMP";
            param=parameter,vertices=[-1,-2],strict=true)
        expected_code,expected_orientation=expected_by_parameter[parameter]
        @test transport.target.gauss==expected_code
        @test transport.target.orientation==expected_orientation
        @test strand_cuts(transport.target)==[true,true,false,false]
        @test _boundary_inventory(transport.target)==_boundary_inventory(source)

        inverse=move_boundary_transport(transport.target,"BMP";
            param=parameter,vertices=[-1,-2,-3],strict=true)
        @test inverse.target.gauss==source.gauss
        @test inverse.target.orientation==source.orientation
    end
end

@testset "MP fragment/capacity recombination over three cuts" begin
    source=_mp_three_cut_fixture()
    transport=move_boundary_transport(source,"MP";
        param="C3",vertices=[1,2],strict=true)
    target=transport.target
    @test target.gauss==[[2,-3,0],[3,1,0],[-1,-2,0]]
    @test target.orientation==[-1,1,1]
    @test strand_cuts(target)==[true,true,true]
    @test _transport_inventory(transport)==sort(_MP_ENDPOINTS)
    @test _boundary_inventory(target)==_boundary_inventory(source)

    inverse=move_boundary_transport(target,"MP";
        param="C3",vertices=[1,2,3],strict=true)
    @test inverse.target.gauss==source.gauss
    @test inverse.target.orientation==source.orientation
    @test _boundary_inventory(inverse.target)==_boundary_inventory(source)
end

@testset "weighted MP carries m+1 portions at every open endpoint" begin
    source=_mp_three_cut_fixture(weighted=true)
    transport=move_boundary_transport(source,"MP";
        param="C3",vertices=[1,2],strict=true)
    target=transport.target
    @test target.gauss==[[2,-3,0],[3,1,0],[-1,-2,0]]
    # Capacity is unchanged: the second fragment's outgoing rule is -1,
    # so its own boundary portion 22 becomes 21; the third keeps 31 and 32.
    @test target.weights==[[11.0,0.0,13.0],[21.0,1.0,21.0],[31.0,0.0,32.0]]
    @test _boundary_inventory(target)==_boundary_inventory(source)

    inverse=move_boundary_transport(target,"MP";
        param="C3",vertices=[1,2,3],strict=true)
    @test inverse.target.gauss==source.gauss
    @test inverse.target.weights==source.weights
    @test inverse.target.orientation==source.orientation
end
