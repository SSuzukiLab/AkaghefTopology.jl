using Test
using Topology

# Independent transcription of Taguchi's MP_list_AB_int.jpg and
# MP_list_CD_int.jpg, compared with MST25 Z_framed_MP_move.pdf. The split lab
# archive has framing.tex in zip-003, lines 188-203 (definition and figure
# labels), and the original figures in zip-004. This oracle does not call
# _reverse_template/_mp_patterns.
# A/B fragment 1 follows lower-left -> lower -> right -> upper -> upper-left;
# fragment 2 is horizontal; fragment 3 is vertical. C/D fragment 1 is horizontal,
# fragment 2 joins the left boundary to the right boundary around the bend,
# and fragment 3 is vertical. The ordered words below specify all passage signs.
# The first/last rule entries are boundary corrections; interior entries are
# required edge weights. In particular these are INTEGER weights, not mod 2.
const MP_SOURCE_ROWS = [
    ("A1", [[-1],[1,2],[-2]], [1,1], [[1,-3,2],[3],[-1,-2]], [1,-1,1]),
    ("A2", [[-1],[1,2],[-2]], [1,-1], [[1,-3,2],[3],[-2,-1]], [-1,1,1]),
    ("A3", [[-1],[2,1],[-2]], [-1,-1], [[1,-3,2],[3],[-1,-2]], [1,-1,-1]),
    ("A4", [[-1],[2,1],[-2]], [-1,1], [[1,-3,2],[3],[-2,-1]], [-1,1,-1]),
    ("B1", [[1],[-1,-2],[2]], [-1,-1], [[-1,3,-2],[-3],[1,2]], [-1,1,-1]),
    ("B2", [[1],[-1,-2],[2]], [-1,1], [[-1,3,-2],[-3],[2,1]], [1,-1,-1]),
    ("B3", [[1],[-2,-1],[2]], [1,1], [[-1,3,-2],[-3],[1,2]], [-1,1,1]),
    ("B4", [[1],[-2,-1],[2]], [1,-1], [[-1,3,-2],[-3],[2,1]], [1,-1,1]),
    ("C1", [[-1,2],[1],[-2]], [-1,1], [[2,-3],[1,3],[-1,-2]], [1,1,-1]),
    ("C2", [[-1,2],[1],[-2]], [-1,-1], [[2,-3],[1,3],[-2,-1]], [-1,-1,-1]),
    ("C3", [[-1,2],[1],[-2]], [1,1], [[2,-3],[3,1],[-1,-2]], [-1,1,1]),
    ("C4", [[-1,2],[1],[-2]], [1,-1], [[2,-3],[3,1],[-2,-1]], [1,-1,1]),
    ("D1", [[2,-1],[1],[-2]], [1,-1], [[-3,1],[3,2],[-1,-2]], [-1,-1,1]),
    ("D2", [[2,-1],[1],[-2]], [1,1], [[-3,1],[3,2],[-2,-1]], [1,1,1]),
    ("D3", [[2,-1],[1],[-2]], [-1,-1], [[-3,1],[2,3],[-1,-2]], [-1,1,-1]),
    ("D4", [[2,-1],[1],[-2]], [-1,1], [[-3,1],[2,3],[-2,-1]], [1,-1,-1]),
]

const MP_SOURCE_NONZERO = Dict(
    "A2" => (1,[0,0,-1,1]), "B2" => (1,[0,0,-1,1]),
    "A4" => (1,[1,-1,0,0]), "B4" => (1,[1,-1,0,0]),
    "C3" => (2,[0,1,-1]), "D3" => (2,[-1,1,0]),
)

@testset "MP all 16 source diagrams and signed weight positions" begin
    data=make_move_data()
    for (name,left,lo,right,ro) in MP_SOURCE_ROWS
        lhs=only(filter(p -> p.parameter==name,data["MP-L"]))
        rhs=only(filter(p -> p.parameter==name,data["MP-R"]))
        expected=[zeros(Int,length(w)+1) for w in right]
        if haskey(MP_SOURCE_NONZERO,name)
            fragment,rule=MP_SOURCE_NONZERO[name]
            expected[fragment]=rule
        end
        @test lhs.gauss==left
        @test lhs.orientation==lo
        @test lhs.weights==[zeros(Int,length(w)+1) for w in left]
        @test rhs.gauss==right
        @test rhs.orientation==ro
        @test rhs.weights==expected
    end
end

@testset "A2/B2 former assignment differs by one explicit H move" begin
    for (name,words,orientation) in (
        ("A2",[[1,-3,2],[3],[-2,-1]],[-1,1,1]),
        ("B2",[[-1,3,-2],[-3],[2,1]],[1,-1,-1]))
        old=VirtualLink()
        set_data!(old;rgauss=words,orientation=orientation)
        # Old +1 on incoming vertical leg, -1 on its internal edge.
        set_weight!(old,[[7,0,0],[-11],[14,-1]])
        move!(old,"H";param=1,vertices=[2])
        # Direct source +1 on outgoing U-strand leg, -1 on upper-right arc.
        @test old.weights==[[8,0,-1],[-11],[13,0]]
    end
end

# Compare decorated oriented edges, avoiding an arbitrary cyclic word origin.
# All labels and components remain significant: no canonicalization is used.
function mp_test_edges(v)
    sort([(e.strand,e.crossing,e.weight) for e in real_edge_table(v)])
end

# Finite seeds from the existing ograph-c2-area2-n5-260907-v3.sqlite3 states.
# These literal INITIAL solutions were found once with calculate_weight; no
# test, MP application, or inverse below solves for replacement weights.
const MP_FRAMED_SEEDS = Dict(
    3 => ([-1,1,-2,2],[1,1],[1,-1,0,-1]),
    4 => ([-1,1,-2,2],[-1,1],[-2,1,-1,1]),
    5 => ([-1,1,-2,2],[-1,-1],[-2,1,-1,1]),
    6 => ([-1,2,-2,1],[1,1],[1,-2,1,-1]),
    7 => ([-1,2,-2,1],[-1,1],[1,-2,1,-1]),
    8 => ([-1,2,-2,1],[-1,-1],[-1,1,-1,0]),
    9 => ([-1,-2,2,1],[1,1],[1,-1,-1,0]),
    11 => ([-1,-2,2,1],[-1,-1],[-1,-1,1,0]),
    26 => ([-1,-2,1,2,-3,3],[-1,1,-1],[-2,0,1,0,-1,1]),
    27 => ([-1,-2,1,2,-3,3],[1,-1,-1],[-2,0,1,0,-1,1]),
)

# name, seed id, ordered local vertices, explicit H vertex, H coefficient.
# H changes the unique LHS internal edge to the required weight zero.
const MP_FRAMED_CASES = [
    ("A1",9,[2,1],2,0), ("A2",27,[1,2],1,0),
    ("A3",11,[1,2],2,0), ("A4",27,[2,1],1,0),
    ("B1",11,[1,2],1,1), ("B2",26,[1,2],1,0),
    ("B3",9,[2,1],1,1), ("B4",26,[2,1],1,0),
    ("C1",7,[1,2],1,2), ("C2",8,[1,2],1,-1),
    ("C3",6,[1,2],1,2), ("C4",7,[2,1],2,1),
    ("D1",4,[2,1],1,1), ("D2",3,[1,2],2,-1),
    ("D3",5,[1,2],2,2), ("D4",4,[1,2],2,2),
]

function mp_test_disk_residual(v)
    # Reconstruct the equation by individual signed boundary occurrences,
    # retaining repeated incidences, instead of calling weight_residual or
    # multiplying the production aggregated delta matrix. The boundary/dot
    # geometry still comes from the existing disk_table implementation.
    disks=disk_table(v)
    weights=BigInt[e.weight for e in real_edge_table(v)]
    [2sum(sign(edge)*weights[abs(edge)] for edge in circuit)+2-dots
     for (circuit,dots) in zip(disks.edges_s,disks.Ndots)]
end

@testset "MP all 16 closed-normal framing representatives" begin
    @test Set(first.(MP_FRAMED_CASES))==Set(first.(MP_SOURCE_ROWS))
    for (name,id,vertices,hvertex,coefficient) in MP_FRAMED_CASES
        word,orientation,weights=MP_FRAMED_SEEDS[id]
        source=VirtualLink()
        set_data!(source;rgauss=[word],orientation=orientation)
        set_weight!(source,weights)
        original=deepcopy(source)
        @test is_closed(source)==(true,(true,true,true))
        @test all(iszero,mp_test_disk_residual(source))
        @test validate_weight(source;framed=true)

        # Independently predict the H correction from its signed incidences,
        # including loops whose incoming/outgoing contributions cancel.
        expected_h=[e.weight+coefficient*((abs(e.crossing[1])==hvertex)-
                    (abs(e.crossing[2])==hvertex)) for e in real_edge_table(source)]
        move!(source,"H";param=coefficient,vertices=[hvertex])
        @test vcat(source.weights...)==expected_h
        @test all(iszero,mp_test_disk_residual(source))
        candidates=movable(source,"MP";param=name,vertices=vertices,strict=true)
        @test length(candidates)==1
        isempty(candidates) && continue
        candidate=only(candidates)
        fragment=only(findall(w -> length(w)==2,candidate.before.gauss))
        component,start=candidate.locations[fragment]
        @test source.weights[component][mod1(start+1,length(source.gauss[component]))]==0

        target=move(source,candidate)
        @test is_closed(target)==(true,(true,true,true))
        @test all(iszero,mp_test_disk_residual(target))
        @test validate_weight(target;framed=true)
        @test all(isinteger,vcat(target.weights...))
        inverse=only(movable(target,"MP";param=name,
            vertices=vcat(vertices,length(orientation)+1),strict=true))
        restored=move(target,inverse)
        @test is_closed(restored)==(true,(true,true,true))
        @test all(iszero,mp_test_disk_residual(restored))
        @test restored.orientation==source.orientation
        @test mp_test_edges(restored)==mp_test_edges(source)
        move!(restored,"H";param=-coefficient,vertices=[hvertex])
        @test all(iszero,mp_test_disk_residual(restored))
        @test mp_test_edges(restored)==mp_test_edges(original)
    end
end

@testset "MP source rules through public candidates in both directions" begin
    for (name,left,lo,right,ro) in MP_SOURCE_ROWS
        rules=[zeros(Int,length(w)+1) for w in right]
        if haskey(MP_SOURCE_NONZERO,name)
            fragment,rule=MP_SOURCE_NONZERO[name]
            rules[fragment]=rule
        end
        # Three independent closures, or one cyclic closure with every possible
        # source origin. These closures test integral transport; they make no
        # claim that the chosen arbitrary boundary weights define a framing.
        for separate in (true,false), root in (separate ? (0:0) : (0:3))
            free=[7,-11,13]
            source_weights=[vcat(free[i],zeros(Int,length(left[i])-1)) for i in 1:3]
            target_weights=[vcat(free[i]+rules[i][1],rules[i][2:end-1]) for i in 1:3]
            for i in 1:3
                next=separate ? i : mod1(i+1,3)
                target_weights[next][1]+=rules[i][end]
            end
            source_code=separate ? deepcopy(left) : [circshift(vcat(left...),root)]
            source_weight=separate ? source_weights : [circshift(vcat(source_weights...),root)]
            target_code=separate ? deepcopy(right) : [vcat(right...)]
            target_weight=separate ? target_weights : [vcat(target_weights...)]
            source=VirtualLink()
            set_data!(source;rgauss=source_code,orientation=lo)
            set_weight!(source,source_weight)
            expected=VirtualLink()
            set_data!(expected;rgauss=target_code,orientation=ro)
            set_weight!(expected,target_weight)
            candidates=movable(source,"MP";param=name,vertices=[1,2])
            @test !isempty(candidates)
            isempty(candidates) && continue
            candidate=first(candidates)
            target=deepcopy(source)
            move!(target,candidate)
            @test target.is_weighted
            @test target.orientation==expected.orientation
            @test mp_test_edges(target)==mp_test_edges(expected)
            @test source.gauss==source_code
            @test source.weights==source_weight
            inverse_candidates=movable(target,"MP";param=name,vertices=[1,2,3])
            @test !isempty(inverse_candidates)
            isempty(inverse_candidates) && continue
            inverse=first(inverse_candidates)
            recovered=deepcopy(target)
            move!(recovered,inverse)
            @test recovered.orientation==source.orientation
            @test mp_test_edges(recovered)==mp_test_edges(source)

            # Perturb every prescribed source internal edge in EACH direction.
            # A stale candidate must reject this rather than discard the weight.
            for (state,witness) in ((source,candidate),(target,inverse))
                for fragment in eachindex(witness.before.gauss),
                    offset in 1:length(witness.before.gauss[fragment])-1
                    component,start=witness.locations[fragment]
                    position=mod1(start+offset,length(state.weights[component]))
                    invalid=deepcopy(state)
                    invalid.weights[component][position]+=2
                    snapshot=deepcopy(invalid)
                    @test isempty(movable(invalid,"MP";param=name,vertices=witness.vertices))
                    @test_throws ArgumentError move!(invalid,witness)
                    @test invalid.gauss==snapshot.gauss
                    @test invalid.weights==snapshot.weights
                    @test invalid.is_weighted
                end
            end
        end
    end
end
