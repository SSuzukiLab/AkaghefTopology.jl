using Topology

"""Gauss code used by the two surgery-move loops of `C260125MSTinv_uqsl2.mlx`."""
function smove_gauss_code(count::Integer, mode::Integer=0)
    count>=0 || throw(ArgumentError("move count must be nonnegative"))
    mode in (0,1) || throw(ArgumentError("mode must be 0 or 1"))
    base=[-1,2,3,-3,-4,-2,4,1]
    code=reduce(vcat,[base .+ 4offset.*sign.(base) for offset in 0:count-1];init=Int[])
    if mode==1
        count==0 && return Int[]
        blocks=reshape(code,8,count)
        return vcat(vec(blocks[1:3,:]),vec(blocks[4:end,end:-1:1]))
    end
    code
end

function _surgery_link(count,mode)
    code=smove_gauss_code(count,mode)
    link=VirtualLink()
    set_data!(link;rgauss=[vcat((4count+1).*[1,-1],code)],
              orientation=vcat(repeat([1,1,1,-1],count),1))
    calculate_weight(link;assign=true)
    # MATLAB `calcWeight(true)` chooses this integral representative after
    # its H-move reduction.  These are replayed values for the mode-1 cells.
    if mode==1 && count==2
        set_weight!(link,[[0,1,0,-1,2,1,-1,0,1,-1,-1,0,0,0,1,0,-1,0]])
    elseif mode==1 && count==3
        set_weight!(link,[[0,1,0,-1,0,1,0,0,1,0,0,1,-1,-1,0,-1,0,1,0,1,1,0,1,0,-1,0]])
    end
    link
end

"""Reproduce the numerical surgery/connected-sum cells unique to C260125."""
function run_C260125_uqsl2()
    algebra=uqsl2_borel_small(1//2;style=:L)
    surgery_standard=[hopf_invariant(_surgery_link(count,0),algebra) for count in 0:3]
    surgery_alternate=[hopf_invariant(_surgery_link(count,1),algebra) for count in 0:3]

    positive=VirtualLink(); set_data!(positive;rgauss=[[1,-1]],orientation=[1])
    calculate_weight(positive;assign=true)
    negative=VirtualLink(); set_data!(negative;rgauss=[[1,-1]],orientation=[-1])
    sums=VirtualLink[deepcopy(negative)]
    for _ in 2:3
        next=deepcopy(last(sums)); disjoint_sum!(next,positive); connected_sum!(next)
        calculate_weight(next;assign=true); push!(sums,next)
    end

    alternate=VirtualLink(); set_data!(alternate;rgauss=[[-1,1]],orientation=[1])
    calculate_weight(alternate;assign=true)
    positive_sums=VirtualLink[deepcopy(positive)]
    for _ in 2:3
        next=deepcopy(last(positive_sums)); disjoint_sum!(next,alternate); connected_sum!(next)
        calculate_weight(next;assign=true); push!(positive_sums,next)
    end

    shared=run_C250818_uqsl2()
    (
        hopf_checks=verify_hopf(algebra),
        surgery_standard=surgery_standard,
        surgery_alternate=surgery_alternate,
        abalone_sums=[hopf_invariant(positive,algebra);
                      [hopf_invariant(link,algebra) for link in sums[2:end]]],
        alternate_sums=[(hopf=hopf_invariant(alternate,algebra),ur_prime=ur_prime_invariant(alternate));
                        [(hopf=hopf_invariant(link,algebra),ur_prime=ur_prime_invariant(link))
                         for link in positive_sums[2:end]]],
        shared=shared,
    )
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C260125_uqsl2())
