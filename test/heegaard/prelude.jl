# Shared prelude of the Heegaard suite: imports and fixtures used by its files.
using Test
using Topology
using JSON

const H = Topology.Heegaard

function s3_fixture()
    code = H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)])
    H.certify_heegaard(H.decorate(code, [[1]]; metadata=Dict("fixture"=>"s3")))
end

function lens_fixture(p)
    code = H.heegaard_code([collect(1:p)], [collect(1:p)],
        [H.Intersection(i, 1, 1, 1) for i in 1:p])
    H.certify_heegaard(H.decorate(code, [[i] for i in 1:p]; metadata=Dict("fixture"=>"lens")))
end
