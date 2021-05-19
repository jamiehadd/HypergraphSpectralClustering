using Revise
using HypergraphDetectability
using SparseArrays
using LinearAlgebra
using Arpack
using Test
##


## 
# get some fake data to play with
n  = 100
c₂ = 10
c₃ = 10
p₂ = 0.9
p₃ = 0.7

H = detectabilityData(n, c₂, c₃, p₂, p₃);


@testset "nonbacktracking matrix" begin
"""
compare functions for computing the entire nonBacktrackingMatrix 
and the version broken down by size

since the edges may be not be consistently indexed, we do this by just comparing the top two eigenvalues. 
"""

    # "ground truth": full nonbacktracking matrix
    B = nonBacktrackingMatrix(H);
    E = eigs(B; nev = 2);

    # partitioned up by edge size
    Bs = nonBacktrackingMatrices(H);
    B_ = sum(Bs)
    E_ = eigs(B_; nev = 2)

    @test E[1] ≈ E_[1]

    # partitioned by edge size and returning edge indices
    Bs, ix = nonBacktrackingMatrices(H; return_indices = true);
    B_ = sum(Bs)
    E_ = eigs(B_; nev = 2)

    @test E[1] ≈ E_[1]

end


@testset "compute clusters" begin
"""
compare functions for computing the entire nonBacktrackingMatrix 
and the version broken down by size

since the edges may be not be consistently indexed, we do this by just comparing the top two eigenvalues. 
"""


    # partitioned by edge size and returning edge indices
    Bs, ix = nonBacktrackingMatrices(H; return_indices = true);

    B = sum(Bs)
    E = eigs(B; nev = 2, ritzvec = true)

    v = E[2][:,2]

    u = aggregateEigenvector(v, ix)

    # sign of u should correspond to clusters, should be 50 in each one
    @test sum(u .> 0) == 50

    # packages up the above computations 
    # starting from the computation of the combined
    # matrix B
    u_ = computeBinaryClusters(B, ix)
end

