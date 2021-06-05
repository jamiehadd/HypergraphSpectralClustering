"""
unit tests
to run all tests:

>>> cd HypergraphDetectability
>>> julia

julia> ]

pkg> activate .
pkg> test

"""

using Revise
using Test
using HypergraphDetectability

# using SparseArrays
# using LinearAlgebra
using Arpack
using Statistics
##


## 
# get some fake data to play with
n  = 100
c₂ = 10
c₃ = 10
p₂ = 0.8
p₃ = 0.8

H = detectabilityData(n, c₂, c₃, p₂, p₃);

z = 1 .+ (1:n .> n/2);


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
    Bs, ix = nonBacktrackingMatrices(H);
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
    Bs, ix = nonBacktrackingMatrices(H);

    B = sum(Bs)
    E = eigs(B; nev = 2, ritzvec = true)

    v = E[2][:,2]

    u = aggregateEigenvector(v, ix)

    # sign of u should correspond to clusters, should be 50 in each one
    # random, so not a great test all things considered
    @test sum(u .> 0) >= n/2 - 10 # should be exactly n/2, but close is ok

    # packages up the above computations 
    # starting from the computation of the combined
    # matrix B
    # z = binaryClusters(B, ix)

end

@testset "degree tensor" begin
    """
    test of a combinatorial identity in the estimated degree tensor
    """

    c, C = degreeTensor(H, z)

    q = 1/n * [sum(z .== i) for i in unique(z)]

    @test mean([q' * ((1/(k-1))*C[k,:,:]*q) ≈ c[k] for k ∈ 2:maximum(keys(H.E))]) == 1
end

@testset "BP linearization matrix" begin

    BP_mat, ix = linearizedBPMatrix(H, z)

    # get the eigenvector corresponding to the overall linearization
    u = aggregateEigenvector(BP_mat, ix)
    
end


@testset "faster nonbacktracking matrix" begin

    # the actual fast version
    Bs, ix = nonBacktrackingMatrices(H);

    # test against deprecated slow version for correctness
    Bs2, ix2 = nonBacktrackingMatrices_(H; return_indices = true);
    
    B = sum(Bs)
    B2 = sum(Bs2)

    E = eigs(B; nev = 2, ritzvec = true)
    E2 = eigs(B2; nev = 2, ritzvec = true)

    @test E[1] ≈ E2[1]

    v = E[2][:,2]

    u = aggregateEigenvector(v, ix)

end


