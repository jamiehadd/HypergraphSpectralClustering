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
    B = nonBacktrackingMatrix(H);
    E = eigs(B; nev = 2);

    Bs = nonBacktrackingMatrices(H);
    B_ = sum(Bs)
    E_ = eigs(B_; nev = 2)

    @test E[1] ≈ E_[1]
end
