"""
unit tests
to run all tests:

>>> cd HypergraphNB
>>> julia

julia> ]

pkg> activate .
pkg> test

"""

using Revise
using Test
using HypergraphNB

# using SparseArrays
# using LinearAlgebra
using Arpack
using Statistics
##


N = [50, 50]
C = [NaN, 5.0, 5.0]
P = [NaN, 0.1, 0.9]

H = plantedPartitionHypergraph(N, C, P; enforce_distinct = true)
z = vcat([repeat([z], N[z]) for z ∈ 1:length(N)]...)

n = length(H.N)

@testset "reduced nonbacktracking matrix" begin
    """
    test of the generalized Ihara-Bass formula for the hypergraph nonbacktracking matrix
    """
        # "ground truth": full nonbacktracking matrix
        B = nonBacktrackingMatrix(H);
        E = eigs(B; nev = 5);

        # reduced version from Ihara-Bass
        B_ = reducedNonBacktrackingMatrix(H)
        E_ = eigs(B_; nev = 5)

        # must agree by Ihara-Bass Theorem
        @test E_[1] ≈ E[1]
end

@testset "degree tensor" begin
    """
    test of a combinatorial identity in the estimated degree tensor
    """

    c, C = degreeTensor(H, z)

    q = 1/n * [sum(z .== i) for i in unique(z)]

    @test mean([q' * ((1/(k-1))*C[k-1,:,:]*q) ≈ c[k-1] for k ∈ 2:maximum(keys(H.E))]) == 1
end

@testset "jacobian matrix" begin
    """
    test of the generalized Ihara-Bass formula for the BP Jacobian
    """

    J = BPJacobian(H, z)
    E = eigs(J; nev = 5);

    J_ = reducedBPJacobian(H, z)
    E_ = eigs(J_; nev = 5);

    # must agree by generalized Ihara-Bass for the Jacobian
    @test E_[1] ≈ E[1]
end


