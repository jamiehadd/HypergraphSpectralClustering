using Revise
using HypergraphDetectability
using Arpack

E = Dict(
    2 => Dict([1, 3] => 1, [2, 3] => 1, [2, 1] => 1)
)

H = hypergraph([1, 2, 3], E, [])

HypergraphDetectability.computeDegrees!(H)


Bs, ix = nonBacktrackingMatrices(H)
Bs[1]

eigvals = Arpack.eigs(Bs[1]; nev = 4)[1]



B = reducedNonBacktrackingMatrices(H)[1]
eigvals = Arpack.eigs(B; nev = 10)[1]








