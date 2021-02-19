using Revise
using HypergraphDetectability

using SparseArrays
using LinearAlgebra
using Arpack

using Plots

# get some fake data to play with
# these settings are enough to see 
# the second-largest eigenvalue separated
# from the bulk

n = 500
c₂ = 3
c₃ = 3
p₂ = 5/6
p₃ = 0.5

H = detectabilityData(n, c₂, c₃, p₂, p₃)

# edge indices

B = nonBacktrackingMatrix(H)
@time E = eigs(B; nev = size(B)[1])

ENV["GKSwstype"] = "100"
scatter(E[1], label = "")
scatter!([E[1][2]], label = "")


e = real.(E[1])

xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (500, 500))

# not obvious that these are correlated with the 
# partition labels, but what one should do is 
# aggregate up to the node level and then check
# plot(real.(E[2][:,2]))