##
using Revise
using HypergraphDetectability
using SparseArrays
using LinearAlgebra
using Arpack
using Plots
##

# get some fake data to play with
# these settings are enough to see 
# the second-largest eigenvalue separated
# from the bulk

## 
n  = 100
c₂ = 5
c₃ = 5
p₂ = 0.9
p₃ = 0.7

@time H = detectabilityData(n, c₂, c₃, p₂, p₃);

##
@time B = nonBacktrackingMatrix(H);
@time E = eigs(B; nev = 10);
## 

@time Bs = nonBacktrackingMatrices(H);
B_ = sum(Bs)
@time E_ = eigs(B_; nev = 10)

## 
ENV["GKSwstype"] = "100"
scatter(E[1], label = "", alpha = 0.2, color = "black")
scatter!([E[1][2]], label = "")

e = real.(E[1])

# xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
# ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (400, 250))
## 

#########################
# Reduced nonbacktracking matrix
#########################

## 
@time B_ = HypergraphDetectability.reducedNonBacktrackingMatrix(H);

@time E_ = eigs(B_; nev = 3);
##

## 
ENV["GKSwstype"] = "100"
scatter(E_[1], label = "")
scatter!([E_[1][2]], label = "")

e = real.(E_[1])

xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (500, 500))
scatter!([10^(0.5)], [0])
##


##
# Comparison
scatter(E[1], label = "B", ms = 4, alpha = 1.0, color = "white")
scatter!(E_[1], label = "B'", ms = 2, alpha = 1.0,
color = "firebrick")
xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
# ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (400, 250))
# png("fig/spectra-comparison.png")

## 