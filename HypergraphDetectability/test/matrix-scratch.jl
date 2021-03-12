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
n = 500
c₂ = 3
c₃ = 3
p₂ = 0.5
p₃ = 0.9

H = detectabilityData(n, c₂, c₃, p₂, p₃)
## 

##
B = nonBacktrackingMatrix(H);
@time E = eigs(B; nev = 500);
## 

##
ENV["GKSwstype"] = "100"
scatter(E[1], label = "")
scatter!([E[1][2]], label = "")

e = real.(E[1])

xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (500, 500))
## 


#########################
# Reduced nonbacktracking matrix
#########################

## 
B_ = HypergraphDetectability.reducedNonBacktrackingMatrix(H);

@time E_ = eigs(B_; nev = 500);
##

## 
ENV["GKSwstype"] = "100"
scatter(E_[1], label = "")
scatter!([E_[1][2]], label = "")

e = real.(E_[1])

xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (500, 500))
##


##
# Comparison
scatter(E[1], label = "B", ms = 4, alpha = 1.0, color = "white")
scatter!(E_[1], label = "B'", ms = 2, alpha = 1.0,
color = "firebrick")
xlims!((minimum(e) - 0.1, maximum(e) + 0.1))
# ylims!((minimum(e) - 0.1, maximum(e) + 0.1))
plot!(size = (400, 250), dpi = 300)
png("fig/spectra-comparison.png")
## 