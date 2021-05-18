using Revise
using HypergraphDetectability
using Statistics # for mean


# create some parameters
n = 20                                         # number of nodes
Z = rand(1:5, n)                               # clusters
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)  # degree parameters
kmax = 4                                       # size of largest hyperedge

# affinity function
# governs the likelihood of an edge between nodes
# based on their cluster labels
# can be any parameterized function of the partition
# of label vector, which counts numbers of nodes in each group 
# and ignores labels. For example, the partition vector
# of [2, 2, 1, 3, 3, 3, 3, 4] is [4, 2, 1, 1]. The second
# argument α is interpreted as a vector of parameters. 

# nothing special about this function
ω(p,α) = 100/factorial(sum(p)) * prod(p.^α)*float(n)^(-sum(p))

# nothing special about this parameter value
α0 = [1.0]

# create an AffinityFunction object using this affinity function
Ω = partitionAffinityFunction(ω,  kmax)

# finally, time to actually sample the hypergraph
H = sampleSBM(Z, ϑ, Ω;α=α0, kmax=kmax, kmin = 1)

# basic stats
for k ∈ 1:4
    println("There are $(length(H.E[k])) edges of size $k")
end
