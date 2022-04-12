using Revise
using HypergraphNB

using DataFrames
using Clustering
using Statistics
using RCall
using Arpack
using MultivariateStats
using TSne
using Random


# make some fake data to play with 

n_groups = 2
n_nodes = 100
N = repeat([n_nodes], n_groups);

# NaN refers to 1 edges, which don't matter
P = [NaN, 0.9, 0.9, .8];
C = [NaN, 5, 5, 5];
H = plantedPartitionHypergraph(N, C, P);
z = vcat([repeat([z], N[z]) for z ∈ 1:length(N)]...); # true labels
B = reducedNonBacktrackingMatrix(H)
E = Arpack.eigs(B; nev = 50);

E[1][1]

D = HypergraphNB.degreeVectors(H)

sum([mean(D[i])*(i-1) for i ∈ 2:length(C)]), E[1][1]

E[1][2]



c, T = degreeTensor(H, z)

0.5*sum(T[k-1,:,:][1] - T[k-1,:,:][2] for k ∈ 2:length(P)), E[1][2]



function estimate_ckin(C, P)
    return [2(k-1)*C[k]*(P[k] + (1-P[k])*(1 - 2.0^(2-k))/(2-2.0^(2-k))) for k ∈ 1:length(P)]
end

estimate_ckin(C, P)

[T[k,:,:][1] for k ∈ 1:(length(P) - 1)]


k = 3
1/(k-1)*(T[k-1,:,:][1] + T[k-1,:,:][2])/2, c[k-1]





