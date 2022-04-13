using HypergraphNB
using Arpack
using StatsBase
using Clustering
using TSne
using LinearAlgebra
using Revise
using Random 
using DelimitedFiles

include("scripts/utils.jl")

function experiment(H, n_groups, n_rounds = 1, nev=10, km_iter = 10, projected = false)
    if projected
        H = HypergraphNB.projectedGraph(H);
    end
    
    n = length(H.N)
    M = zeros(Int64, n, n)
    ẑ = rand(1:n_groups, n);
    
    for i ∈ 1:n_rounds
        println("Round $i")
        B = reducedBPJacobian(H, ẑ);
        E = Arpack.eigs(B; nev = nev);
        V = hcat([HypergraphNB.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:nev]...);
        V = real.(V)'
        V = 1.0*(V .> 0 )
            
        for j ∈ 1:km_iter
            clus = Clustering.kmeans(V, n_groups)
            if i >= 3
                ẑ = Clustering.assignments(clus);
                M += (ẑ .== ẑ')
            end    
        end
    end
    return M
end



Random.seed!(54321)

base_dir = "throughput/math-sx-cc"

if !isdir(base_dir)
    Base.Filesystem.mkpath(base_dir)
end


H, labels = HypergraphNB.read_unlabeled_data("tags-math-sx");

H, node_map = HypergraphNB.degreeFilter(H, 20);

n = length(H.N)

node_map_reverse = Dict(node_map[i] => i for i in keys(node_map))

label_vec = labels[[node_map_reverse[i] for i in H.N]]


n_groups = 4
nev = 20
n_epochs = 10
n_rounds = 10
km_iter = 100


for projected ∈ [true, false]

    structure = projected ? "graph" : "hypergraph"

    
    path = base_dir*"/"*structure

    clearDir!(path)

    M = sum(experiment(H, n_groups, n_rounds, nev, km_iter, projected) for i ∈ 1:n_epochs)

    writedlm("$base_dir/$structure/M.txt", M)
    writedlm("$base_dir/$structure/labels.txt", label_vec)
end

##############################################################################
# VANILLA MATRIX
##############################################################################

B = HypergraphNB.reducedNonBacktrackingMatrix(H)

for nev ∈ [5, 10]

    E = Arpack.eigs(B; nev = nev)
    V = real.(E[2])
    κ = length(keys(H.E))

    structure = "hypergraph-vanilla-$nev"
    path = base_dir*"/"*structure
    clearDir!(path)

    V̄ = sum([V[((i)*n + 1):(i+1)*n,:] for i ∈ 0:(κ - 1)])

    V̄ = 1.0*(V̄ .> 0)

    M = zeros(Int64, n, n)

    ẑ = zeros(n)

    for j ∈ 1:(km_iter*n_rounds*n_epochs)
        clus = Clustering.kmeans(V̄', n_groups)
        if j >= 3
            ẑ = Clustering.assignments(clus);
            M += (ẑ .== ẑ')
        end    
    end

    writedlm("$path/M.txt", M)
    writedlm("$path/labels.txt", label_vec)
end