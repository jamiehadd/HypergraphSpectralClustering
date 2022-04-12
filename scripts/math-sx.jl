using HypergraphNB
using Arpack
using StatsBase
using Clustering
using TSne
using LinearAlgebra
using Revise
using Random 
using DelimitedFiles


function experiment(H, n_groups, n_rounds = 1, nev=10, km_iter = 10, projected = false)
    if projected
        H = HypergraphNB.projectedGraph(H);
    end
    
    ẑ = rand(1:n_groups, length(H.N));
    obj = 1
    V̂ = 0
    
    for i ∈ 1:n_rounds
        println("Round $i")
        B = reducedBPJacobian(H, ẑ);
        E = Arpack.eigs(B; nev = nev);
        V = hcat([HypergraphNB.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:nev]...);
        V = real.(V)'
        V = 1.0*(V .> 0 )
        tot_ss = (V .- mean(V, dims = 1)).^2 |> sum
            
        for j ∈ 1:km_iter
            clus = Clustering.kmeans(V, n_groups)
            cost_ss = clus.totalcost
            if cost_ss / tot_ss < obj
                ẑ = Clustering.assignments(clus);
                obj = cost_ss / tot_ss
                V̂ = V
            end    
        end
    end
    return ẑ, V̂, obj
end

Random.seed!(54321)

base_dir = "throughput/math-sx"

if !isdir(base_dir)
    Base.Filesystem.mkpath(base_dir)
end

for suffix ∈ ["graph", "hypergraph"]
    path = base_dir*"/"*suffix
    if isdir(path)
        Base.Filesystem.rm(path, recursive = true)
    end
    Base.Filesystem.mkpath(path)
end

H, labels = HypergraphNB.read_unlabeled_data("tags-math-sx");

H, node_map = HypergraphNB.degreeFilter(H, 20);

n = length(H.N)

node_map_reverse = Dict(node_map[i] => i for i in keys(node_map))

label_vec = labels[[node_map_reverse[i] for i in H.N]]


n_epochs = 50
nev = 15
n_rounds = 50
km_iter = 20
n_groups = 4

for projected ∈ [true, false]

    structure = projected ? "graph" : "hypergraph"

    best_ẑ_ = rand(1:n_groups, n);
    best_V̂_ = nothing
    best_obj_ = 1.0

    Ẑs   = zeros(n_epochs, n)
    Objs = ones(n_epochs)
    V̂s   = zeros(n_epochs, n_groups*nev, n)

    for i ∈ 1:n_epochs
        println("starting epoch $i")
        ẑ, V̂, obj = experiment(H, n_groups, n_rounds, nev, km_iter, projected)
        Ẑs[i,:] = ẑ
        Objs[i] = obj
        V̂s[i,:,:] = V̂
        println("Current best objective is $(minimum(Objs))")
    end

    i = argmin(Objs)

    best_obj_ = Objs[i]
    best_ẑ_ = Ẑs[i,:]
    best_V̂_ = V̂s[i,:,:]

    writedlm("$base_dir/$structure/obj.txt", best_obj_)
    writedlm("$base_dir/$structure/z.txt", best_ẑ_)
    writedlm("$base_dir/$structure/V.txt", best_V̂_)
    writedlm("$base_dir/$structure/labels.txt", label_vec)
end