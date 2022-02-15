using Revise
using HypergraphDetectability
using Arpack
using Clustering
using StatsBase
using LinearAlgebra
using DataFrames
using CSV

function experiment(H, z; N_ev, N_groups, km_iter, n_rounds, data_name, n_reps = 1)

    # random group vector
    # this is used for estimating the parameters in the linearized BP matrix, and is progressively refined throughout the rounds. 
    n = length(H.N)
    
    NGROUPS = zeros(0)
    NEV = zeros(0)
    round = zeros(0)
    tot_SS = zeros(0)
    cost_SS = zeros(0)
    ARI = zeros(0)
    data = []

    for n_reps ∈ 1:n_reps, n_ev ∈ N_ev, n_groups ∈ N_groups
        ẑ = rand(1:n_groups, n);
        for i ∈ 1:n_rounds
            println("round ", i)

            B = reducedBPJacobian(H, ẑ)
            E = Arpack.eigs(B; nev = n_ev, ritzvec = true)            

            V = hcat([HypergraphDetectability.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:n_ev]...)
            V = real.(V)'
            V = V .> 0 
            
            tot_ss = (V .- mean(V, dims = 1)).^2 |> sum
            
            best_kmeans_val = 1.0
            for j ∈ 1:km_iter        
                clus = Clustering.kmeans(V, n_groups)
                
                ẑ_ = Clustering.assignments(clus)
                
                obj = clus.totalcost/tot_ss

                if obj < best_kmeans_val
                    best_kmeans_val = obj
                    ẑ = ẑ_
                end

                append!(NGROUPS, n_groups)
                append!(NEV, n_ev)
                append!(tot_SS, tot_ss)
                append!(cost_SS, clus.totalcost)
                append!(ARI, randindex(ẑ_, z)[1])
                append!(round, i)
                push!(data, data_name)
            end
        end
    end
    DF = DataFrame(
        ngroups = NGROUPS, 
        nev = NEV, 
        tot_SS = tot_SS,
        cost_SS = cost_SS, 
        ari = ARI, 
        round = round,
        data = data
    )
end

data_name = "SN-congress-bills"
H, z = HypergraphDetectability.read_hypergraph_data(data_name);

data_name = "SN-congress-bills"
n_ev = 5
n_groups = 2
kmax = 10

H = hypergraph(H.N, Dict(k => v for (k, v) ∈ H.E if k <= kmax))

experiment(H, z; N_ev = [5], N_groups = [2], km_iter = 100, n_rounds = 5, data_name = data_name, n_reps = 1) |> 
    CSV.write("throughput/data-throughput/"*data_name*".csv")



G = HypergraphDetectability.projectedGraph(H)

experiment(G, z; N_ev = [5], N_groups = [2], km_iter = 100, n_rounds = 5, data_name = data_name, n_reps = 1) |> 
    CSV.write("throughput/data-throughput/"*data_name*"-projected.csv")

