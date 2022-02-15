using Revise
using HypergraphDetectability
using Arpack
using Clustering
using StatsBase
using LinearAlgebra
using DataFrames
using CSV

##########################################
# FIRST EXPERIMENT: CONTACT HIGH SCHOOL
##########################################

function experiment(H, z; N_ev, N_groups, km_iter, n_rounds, data_name, path = "throughput/data-throughput/"*data_name*".csv")

    # random group vector
    # this is used for estimating the parameters in the linearized BP matrix, and is progressively refined throughout the rounds. 
    n = length(H.N)
    
    for n_ev ∈ N_ev, n_groups ∈ N_groups
        NGROUPS = zeros(0)
        NEV = zeros(0)
        round = zeros(0)
        tot_SS = zeros(0)
        cost_SS = zeros(0)
        ARI = zeros(0)
        data = []

        ẑ = rand(1:n_groups, n);
        for i ∈ 1:n_rounds
            B = reducedBPJacobian(H, ẑ)
            E = Arpack.eigs(B; nev = n_ev, ritzvec = true)            
            V = hcat([HypergraphDetectability.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:n_ev]...)
            V = real.(V)'
            V = V .> 0 # experiment
            
            tot_ss = (V .- mean(V, dims = 1)).^2 |> sum
            
            for j ∈ 1:km_iter        
                clus = Clustering.kmeans(V, n_groups)
                
                ẑ = Clustering.assignments(clus)
                
                append!(NGROUPS, n_groups)
                append!(NEV, n_ev)
                append!(tot_SS, tot_ss)
                append!(cost_SS, clus.totalcost)
                append!(ARI, randindex(ẑ, z)[1])
                append!(round, i)
                push!(data, data_name)
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

        CSV.write(path, DF, append = isfile(path))
    end
end

n_reps = 20

for rep ∈ 1:n_reps

    println("Rep ", rep)

    # primary school classes

    data_name = "contact-primary-school-classes"
    H, z = HypergraphDetectability.read_hypergraph_data(data_name);

    println("Beginning "*data_name)

    experiment(H, z; N_ev = [2, 5, 10, 20, 30, 40, 50], N_groups = [2, 3, 5, 7, 9, 11, 13, 15], km_iter = 20,  n_rounds = 10, data_name = data_name) 

    println("Beginning "*data_name*", projected graph")
    G = HypergraphDetectability.projectedGraph(H)

    experiment(G, z; N_ev = [2, 5, 10, 20, 30, 40, 50], N_groups = [2, 3, 5, 7, 9, 11, 13, 15], km_iter = 20,  n_rounds = 10, data_name = data_name*"-projected") 

    # high school classes

    data_name = "contact-high-school-classes"
    H, z = HypergraphDetectability.read_hypergraph_data(data_name);

    println("Beginning "*data_name)
    experiment(H, z; N_ev = [2, 5, 10, 20, 30, 40, 50], N_groups = [2, 3, 5, 7, 9, 11, 13, 15], km_iter = 20,  n_rounds = 10, data_name = data_name) 

    println("Beginning "*data_name*", projected graph")
    G = HypergraphDetectability.projectedGraph(H)

    experiment(G, z; N_ev = [2, 5, 10, 20, 30, 40, 50], N_groups = [2, 3, 5, 7, 9, 11, 13, 15], km_iter = 20,  n_rounds = 10, data_name = data_name*"-projected") 

    # senate bills

    data_name = "SN-congress-bills"
    H, z = HypergraphDetectability.read_hypergraph_data(data_name);

    println("Beginning "*data_name)
    experiment(H, z; N_ev = [2, 3, 4, 5, 10, 15, 20], N_groups = [2, 3, 4, 5], km_iter = 20,  n_rounds = 10, data_name = data_name) 

    println("Beginning "*data_name*", projected graph")
    G = HypergraphDetectability.projectedGraph(H)

    experiment(G, z; N_ev = [2, 3, 4, 5, 10, 15, 20], N_groups = [2, 3, 4, 5], km_iter = 20,  n_rounds = 10, data_name = data_name*"-projected") 
end

