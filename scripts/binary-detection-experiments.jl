using Revise
using HypergraphDetectability
using Arpack
using Clustering
using StatsBase
using LinearAlgebra
using DataFrames
using CSV
using MultivariateStats


function experiment_inner(N, C, P₂, P₃, P₄, projected = false, kmeans_reps = 100, nev = length(N))
    
    P2, P3, P4 = zeros(0), zeros(0), zeros(0)
    C2, C3, C4 = zeros(0), zeros(0), zeros(0)
    N1, N2 = zeros(0), zeros(0)
    ARI = zeros(0)

    for p₂ ∈ P₂, p₃ ∈ P₃, p₄ ∈ P₄
        P = [NaN, p₂, p₃, p₄]
        
        H = plantedPartitionHypergraph(N, C, P; enforce_distinct = true)
        if projected
            H = HypergraphDetectability.projectedGraph(H)
        end    
        # only used for estimating the degree tensor
        z = vcat([repeat([z], N[z]) for z ∈ 1:length(N)]...)
        B = reducedBPJacobian(H, z)
                
        try 
            E = Arpack.eigs(B; nev = 2, ritzvec = true)
            
            # λ₁ = maximum(real.(E[1]))
            # ix = (imag.(E[1]) .== 0 ) .& (abs.(real.(E[1])) .> sqrt(λ₁))
            # ix = findall(ix)

            V = hcat([HypergraphDetectability.transform_eigenvector(real.(E[2][:,i]), H) for i ∈ 1:2]...)
            V = real.(V)
            V = V .> 0 # experiment

            # normalize
            # V = V ./ sqrt.(sum(1.0*V.^2, dims = 1))
            
            # # try with PCA
            # m = fit(PCA, 1.0*V; maxoutdim=2)
            # V = m.proj

            # do kmeans a bunch of times
            best_cost = Inf
            best_clus = nothing
            for i ∈ 1:kmeans_reps
                clus = Clustering.kmeans(V', 2)
                if clus.totalcost < best_cost 
                    best_cost = clus.totalcost
                    best_clus = clus
                end
            end

            clusters = Clustering.assignments(best_clus)
            
            ari = randindex(clusters, z)[1]
            append!(P2, p₂)
            append!(P3, p₃)
            append!(P4, p₄)

            append!(C2, C[2])
            append!(C3, C[3])
            append!(C4, C[4])

            append!(N1, N[1])
            append!(N2, N[2])

            append!(ARI, ari)
        catch e
            nothing
        end
    end
    DF = DataFrame(
        P_2 = P2, 
        P_3 = P3,
        P_4 = P4,
        C_2 = C2, 
        C_3 = C3, 
        C_4 = C4,
        N_1 = N1, 
        N_2 = N2,
        ARI = ARI
    )

    return DF
end

function experiment(N, C, P₂, P₃, P₄, projected, n_reps, fname; kmeans_reps = 100, nev = length(N))
    fname = projected ? fname*"-projected" : fname
    path = "throughput/bulk-throughput/"*fname*".csv"

    if ! isfile(path)
        println("beginning "*fname)
        for rep ∈ 1:n_reps
            println("round ", rep)
            df = experiment_inner(N, C, P₂, P₃, P₄, projected, kmeans_reps, nev)
            CSV.write(path, df, append = isfile(path))
        end
    end
end

##########################################
# GLOBAL PARAMETERS
##########################################

# spacing between values of p₂, p₃, and p₄
δ = 0.01
# number of repetitions, used for averaging
n_reps = 20


##########################################
# FIRST EXPERIMENT: no 4-edges
##########################################

N = [100, 100]
C = [NaN, 5, 5, 0]
P₂ = 0.0:δ:1.0
P₃ = 0.0:δ:1.0
P₄ = [0]

experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-1"; kmeans_reps = 100, nev = 5)
experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-1"; kmeans_reps = 100, nev = 5) 

##########################################
# SECOND EXPERIMENT: no 3-edges
##########################################

N = [100, 100]
C = [NaN, 3, 0, 3]
P₂ = 0.0:δ:1.0
P₃ = [0]
P₄ = 0.0:δ:1.0

experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-2")
experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-2") 

##########################################
# THIRD EXPERIMENT: no 2-edges
##########################################

N = [100, 100]
C = [NaN, 0, 3, 3]
P₂ = [0]
P₃ = 0.0:δ:1.0
P₄ = 0.0:δ:1.0

experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-3") 
experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-3")

##########################################
# FOURTH EXPERIMENT: balanced clusters, no 4-edges, 
# very different degrees
##########################################

N = [100, 100]
C = [NaN, 5, 50, 0]
P₂ = 0.0:δ:1.0
P₃ = 0.0:δ:1.0
P₄ = [0]

experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-4")
experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-4") 

##########################################
# FIFTH EXPERIMENT: imbalanced clusters, no 3-edges, 
# very different degrees
##########################################

N = [100, 100]
C = [NaN, 50, 0, 5]
P₂ = 0.0:δ:1.0
P₃ = [0]
P₄ = 0.0:δ:1.0

experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-5")
experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-5") 


##########################################
# SIXTH EXPERIMENT: balanced clusters, some 2-edges, 2-edges not very informative
##########################################

N = [100, 100]
C = [NaN, 2, 3, 3]
P₂ = [.7]
P₃ = 0.0:δ:1.0
P₄ = 0.0:δ:1.0

experiment(N, C, P₂, P₃, P₄, true, n_reps, "exp-6")
experiment(N, C, P₂, P₃, P₄, false, n_reps, "exp-6")