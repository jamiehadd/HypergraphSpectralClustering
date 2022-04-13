using Revise
using HypergraphNB
using Arpack
using Clustering
using DataFrames
using CSV

include("utils.jl")



function experiment_inner(N, C, P₂, P₃)
    P2, P3,  ARI, EV = zeros(0), zeros(0), zeros(0), zeros(0)
    for p₂ ∈ P₂, p₃ ∈ P₃
        
        P = [NaN, p₂, p₃]

        H = plantedPartitionHypergraph(N, C, P; enforce_distinct = true)
        n = length(H.N)

        B = reducedNonBacktrackingMatrix(H)

        z = vcat([repeat([z], N[z]) for z ∈ 1:length(N)]...)

        try
            E = Arpack.eigs(B; nev = 2, ritzvec = true)
            
            u = E[2][:,2][1:n] + E[2][:,2][(n+1):2n]
            if !(mean(abs.(imag.(u))) ≈ 0)
                ari = 0.0
                append!(P2, p₂)
                append!(P3, p₃)
                append!(ARI, ari)
                println("whoops")
            else
                clusters = 1 .+ (real.(u) .> 0)
                ari = randindex(clusters, z)[1]
                append!(P2, p₂)
                append!(P3, p₃)
                append!(ARI, ari)
                println("yay")
            end
        catch e
            nothing
        end        
    end
    DF = DataFrame(
        P_2 = P2, 
        P_3 = P3, 
        ARI = ARI
    )
    return DF
end

function experiment(N, C, P₂, P₃, n_reps, path)

    clearFile!(path)
    
    for rep ∈ 1:n_reps
        println("round ", rep)
        df = experiment_inner(N, C, P₂, P₃)
        CSV.write(path, df, append = isfile(path))
    end
end


path = "throughput/bulk-throughput/exp-vanilla.csv"

# spacing between values of p₂, p₃, and p₄
δ = 0.01
# number of repetitions, used for averaging
n_reps = 20

N = [100, 100]
C = [NaN, 5, 5, 0]
P₂ = 0.0:δ:1.0
P₃ = 0.0:δ:1.0


experiment(N, C, P₂, P₃, n_reps, path) 
