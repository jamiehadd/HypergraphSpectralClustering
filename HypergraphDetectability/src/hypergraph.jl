# largely copy-pasta'd from 
# https://github.com/PhilChodrow/HypergraphModularity/blob/master/src/HSBM.jl

mutable struct hypergraph
    """
    A very simple hypergraph composite type, designed to hold a node list N, an edge list E, a degree sequence D,
    TODO: remove H.D
    """

    N::Vector{Int64}
    E::Dict{Int64, Dict}
    mat::Dict{String, Dict{Int64, SparseArrays.SparseMatrixCSC}}
    function hypergraph(N, E)
        delete!(E, 1)             # ignore size 1 edges
        for k ∈ keys(E)
            if length(values(E[k])) == 0
                delete!(E,k)
            end
        end
        mat = cacheMatrices(N, E)
        new(N, E, mat)
    end
end

function cacheMatrices(N, E)
    n = length(N)
    K = sort(collect(keys(E)))

    D = Dict(k => degreeMatrix(N, E, k) for k ∈ K)
    A = Dict(k => adjacencyMatrix(N, E, k) for k ∈ K)

    return Dict("adj" => A, "deg" => D)
end

function degreeVectors(H)
    return Dict(k => collect(LinearAlgebra.diag(D)) for (k, D) ∈ H.mat["deg"])
end

Base.copy(H::hypergraph) = hypergraph(H.N, H.E, H.mat)

