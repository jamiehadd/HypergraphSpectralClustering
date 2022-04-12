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

function subhypergraph(h::hypergraph, in_subhypergraph::Vector{Bool})
    # Get new set of edges
    new_edges = []
    for (sz, edges) in h.E
        for (edge, val) in edges
            new_edge = filter(v -> in_subhypergraph[v], edge)
            if length(new_edge) > 1
               push!(new_edges, new_edge)
            end
        end
    end
    
    # renumbering
    node_map = Dict{Int64,Int64}()
    for (i, val) in enumerate(in_subhypergraph)
        if val
            node_map[i] = length(node_map) + 1
        end
    end

    renumber_edge(e) = [node_map[v] for v in e]
    renumbered_new_edges = [renumber_edge(e) for e in new_edges]

    # New edges
    subE = Dict{Integer, Dict}()
    for edge in renumbered_new_edges
        sz = length(edge)
        if !haskey(subE, sz)
            subE[sz] = Dict{}()
        end
        subE[sz][edge] = 1
    end

    # New degrees
    n = length(node_map)
    subD = zeros(Int64, n)
    for (sz, edges) in subE
        for (e, _) in edges
            subD[e] .+= 1
        end
    end
    
    return hypergraph(1:n, subE), node_map
end

function degreeFilter(h, minDeg)
    DVecs = degreeVectors(h)
    d = [sum(D[i] for (k, D) in DVecs) for i in h.N]
    return subhypergraph(h, collect(d .>= minDeg))
end



Base.copy(H::hypergraph) = hypergraph(H.N, H.E, H.mat)

