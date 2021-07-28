# largely copy-pasta'd from 
# https://github.com/PhilChodrow/HypergraphModularity/blob/master/src/HSBM.jl

Parameters.@with_kw mutable struct hypergraph
    """
    A very simple hypergraph composite type, designed to hold a node list N, an edge list E, a degree sequence D,
    """

    N::Vector{Int64}
    E::Dict{Int64, Dict}
    D::Array{Int64, 1} = Array{Int64, 1}()
    mat::Dict{String, Array{Int64, 3}}
    function hypergraph(N, E, D)
        mat = Dict{String, Array{Int64, 3}}()
        new(N, E, D, mat)
    end
end


function computeDegrees(E::Dict{Int64, Dict}, N::Vector{Int64})
    """
    Compute the degree sequence of an edge list.
    """

    d = zeros(length(N))

    for k in keys(E)
        Ek = E[k]
        for e in keys(Ek)
            for i in e
                d[i] += 1
            end
        end
    end
    return(d)
end

function computeDegrees(H::hypergraph)
    return computeDegrees(H.E, H.N)
end

function computeDegrees!(H::hypergraph)
    """
    Compute the degree sequence of a hypergraph and store it as a field of the hypergraph.
    """
    H.D = computeDegrees(H)
end

function countEdges(H::hypergraph, pointed = false)
    """
    count the number of edges in H
    if pointed, a k-edge is counted k times
    """
    if pointed
        return sum([k*length(H.E[k]) for k in keys(H.E)])    
    else
        return sum([length(H.E[k]) for k in keys(H.E)])
    end
end

Base.copy(H::hypergraph) = hypergraph(H.N, H.E, H.D)

