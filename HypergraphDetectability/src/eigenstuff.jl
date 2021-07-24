function sumEigenvector(v::Vector{Complex{Float64}}, ix, mode = "out")
    nodes = unique(pe.point for pe ∈ values(ix))
    
    n = length(nodes) # number of nodes 
    M = length(ix)    # number of pointed edges
    
    # number of distinct edge sizes
    K = unique(length(pe.nodes) for pe ∈ values(ix))
    k̄ = length(K)
    Kix = Dict(K[i] => i for i ∈ 1:k̄)
    # number of communities
    ℓ = length(v) ÷ M

    # outline of aggregated vector
    u = zeros(ℓ * k̄ * n)

    # populate u by looping over v

    for l ∈ 1:ℓ
        v̄ = v[(1+(l-1)*M):(l*M)]
        for edge ∈ keys(ix)
                k = ix[edge].nodes |> length
            if mode == "out"
                i = ix[edge].point
                u[i + (Kix[k]-1)*n + (l-1)*k̄*n] += v̄[edge]            
            elseif mode == "in"
                for i ∈ ix[edge].nodes
                    if i != ix[edge].point
                        u[i + (Kix[k]-1)*n + (l-1)*k̄*n] += v̄[edge]
                    end
                end
            end
        end
    end

    return u
end

function aggregateEigenvector(v::Vector{Complex{Float64}}, ix)
    """
    given an eigenvector v of a nonbacktracking operator and edge list ix corresponding indices of v to pointed edges, compute an aggregate eigenvector of length n. 

    The entry of the new eigenvector corresponding to i is the sum of all entries in v corresponding to edges for which i is the point. 

    this is the eigenvector actually used for community detection applications. 

    corresponds to v_in (eq. 22)  in Angelini et all 

    https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7446987&casa_token=JyWZL63o7IkAAAAA:g5O4PX6MDgIQZp-m-kXcaL1fC663brk3DHbv6oaOe_WwjHNiRV1ZbnbBBdIEUSSVo_dKrl29eA&tag=1
    """

    nodes = unique(pe.point for pe ∈ values(ix))
    n = length(nodes)
    u = zeros(Complex{Float64}, n)

    for edge in keys(ix)
        u[ix[edge].point] += v[edge]
    end

    return u
end

function aggregateEigenvector(B, ix; mode = "single")
    E = Arpack.eigs(B; nev = 2, ritzvec = true)
    
    if !(imag(E[1][2]) ≈ 0)
        println("Warning: 2nd eigenvalue complex")
    end

    v = E[2][:,2]

    if mode == "BP"
        M = length(ix)
        distinct_k = size(B)[1] ÷ M
        u = sum([aggregateEigenvector(v[((i-1)*M + 1):i*M], ix) for i ∈ distinct_k])
        return u
    end

    return aggregateEigenvector(E[2][:,2], ix)
end


function binaryClusters(B, ix, ϵ = 0)
    u = aggregateEigenvector(B, ix)
    return (u .> ϵ) .+ 1
end

