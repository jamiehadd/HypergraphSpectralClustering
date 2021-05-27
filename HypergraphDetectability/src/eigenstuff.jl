function aggregateEigenvector(v::Vector{Complex{Float64}}, ix)
    """
    given an eigenvector v of a nonbacktracking operator and edge list ix corresponding indices of v to pointed edges, compute an aggregate eigenvector of length n. 

    The entry of the new eigenvector corresponding to i is the sum of all entries in v corresponding to edges for which i is the point. 

    this is the eigenvector actually used for community detection applications. 

    corresponds to v_in (eq. 22)  in Angelini et all 

    https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7446987&casa_token=JyWZL63o7IkAAAAA:g5O4PX6MDgIQZp-m-kXcaL1fC663brk3DHbv6oaOe_WwjHNiRV1ZbnbBBdIEUSSVo_dKrl29eA&tag=1
    """

    n = maximum(i for (i, e) in values(ix))

    u = zeros(n)

    for edge in keys(ix)
        u[ix[edge][1]] += v[edge]
    end

    return u
end

function aggregateEigenvector(B, ix)
    E = Arpack.eigs(B; nev = 2, ritzvec = true)

    if !(imag(E[1][2]) ≈ 0)
        print("Warning: 2nd eigenvalue complex")
    end

    return aggregateEigenvector(E[2][:,2], ix)
end

function binaryClusters(B, ix, ϵ = 0)
    u = aggregateEigenvector(B, ix)
    return (u .> ϵ) .+ 1
end

function degreeTensor(H, z)
   """
   Intended to give empirical estimates of the two-point correlations (eq. (9) in Angelini). The rs-th entry of C_k should be the average number of edges between a given node in cluster r and all nodes in cluster s. 

   Probably not correct yet! Runs and has some of the correct ideas. 
   """
    k̄ = maximum(keys(H.E))
    ℓ = length(unique(z))
    C = zeros(k̄, ℓ, ℓ)

    for k ∈ 1:k̄, e ∈ keys(H.E[k]), (i, j) ∈ Combinatorics.combinations(e, 2)
        C[k, z[i], z[j]] += 1
        C[k, z[j], z[i]] += 1
    end

    q = 1/length(z) * [sum(z .== 1), sum(z .== 2)]

    for k ∈ keys(H.E)
        C[k,:,:] = C[k,:,:] ./  (length(z) * factorial(maximum([k-2, 0])))
        C[k,:,:] = C[k,:,:] ./ (q * q')
        # C[k,:,:] = diagm(1.0 ./ [sum(z .== 1), sum(z .== 2)]) * C[k,:,:]
    end


    return C/2
end
