function degreeTensor(H, z)
    """
    Intended to give empirical estimates of the two-point correlations (eq. (9) in Angelini). The rs-th entry of C_k should be the average number of edges between a given node in cluster r and all nodes in cluster s. 
 
    Probably not correct yet! Runs and has some of the correct ideas. 
    """
    k̄ = maximum(keys(H.E))
    ℓ = length(unique(z))
    n = length(H.D)

    c = zeros(k̄)

    for k ∈ 1:k̄, e ∈ keys(H.E[k]), i ∈ e
        c[k] += 1
    end
    c = c ./ n

    C = zeros(k̄, ℓ, ℓ)

    for k ∈ 1:k̄, e ∈ keys(H.E[k]), (i, j) ∈ Combinatorics.combinations(e, 2)
        C[k, z[i], z[j]] += 1
        C[k, z[j], z[i]] += 1
    end

    counts = [sum(z .== i) for i in unique(z)]

    for k ∈ keys(H.E)
        C[k,:,:] = n*C[k,:,:] ./  (counts * counts')
    end

    return c, C
 end
 