function degreeTensor(H, z; normalized = false)
    """
    Intended to give empirical estimates of the two-point correlations (eq. (9) in Angelini). The rs-th entry of C_k should be the average number of edges between a given node in cluster r and all nodes in cluster s. 
    """
    K = sort(collect(keys(H.E)))

    ℓ = length(unique(z))
    n = length(H.N)

    c = zero(K)

    for ix ∈ 1:length(K), e ∈ keys(H.E[K[ix]]), i ∈ e
        c[ix] += 1
    end

    c = c ./ n

    C = zeros(length(K), ℓ, ℓ)

    for ix ∈ 1:length(K), e ∈ keys(H.E[K[ix]]), (i, j) ∈ Combinatorics.combinations(e, 2)
        C[ix, z[i], z[j]] += 1
        C[ix, z[j], z[i]] += 1
    end

    counts = [sum(z .== i) for i in unique(z)]

    for k ∈ 1:length(K)
        C[k,:,:] = n*C[k,:,:] ./  (counts * counts')
    end

    if normalized
        q = 1/n * StatsBase.counts(z)
        T = zero(C)

        for i ∈ 1:length(K)
            k = K[i]
            T[i,:,:] = (C[i,:,:] / ((k - 1) * c[i]) .- 1) .* q
        end
        return c, T
    end

    return c, C
 end
