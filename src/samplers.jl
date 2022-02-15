function detectabilityData(n, c₂, c₃, p₂, p₃)
    """
    construct a hypergraph of two clusters in which: 
    - n is the number of nodes (should be even)
    - c₂ is the mean number of 2-edges attached to a node
    - c₃ is the mean number of 3-edges attached to a node
    - p₂ is the proportion of 2-edges within clusters
    - p₃ is the proportion of 3-edges within clusters
    
    p₂ = 0 and p₃ = 1 for example, means that all 2-edges are between 
    cluster and all 3-edges are within cluster. 
    """
    n_group = n ÷ 2
    
    # numbers of edges of sizes 2 and 3
    m₂ = c₂*n/2
    m₃ = c₃*n/3
    
    # initialize edge dict
    E = Dict(k => Dict() for k ∈ 2:3)
    
    # assign 2-edges
    for i ∈ 1:m₂
        # within-cluster 2-edges
        if rand() < p₂
            s, t = 0, 0
            while s == t
                s, t = rand(1:n_group), rand(1:n_group)
            end
            if rand() < 0.5
                s += n_group
                t += n_group
            end
        # between-cluster 2-edges
        else
            s = rand(1:n_group)
            t = rand(1:n_group) + n_group
        end
        s, t = sort([s, t])
        E[2][[s,t]] = get(E[2], [s,t], 0) + 1
    end
    
    # assign 3-edges
    
    for i ∈ 1:m₃
        # within-cluster 3-edges
        if rand() < p₃
            s, t, v = 0, 0, 0
            while length(unique([s, t, v])) < 3
                s, t, v = rand(1:n_group), rand(1:n_group), rand(1:n_group)
            end
            if rand() < 0.5
                s += n_group
                t += n_group
                v += n_group
            end
        # between-cluster 2-edges
        else
            s, t, v = 0, 0, 0
            while length(unique([s, t, v])) < 3
                s = rand(1:n_group)
                t = rand(1:n_group)
                v = rand(1:n_group) + n_group
            end
            if rand() < 0.5
                s += n_group
                t += n_group
                v -= n_group
            end
            s, t, v = sort([s, t, v])
        end
        E[3][[s, t, v]] = get(E[3], [s, t, v], 0) + 1
    end

    H = hypergraph(1:n, E)
    return(H)
end


function plantedPartitionHypergraph(N::Vector,  C::Vector, P::Vector; enforce_distinct = true)

    E = Dict(k => Dict() for k ∈ 1:length(C))

    tot_N = sum(N)
    nodes = collect(1:tot_N)

    Q = StatsBase.Weights(N ./ tot_N)

    clusters = Dict()

    ℓ = 1
    for i ∈ 1:length(N)
        clusters[i] = ℓ:(ℓ - 1 + N[i])
        ℓ += N[i]
    end

    z = vcat([repeat([y], N[y]) for y ∈ 1:length(N)]...)
    # for each hyperedge size: 
    for k ∈ 1:length(C)

        if isnan(C[k])
            continue
        end

        # generate number of hyperedges of the appropriate size
        m = tot_N*C[k]/k

        # number of edges to generate
        M = Poisson(m) |> rand 

        for i ∈ 1:M
            # random case
            if rand() > P[k]
                if enforce_distinct
                    distinct = false
                    while !distinct
                        edge = sample(nodes, k; replace = false) |> sort
                        distinct = length(unique(z[edge])) > 1 
                    end
                else
                    edge = sample(nodes, k; replace = false) |> sort
                end
            # clustered case
            else
                y = sample(1:length(N), Q)
                choices = collect(clusters[y])
                edge = sample(choices, k; replace = false) |> sort
            end
            E[k][edge] = get(E[k], edge, 0) + 1
        end
    end

    H = hypergraph(1:tot_N, E)
    return(H)
end