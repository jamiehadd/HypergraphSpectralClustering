function nonBacktrackingMatrix(H)
    return nonBacktrackingMatrices(H)[1] |> sum
end

function nonBacktrackingMatrices(H; K = sort(collect(keys(H.E))), return_ix = false)
    """
    compute the size-specific nonbacktracking matrices for 
    a nonuniform hypergraph. 

    Return: Bs, an Array, whose kth entry gives B_k, the Kth nonbacktracking operator as described in the notes. The argument K can be used to restrict the values of k for which B_k is computed, although ATM Phil can't think of any reason why one would want to do this.  

    Aggregated nonbacktracking matrix can be obtained by computing 
    sum(Bs), where Bs is the return value of this function (notation needs improvement)
    """
    
    edgeList = [[edge, k] for k ∈ keys(H.E) for (edge, m) ∈ H.E[k] for j ∈ 1:m]
    push!.(edgeList, 1:length(edgeList)) 
    pointedEdges = [pointedEdge(edge, v, eid) for (edge, k, eid) ∈ edgeList for v ∈ edge]

    edgeIDs = Dict(i => pe for (i, pe) ∈ enumerate(pointedEdges))

    # i => ids of edges incident to i
    incidence = Dict{Int64, Vector{Int64}}()
    for (id, pe) ∈ edgeIDs, i ∈ pe.nodes
        val = get(incidence, i, [])
        append!(val, id)
        incidence[i] = val
    end

    Ix = Vector{Int64}()
    Jx = Vector{Int64}()
    V  = Vector{Int64}()

    # nonbacktracking logic -- check this?
    for (edgeID, pe) ∈ edgeIDs, edgeID2 ∈ incidence[pe.point]
        pe2 = edgeIDs[edgeID2]
        
        if (pe2.point != pe.point) && (pe.eid != pe2.eid)
            # print(edgeID)
            # print(" ")
            # println(edgeID2)
            push!(Ix, edgeID)
            push!(Jx, edgeID2)
            push!(V, 1)
        end
    end

    B = []

    M = length(edgeIDs)
    for k ∈ K
        size_matches = [length(edgeIDs[j].nodes) == k for j ∈ Jx]
        b = SparseArrays.sparse(Ix[size_matches], 
                                Jx[size_matches], 
                                V[size_matches],
                                M, M)
        push!(B, b)
    end

    return (B, edgeIDs)
end

#########################
# Reduced nonbacktracking matrix
#########################


function degreeMatrix(N, E, k)
    """
    kth degree matrix
    """

    n = length(N)
    Ix = 1:n
    Jx = 1:n
    Vx = zeros(Int64, n)
    
    for (e, m) ∈ E[k], i ∈ e
        Vx[i] += m         
    end

    D = SparseArrays.sparse(Ix, Jx, Vx, n, n)
    return D
end

degreeMatrix(H::hypergraph, k) = H.mat["deg"][k]

function adjacencyMatrix(N, E, k)
    """
    kth adjacency matrix
    """
    n = length(N)
    
    Ix, Jx, Vx = Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
    
    for (e, m) ∈ E[k], (i, j) ∈ Combinatorics.combinations(e, 2)
        push!(Ix, i)
        push!(Jx, j)
        push!(Vx, m)
        push!(Ix, j)
        push!(Jx, i)
        push!(Vx, m)
    end

    A = SparseArrays.sparse(Ix, Jx, Vx, n, n)
    
    return A
end

adjacencyMatrix(H, k) = H.mat["adj"][k]

function adjacencyBlockMatrix(H)
    K = sort(collect(keys(H.E)))
    a = vcat([adjacencyMatrix(H, k) for k ∈ K]...)
    A = hcat([a for k ∈ K]...)
    return A
end

function degreeBlockMatrix(H)
    K = sort(collect(keys(H.E)))
    d = vcat([degreeMatrix(H, k) for k ∈ K]...)
    D = hcat([d for k ∈ K]...)
    return D
end

# this seems to work, so now we know what we should prove. How to get there...
function reducedNonBacktrackingMatrix(H)
    K = sort(collect(keys(H.E)))
    K = diagm(K)

    n = length(H.N)

    D = degreeBlockMatrix(H)

    A = adjacencyBlockMatrix(H)

    upperLeft = zero(D)
    upperRight = D-I

    lowerLeft = (I-K)⊗I(n)
    lowerRight = A+(2I-K)⊗I(n)

    B_ = hcat(upperLeft, upperRight)
    B_ = vcat(B_, hcat(lowerLeft, lowerRight))
    return B_
end

function BPJacobian(H, ẑ; return_ix = false)
    """
    ideally, should give indices as well as the relevant matrices
    """    
    Bs, ix = nonBacktrackingMatrices(H; return_ix = true)

    K = sort(collect(keys(H.E)))

    c, C = degreeTensor(H, ẑ)

    n = length(ẑ)
    q = 1/n * StatsBase.counts(ẑ)

    T = zero(C)

    for i ∈ 1:length(K)
        k = K[i]
        T[i,:,:] = (C[i,:,:] / ((k - 1) * c[i]) .- 1) .* q
    end

    BP_mat = sum(Kronecker.kronecker(T[i,:,:], Bs[i]) for i ∈ 1:length(K))

    if return_ix
        return BP_mat, ix
    end

    return BP_mat
end

function reducedBPJacobian(H, ẑ)
    """
    Needs writing eventually, and a proof about why the matrix looks *exactly like this* would certainly be helpful. 
    Also, figuring out a bit about the kernel of the transformation would be good. 

    Finally, maybe we can make this faster? Pretty slow ATM e.g. on 1K nodes just to construct the matrix. 
    """

    # println("Warning: although this function is written in code that appears generalizable, it has only been used and tested for hypergraphs with edge sizes 2 and 3, and with only two clusters.")
    # edge sizes
    K_ = sort(collect(keys(H.E))) # list of sizes
    k̄ = length(K_)                # number of distinct sizes
    K = diagm(K_)                 # diagonal matrix of sizes

    ℓ = length(unique(ẑ))         # number of clusters
    n = length(H.N)               # number of nodes

    # graph structure matrices 
    # degree diagonal matrix
    d = [HypergraphNB.degreeMatrix(H, k) for k ∈ K_]
    D = cat(d..., dims = (1, 2))

    # adjacency block diagonal matrix
    a = [HypergraphNB.adjacencyMatrix(H, k) for k ∈ K_]
    A = cat(a..., dims = (1, 2));

    # parameter array: basic version
    c, G = degreeTensor(H, ẑ);
    q = 1/n * StatsBase.counts(ẑ)
    G_ = zero(G)
    for i ∈ 1:length(K_)
        G_[i,:,:] = (G[i,:,:] / ((K_[i] - 1) * c[i]) .- 1) .* q
    end
    

    # expanded parameter array #1
    dC = zeros(k̄*ℓ, k̄*ℓ)
    for k ∈ 1:k̄, s ∈ 1:ℓ, t ∈ 1:ℓ
        dC[k + (s-1)*k̄, k + (t-1)*k̄] = G_[k, s, t]
    end
    dC = sparse(dC)

    # expanded parameter array #2
    C = zeros(k̄*ℓ, k̄*ℓ)
    for k ∈ 1:k̄, s ∈ 1:ℓ, t ∈ 1:ℓ, k_ ∈ 1:k̄
        C[k + (s-1)*k̄, k_ + (t-1)*k̄] = G_[k, s, t]
    end
    C = sparse(C)

    # construct main blocks 

    upperRight = (sparse(C ⊗ I(n)) * sparse(I(ℓ) ⊗ D) - sparse(dC ⊗ I(n)))'
    lowerLeft  = sparse((dC * (I(ℓ)⊗(I - K))) ⊗ I(n))'
    lowerRight = (sparse(C ⊗ I(n)) * sparse((I(ℓ) ⊗ A)) - sparse((dC*(I(ℓ)⊗(K - 2I)))⊗I(n)))'

    B_ = hcat(zero(upperRight), upperRight);
    B_ = vcat(B_, hcat(lowerLeft, lowerRight));

    # M = n*ℓ*k̄
    # B_ = spzeros(2M, 2M)

    # B_[1:M, (M+1):end] = (sparse(C ⊗ I(n)) * sparse(I(ℓ) ⊗ D) - sparse(dC ⊗ I(n)))'
    # B_[(M+1):end, 1:M] = sparse((dC * (I(ℓ)⊗(I - K))) ⊗ I(n))'
    # B_[(M+1):end, (M+1):end] = (sparse(C ⊗ I(n)) * sparse((I(ℓ) ⊗ A)) - sparse((dC*(I(ℓ)⊗(K - 2I)))⊗I(n)))'
    return B_
end


function reducedBPJacobian_(H, ẑ)
    """
    Needs writing eventually, and a proof about why the matrix looks *exactly like this* would certainly be helpful. 
    Also, figuring out a bit about the kernel of the transformation would be good. 

    Finally, maybe we can make this faster? Pretty slow ATM e.g. on 1K nodes just to construct the matrix. 
    """

    # println("Warning: although this function is written in code that appears generalizable, it has only been used and tested for hypergraphs with edge sizes 2 and 3, and with only two clusters.")
    # edge sizes
    K_ = sort(collect(keys(H.E))) # list of sizes
    k̄ = length(K_)                # number of distinct sizes
    K = diagm(K_)                 # diagonal matrix of sizes

    ℓ = length(unique(ẑ))         # number of clusters
    n = length(H.N)               # number of nodes

    # graph structure matrices 
    # degree diagonal matrix
    d = [HypergraphNB.degreeMatrix(H, k) for k ∈ K_]
    D = cat(d..., dims = (1, 2))

    # adjacency block diagonal matrix
    a = [HypergraphNB.adjacencyMatrix(H, k) for k ∈ K_]
    A = cat(a..., dims = (1, 2));

    # parameter array: basic version
    c, G = degreeTensor(H, ẑ);
    q = 1/n * StatsBase.counts(ẑ)
    G_ = zero(G)
    for i ∈ 1:length(K_)
        G_[i,:,:] = (G[i,:,:] / ((K_[i] - 1) * c[i]) .- 1) .* q
    end
    

    # expanded parameter array #1
    dC = zeros(k̄*ℓ, k̄*ℓ)
    for k ∈ 1:k̄, s ∈ 1:ℓ, t ∈ 1:ℓ
        dC[k + (s-1)*k̄, k + (t-1)*k̄] = G_[k, s, t]
    end
    dC = sparse(dC)

    # expanded parameter array #2
    C = zeros(k̄*ℓ, k̄*ℓ)
    for k ∈ 1:k̄, s ∈ 1:ℓ, t ∈ 1:ℓ, k_ ∈ 1:k̄
        C[k + (s-1)*k̄, k_ + (t-1)*k̄] = G_[k, s, t]
    end
    C = sparse(C)

    # construct main blocks 


    M = n*ℓ*k̄
    B_ = spzeros(2M, 2M)

    B_[1:M, (M+1):end] = (sparse((C ⊗ I(n)) * (I(ℓ) ⊗ D)) - sparse(dC ⊗ I(n)))'
    B_[(M+1):end, 1:M] = ((dC * (I(ℓ)⊗(I - K))) ⊗ I(n))'
    B_[(M+1):end, (M+1):end] = (sparse(C ⊗ I(n)) * sparse((I(ℓ) ⊗ A)) - sparse((dC*(I(ℓ)⊗(K - 2I)))⊗I(n)))'
    return B_
end