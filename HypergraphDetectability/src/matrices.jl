function nonBacktrackingMatrix(H; k = "all")

    # these are likely memory-intensive and should be reimplemented
    # as generators
    # TODO: modify this so that we do multiple edges e if edge appears multiple times
    edgeList = [e for k in keys(H.E) for (e, m) in H.E[k] for j ∈ 1:m]
    pointedEdges = [[v, e] for e in edgeList for v in e]


    M = length(pointedEdges)
    push!.(pointedEdges, 1:M)

    Ix = Vector{Int64}()
    Jx = Vector{Int64}()
    V = Vector{Int64}()

    for (v₁, e₁, i) ∈ pointedEdges, (v₂, e₂, j) ∈ pointedEdges
        if (v₂ ∈ e₁) && (e₂ != e₁) && (v₂ != v₁) 
            push!(Ix, i)
            push!(Jx, j)
            push!( V, 1)
        end
    end

    B = SparseArrays.sparse(Ix, Jx, V, M, M)

    return B
end

function nonBacktrackingMatrices_(H; K = sort(collect(keys(H.E))), return_indices = false)
    """
    compute the size-specific nonbacktracking matrices for 
    a nonuniform hypergraph. 

    Return: Bs, an Array, whose kth entry gives B_k, the Kth nonbacktracking operator as described in the notes. The argument K can be used to restrict the values of k for which B_k is computed, although ATM Phil can't think of any reason why one would want to do this.  

    Aggregated nonbacktracking matrix can be obtained by computing 
    sum(Bs), where Bs is the return value of this function (notation needs improvement)

    VERY SLOW, definitely in need of performance improvements.
    Likely there are many unnecessary loops in here. 
    The main issue is the necessity of assigning to each pointed edge an index (to locate it in the nonbacktracking matrix)
    """
    # these are likely memory-intensive and should be reimplemented
    # as generators or something like that
    edgeList = [[e, k] for k in keys(H.E) for e in keys(H.E[k])]
    pointedEdges = [[v, e, k] for (e, k) in edgeList for v in e]


    M = length(pointedEdges)
    # index each edge by number
    push!.(pointedEdges, 1:M)

    Ix = Vector{Tuple{Int64, Int64}}()
    Jx = Vector{Tuple{Int64, Int64}}()
    V = Vector{Int64}()

    for (v₁, e₁, k₁, i) ∈ pointedEdges, (v₂, e₂, k₂, j) ∈ pointedEdges
        if (v₂ != v₁) && (v₂ ∈ e₁) && (e₂ != e₁) 
            push!(Ix, (i, k₁))
            push!(Jx, (j, k₂))
            push!( V, 1)
        end
    end

    B = []
    
    for k ∈ K
        size_matches = [k₂ == k for (j, k₂) ∈ Jx]
        ix = [i for (i, k₁) ∈ Ix]
        jx = [j for (j, k₂) ∈ Jx]
        b = SparseArrays.sparse(ix[size_matches],
                                jx[size_matches], 
                                V[size_matches], 
                                M, M)
        push!(B, b)
    end

    edgeIndices = Dict(ix => (v, e) for (v, e, k, ix) ∈ pointedEdges)

    return return_indices ? (B, edgeIndices) : B

end

function nonBacktrackingMatrices(H; K = sort(collect(keys(H.E))), return_indices = false)
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
    # println(edgeIDs)

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
            push!( V, 1)
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

# Indices: run through all nodes, and then increment edge sizes
# can start with kmin and kmax for simplicity, although this should
# later be improved for computational savings

function diagonalMatrixSkeleton(N, E)
    K = sort(collect(keys(E)))
    k̄ = length(K)
    n = length(N)
    M = n*k̄

    Ix = zeros(Int64, M)
    Jx = zeros(Int64, M)
    Vx = zeros(Int64, M)

    return K, n, Ix, Jx, Vx
end

diagonalMatrixSkeleton(H::hypergraph) = diagonalMatrixSkeleton(H.N, H.E)     


function degreeMatrix(H,k = nothing)
    """
    TODO: redo the handling of K so that it can handle missing edge sizes. 
    """

    K, n, Ix, Jx, Vx = diagonalMatrixSkeleton(H)

    if !isnothing(k)
        K = [k]
        Ix = zeros(Int64, n)
        Jx = zeros(Int64, n)
        Vx = zeros(Int64, n)
    end

    k₀ = K[1]
    for i ∈ 1:n, k ∈ K
        Ix[i + (k-k₀)*n] = i + (k-k₀)*n
        Jx[i + (k-k₀)*n] = i + (k-k₀)*n

        edges = [val for (key, val) in H.E[k] if i ∈ key]

        Vx[i + (k-k₀)*n] = length(edges) == 0 ? 0 : sum(edges) 
    end

    N = n*length(K)
    D = SparseArrays.sparse(Ix, Jx, Vx, N, N)
    return D
end

function adjacencyMatrix(H, k = nothing)
    n = length(H.D)
    Ix, Jx, Vx = Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
    K = k === nothing ? keys(H.E) : [k]

    for k ∈ K, (e, m) ∈ H.E[k], (i, j) ∈ Combinatorics.combinations(e, 2)
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

    n = length(H.D)

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

function linearizedBPMatrix(H, ẑ)
    """
    ideally, should give indices as well as the relevant matrices
    """    
    Bs, ix = nonBacktrackingMatrices(H; return_indices = true)

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

    return BP_mat, ix
end

function reducedBPJacobian(H, ẑ)
    """
    Needs writing eventually, and a proof about why the matrix looks *exactly like this* would certainly be helpful. 
    Also, figuring out a bit about the kernel of the transformation would be good. 

    Finally, maybe we can make this faster? Pretty slow ATM e.g. on 1K nodes just to construct the matrix. 
    """

    println("Warning: although this function is written in code that appears generalizable, it has only been used and tested for hypergraphs with edge sizes 2 and 3, and with only two clusters.")
    # edge sizes
    K_ = sort(collect(keys(H.E))) # list of sizes
    k̄ = length(K_)                # number of distinct sizes
    K = diagm(K_)                 # diagonal matrix of sizes

    ℓ = length(unique(ẑ))         # number of clusters
    n = length(H.D)               # number of nodes

    # graph structure matrices 
    # degree diagonal matrix
    d = [HypergraphDetectability.degreeMatrix(H, k) for k ∈ K_]
    D = cat(d..., dims = (1, 2))

    # adjacency block diagonal matrix
    a = [HypergraphDetectability.adjacencyMatrix(H, k) for k ∈ K_]
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

    println("forming main blocks")
    println("d")
    upperRight = (sparse(C ⊗ I(n)) * sparse(I(ℓ) ⊗ D) - sparse(dC ⊗ I(n)))'
    lowerLeft  = sparse((dC * (I(ℓ)⊗(I - K))) ⊗ I(n))'
    lowerRight = (sparse(C ⊗ I(n)) * sparse((I(ℓ) ⊗ A))   - sparse((dC*(I(ℓ)⊗(K - 2I)))⊗I(n)))'

    B_ = hcat(zero(upperRight), upperRight);
    B_ = vcat(B_, hcat(lowerLeft, lowerRight));

    return B_
end