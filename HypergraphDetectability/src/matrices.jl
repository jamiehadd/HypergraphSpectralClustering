function nonBacktrackingMatrix(H; k = "all")

    # these are likely memory-intensive and should be reimplemented
    # as generators
    edgeList = [e for k in keys(H.E) for e in keys(H.E[k])]
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

function nonBacktrackingMatrices(H; K = sort(collect(keys(H.E))))
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
        if (v₂ ∈ e₁) && (e₂ != e₁) && (v₂ != v₁) 
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

    return B
end






#########################
# Reduced nonbacktracking matrix
#########################

# Indices: run through all nodes, and then increment edge sizes
# can start with kmin and kmax for simplicity, although this should
# later be improved for computational savings


function diagonalMatrixSkeleton(H)
    """
    Handy outline for constructing these matrices
    """
    k₀, k₁ = maximum([minimum(keys(H.E)), 2]), maximum(keys(H.E))   

    n  = length(H.D)
    N  = n*(k₁ - k₀ + 1)

    Ix = zeros(Int64, N)
    Jx = zeros(Int64, N)
    Vx = zeros(Int64, N)

    return k₀, k₁, n, Ix, Jx, Vx
end

function degreeMatrix(H)

    k₀, k₁, n, Ix, Jx, Vx = diagonalMatrixSkeleton(H)

    for i ∈ 1:n, k ∈ k₀:k₁
        Ix[i + (k-k₀)*n] = i + (k-k₀)*n
        Jx[i + (k-k₀)*n] = i + (k-k₀)*n

        edges = [val for (key, val) in H.E[k] if i ∈ key]

        Vx[i + (k-k₀)*n] = length(edges) == 0 ? 0 : sum(edges) 
    end

    N = n*(k₁ - k₀ + 1)
    D = SparseArrays.sparse(Ix, Jx, Vx, N, N)
    return D
end

function sizeScalingMatrix(H)
    k₀, k₁, n, Ix, Jx, Vx = diagonalMatrixSkeleton(H)

    for i ∈ 1:n, k ∈ k₀:k₁
        Ix[i + (k-k₀)*n] = i + (k-k₀)*n
        Jx[i + (k-k₀)*n] = i + (k-k₀)*n
        Vx[i + (k-k₀)*n] = k
    end

    N = n*(k₁ - k₀ + 1)
    J = SparseArrays.sparse(Ix, Jx, Vx, N, N)
    return J
end

function sizeAggregationMatrix(H)
    """
    Think this is also E⊗I in the case that we 
    fully enumerate all edge sizes. 
    """
    k₀, k₁ = maximum([minimum(keys(H.E)), 2]), maximum(keys(H.E))    
    n = length(H.D)

    N = n*(k₁ - k₀ + 1)^2

    Ix = zeros(Int64, N)
    Jx = zeros(Int64, N)
    Vx = zeros(Int64, N)
    
    ix = 1
    for i ∈ 1:n, k ∈ k₀:k₁, k_ ∈ k₀:k₁

        row = i + (k  - k₀)*n
        col = i + (k_ - k₀)*n

        Ix[ix] = row
        Jx[ix] = col
        Vx[ix] = 1

        ix += 1
    end

    N = n*(k₁ - k₀ + 1)
    S = SparseArrays.sparse(Ix, Jx, Vx, N, N)
    return S
end

function adjacencyAggregationMatrix(H)
    Ix = Vector{Int64}()
    Jx = Vector{Int64}()
    Vx = Vector{Int64}()

    k₀, k₁ = maximum([minimum(keys(H.E)), 2]), maximum(keys(H.E))    
    n = length(H.D)

    for k ∈ k₀:k₁, e in keys(H.E[k]), i ∈ e, j ∈ e
        if i != j
            push!(Ix, i + (k - k₀)*n)
            push!(Jx, j + (k - k₀)*n)
            push!(Vx, H.E[k][e])
        end
    end

    N = n*(k₁ - k₀ + 1)
    S = SparseArrays.sparse(Ix, Jx, Vx, N, N)
end


function reducedNonBacktrackingMatrix(H)
    D = HypergraphDetectability.degreeMatrix(H);
    J = HypergraphDetectability.sizeScalingMatrix(H);
    S = HypergraphDetectability.sizeAggregationMatrix(H);
    A = HypergraphDetectability.adjacencyAggregationMatrix(H);
    
    I = LinearAlgebra.UniformScaling(1)

    B_ = hcat(vcat(zero(A), I - J), vcat(D*S - I, A*S + 2I - J))

    return B_
end

