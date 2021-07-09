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

    Hopefully a faster version, we'll see
    """
    
    # given a pointed edge, I would like to know what OTHER pointed edges
    # are adjacent through the point. So, I would like to lookup from an edge
    # to its point, and then from the point to the edges incident on that point. 

    edgeList = [[edge, k] for k ∈ keys(H.E) for edge ∈ keys(H.E[k])]
    pointedEdges = [pointedEdge(edge, v) for (edge, k) ∈ edgeList for v ∈ edge]

    edgeIDs = Dict(i => pe for (i, pe) ∈ enumerate(pointedEdges))

    incidence = Dict{Int64, Vector{Int64}}()
    for (id, pe) ∈ edgeIDs, i ∈ pe.nodes
        val = get(incidence, i, [])
        append!(val, id)
        incidence[i] = val
    end

    Ix = Vector{Int64}()
    Jx = Vector{Int64}()
    V  = Vector{Int64}()

    for (edgeID, pe) ∈ edgeIDs, edgeID2 ∈ incidence[pe.point]
        pe2 = edgeIDs[edgeID2]
        if (pe2.point != pe.point) && (pe.nodes != pe2.nodes)
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

function adjacencyMatrix(H, k = nothing)
    n = length(H.D)
    Ix, Jx, Vx = Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
    K = k === nothing ? keys(H.E) : [k]

    for k ∈ K, e ∈ keys(H.E[k]), (i, j) ∈ Combinatorics.combinations(e, 2)
        push!(Ix, i)
        push!(Jx, j)
        push!(Vx, 1)
        push!(Ix, j)
        push!(Jx, i)
        push!(Vx, 1)
    end

    A = SparseArrays.sparse(Ix, Jx, Vx, n, n)
    
    return A
end

function degreeDiagonalMatrix(H, k = nothing)

    n = length(H.D)

    Ix, Jx, Vx = Vector{Int64}(), Vector{Int64}(), Vector{Int64}()
    K = k === nothing ? keys(H.E) : [k]
    
    for k ∈ K, e ∈ keys(H.E[k]), i ∈ e
        push!(Ix, i)
        push!(Jx, i)
        push!(Vx, 1)
    end

    D = SparseArrays.sparse(Ix, Jx, Vx, n, n)
    return D
end

function reducedNonBacktrackingMatrices(H, K = sort(collect(keys(H.E))))
    Bs = []
    for k ∈ K
        if k == 1
            n = length(H.D)
            B̂ = SparseArrays.sparse([], [], [], n, n)
            push!(Bs, B̂)
        else
            A = adjacencyMatrix(H, k)
            D = degreeDiagonalMatrix(H, k)
            I = LinearAlgebra.UniformScaling(1)
            B̂ = hcat(vcat(zero(A), (1-k)*I), vcat(D - I, A - (k-2)*I))
            push!(Bs, B̂)        
        end
    end
    return Bs
end



function linearizedBPMatrix(H, ẑ; reduced = false)
    """
    ideally, should give indices as well as the relevant matrices
    """
    if reduced
        Bs, ix = reducedNonBacktrackingMatrices(H), nothing
    else
        Bs, ix = nonBacktrackingMatrices(H; return_indices = true)
    end

    c, C = degreeTensor(H, ẑ)

    n = length(ẑ)
    q = 1/n * StatsBase.counts(ẑ)

    T = zero(C)
    for k ∈ 1:size(C)[1]
        T[k,:,:] = (C[k,:,:] / ((k - 1) * c[k]) .- 1) .* q
    end

    not_nan = [k for k ∈ 1:size(T)[1] if !isnan(T[k,1,1])]

    BP_mat = sum(Kronecker.kronecker(T[k,:,:], Bs[k]) for k ∈ not_nan)

    return reduced ? BP_mat : BP_mat, ix
end