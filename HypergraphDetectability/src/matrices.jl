function nonBacktrackingMatrix(H)
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
            push!(V, 1)
        end
    end

    B = SparseArrays.sparse(Ix, Jx, V)

    return B
end