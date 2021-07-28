
function NBSC(H; nev, n_clusters)
    # ẑ = 
    println("not implemented")
end

function NBSC_round(H, ẑ; nev, ngroups)
    J = reducedBPJacobian(H, ẑ)

    E = Arpack.eigs(J; nev = nev)

    # figure out how many eigenvalues/vectors are real
    # always ≤ nev
    nreal = sum(imag.(E[1]) .≈ 0 )
    evals = real.(E[1][1:nreal])
    evecs = real.(E[2][:,1:nreal])

    # form consensus matrix
    ẑ = consensusSpectralClustering(H, evecs, ngroups)
    return ẑ
end

function consensusMatrix(H, evecs)
    
    n = length(H.N)
    C = zeros(n, n)
    k̄ = length(keys(H.E))
    ℓ = size(evecs)[1] ÷ (2*n*k̄)

    for h ∈ 1:size(evecs)[2]
        v = evecs[:, h]
        U = transform_eigenvector(v, H);
        for s ∈ 1:ℓ, i ∈ 1:n, j ∈ 1:n
            C[i,j] += (sign(U[i, s]) == sign(U[j, s]))
        end
    end
    return C
end

function consensusSpectralClustering(H, evecs, ngroups)
    C = consensusMatrix(H, evecs)
    D = diagm(vcat(sum(C, dims = 1)...))
    D_inv = diagm(vcat(1 ./ sum(C, dims = 1)...))
    
    # could consider other laplacians
    L = D_inv*(D - C)

    E = eigen(L)
    V = real.(E.vectors[:,2:(ngroups+1)])

    # kmeans on this space (other algs might be better)
    clus = Clustering.kmeans(V', ngroups)
    ẑ = Clustering.assignments(clus)
    return ẑ
end


function transform_eigenvector(v, H)
    n = length(H.N)
    k̄ = length(keys(H.E))
    α = v[((length(v)÷2)+1):end]
    ℓ = length(α) ÷ (n*k̄)
    U = zeros(n, ℓ)
    for i ∈ 1:n, k ∈ 1:k̄, s ∈ 1:ℓ
        U[i, s] += α[i + (k-1)*n + (s-1)*(n*k̄)]
    end
    return U
end
