function partitionize(a::Vector{<:Integer})
    """
    For a given integer vector a, return the partition corresponding to that
    vector. Useful for both counting corrections when sampling and computing
    likelihoods, and when computing partition-based values of Ω.
    This is the fastest version I could come up with.
    """
    a = sort(a)
    k = length(a)
    v = zero(a)

    v[1] = 1
    current = 1

    for i = 2:k
        if a[i] == a[i-1]
            v[current] += 1
        else
            current += 1
            v[current] = 1
        end
    end
    # v = v[v.>0]
    # return sort(v, rev = true)
    sort!(v, rev = true)
    return sortedremovezeros(v)
end

function sortedremovezeros(p::Vector{<:Integer})
    for i = 2:length(p)
        if p[i] == 0
            return p[1:i-1]
        end
    end
    return p
end

function counting_coefficient(z::Array{T, 1}) where {T<:Integer}
    p = partitionize(z)
    return Combinatorics.multinomial(p...)
end

function poisson_pdf(x::Integer, λ::Float64)
    exp(-λ)*λ^x/factorial(big(x))
end

function projectedGraph(H)
    E_ = Dict(2 => Dict())
    for k ∈ keys(H.E)
        for (e, m) ∈ H.E[k], (i, j) ∈ Combinatorics.combinations(e, 2)
            E_[2][e] = get(E_[2], e, 0) + m
        end
    end
    return hypergraph(H.N, E_)
end

function read_unlabeled_data(dataname)
    ix = Int64[]
    open("data/$dataname/hyperedges-$dataname.txt") do f 
        for line in eachline(f)
            push!(ix, parse(Int64, line))
        end
    end
    
    E = Dict()
    i = 1
    open("data/$dataname/nverts-$dataname.txt") do f
        for line in eachline(f)
            m = parse(Int64, line)
                
            try
                e = ix[i:(m+i-1)]
                i += m

                e = sort(e)
                k = length(e)
                if k ∈ keys(E)
                    E[k][e] = get(E[k], e, 0) + 1
                else
                    E[k] = Dict()
                end
            catch er
            end
        end
    end
    
    labels = String[]
    open("data/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(labels, line)
        end
    end
    
    return hypergraph(unique(ix), E), labels
end