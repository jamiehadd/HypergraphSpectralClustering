function read_hypergraph_label_names(dataname::String)
    names = String[]
    open("data/$dataname/label-names-$dataname.txt") do f
        for line in eachline(f)
            push!(names, line)
        end
    end
    return names
end

function read_hypergraph_labels(dataname::String)
    labels = Int64[]
    open("data/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    return labels
end

function read_hypergraph_edges(dataname::String, maxsize::Int64=25, minsize::Int64=2)
    E = Dict{Integer, Dict}()
    open("data/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if minsize <= length(edge) <= maxsize
                sz = length(edge)
                if !haskey(E, sz)
                    E[sz] = Dict{}()
                end
                E[sz][edge] = 1
            end
        end
    end
    return E
end


function read_hypergraph_data(dataname::String, maxsize::Int64=25, minsize::Int64=2, return_labels=true)
    E = read_hypergraph_edges(dataname, maxsize, minsize)
    
    n = maximum([maximum(e) for k in keys(E) for e in keys(E[k])])
    
    maxedges = maximum(keys(E))
    for k in 1:maxedges
        if !haskey(E, k)
            E[k] = Dict{}()
        end
    end
    
    N = 1:n
    
    if return_labels
        labels = read_hypergraph_labels(dataname)
        return hypergraph(N, E), labels[N]
    end
    
    return hypergraph(N, E)
end
;