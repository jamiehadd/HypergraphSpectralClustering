mutable struct pointedEdge
    nodes::Vector{Int64}
    point::Int64
    eid::Int64 # id number of unpointed edge
    pointedEdge(nodes, point, eid) = (point ∈ nodes) ? new(sort(nodes), point, eid) : error("point not in edge")
end