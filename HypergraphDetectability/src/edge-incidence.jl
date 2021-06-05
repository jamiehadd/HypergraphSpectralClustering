mutable struct pointedEdge
    nodes::Vector{Int64}
    point::Int64
    pointedEdge(nodes, point) = (point ∈ nodes) ? new(nodes, point) : error("point not in edge")
end