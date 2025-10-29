"""
An un-directed graph for a simple network specified by an adjacency list (no multi-links)
"""
mutable struct Graph
    adjacency_list::Vector{Vector{Int64}}
    degrees::Vector{Int64}
    N::Int64

    function Graph(adjacency_list::Vector{Vector{Int64}})
        N = length(adjacency_list)
        degrees = zeros(Int64, N)

        for (i, list) in enumerate(adjacency_list)
            degrees[i] = length(list)
        end

        new(adjacency_list, degrees, N)
    end
end

"""
Returns a vector of neighbors for `node` in network `g`
"""
function neighbors(g::Graph, node::Int64)
    return g.adjacency_list[node]
end

"""
Average degree of a network `g`
"""
function k_average(g::Graph)
	return sum(g.degrees) / g.N
end

"""
Compute the degree distribution of a network `g`.
Returns two objects - the first is the list of degrees
and the second is the distribution
"""
function degree_distribution(g::Graph)
    k_max = maximum(g.degrees)
    k_min = minimum(g.degrees)
    p = zeros(k_max)
    for k ∈ g.degrees
        if k ≠ 0
            p[k] += 1 / g.N
        end
    end
    
    k_min == 0 ? k_min = 1 : nothing
    
    return collect(k_min:k_max), p[k_min:k_max]
end
