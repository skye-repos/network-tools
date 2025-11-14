abstract type AbstractGraph end
"""
Un-weighted Graph object specified a dictionary where keys (node labels) have values that are vectors of node labels. Note that if the graph is directed, the key stores the out-links in the values
"""
mutable struct Graph{T} <: AbstractGraph where {T<:Union{String,Integer}}
    adjacency_list::Dict{T,Vector{T}}
end

"""
Default Graph constructor to give integer valued node labels
"""
function Graph()
    return Graph(Dict{Int,Vector{Int}}())
end

"""
Converts an vector of adjacency lists to a dictionary where the value attached to a key (a node label) are the (out)neighbors of the node.
"""
function adjList_to_dict(adjacency::Vector{Vector{T}}) where {T<:Integer}
    dict = Dict{T,Vector{T}}()

    for i ∈ eachindex(adjacency)
        if adjacency[i] ≠ []
            get!(dict, i, adjacency[i])
        end
    end

    return dict
end

"""
Create a copy of a mutable graph
"""
function copy(g::Graph)
    a = Base.deepcopy(g.adjacency_list)
    return Graph(a)
end

"""
Gets the number of nodes in the network `g`
"""
function size(g::Graph{T}) where {T<:Union{String,Integer}}
    return length(keys(g.adjacency_list))
end

"""
Returns a vector of neighbors for `node` in graph `g`. If the graph was directed, the vector is a list of out-neighbors (i.e nodes who have a link pointing to them from the node label)
"""
function neighbors(g::Graph{T}, node::T) where {T<:Union{String,Integer}}
    return get(g.adjacency_list, node, nothing)
end

"""
Returns a dictionary where the keys are node labels, and the values are the degrees. If the graph `g` is directed, then the values represent out-degrees.
"""
function degrees(g::Graph{T}) where {T<:Union{String,Integer}}
    d = Dict{T,Int}()
    for key ∈ keys(g.adjacency_list)
        get!(d, key, length(g.adjacency_list[key]))
    end

    return d
end

"""
Average degree of a network `g`
"""
function degree_average(g::Graph)
    return sum(values(degrees(g))) / size(g)
end

"""
Compute the degree distribution of a network `g`.
Returns two objects - the first is the list of degrees
and the second is the distribution
"""
function degree_distribution(g::Graph)
    k = values(degrees(g))
    N = size(g)
    k_max = maximum(k)
    k_min = minimum(k)
    p = zeros(k_max)
    for value ∈ k
        if k ≠ 0
            p[value] += 1 / N
        end
    end

    k_min == 0 ? k_min = 1 : nothing

    return collect(k_min:k_max), p[k_min:k_max]
end

"""
Add a link to graph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function add_link!(
    g::Graph{T},
    src::T,
    dst::T;
    directed=false) where {T<:Union{String,Integer}}

    if directed
        haskey(g.adjacency_list, src) ? push!(g.adjacency_list[src], dst) : g.adjacency_list[src] = [dst]
    else
        haskey(g.adjacency_list, src) ? push!(g.adjacency_list[src], dst) : g.adjacency_list[src] = [dst]
        haskey(g.adjacency_list, dst) ? push!(g.adjacency_list[dst], src) : g.adjacency_list[dst] = [src]
    end

    return nothing
end

"""
Add a link to graph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function add_link!(
    g::Graph{T},
    src::T,
    dsts::Vector{T};
    directed=false) where {T<:Union{String,Integer}}

    for dst ∈ dsts
        add_link!(g, src, dst, directed=directed)
    end

    return nothing
end

"""
Remove a link from graph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function delete_link!(
    g::Graph{T},
    src::T,
    dst::T;
    directed=false) where {T<:Union{String,Integer}}

    if !haskey(g.adjacency_list, src)
        error("Source $(src) node label does not exist as a key in the graph")
    end

    if !haskey(g.adjacency_list, dst)
        error("Destination $(dst) node label does not exist as a key in the graph")
    end

    if directed
        filter!(e -> e ≠ dst, g.adjacency_list[src])
    else
        filter!(e -> e ≠ dst, g.adjacency_list[src])
        filter!(e -> e ≠ src, g.adjacency_list[dst])
    end

    if g.adjacency_list[src] == []
        delete!(g.adjacency_list, src)
    end

    if g.adjacency_list[dst] == []
        delete!(g.adjacency_list, dst)
    end

    return nothing
end

"""
Remove a link from graph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function delete_link!(
    g::Graph{T},
    src::T,
    dsts::Vector{T};
    directed=false) where {T<:Union{String,Integer}}

    for dst ∈ dsts
        delete_link!(g, src, dst, directed=directed)
    end

    return nothing
end

"""
Return an edge-list for a graph `g`. If the network is directed the values are of the form (src, dst). If the network is un-directed, (src, dst) and (dst, src) are treated as the same entry and only one will be present in the edge list
"""
function edges(
    g::Graph{T};
    directed=false) where {T<:Union{String,Integer}}
    el::Vector{Tuple{T,T}} = []

    for key ∈ keys(g.adjacency_list), val ∈ g.adjacency_list[key]
        push!(el, (key, val))
    end

    if !directed
        for (src, dst) ∈ el
            if (dst, src) ∈ el
                filter!(e -> e ≠ (dst, src), el)
            end
        end
    end

    return el
end

"""
Count the number of edges in the graph `g`
"""
function edge_count(g::Graph; directed=false)
    return length(edges(g, directed=directed))
end

"""
Returns a list of nodes in the graph `g`
"""
function nodes(g::Graph)
    return keys(g.adjacency_list)
end

g = Graph()
add_link!(g, 3, 55)
add_link!(g, 3, [33, 25, 5, 1, 2])
add_link!(g, 2, [5, 33, 25])
add_link!(g, 25, [55, 33])
deg = degrees(g)
degree_average(g)
delete_link!(g, 3, 55)
delete_link!(g, 3, [33, 25, 5, 1, 2])
nodes(g)


