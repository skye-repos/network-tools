using Base: delete
using SparseArrays
include("graph.jl")

"""
A weighted graph stored as a dictionary of dictionaries where adjacency_list[i][j] represents the weight of the link from node i to node j
"""
mutable struct WeightedGraph{T,U} <: AbstractGraph where {T<:Union{String,Integer}} where {U<:Number}
    adjacency_list::Dict{T,Vector{Tuple{T,U}}}
end

"""
Default WeightedGraph constructor to initialzie with integer valued node labels and float valued weights
"""
function WeightedGraph()
    return WeightedGraph(Dict{Int,Vector{Tuple{Int,Float64}}}())
end

"""
Convert a matrix of `weights` to a WeightedGraph
"""
function weights_to_weightedGraph(weights::AbstractMatrix{U}) where {U<:Number}
    adj = Dict{Int,Vector{Tuple{Int,U}}}()

    for j ∈ axes(weights, 2), i ∈ axes(weights, 1)
        w = weights[i, j]
        if !iszero(w)
            inner = get!(Vector{Tuple{Int,U}}, adj, i)
            push!(inner, (j, w))
        end
    end

    return WeightedGraph(adj)
end

"""
Convert a WeightedGraph `g` into a (possible sparse) matrix of weights
"""
function weightedGraph_to_weights(g::WeightedGraph{Int,U}) where {U<:Number}
    N = maximum(keys(g.adjacency_list))
    w = zeros(U, N, N)

    for i ∈ nodes(g), j ∈ neighbors(g, i)
        idx = findfirst(e -> e == j, neighbors(g, i))
        w[i, j] = last(neighbors(g, i)[idx])
    end

    if issparse(w)
        return sparse(w)
    else
        return w
    end
end

"""
Create a copy of a WeightedGraph `g`
"""
function copy(g::WeightedGraph)
    a = Base.deepcopy(g.adjacency_list)
    return WeightedGraph(a)
end

"""
List the neighbors of `node` in WeightedGraph `g`
"""
function neighbors(g::WeightedGraph{T,U}, node::T) where {T<:Union{String,Integer}} where {U<:Number}
    if !isnothing(get(g.adjacency_list, node, nothing))
        return first.(g.adjacency_list[node])
    else
        return nothing
    end
end

"""
Gets the weight for the link from `src` to `dst` in WeightedGraph `g`
"""
function get_weight(g::WeightedGraph, src, dst)
    idx = findfirst(e -> e == dst, neighbors(g, src))
    return g.adjacency_list[src][idx][2]
end

"""
Gets the weight for the links from `src` to `dsts` in WeightedGraph `g`
"""
function get_weight(g::WeightedGraph, src, dsts::AbstractVector)
    wts = []

    for i ∈ eachindex(dsts)
        push!(wts, get_weight(g, src, dsts[i]))
    end

    return wts
end

"""
Add a link to WeightedGraph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function add_link!(
    g::WeightedGraph{T,U},
    src::T,
    dst::T,
    weight::U;
    directed=false) where {T<:Union{String,Integer}} where {U<:Number}

    haskey(g.adjacency_list, src) ? push!(g.adjacency_list[src], (dst, weight)) : g.adjacency_list[src] = [(dst, weight)]

    if !directed
        haskey(g.adjacency_list, dst) ? push!(g.adjacency_list[dst], (src, weight)) : g.adjacency_list[dst] = [(src, weight)]
    end

    return nothing
end

"""
Remove a link from WeightedGraph `g` who's directedness is `directed` between `src` and `dst` or between `src` and a list of nodes specified by `dsts`
"""
function delete_link!(
    g::WeightedGraph{T,U},
    src::T,
    dst::T;
    directed=false) where {T<:Union{String,Integer}} where {U<:Number}

    if !haskey(g.adjacency_list, src)
        error("Source $(src) node label does not exist as a key in the graph")
    end

    if !haskey(g.adjacency_list, dst) && directed
        error("Destination $(dst) node label does not exist as a key in the graph")
    end

    idx_dst = findfirst(e -> e == dst, first.(g.adjacency_list[src]))
    idx_src = findfirst(e -> e == src, first.(g.adjacency_list[dst]))

    popat!(g.adjacency_list[src], idx_dst, nothing)

    if !directed
        popat!(g.adjacency_list[dst], idx_src, nothing)
    end

    if g.adjacency_list[src] == []
        delete!(g.adjacency_list, src)
    end

    if g.adjacency_list[dst] == []
        delete!(g.adjacency_list, dst)
    end

    return nothing
end
