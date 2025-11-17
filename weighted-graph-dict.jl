using Base: delete
using SparseArrays
include("graph-dict.jl")

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
function matrix_from_weightedGraph(g::WeightedGraph{Int,U}) where {U<:Number}
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
