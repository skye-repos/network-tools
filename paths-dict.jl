include("graph-dict.jl")
include("weighted-graph-dict.jl")
using DataStructures

struct DijkstraState{T<:Integer,U<:Number}
    parents::Dict{T,U}
    dists::Dict{T,U}
    predecessors::Dict{T,Vector{U}}
    path_counts::Dict{T,Float64}
    closest_vertices::Vector{U}
end

function dijkstra(
    g::WeightedGraph{T,U},
    srcs::Vector{T};
    allpaths=false,
    trackvertices=false,
    maxdist=typemax{U}
) where {U<:Number} where {T<:Union{String,Integer}}

    N = size(g)
    dists = DefaultDict{T,U}(maxdist)
    parents = DefaultDict{T,T}(zero(T))
    visited = DefaultDict{T,Bool}(false)

    pathcounts = DefaultDict{T,Float64}(zero(Float64))
    preds = DefaultDict{T,Vector{T}}(Vector{T}())
    H = PriorityQueue{T,U}()

    for src ∈ srcs
        dists[src] = zero(U)
        visited[src] = true
        pathcounts[src] = one(Float64)
        H[src] = zero(U)
    end

    closest_vertices = Vector{T}()
    sizehint!(closest_vertices, N)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u]

        for v ∈ neighbors(g, u)
            alt = d + g.adjacency_list[u][v]

            alt > maxdist && continue


            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]

                if allpaths
                    preds[v] = [u;]
                end

                H[v] = alt

            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                pathcounts[v] = pathcounts[u]

                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end

                H[v] = alt

            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]

                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for i ∈ nodes(g)
            if !visited[i]
                push!(closest_vertices, i)
            end
        end
    end

    for src ∈ srcs
        pathcounts[src] = one(Float64)
        parents[src] = 0
        empty!(preds[src])
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

function dijkstra(
    g::WeightedGraph{T,U},
    src::Integer,
    allpaths=false,
    trackvertices=false,
    maxdist=typemax(U)
) where {T<:Integer} where {U<:Number}

    return dijkstra(
        g, [src;]; allpaths=allpaths, trackvertices=trackvertices, maxdist=maxdist)
end
