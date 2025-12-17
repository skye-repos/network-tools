include("graph.jl")
include("weighted-graph.jl")
using DataStructures

struct DijkstraState{T<:Union{String, Integer}, U<:Number}
    parents::Dict{T,T}
    dists::Dict{T,U}
    predecessors::Dict{T,Vector{T}}
    path_counts::Dict{T,Float64}
    closest_vertices::Vector{T}
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
            alt = d + get_weight(g, u, v)

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

    for src ∈ srcs
        pathcounts[src] = one(Float64)
        parents[src] = 0
        empty!(preds[src])
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

function dijkstra(
    g::WeightedGraph{T,U},
    src::Integer;
    allpaths=false,
    trackvertices=false,
    maxdist=typemax(U)
) where {T<:Integer} where {U<:Number}

    return dijkstra(
        g, [src;]; allpaths=allpaths, trackvertices=trackvertices, maxdist=maxdist)
end

function betweenness_centrality(
    g::WeightedGraph{T,U};
    nodes::Vector{T}=collect(nodes(g)),
    normalize=true,
    endpoints=false,
) where {T<:Union{String,Integer}} where {U<:Number}

    N = size(g)
    k = length(nodes)

    betweenness = DefaultDict{T,Float64}(0)

    for src ∈ nodes
        if degrees(g)[src] > 0
            state = dijkstra(g, src, allpaths=true, trackvertices=true)
            if endpoints
                _accumulate_endpoints!(betweenness, state, src)
            else
                _accumulate_basic!(betweenness, state, src)
            end
        end
    end

    _rescale!(betweenness, N, normalize, false, k)

    return betweenness
end

function _accumulate_basic!(
    betweenness::AbstractDict{T,Float64},
    state::DijkstraState,
    si::T
) where {T<:Union{String,Integer}}

    δ = DefaultDict{T,Float64}(0)
    σ = state.path_counts
    P = state.predecessors

    P[si] = []

    S = reverse(state.closest_vertices)
    for w ∈ S
        coeff = (1.0 + δ[w]) / σ[w]
        for v ∈ P[w]
            if v > 0
                δ[v] += (σ[v] * coeff)
            end
        end

        if w ≠ si
            betweenness[w] += δ[w]
        end
    end

    return nothing
end

function _accumulate_endpoints!(
    betweenness::AbstractDict{T,Float64},
    state::DijkstraState,
    si::T
) where {T<:Union{String,Integer}}

    δ = DefaultDict{T,Float64}(0)
    σ = state.path_counts
    P = state.predecessors

    S = reverse(state.closest_vertices)
    betweenness[si] += length(S) - 1

    for w ∈ S
        coeff = (1.0 + δ[w]) / σ[w]
        for v ∈ P[w]
            δ[v] += σ[v] * coeff
        end

        if w ≠ si
            betweenness[w] += (δ[w] + 1)
        end
    end


    return nothing
end

function _rescale!(
    betweenness::AbstractDict{T,Float64},
    N::Integer,
    normalize::Bool,
    directed::Bool,
    k::Integer
) where {T<:Union{String,Integer}}

    scale = nothing
    if normalize
        if 2 < N
            scale = 1.0 / ((N - 1) * (N - 2))
        end
    elseif !directed
        scale = 0.5
    end

    if !isnothing(scale)
        if k > 0
            scale = scale * N / k
        end
        map!(x -> x * scale, values(betweenness))
    end

    return nothing
end

function components(label::Vector{T}) where {T<:Integer}
    d = Dict{T,T}()
    c = Vector{Vector{T}}()
    i = one(T)

    for (v, l) in enumerate(label)
        index = get!(d, l, i)
        if length(c) ≥ index
            push!(c[index], v)
        else
            push!(c, [v])
            i += 1
        end
    end

    return c, d
end

function connected_componenents(g::AbstractGraph)
    N = maximum(nodes(g))
    label = zeros(Int, N)

    for u ∈ nodes(g)
        label[u] != 0 && continue
        label[u] = u
        Q = Vector{Int}()
        push!(Q, u)

        while !isempty(Q)
            src = popfirst!(Q)
            for v in neighbors(g, src)
                if label[v] == 0
                    push!(Q, v)
                    label[v] = u
                end
            end
        end
    end

    c, _ = components(label)

    return c
end

function giant_component(g::AbstractGraph)
    components = connected_componenents(g)
    _, i = findmax(length.(components))

    return components[i]
end
