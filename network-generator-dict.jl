using SparseArrays: SparseMatrixCSCSymmHerm
include("graph-dict.jl")
include("weighted-graph-dict.jl")
using StatsBase: sample, Weights

"""
Creates a completely connected network with `N0` nodes
"""
function connected_graph(N0::Int)
    g = Graph()
    for i ∈ 1:N0, j ∈ i+1:N0
        add_link!(g, i, j)
    end

    return g
end

"""
Computes the next Barabasi-Albert iteration of a network `g` where each new node brings `m` half-links
"""
function ba_step!(g::Graph{T}, m::Int64) where {T<:Union{String,Integer}}
    # Create a pool of nodes weighted by their degree
    old_nodes = collect(nodes(g))
    w = collect(values(degrees(g)))
    new_id = convert(T, size(g) + 1)
    targets = sample(old_nodes, Weights(w), m, replace=false)

    for v ∈ targets
        add_link!(g, new_id, v)
    end

    return nothing
end

"""
Given an initial network `g`, compute the Barabasi-Albert evolution of `g` after `n_steps` nodes are added each of which brings `m` half-links
"""
function barabasi_albert!(g::Graph, m::Integer, n_steps::Integer)
    for _ ∈ 1:n_steps
        ba_step!(g, m)
    end

    return nothing
end

"""
Generate Erdos-Renyi random network G(N, p)
"""
function erdos_renyi(N, p)
    g = Graph()

    n_edges = Int(ceil(N * (N - 1) * p / 2))
    adj = g.adjacency_list

    if p < 0.1
        for _ = 1:n_edges
            i = Int(rand(1:N))
            j = Int(rand(1:N))

            if i ∉ nodes(g)
                add_link!(g, i, j)
            elseif j ∉ nodes(g)
                add_link!(g, i, j)
            else
                while j ∈ adj[i]
                    j = Int(rand(1:N))
                end
                add_link!(g, i, j)
            end
        end
    else
        for i = 1:N
            for j = 1:N
                u = rand()
                if u < p
                    add_link!(g, Int(i), Int(j))
                end
            end
        end
    end
    return g
end

"""
Generate a ζ-model network of side length `L`, characteristic link-length of
`ζ`, and average degree `k`
"""

function ζ_model(L::Integer, ζ::Number, k_avg::Number)
    local N = L^2

    coords = Vector{Tuple{Int, Int}}()

    for i ∈ 1:L, j ∈ 1:L
        push!(coords, (i, j))
    end

    @inline function distance(a::Tuple{Int, Int}, b::Tuple{Int, Int})
        x = min(abs(a[1] - b[1]), L - abs(a[1] - b[1]))
        y = min(abs(a[2] - b[2]), L - abs(a[2] - b[2]))
        return sqrt(x^2 + y^2)
    end

    P = zeros(N, N)
    @inbounds for j ∈ 1:N, i ∈ 1:j-1
        P[i, j] = exp(-distance(coords[i], coords[j]) / ζ)
    end

    P .*= (N * k_avg / 2) / sum(P)

    g = Graph()

    @inbounds for j ∈ 1:N, i ∈ 1:j-1
        rand() < P[i, j] && add_link!(g, Int(i), Int(j))
    end

    return g
end

# @time grph = ζ_model(150, 15, 5);
# @time wts = θ_weights(grph, θ=0.5);
# @time sprs_wts = sparse(wts)
# @time wtd_grph = weights_to_weightedGraph(sprs_wts)

# varinfo(r"grph|wtd_grph|wts|sprs_wts")
