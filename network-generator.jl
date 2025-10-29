include("graph.jl")

"""
Computes the edge list for a connected graph with N0 nodes
"""
function connected_graph(N0::Int)
    adjacency_list::Vector{Vector{Int}} = []
    nodes = collect(1:N0)
    for i = 1:N0
        push!(adjacency_list, filter(e->e≠i, nodes))
    end
    return Graph(adjacency_list)
end

"""
Computes the next barabasi-albert iteration of a network `g` where a new node brings with it `m` half links.
"""
function ba_step!(g::Graph, m::Int64)
    # Create a pool of nodes weighted by their degree
    pool = vcat([fill(i, g.degrees[i]) for i in 1:g.N]...)
    targets = Int[]
    while length(targets) < m
        v = rand(pool)
        v ∉ targets && push!(targets, v)
    end

    for v in targets
        push!(g.adjacency_list[v], g.N+1)
        g.degrees[v] += 1
    end

    g.N += 1
    push!(g.degrees, m)
    push!(g.adjacency_list, targets)

    return g
end

"""
Given an initial network `g`, computes the Barabasi-Albert graph after `n_steps` nodes are added each of which brings with it `m` half-links
"""
function barabasi_albert!(g::Graph, m::Int64, n_steps::Int64)
    for _ = 1:n_steps
        ba_step!(g, m)
    end
    return g
end

"""
Generate Erdos-Renyi random network G(N, p)
""" 
function erdos_renyi(N, p)
    adjacency_list = [Vector{Int64}() for _ = 1:N]

    n_edges = Int(ceil(N*(N-1)*p/2))

    println(n_edges)

    if p < 0.1
        for _ = 1:n_edges
            i = Int.(rand(1:N))
            j = Int.(rand(1:N))
            while j ∈ adjacency_list[i]
                j = Int.(rand(1:N))
            end
            push!(adjacency_list[i], j)
            push!(adjacency_list[j], i)
        end
    else
        for i = 1:N
            for j = 1:N
                u = rand()
                if u < p
                    push!(adjacency_list[Int(i)], Int(j))
                    push!(adjacency_list[Int(j)], Int(i))
                end
            end
        end
    end
    return Graph(adjacency_list)
end
