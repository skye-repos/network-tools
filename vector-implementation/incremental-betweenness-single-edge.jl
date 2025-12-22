include("../lib/graph.jl")
include("../lib/paths.jl")
using DataStructures

function find_affected_sources(
    g::Graph,
    u::Int64,
    v::Int64,
    ω::AbstractMatrix,
    d::AbstractMatrix,
    )

    N = size(g)
    affected = Set{Int64}()
    visited = falses(N)

    Q = Queue{Int64}()
    enqueue!(Q, v)

    visited[v] = true

    while !isempty(Q)
        t = dequeue!(Q)

        for s ∈ neighbors(g, t)
            if !visited[s]
                if d[s, v] ≥ d[s, u] + ω[u, v]
                    visited[s] = true
                    push!(affected, s)
                    enqueue!(Q, s)
                end
            end
        end
    end
    
    return affected
end

function iBet_betweenness_update!(
    g::Graph,
    u::Int64,
    v::Int64,
    ω::AbstractMatrix,
    d::AbstractMatrix,
    σ::AbstractMatrix,
    betweenness::Vector{Float64},
    )

    N = size(g)
    visited = falses(N)
    S = [Set{Int64}() for _ ∈ 1:N]

    new_d = deepcopy(d)
    new_σ = deepcopy(σ)

    if ω[u, v] ≤ d[u, v]
        S[v] = find_affected_sources(g, u, v, ω, d)
        d[u, v] = ω[u, v]

        Q = Queue{Int64}()
        p = Vector{Int64}(undef, N)

        enqueue!(Q, v)
        visited[v] = true
        p[v] = v

        while !isempty(Q)
            t = dequeue!(Q)

            for s ∈ S[p[t]]
                d_star = d[s, u] + ω[u, v] + d[v, t]
                if d[s, t] ≥ d_star
                    if d[s, t] > d_star
                        new_d[s, t] = d_star
                        new_σ[s, t] = 0
                    end

                    new_σ[s, t] += σ[s, u] * σ[v, t]
                    
                    if t != v
                        push!(S[t], s)
                    end
                end
            end

            for w ∈ neighbors(g, t)
                if !visited[w] && d[u, w] ≥ d[v, w] + ω[u, v]
                    enqueue!(Q, w)
                    visited[w] = true
                    p[w] = t
                end
            end
        end

        for s ∈ S[v]
            PQ = PriorityQueue{Int64, Float64}(Base.Order.Reverse)
            new_PQ = PriorityQueue{Int64, Float64}(Base.Order.Reverse)

            for q ∈ S[s]
                PQ[q] = d[s, q]
                new_PQ[q] = new_d[s, q]
            end
            
            iBet_dependency_accumulation!(g, PQ, S[s], betweenness, s, d, σ, ω, :decrease)
            iBet_dependency_accumulation!(g, new_PQ, S[s], betweenness, s, new_d, new_σ, ω, :increase)
        end
    end

    d = new_d
    σ = new_σ
    
    return betweenness, d, σ
end

function iBet_dependency_accumulation!(
    g::Graph,
    PQ::PriorityQueue{Int64, Float64},
    affected::Set{Int64},
    betweenness::Vector{Float64},
    s::Int64,
    d::AbstractMatrix,
    σ::AbstractMatrix,
    ω::AbstractMatrix,
    mode::Symbol
    )
    
    N = size(g)
    Δ = zeros(N)

    while !isempty(PQ)
        w = dequeue!(PQ)

        coeff = (mode == :decrease) ? -1.0 : 1.0
        betweenness[w] += coeff * Δ[w]

        for y ∈ neighbors(g, w)
            if y ≠ s && d[s, w] == d[s, y] + ω[y, w]
                if w ∈ affected
                    c = (σ[s, y] / σ[s, w]) * (1 + Δ[w])
                else
                    c = (σ[s, y] / σ[s, w]) * Δ[w]
                end
                if !haskey(PQ, y)
                    PQ[y] = d[s, y]
                end
                Δ[y] += c
            end
        end    
    end
end
