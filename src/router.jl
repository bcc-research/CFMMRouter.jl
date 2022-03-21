export Router, route!
export netflows!, netflows, update_reserves!

struct Router{T, V, VT, O <: Objective}
    cfmms::Vector{CFMM}
    objective::O
    Δs::VT
    Λs::VT
    v::V
end

"""
    Router(objective, cfmms, n_tokens)

Constructs a router that finds a trade (router.Δs, router.Λs) through `cfmms` 
which maximizes `objective`. The number of tokens `n_tokens` must be specified.
"""
function Router(objective::O, cfmms::Vector{C}, n_tokens) where {T, O <: Objective, C <: CFMM{T}}
    V = Vector{T}
    VT = Vector{AbstractVector{T}}
    Δs = VT()
    Λs = VT()

    for c in cfmms
        push!(Δs, zero(c.R))
        push!(Λs, zero(c.R))
    end

    return Router{T, V, VT, O}(
        cfmms,
        objective,
        Δs,
        Λs,
        zeros(T, n_tokens)
    )
end

Router(objective, n_tokens) = Router(objective, Vector{CFMM{Float64}}(), n_tokens)

function find_arb!(r::Router, v)
    Threads.@threads for i in 1:length(r.Δs)
        find_arb!(r.Δs[i], r.Λs[i], r.cfmms[i], v[r.cfmms[i].Ai])
    end
end

"""
    find_arb!(Δ, Λ, cfmm, v)

Solves the arbitrage problem for `cfmm` given price vector `v`,
```math
\begin{array}{ll}
\text{minimize} & \nu^T(\Lambda - \Delta) \\
\text{subject to} & \varphi(R + \gamma\Delta - \Lambda) = \varphi(R) \\
& \Delta, \Lambda \geq 0.
\end{array}
Overwrites the variables `Δ` and `Λ`.
```
"""
function route!(r::R; verbose=false) where {R <: Router}
    # Optimizer set up
    optimizer = L_BFGS_B(length(r.v), 17)
    r.v .= 2*ones(length(r.v)) # We should use the initial marginal price here

    bounds = zeros(3, length(r.v))
    bounds[1, :] .= 2
    bounds[2, :] .= lower_limit(r.objective)
    bounds[3, :] .= upper_limit(r.objective)

    # Objective function
    function fn(v)
        if !all(v .== r.v)
            find_arb!(r, v)
            r.v .= v
        end

        acc = 0.0

        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
            acc += @views dot(Λ, v[c.Ai]) - dot(Δ, v[c.Ai])
        end

        return f(r.objective, v) + acc
    end

    # Derivative of objective function
    function g!(G, v)
        G .= 0

        if !all(v .== r.v)
            find_arb!(r, v)
            r.v .= v
        end
        grad!(G, r.objective, v)

        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
            @views G[c.Ai] .+= Λ - Δ 
        end

    end

    find_arb!(r, r.v)
    _, v = optimizer(fn, g!, r.v, bounds, m=5, factr=1e1, pgtol=1e-5, iprint=verbose ? 1 : -1, maxfun=15000, maxiter=15000)
    r.v .= v
    find_arb!(r, v)
end

# ----- Convenience functions
function netflows!(ψ, r::Router)
    fill!(ψ, 0)

    for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
        ψ[c.Ai] += Λ - Δ
    end

    return nothing
end

function netflows(r::Router)
    ψ = zero(v)
    netflows!(ψ, r)
    return ψ
end

function update_reserves!(r::Router)
    for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
        c.R .+= Δ - Λ
    end

    return nothing
end