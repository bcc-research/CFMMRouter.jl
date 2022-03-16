export Router, route!

# This is unfortunately complicated; is there a better way
struct Router{T, V, VT, O <: Objective}
    cfmms::Vector{CFMM}
    objective::O
    Δs::VT
    Λs::VT
    v::V
end

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
    for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
        find_arb!(Δ, Λ, c, v[c.Ai])
    end
end

function route!(r::R) where {R <: Router}
    # function fg!(func, g, v)
    #     # Update all trade parameters given prices v
    #     for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
    #         find_arb!(Δ, Λ, c, v)
    #     end

    #     # Return gradient, if needed
    #     if g !== nothing
    #         g .= 0

    #         for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
    #             @views G[c.Ai] .+= Δ - Λ
    #         end

    #         return -f(r.objective, v) + dot(v, G)
    #     end
    #     # Otherwise just return the objective
    #     if func !== nothing
    #         acc = 0.0

    #         for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
    #             acc += @views dot(v[c.Ai], Δ) - dot(v[c.Ai], Λ)
    #         end

    #         return acc - f(r.objective, v)
    #     end
    # end

    # Optimizer set up
    optimizer = L_BFGS_B(length(r.v), 17)
    v_0 = ones(length(r.v)) # We should use the marginal price here
    bounds = zeros(3, length(r.v))
    bounds[1, :] .= 2
    bounds[2, :] .= lower_limit(r.objective)
    bounds[3, :] .= upper_limit(r.objective)

    function fn(v)
        find_arb!(r, v)

        acc = 0.0

        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
            acc += @views dot(Λ, v[c.Ai]) - dot(Δ, v[c.Ai])
        end

        return f(r.objective, v) + acc
    end

    function g!(G, v)
        G .= 0

        find_arb!(r, v)

        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
            @views G[c.Ai] .+= Λ - Δ 
        end

        return -grad!(zero(v), r.objective, v) + G
    end

    _, v = optimizer(fn, g!, v_0, bounds, m=5, factr=1e7, pgtol=1e-5, iprint=1, maxfun=15000, maxiter=15000)
    find_arb!(r, v)

    return v
end