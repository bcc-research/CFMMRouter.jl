export Router, route!

# This is unfortunately complicated; is there a better way
struct Router{T, V, VT, O <: Objective}
    cfmms::Vector{CFMM}
    objective::O
    Δs::VT
    Λs::VT
    v::V
end

function Router(T::Type, objective::O, cfmms::Vector{CFMM}, n_tokens) where {O <: Objective}
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

Router(T, objective, n_tokens) = Router(T, objective, Vector{CFMM}(), n_tokens)
Router(objective, n_tokens) = Router(Float64, objective, n_tokens)

function route!(r::R) where {R <: Router}
    function fg!(func, g, v)
        # Update all trade parameters given prices v
        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
            find_arb!(Δ, Λ, c, v)
        end

        # Return gradient, if needed
        if g !== nothing
            g .= 0

            for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
                @views G[c.Ai] .+= Δ - Λ
            end

            return -f(r.objective, v) + dot(v, G)
        end
        # Otherwise just return the objective
        if func !== nothing
            acc = 0.0
        end

        for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.c)
            acc += @views dot(v[c.Ai], Δ) - dot(v[c.Ai], Λ)
        end

        return acc - f(r.objective, v)
    end

    return optimize(Optim.only_fg!(fg!), lower_limit(r.objective), upper_limit(r.objective), v, LBFGS())
end