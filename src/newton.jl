# We solve problems of the form
# max νᵀ(Λ - Δ) s.t. ϕ(R + γΔ - Λ) - ϕ(R) ≥ 0; Δ, Λ ≥ 0
# ⟺ min νᵀ(Λ - Δ) s.t. ϕ(R) - ϕ(R + γΔ - Λ) ≤ 0, -Δ ≤ 0, -Λ ≤ 0, 
# x ⧋ [Δ; Λ] ⟹ ϕ(R + γΔ - Λ) = ϕ(R + [γI -I]x)

mutable struct NewtonSolver{T}
    trade::Trade{T}
    rdual::Vector{T}
    rcent::Vector{T}
    x::Vector{T}
    dx::Vector{T}
    λ::Vector{T}
    dλ::Vector{T}
    μ::T,
    η̂::T,
    t::T,
    tol_feas::T,
    tol::T,
    ni::Int,
    m::Int,
    cache
    function NewtonSolver(trade::Trade{T}) where {T}
        ni = length(trade.cfmm, μ=10, tol_feas=1e-10, tol=1e-10)
        m = 2ni + 1
        cache = (
            Hpd = zeros(T, 2ni, 2ni),
            bpd = zeros(T, 2ni),
            fx = zeros(T, m),
            # TODO: make this its own efficient thing
            DfxT = sparse(1:2ni, 1:2ni, -ones(T, 2ni), 2ni, m),
            # ∇ϕ = zeros(2ni),
            Rnew = zeros(T, ni),
            x⁺ = zeros(T, 2ni),
            λ⁺ = zeros(T, m),
        )
        η̂ = Inf
        return new{T}(trade, zeros(T, 2ni, 2ni), zeros(T, 2ni), ones(T, 2ni),
                      zeros(T, m), zeros(T, m),
                      μ, η̂, zero(T), tol_feas, tol, ni, m, cache
        )
    end
end


# Updates f(x) = [f₁(x) ... fₘ(x)]ᵀ, Df(x) = [∇f₁(x) ... ∇fₘ(x)]ᵀ, rcent, rdual
function update_state!(ns::NewtonSolver)
    # 1. update cache: f(x) and Df(x)
    # Update f(x)
    @. @views ns.cache.Rnew = ns.trade.cfmm.R + ns.trade.cfmm.γ * ns.x[1:ni] - ns.x[ni+1:end]
    @. ns.cache.fx[1:2ni] = -x
    ns.cache.fx[end] = ns.trade.cfmm.ϕ(ns.trade.cfmm.R) - ns.trade.cfmm.ϕ(ns.cache.Rnew)
    
    # Update Df(x) [Note that the top 2ni x 2ni block is -I]: [-I ; ∇ϕᵀ]
    @views ns.trade.cfmm.∇ϕ!(ns.cache.∇ϕ[1:ni], ns.cache.Rnew)
    ns.cache.∇ϕ[1:ni] .*= γ
    ns.cache.∇ϕ[ni+1:end] .= -ns.cache.∇ϕ[1:ni]

    # Compute η̂ = -f(Δ, Λ)ᵀλ → update t
    η̂ = -dot(λ, ns.cache.fx)
    ns.t = ns.μ * ns.m / η̂


    # 2. Compute residuals
    @. ns.rcent = -ns.λ * ns.cache.fx - 1/t
    # rdual = ∇²f₀(x) + Df(x)ᵀλ = 
    @. @views ns.rdual = ns.λ[1:2ni] + ns.λ[end] * ns.cache.∇ϕ

    return nothing
end 


# Newtown system solve
#   [∇²f₀(x) + ∑λᵢ∇²fᵢ(x)   Df(x)ᵀ     ] [dx]   =   [∇f₀(x) + Df(x)ᵀλ]
#   [-diag(λ)Df(x)          -diag(f(x))] [dλ]       [-diag(λ)f(x) - (1/t)𝟏]
# Uses block elimination:
#   dλ = -diag(f(x))⁻¹ * (diag(λ)*Df(x)*dx - rcent)
#   ⟹ [∇²f₀(x) + ∑λᵢ∇²fᵢ(x) + ∑(λᵢ / -fᵢ(x))*∇fᵢ(x)∇fᵢ(x)ᵀ] dx = -[rdual + Df(x)ᵀ*diag(f(x))⁻¹rcent]
# where
# f₀(x) = -νᵀ[-I I]x    ⟹ ∇f₀(x) = -νᵀ[-I I],   ∇²f₀(x) = 0
# fᵢ(x) = -eᵢᵀx         ⟹ ∇fᵢ(x) = -eᵢ,         ∇²fᵢ(x) = 0     i = 1, ..., 2ni
# fᵢ(x) = ϕ(R) - ϕ(R + [γI -I]x) 
#                       ⟹ ∇fᵢ(x) = [γI -I]ᵀ*∇ϕ(R + [γI -I]x) = [γ∇ϕ(Rnew); - ∇ϕ(Rnew)]
#                       ⟹ ∇²fᵢ(x) = [γI -I]ᵀ*∇²ϕ(R + [γI -I]x)*[γI -I]
# We assume access to the gradient and hessian oracles of the CFMM
# We assume state has been updated (f(x), Df(x), rcent, rdual)
function compute_search_direction!(ns::NewtonSolver)    
    Hpd = ns.cache.Hpd
    bpd = ns.cache.bpd

    # 1. Compute dx = [∇²f₀(x) + ∑λᵢ∇²fᵢ(x) + ∑(λᵢ / -fᵢ(x))*∇fᵢ(x)∇fᵢ(x)ᵀ] \ -[rdual + Df(x)ᵀ*diag(f(x))⁻¹rcent]
    Hpd .= ns.λ[end]*[γ*I -I]'*ns.trade.cfmm.∇²ϕ(ns.cache.Rnew)*[γ*I -I]
    Hpd[diagind(Hpd)[1:2ni]] .+= -ns.λ[1:2ni] ./ ns.cache.fx[1:2ni]
    ∇fₘ = [γ*ns.trade.cfmm.∇ϕ(ns.cache.Rnew); -ns.trade.cfmm.∇ϕ(ns.cache.Rnew)]
    Hpd .+= -ns.λ[end] / ns.cache.fx[end] * ∇fₘ*∇fₘ'
    bpd .= ns.cache.DfxT * (rcent ./ ns.cache.fx)
    bpd .+= -rdual
    ldiv!(ns.dx, Hpd, bpd)

    # 2. Compute dλ = -diag(f(x))⁻¹ * (diag(λ)*Df(x)*dx - rcent)
    ns.dλ = -Diagonal(1 ./ ns.cache.fx) * (Diagonal(ns.λ) * ns.cache.DfxT' * ns.dx - ns.rcent)

    return nothing
end


# Backtracking line search
function take_step!(ns::NewtonSolver{T}; α=0.05, β=0.5) where {T}
    x = ns.x
    λ = ns.λ
    dx = ns.dx
    dλ = ns.dλ
    x⁺ = cache.x⁺
    λ⁺ = cache.λ⁺

    # Largest positive step length ≤ 1 that gives λ⁺ ≥ 0 (BV pg 613)
    smax = min(one(T), minimum(i -> dλ[i] < 0 ? -λ[i]/dλ[i] : Inf , 1:length(λ)))
    s = 0.99smax

    # Line Search
    @. x⁺ = x + s*dx
    @. λ⁺ = λ + s*dλ
    while residual(x⁺, λ⁺) > (1 - α * s) * residual(x, λ)
        s *= β
        @. x⁺ = x + s*dx
        @. λ⁺ = λ + s*dλ
    end

    # Update variables
    x .= x⁺
    λ .= λ⁺

    return nothing
end

function residual(ns, x, λ)
    Rnew = ns.trade.cfmm.R + ns.trade.cfmm.γ * x[1:ni] - x[ni+1:end]
    fx = vcat(-x, s.trade.cfmm.ϕ(ns.trade.cfmm.R) - ns.trade.cfmm.ϕ(Rnew))
    ∇fₘ = [γI -I]'*ns.trade.cfmm.∇ϕ(Rnew)
    rcent = -λ * fx - 1/ns.t
    rdual = λ[1:2ni] + λ[end] * ∇fₘ

    return sqrt(sum(x->x^2, rcent) + sum(x->x^2, rdual))
end

function solve!(ns::NewtonSolver; max_iters=100)
    iters = 0
    while (ns.η̂ > ns.tol || norm(ns.rdual) > ns.tol_feas) && iters < max_iters
        update_state!(ns)
        compute_search_direction!(ns)
        take_step!(ns)

        iters += 1
    end
end