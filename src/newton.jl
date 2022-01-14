# We solve problems of the form
# max νᵀ(Λ - Δ) s.t. ϕ(R + γΔ - Λ) - ϕ(R) ≥ 0; Δ, Λ ≥ 0
# ⟺ min νᵀ(Λ - Δ) s.t. ϕ(R) - ϕ(R + γΔ - Λ) ≤ 0, -Δ ≤ 0, -Λ ≤ 0, 
# x ⧋ [Δ; Λ] ⟹ ϕ(R + γΔ - Λ) = ϕ(R + [γI -I]x)

mutable struct NewtonSolver{T}
    cfmm::CFMM{T}
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
    function NewtonSolver(cfmm::CFMM{T}; μ=10, tol_feas=1e-10, tol=1e-10) where {T}
        ni = length(cfmm)
        m = 2ni + 1
        cache = (
            Hpd = zeros(T, 2ni, 2ni),
            bpd = zeros(T, 2ni),
            fx = zeros(T, m),
            # TODO: make this its own efficient thing with custom mul! & mul!(_, adjoint(Dfx), x)
            DfxT = vcat(Matrix(I, 2ni, 2ni), ones(2ni)'), #sparse(1:2ni, 1:2ni, -ones(T, 2ni), 2ni, m),
            ∇fₘ = zeros(2ni),
            Rnew = zeros(T, ni),
            x⁺ = zeros(T, 2ni),
            λ⁺ = zeros(T, m),
        )
        η̂ = Inf
        return new{T}(cfmm, zeros(T, 2ni, 2ni), zeros(T, 2ni), ones(T, 2ni),
                      zeros(T, m), zeros(T, m),
                      μ, η̂, zero(T), tol_feas, tol, ni, m, cache
        )
    end
end


# Updates f(x) = [f₁(x) ... fₘ(x)]ᵀ, Df(x) = [∇f₁(x) ... ∇fₘ(x)]ᵀ, rcent, rdual
function update_state!(ns::NewtonSolver)
    R = ns.cfmm.R
    γ = ns.cfmm.γ
    ϕ = ns.cfmm.ϕ
    ∇ϕ! = ns.cfmm.∇ϕ!
    Rnew = ns.cache.Rnew

    x, λ = ns.x, ns.λ
    fx = ns.cache.fx
    ∇fₘ = ns.cache.∇fₘ

    t, μ = ns.t, ns.μ
    m, ni = ns.m, ns.ni
    rcent, rdual = ns.rcent, ns.rdual

    # 1. update cache: f(x) and Df(x)
    # Update f(x)
    @. @views Rnew = R + γ * x[1:ni] - x[ni+1:end]
    @. fx[1:2ni] = -x
    fx[end] = ϕ(R) - ϕ(Rnew)
    
    # Update Df(x) [Note that the top 2ni x 2ni block is -I]: [-I ; ∇ϕᵀ]
    @views ∇ϕ!(∇fₘ[1:ni], Rnew)
    ∇fₘ[1:ni] .*= γ
    ∇fₘ[ni+1:end] .= -∇fₘ[1:ni]

    # Compute η̂ = -f(Δ, Λ)ᵀλ → update t
    η̂ = -dot(λ, fx)
    t = μ * m / η̂


    # 2. Compute residuals
    @. rcent = -λ * fx - 1/t
    # rdual = ∇²f₀(x) + Df(x)ᵀλ = 
    @. @views rdual = λ[1:2ni] + λ[end] * ∇fₘ

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
    dx = ns.dx
    λ = ns.λ
    ∇²ϕ = ns.cfmm.∇²ϕ
    ∇ϕ = ns.cfmm.∇ϕ
    fx = ns.cache.fx
    Rnew = ns.cache.Rnew
    DfxT = ns.cache.DfxT
    rcent = ns.rcent

    # 1. Compute dx = [∇²f₀(x) + ∑λᵢ∇²fᵢ(x) + ∑(λᵢ / -fᵢ(x))*∇fᵢ(x)∇fᵢ(x)ᵀ] \ -[rdual + Df(x)ᵀ*diag(f(x))⁻¹rcent]
    Hpd .= λ[end]*[γ*I -I]'*∇²ϕ(Rnew)*[γ*I -I]
    Hpd[diagind(Hpd)[1:2ni]] .+= -λ[1:2ni] ./ fx[1:2ni]
    ∇fₘ = [γ*∇ϕ(Rnew); -∇ϕ(Rnew)]
    Hpd .+= -λ[end] / fx[end] * ∇fₘ*∇fₘ'
    bpd .= DfxT * (rcent ./ fx)
    bpd .+= -rdual
    ldiv!(dx, Hpd, bpd)

    # 2. Compute dλ = -diag(f(x))⁻¹ * (diag(λ)*Df(x)*dx - rcent)
    ns.dλ = -Diagonal(1 ./ fx) * (Diagonal(λ) * DfxT' * dx - rcent)

    return nothing
end


# Backtracking line search
function take_step!(ns::NewtonSolver{T}; α=0.05, β=0.5) where {T}
    x, λ = ns.x, ns.λ
    dx, dλ = ns.dx, ns.dλ
    x⁺, λ⁺ = cache.x⁺, cache.λ⁺
    rcent, rdual = ns.rcent, ns.rdual

    # Largest positive step length ≤ 1 that gives λ⁺ ≥ 0 (BV pg 613)
    smax = min(one(T), minimum(i -> dλ[i] < 0 ? -λ[i]/dλ[i] : Inf , 1:length(λ)))
    s = 0.99smax

    # Compute current residual
    residual_current = sqrt(sum(x->x^2, rcent) + sum(x->x^2, rdual))

    # Line Search
    @. x⁺ = x + s*dx
    @. λ⁺ = λ + s*dλ
    while residual(x⁺, λ⁺) > (1 - α * s) * residual_current
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
    R, γ = ns.cfmm.R, ns.cfmm.γ
    ϕ, ∇ϕ = ns.cfmm.ϕ, ns.cfmm.∇ϕ
    x, λ = ns.x, ns.λ
    ni = ns.ni
    t = ns.t
    Rnew, fx, ∇fₘ = ns.cache.Rnew, ns.cache.fx, ns.cache.∇fₘ
    rcent⁺, rdual⁺ = ns.cache.rcent⁺, ns.cache.rdual⁺

    @views @. Rnew = R + γ * x[1:ni] - x[ni+1:end]
    fx = vcat(-x, ϕ(R) - ϕ(Rnew))
    ∇fₘ = [γI -I]'*∇ϕ(Rnew)
    rcent⁺ = -λ * fx - 1/t
    @views rdual⁺ = λ[1:2ni] + λ[end] * ∇fₘ

    return sqrt(sum(x->x^2, rcent⁺) + sum(x->x^2, rdual⁺))
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