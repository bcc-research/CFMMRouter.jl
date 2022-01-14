# We solve problems of the form
# max νᵀ(Λ - Δ) s.t. ϕ(R + γΔ - Λ) - ϕ(R) ≥ 0; Δ, Λ ≥ 0
# ⟺ min νᵀ(Λ - Δ) s.t. ϕ(R) - ϕ(R + γΔ - Λ) ≤ 0, -Δ ≤ 0, -Λ ≤ 0, 
# x ⧋ [Δ; Λ] ⟹ ϕ(R + γΔ - Λ) = ϕ(R + [γI -I]x)

mutable struct NewtonSolver{T}
    cfmm::CFMM{T}
    ν::Vector{T}
    rdual::Vector{T}
    rcent::Vector{T}
    x::Vector{T}
    dx::Vector{T}
    λ::Vector{T}
    dλ::Vector{T}
    μ::T
    η̂::T
    t::T
    tol_feas::T
    tol::T
    cache
    function NewtonSolver(cfmm::CFMM{T}, ν::Vector{T}; μ=10, tol_feas=1e-10, tol=1e-10) where {T}
        ni = length(cfmm)
        m = 2ni + 1
        cache = (
            Hpd = zeros(T, 2ni, 2ni),
            bpd = zeros(T, 2ni),
            fx = zeros(T, m),
            # TODO: make this its own efficient thing with custom mul! & mul!(_, adjoint(Dfx), x)
            Dfx = vcat(Matrix(I, 2ni, 2ni), ones(2ni)'), #sparse(1:2ni, 1:2ni, -ones(T, 2ni), 2ni, m),
            ∇fₘ = zeros(2ni),
            Rnew = zeros(T, ni),
            x⁺ = zeros(T, 2ni),
            λ⁺ = zeros(T, m),
            rdual⁺ = zeros(T, 2ni),
            rcent⁺ = zeros(T, m)
        )
        η̂ = T(Inf)
        return new{T}(cfmm, ν,
                      zeros(T, 2ni), zeros(T, m),
                      vcat(ones(T, ni), T(0.5)*ones(T, ni)), zeros(T, 2ni),
                      ones(T, m), zeros(T, m),
                      T(μ), η̂, zero(T), T(tol_feas), T(tol),
                      cache
        )
    end
end

n_constraints(ns::NewtonSolver) = length(ns.λ)
n_coins(ns::NewtonSolver) = length(ns.cfmm)


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
    Dfx = ns.cache.Dfx

    m, ni = n_constraints(ns), n_coins(ns)
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
    Dfx[end, :] .= ∇fₘ

    # Compute η̂ = -f(Δ, Λ)ᵀλ → update t
    ns.η̂ = -dot(λ, fx)
    ns.t = ns.μ * m / ns.η̂


    # 2. Compute residuals
    @. rcent = -λ * fx - 1/ns.t
    # rdual = ∇f₀(x) + Df(x)ᵀλ = 
    rdual[1:ni] .= ns.ν
    rdual[ni+1:2ni] .= -ns.ν
    @. @views rdual .+= λ[1:2ni] + λ[end] * ∇fₘ

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
function compute_search_direction!(ns::NewtonSolver{T}) where {T}    
    Hpd = ns.cache.Hpd
    bpd = ns.cache.bpd
    γ = ns.cfmm.γ
    dx = ns.dx
    λ = ns.λ
    ∇²ϕ! = ns.cfmm.∇²ϕ!
    ∇ϕ! = ns.cfmm.∇ϕ!
    fx = ns.cache.fx
    Rnew = ns.cache.Rnew
    Dfx = ns.cache.Dfx
    rcent, rdual = ns.rcent, ns.rdual

    # 1. Compute dx = [∇²f₀(x) + ∑λᵢ∇²fᵢ(x) + ∑(λᵢ / -fᵢ(x))*∇fᵢ(x)∇fᵢ(x)ᵀ] \ -[rdual + Df(x)ᵀ*diag(f(x))⁻¹rcent]
    ni = length(ns.cfmm)
    Hmid = zeros(T, ni, ni)
    ∇²ϕ!(Hmid, Rnew)
    Hpd .= λ[end]*[γ^2*Hmid -γ*Hmid; -γ*Hmid Hmid]
    Hpd[diagind(Hpd)[1:2ni]] .+= -λ[1:2ni] ./ fx[1:2ni]
    grad_cache = zeros(T, ni)
    ∇ϕ!(grad_cache, Rnew)
    ∇fₘ = [γ*grad_cache; -grad_cache]
    Hpd .+= -λ[end] / fx[end] * ∇fₘ*∇fₘ'
    bpd .= Dfx' * (rcent ./ fx)
    bpd .+= -rdual
    # bunchkaufman!(Symmetric(Hpd))
    dx .= Hpd \ bpd

    # 2. Compute dλ = -diag(f(x))⁻¹ * (diag(λ)*Df(x)*dx - rcent)
    ns.dλ = -Diagonal(1 ./ fx) * (Diagonal(λ) * Dfx * dx - rcent)

    return nothing
end


# Backtracking line search
function take_step!(ns::NewtonSolver{T}; α=0.05, β=0.5) where {T}
    x, λ = ns.x, ns.λ
    dx, dλ = ns.dx, ns.dλ
    x⁺, λ⁺ = ns.cache.x⁺, ns.cache.λ⁺
    rcent, rdual = ns.rcent, ns.rdual

    # Largest positive step length ≤ 1 that gives λ⁺ ≥ 0 (BV pg 613)
    smax = min(one(T), minimum(i -> dλ[i] < 0 ? -λ[i]/dλ[i] : Inf , 1:length(λ)))
    s = 0.99smax

    # Compute current residual
    residual_current = sqrt(sum(x->x^2, rcent) + sum(x->x^2, rdual))

    # Line Search
    @. x⁺ = x + s*dx
    @. λ⁺ = λ + s*dλ
    iter = 0
    #TODO: caching efficiency
    # x in particular
    while f_x⁺_geq_0(ns, x⁺) || residual(ns, x⁺, λ⁺) > (1 - α * s) * residual_current
        s *= β
        @. x⁺ = x + s*dx
        @. λ⁺ = λ + s*dλ
        iter += 1
        iter > 10 && (Main.x[] = smax; error("Maximum iterations hit in line search"))
    end

    # Update variables
    x .= x⁺
    λ .= λ⁺

    return nothing
end

@inline function f_x⁺_geq_0(ns, x⁺)
    return any(xi -> -xi ≥ 0, x⁺) || ns.cfmm.ϕ(ns.cfmm.R) - ns.cfmm.ϕ(ns.cfmm.R + ns.cfmm.γ * x⁺[1:n_coins(ns)] - x⁺[n_coins(ns)+1:end]) ≥ 0
end
    

function residual(ns::NewtonSolver{T}, x, λ) where {T}
    R, γ = ns.cfmm.R, ns.cfmm.γ
    ϕ, ∇ϕ! = ns.cfmm.ϕ, ns.cfmm.∇ϕ!
    ni = n_coins(ns)
    t = ns.t
    Rnew, fx, ∇fₘ = ns.cache.Rnew, ns.cache.fx, ns.cache.∇fₘ
    rcent⁺, rdual⁺ = ns.cache.rcent⁺, ns.cache.rdual⁺

    @views @. Rnew = R + γ * x[1:ni] - x[ni+1:end]
    fx = vcat(-x, ϕ(R) - ϕ(Rnew))
    grad_cache = zeros(T, ni)
    ∇ϕ!(grad_cache, Rnew)
    ∇fₘ .= [γ*grad_cache; -grad_cache]
    @. rcent⁺ = -λ * fx - 1/t
    Main.x[] = λ
    @views @. rdual⁺ = λ[1:2ni] + λ[end] * ∇fₘ

    return sqrt(sum(x->x^2, rcent⁺) + sum(x->x^2, rdual⁺))
end

# TODO: looks like rdual and rcent have an issue
function solve!(ns::NewtonSolver; max_iters=100, verbose=false)
    if verbose
        headers = ["iter", "rdual", "rcent", "dual gap", "obj", "time"]
        solve_time_start =  time_ns()
        print_header(headers)
    end

    iters = zero(UInt)
    while (ns.η̂ > ns.tol || norm(ns.rdual) > ns.tol_feas) && iters < max_iters
        update_state!(ns)
        compute_search_direction!(ns)
        take_step!(ns)

        iters += 1
        if verbose
            ni = n_coins(ns)
            arbed = ns.cache.x⁺
            @views @. arbed[1:ni] = ns.x[1:ni] - ns.x[ni+1:2ni]
            print_iter_func((
                string(iters),
                norm(ns.rdual),
                norm(ns.rcent),
                ns.η̂,
                dot(ns.ν, arbed[1:ni]),
                (time_ns() - solve_time_start) / 1e9
            ))
        end
    end
    verbose && print_footer()

    return nothing
end