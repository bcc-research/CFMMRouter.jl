export CFMM, ProductTwoCoin, GeometricMeanTwoCoin, UniV3
export find_arb!
export update_reserves!

abstract type CFMM{T} end

@def add_generic_fields begin
    R::Vector{T}
    γ::T
    Ai::Vector{Int}                 # idx vector: jth coin in CFMM is Ai[j]
end

@def add_two_coin_fields begin
    R::Vector{T}
    γ::T
    Ai::Vector{Int}
end

Base.length(c::CFMM) = length(c.Ai)

# This def is for the docstring
@doc raw"""
    find_arb!(Δ, Λ, cfmm, v)

Solves the arbitrage problem for `cfmm` given price vector `v`,
```math
\begin{array}{ll}
\text{minimize} & \nu^T(\Lambda - \Delta) \\
\text{subject to} & \varphi(R + \gamma\Delta - \Lambda) = \varphi(R) \\
& \Delta, \Lambda \geq 0.
\end{array}
```
Overwrites the variables `Δ` and `Λ`.
"""
function find_arb! end

@doc raw"""
    ϕ(c::CFMM)

Computes the trading function for CFMM `c`.
"""
function ϕ end

@doc raw"""
    ∇ϕ!(x, c::CFMM)

Computes the gradient of the trading function for CFMM `c`.
The result is stored in `x`.
"""
function ∇ϕ! end


# ------------------------------------------------------------------------------
# |                              CFMM Definitions                              |
# ------------------------------------------------------------------------------
# Each CFMM needs to implement its find_arb! function
struct Product{T} <: CFMM{T}
    @add_generic_fields
end

struct GeometricMean{T} <: CFMM{T}
    @add_generic_fields
    w::Vector{T}
end

struct Curve{T} <: CFMM{T}
    @add_generic_fields
    α::T
    β::T
end

# Two coin specific cases
function two_coin_check_cast(R, γ, idx)
    length(R) != 2 && throw(ArgumentError("length of R must be 2 for *TwoCoin constructors"))
    length(idx) != 2 && throw(ArgumentError("length of idx must be 2 for *TwoCoin constructors"))

    T = eltype(R)

    if T <: Integer
        T = Float64
    end

    γ_T = convert(T, γ)
    idx_uint = convert.(UInt, idx)

    return γ_T, idx_uint, T
end

@doc raw"""
    ProductTwoCoin(R, γ, idx)

Creates a two coin product CFMM with coins `idx[1]` and `idx[2]`, reserves `R`,
and fee `γ`. Specifically, the invariant is
```math
\varphi(R) = R_1R_2.
```
"""
struct ProductTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
    function ProductTwoCoin(R, γ, idx)
        γ_T, idx_uint, T = two_coin_check_cast(R, γ, idx)
        return new{T}(
            MVector{2,T}(R),
            γ_T,
            MVector{2,UInt}(idx_uint)
        )
    end
end

function ϕ(cfmm::ProductTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    return R[1] * R[2]
end
function ∇ϕ!(R⁺, cfmm::ProductTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    R⁺[1] = R[2]
    R⁺[2] = R[1]
    return nothing
end

# See App. A of "An Analysis of Uniswap Markets"
@inline prod_arb_δ(m, r, k, γ) = max(sqrt(γ*m*k) - r, 0)/γ
@inline prod_arb_λ(m, r, k, γ) = max(r - sqrt(k/(m*γ)), 0)

# Solves the maximum arbitrage problem for the two-coin constant product case.
# Assumes that v > 0 and γ > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::ProductTwoCoin{T}, v::VT) where {T, VT<:AbstractVector{T}}
    R, γ = cfmm.R, cfmm.γ
    k = R[1]*R[2]

    Δ[1] = prod_arb_δ(v[2]/v[1], R[1], k, γ)
    Δ[2] = prod_arb_δ(v[1]/v[2], R[2], k, γ)

    Λ[1] = prod_arb_λ(v[1]/v[2], R[1], k, γ)
    Λ[2] = prod_arb_λ(v[2]/v[1], R[2], k, γ)
    return nothing
end

@doc raw"""
    GeometricMeanTwoCoin(R, w, γ, idx)

Creates a two coin geometric mean CFMM with coins `idx[1]` and `idx[2]`, 
reserves `R`, fee `γ`, and weights `w` such that `w[1] + w[2] == 1.0`.
Specifically, the invariant is
```math
\varphi(R) = R_1^{w_1}R_2^{w_2}.
```
"""
struct GeometricMeanTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
    w::SVector{2,T}
    function GeometricMeanTwoCoin(R, w, γ, idx)
        γ_T, idx_uint, T = two_coin_check_cast(R, γ, idx)
    
        return new{T}(
            MVector{2,T}(R),
            γ_T,
            MVector{2,UInt}(idx_uint),
            SVector{2,T}(w),
        )
    end
end

function ϕ(cfmm::GeometricMeanTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    w = cfmm.w
    return R[1]^w[1] * R[2]^w[2]
end
function ∇ϕ!(R⁺, cfmm::GeometricMeanTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    w = cfmm.w
    R⁺[1] = w[1] * (R[2]/R[1])^w[2]
    R⁺[2] = w[2] * (R[1]/R[2])^w[1]
    return nothing
end

@inline geom_arb_δ(m, r1, r2, η, γ) = max((γ*m*η*r1*r2^η)^(1/(η+1)) - r2, 0)/γ
@inline geom_arb_λ(m, r1, r2, η, γ) = max(r1 - ((r2*r1^(1/η))/(η*γ*m))^(η/(1+η)), 0)

# Solves the maximum arbitrage problem for the two-coin geometric mean case.
# Assumes that v > 0 and w > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::GeometricMeanTwoCoin{T}, v::VT) where {T, VT<:AbstractVector{T}}
    R, γ, w = cfmm.R, cfmm.γ, cfmm.w

    η = w[1]/w[2]

    Δ[1] = geom_arb_δ(v[2]/v[1], R[2], R[1], η, γ)
    Δ[2] = geom_arb_δ(v[1]/v[2], R[1], R[2], 1/η, γ)

    Λ[1] = geom_arb_λ(v[1]/v[2], R[1], R[2], 1/η, γ)
    Λ[2] = geom_arb_λ(v[2]/v[1], R[2], R[1], η, γ)
    return nothing
end


#The price in each tick refers to the lower bound on the interval
#so the intervals are [t_i,t_i+1] with liquidity L_i
#         p2---------
# p1------     |     p3-------
#   |          |         |
#   L1         L2        L3
#   |          |         | 

mutable struct UniV3{T} <: CFMM{T}
    current_price::T
    current_tick::Int
    lower_ticks::Vector{T}
    liquidity::Vector{T}
    γ::T
    Ai::Vector{Int}
end

tick_high_price(cfmm::UniV3{T}, idx) where T = cfmm.lower_ticks[idx]
tick_high_price(cfmm) = tick_high_price(cfmm, cfmm.current_tick)

function tick_low_price(cfmm::UniV3{T}, idx) where T
    if idx < length(cfmm.lower_ticks) 
        return cfmm.lower_ticks[idx + 1]
    end
    return zero(T)
end
tick_low_price(cfmm) = tick_low_price(cfmm, cfmm.current_tick)

struct BoundedProduct{T}
    k::T
    α::T
    β::T
    R_1::T
    R_2::T
end

max_price(t::BoundedProduct{T}) where T = t.k/(t.α^2)
min_price(t::BoundedProduct{T}) where T = (t.β^2)/t.k
curr_price(t::BoundedProduct{T}) where T = (t.R_2 + t.β)/(t.R_1 + t.α)

# Computes the properties of a specific tick (which is just a bounded product
# CFMM) at the specified index
function compute_at_tick(cfmm::UniV3{T}, idx) where T
    k = cfmm.liquidity[idx]
    pminus = tick_low_price(cfmm, idx)
    pplus = tick_high_price(cfmm, idx)
    α = sqrt(k/pplus)
    β = sqrt(k*pminus)

    if idx > cfmm.current_tick
        p = pplus
    elseif idx < cfmm.current_tick
        p = pminus
    else
        p = cfmm.current_price
    end

    R_1 = sqrt(k/p) - α
    if iszero(R_1 + α)
        R_2 = 0.0
    else
        R_2 = k/(R_1 + α) - β
    end

    return BoundedProduct{T}(k, α, β, R_1, R_2)
end

# Returns ([if the tick was saturated], δ, λ)
function find_arb_pos(t::BoundedProduct{T}, price) where T
    if iszero(t.k)
        return true, 0.0, 0.0
    end
    δ = sqrt(t.k/price) - (t.R_1 + t.α)

    if δ <= 0
        return true, 0.0, 0.0
    end

    δ_max = (t.R_1 + t.α)*(t.R_2/t.β) 
    if δ >= δ_max
        return true, δ_max, t.R_2
    end
    
    λ = (t.R_2 + t.β) - sqrt(price*t.k)

    return false, δ, λ
end

# Flips the sides of a bounded product CFMM
flip_sides(t::BoundedProduct{T}) where T = BoundedProduct{T}(t.k, t.β, t.α, t.R_2, t.R_1)

function find_arb!(Δ::VT, Λ::VT, cfmm::UniV3, v::VT) where {T, VT<:AbstractVector{T}}
    p = v[1]/v[2]
    γ = cfmm.γ
    if γ*cfmm.current_price <= p <= cfmm.current_price/γ
        fill!(Δ, 0)
        fill!(Λ, 0)
        return nothing
    end

    if p < γ*cfmm.current_price
        p /= γ
        δ, λ = 0.0, 0.0
        for idx in cfmm.current_tick:length(cfmm.lower_ticks)
            t = compute_at_tick(cfmm, idx)
            issat, δ_, λ_ = find_arb_pos(t, p)
            δ += δ_
            λ += λ_

            if !issat
                break
            end
        end
        Δ .= δ/γ, 0
        Λ .= 0, λ
    else
        p = 1/(γ*p)
        δ, λ = 0.0, 0.0
        for idx in cfmm.current_tick:-1:1
            t = flip_sides(compute_at_tick(cfmm, idx))
            issat, δ_, λ_ = find_arb_pos(t, p)
            δ += δ_
            λ += λ_

            if !issat
                break
            end
        end
        Δ .= 0, δ/γ
        Λ .= λ, 0
    end

    return nothing
end

# --- Testing helper functions below

# Compute max amount that can be traded at current tick
function max_amount_pos(t::BoundedProduct{T}) where T 
    if t.β > 0
        return (t.R_1 + t.α)*t.R_2/t.β
    elseif t.α > 0
        return typemax(T)
    end
    return 0.0
end

function forward_trade(Δ, cfmm::UniV3{T}) where {T}
    # Construct reserves at current tick
    γ = cfmm.γ

    if iszero(Δ)
        return 0.0
    end

    if Δ[1] > 0
        δ = γ*Δ[1]
        λ = 0.0
        for idx in cfmm.current_tick:length(cfmm.lower_ticks)
            t = compute_at_tick(cfmm, idx)
            max_amount = max_amount_pos(t)

            if max_amount > δ
                λ += δ*(t.R_2 + t.β)/(t.R_1 + t.α + δ)
                return λ
            end
            # If not, add all reserves
            λ += t.R_2

            δ -= max_amount
        end

        # We've exhausted all liquidity
        return λ
    else
        δ = γ*Δ[2]
        λ = 0.0
        for idx in cfmm.current_tick:-1:1
            t = compute_at_tick(cfmm, idx)
            max_amount = max_amount_pos(flip_sides(t))

            if max_amount > δ
                λ += δ*(t.R_2 + t.β)/(t.R_1 + t.α + δ)
                return λ
            end
            # If not, add all reserves
            λ += t.R_2

            δ -= max_amount
        end

        # We've exhausted all liquidity
        return λ
    end
end