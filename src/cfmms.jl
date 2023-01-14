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

function forward_trade(Δ::VT, cfmm::UniV3{T}) where {T, VT<:AbstractVector{T}}
    # Construct reserves at current tick
    k = cfmm.liquidity[cfmm.current_tick]
    α = sqrt(k/tick_high_price(cfmm))
    β = sqrt(k*tick_low_price(cfmm))
    R_1 = sqrt(k/cfmm.current_price) - α
    R_2 = k/(R_1 + α) - β

    γ = cfmm.γ

    if Δ[1] > 0
        δ = γ*Δ[1]
        idx = cfmm.current_tick
        λ = 0.0
        while idx <= length(cfmm.lower_ticks)
            # Compute max amount that can be traded at current tick
            max_amount = (R_1 + α)*R_2/β

            if max_amount > δ
                λ += γ*δ*(R_2 + β)/(R_1 + α + δ)
                return λ
            end
            # If not, add all reserves
            λ += R_2

            # Update all new values
            α = sqrt(k/tick_high_price(cfmm, idx))
            β = sqrt(k*tick_low_price(cfmm, idx))
            R_1 = 0
            R_2 = k/α - β

            δ -= max_amount
        end

        # We've exhausted all liquidity
        return λ
    else
        @error "Not implemented yet"
    end
end
