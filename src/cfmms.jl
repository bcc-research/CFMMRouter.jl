export CFMM, ProductTwoCoin, GeometricMeanTwoCoin, UniV3
export find_arb!

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
    GeometricMeanTwoCoin(R, γ, idx, w)

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


struct UniV3{T} <: CFMM{T}
    @add_two_coin_fields
    current_price :: T
    current_tick_index :: Int #This is the index of the maximal tick with price lower than current_price in the ticks dictionary
    ticks #This is the tick mapping sorted by price
    function UniV3(R,γ,idx,current_price,current_tick_index,ticks)
        γ_T, idx_uint, T = two_coin_check_cast(R, γ, idx)
        return new{T}(
            MVector{2,T}(R),
            γ_T,
            MVector{2,UInt}(idx_uint),
            current_price,
            current_tick_index,
            ticks
        )
    end
end

## See univ3 whitepaper
function virtual_reserves(P,L)
    sP = sqrt(P)
    x = L/sP
    y = L*sP
    return([x, y])
end

function find_arb!(Δ::VT, Λ::VT, cfmm::UniV3{T}, v::VT) where {T, VT<:AbstractVector{T}} 
    current_price, current_tick_index, γ, ticks = cfmm.current_price, cfmm.current_tick_index, cfmm.γ, cfmm.ticks
    Δ[1] = 0
    Δ[2] = 0

    Λ[1] = 0
    Λ[2] = 0

    target_price = v[2]/v[1]

    if target_price >= current_price #iterate forwards in the tick mapping
        i = 1
        while true
            next_tick_price = ticks[current_tick_index + i]["price"]
            if next_tick_price > target_price ## so now we know that current_price <= target_price < next_tick_price
                R = virtual_reserves(max(current_price,ticks[current_tick_index + i - 1]["price"]),ticks[current_tick_index + i - 1]["liquidity"])
                k = R[1]*R[2]

                Δ[1] += prod_arb_δ(target_price, R[1], k, γ)
                Δ[2] += prod_arb_δ(1/target_price, R[2], k, γ)

                Λ[1] += prod_arb_λ(1/target_price, R[1], k, γ)
                Λ[2] += prod_arb_λ(target_price, R[2], k, γ)

                break

            elseif next_tick_price <= target_price ## so now we know that current_price <= next_tick_price <= target_price
                R = virtual_reserves(max(current_price,ticks[current_tick_index + i - 1]["price"]),ticks[current_tick_index + i - 1]["liquidity"])
                k = R[1]*R[2]

                Δ[1] += prod_arb_δ(next_tick_price, R[1], k, γ)
                Δ[2] += prod_arb_δ(1/next_tick_price, R[2], k, γ)

                Λ[1] += prod_arb_λ(1/next_tick_price, R[1], k, γ)
                Λ[2] += prod_arb_λ(next_tick_price, R[2], k, γ)
            end
            i += 1
        end

    elseif target_price < current_price #iterate backwards in the tick mapping
        i = 1
        while true
            prev_tick_price = ticks[current_tick_index - i]["price"]
            if prev_tick_price < target_price ## so now we know that prev_tick_price < target_price <=  current_price

                R = virtual_reserves(min(current_price,ticks[current_tick_index - i + 1]["price"]),ticks[current_tick_index - i + 1]["liquidity"])
                k = R[1]*R[2]
                
                Δ[1] += prod_arb_δ(target_price, R[1], k, γ)
                Δ[2] += prod_arb_δ(1/target_price, R[2], k, γ)

                Λ[1] += prod_arb_λ(1/target_price, R[1], k, γ)
                Λ[2] += prod_arb_λ(target_price, R[2], k, γ)

                break

            elseif prev_tick_price <= target_price ## so now we know that current_price <= next_tick_price <= target_price
                R = virtual_reserves(min(current_price,ticks[current_tick_index - i + 1]["price"]),ticks[current_tick_index - i + 1]["liquidity"])
                k = R[1]*R[2]

                Δ[1] += prod_arb_δ(prev_tick_price, R[1], k, γ)
                Δ[2] += prod_arb_δ(1/prev_tick_price, R[2], k, γ)

                Λ[1] += prod_arb_λ(1/prev_tick_price, R[1], k, γ)
                Λ[2] += prod_arb_λ(prev_tick_price, R[2], k, γ)
            end
            i += 1
        end
    end
    return nothing
end

