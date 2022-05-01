export CFMM, ProductTwoCoin, GeometricMeanTwoCoin, StableswapTwoCoin
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
\text{subject to} & \varphi(R + \γ\Delta - \Lambda) = \varphi(R) \\
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
# struct Product{T} <: CFMM{T}
#     @add_generic_fields
# end

# struct GeometricMean{T} <: CFMM{T}
#     @add_generic_fields
#     w::Vector{T}
# end

# struct Curve{T} <: CFMM{T}
#     @add_generic_fields
#     α::T
#     β::T
# end

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
@inline prod_arb_δ(m, r, k, γ) = max(sqrt(γ * m * k) - r, 0) / γ
@inline prod_arb_λ(m, r, k, γ) = max(r - sqrt(k / (m * γ)), 0)

# Solves the maximum arbitrage problem for the two-coin constant product case.
# Assumes that v > 0 and γ > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::ProductTwoCoin{T}, v::VT) where {T,VT<:AbstractVector{T}}
    R, γ = cfmm.R, cfmm.γ
    k = R[1] * R[2]

    Δ[1] = prod_arb_δ(v[2] / v[1], R[1], k, γ)
    Δ[2] = prod_arb_δ(v[1] / v[2], R[2], k, γ)

    Λ[1] = prod_arb_λ(v[1] / v[2], R[1], k, γ)
    Λ[2] = prod_arb_λ(v[2] / v[1], R[2], k, γ)
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
    R⁺[1] = w[1] * (R[2] / R[1])^w[2]
    R⁺[2] = w[2] * (R[1] / R[2])^w[1]
    return nothing
end

@inline geom_arb_δ(m, r1, r2, η, γ) = max((γ * m * η * r1 * r2^η)^(1 / (η + 1)) - r2, 0) / γ
@inline geom_arb_λ(m, r1, r2, η, γ) = max(r1 - ((r2 * r1^(1 / η)) / (η * γ * m))^(η / (1 + η)), 0)

# Solves the maximum arbitrage problem for the two-coin geometric mean case.
# Assumes that v > 0 and w > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::GeometricMeanTwoCoin{T}, v::VT) where {T,VT<:AbstractVector{T}}
    R, γ, w = cfmm.R, cfmm.γ, cfmm.w

    η = w[1] / w[2]

    Δ[1] = geom_arb_δ(v[2] / v[1], R[2], R[1], η, γ)
    Δ[2] = geom_arb_δ(v[1] / v[2], R[1], R[2], 1 / η, γ)

    Λ[1] = geom_arb_λ(v[1] / v[2], R[1], R[2], 1 / η, γ)
    Λ[2] = geom_arb_λ(v[2] / v[1], R[2], R[1], η, γ)
    return nothing
end


struct Stableswap{T} <: CFMM{T}
    @add_generic_fields
end

@doc raw"""
    StableswapTwoCoin(R, γ, idx)

Creates a two coin Stableswap CFMM with coins `idx[1]` and `idx[2]`, reserves `R`,
and fee `γ`. For 2 coins, it looks like the following:
```math
\varphi(x, y) &= \sqrt[3]{d}+\frac{ab}{3\sqrt[3]{d}}\quad\text{with} \\
a &= 4xy \\
b &= 1 - 4A \\
c &= 4A(x+y) \\
d &= \sqrt{\left(\frac{ac}{2}\right)^ 2 - \left(\frac{ab}{3}\right)^3} + \frac{ac}{2}
```
where `A` is the amplification coefficient.
Refer to https://hackmd.io/@prism0x/stableswap-optimal-routing for more details.
"""
struct StableswapTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
    A::T
    function StableswapTwoCoin(R, γ, idx, A)
        γ_T, idx_uint, T = two_coin_check_cast(R, γ, idx)
        return new{T}(
            MVector{2,T}(R),
            γ_T,
            MVector{2,UInt}(idx_uint),
            A
        )
    end
end

function ϕ(cfmm::StableswapTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    x = R[1]
    y = R[2]
    A = cfmm.A

    a = 4 * x * y
    b = 1 - 4 * A
    c = 4 * A * (x + y)
    d = sqrt((a * c / 2)^2 - (a * b / 3)^3) + a * c / 2
    return d^(1 / 3) + a * b / (3 * d^(1 / 3))
end

function ∇ϕ!(R⁺, cfmm::StableswapTwoCoin; R=nothing)
    R = isnothing(R) ? cfmm.R : R
    x = R[1]
    y = R[2]
    A = cfmm.A
    # I used SymPy to get the derivatives, the expressions can probably
    # be further simplified:
    R⁺[1] = 4 * x * y * (1 - 4 * A) * (-8 * A * x * y / 3 - 8 * A * y * (x + y) / 3 - (32 * A^2 * x^2 * y^2 * (2 * x + 2 * y) + 64 * A^2 * x * y^2 * (x + y)^2 - 32 * x^2 * y^3 * (1 - 4 * A)^3 / 9) / (3 * sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))) / (3 * (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(4 / 3)) + 4 * y * (1 - 4 * A) / (3 * (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(1 / 3)) + (8 * A * x * y / 3 + 8 * A * y * (x + y) / 3 + (32 * A^2 * x^2 * y^2 * (2 * x + 2 * y) + 64 * A^2 * x * y^2 * (x + y)^2 - 32 * x^2 * y^3 * (1 - 4 * A)^3 / 9) / (3 * sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))) / (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(2 / 3)
    R⁺[2] = 4 * x * y * (1 - 4 * A) * (-8 * A * x * y / 3 - 8 * A * x * (x + y) / 3 - (32 * A^2 * x^2 * y^2 * (2 * x + 2 * y) + 64 * A^2 * x^2 * y * (x + y)^2 - 32 * x^3 * y^2 * (1 - 4 * A)^3 / 9) / (3 * sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))) / (3 * (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(4 / 3)) + 4 * x * (1 - 4 * A) / (3 * (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(1 / 3)) + (8 * A * x * y / 3 + 8 * A * x * (x + y) / 3 + (32 * A^2 * x^2 * y^2 * (2 * x + 2 * y) + 64 * A^2 * x^2 * y * (x + y)^2 - 32 * x^3 * y^2 * (1 - 4 * A)^3 / 9) / (3 * sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))) / (8 * A * x * y * (x + y) + sqrt(64 * A^2 * x^2 * y^2 * (x + y)^2 - 64 * x^3 * y^3 * (1 - 4 * A)^3 / 27))^(2 / 3)
    return nothing
end

function stableswap_arb_λ(m_p, r, k, γ, A, tolerance=1e-9, epsilon=1e-5, max_iter=256)
    @inline x(Δ) = r - Δ / γ
    @inline α(Δ) = x(Δ) - k + k / (4 * A)
    @inline β(Δ) = (k^3) / (4 * A * x(Δ))

    @inline Π′(Δ) = (
        m_p
        -
        (-α(Δ) / γ + β(Δ) / (2 * γ * x(Δ)))
        /
        (2 * sqrt(α(Δ)^2 + β(Δ)))
        -
        1 / (2 * γ)
    )
    @inline Π′′(Δ) = (
        (2 * α(Δ) * x(Δ) - β(Δ))^2
        -
        4 * (α(Δ)^2 + β(Δ)) * (β(Δ) + x(Δ)^2)
    ) / (8 * γ^2 * x(Δ)^2 * (α(Δ)^2 + β(Δ))^(3 / 2))

    Δ_α = 0
    for i = 1:max_iter
        if α(Δ_α)^2 + β(Δ_α) < 0
            return 0
        end

        val = Π′(Δ_α)
        Δ_α = min(Δ_α - val / Π′′(Δ_α), r * γ * (1 - epsilon))
        println("Lambda ", i, " abs_err:", val, " sln:", Δ_α)
        if abs(val) < tolerance
            println()
            return max(Δ_α, 0)
        end
        if isnan(Δ_α)
            return 0
        end
    end
    # If the solver doesn't converge by this point, return 0
    return 0
end

function stableswap_arb_δ(m_p, r, k, γ, A, tolerance=1e-9, epsilon=1e-5, max_iter=256)
    @inline y(Δ) = r + Δ
    @inline α(Δ) = y(Δ) - k + k / (4 * A)
    @inline β(Δ) = (k^3) / (4 * A * y(Δ))

    @inline Π′(Δ) = (
        γ
        * m_p
        * (1 / 2 - (α(Δ) - k^3 / (8 * A * y(Δ)^2)) / (2 * sqrt((α(Δ))^2 + β(Δ))))
        -
        1
    )
    @inline Π′′(Δ) = (
        γ
        * m_p
        * (
            +((α(Δ) - β(Δ) / (2 * y(Δ)))^2 - (1 + β(Δ) / y(Δ)^2) * ((α(Δ))^2 + β(Δ)))
            /
            (2 * ((α(Δ))^2 + β(Δ))^(3 / 2))
        )
    )

    Δ_β = 0
    for i = 1:max_iter
        if α(Δ_β)^2 + β(Δ_β) < 0
            return 0
        end

        val = Π′(Δ_β)
        Δ_β = max(Δ_β - val / Π′′(Δ_β), -r * (1 - epsilon))
        println("Delta ", i, " abs_err:", val, " sln:", Δ_β)
        if abs(val) < tolerance
            return max(Δ_β, 0)
        end
        if isnan(Δ_β)
            return 0
        end
    end
    # If the solver doesn't converge by this point, return 0
    return 0
end

# Solves the maximum arbitrage problem for the two-coin Stableswap case.
# Assumes that v > 0 and γ > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::StableswapTwoCoin{T}, v::VT) where {T,VT<:AbstractVector{T}}
    R, γ, A = cfmm.R, cfmm.γ, cfmm.A
    k = ϕ(cfmm, R=cfmm.R)

    Δ[1] = stableswap_arb_δ(v[2] / v[1], R[1], k, γ, A)
    Δ[2] = stableswap_arb_δ(v[1] / v[2], R[2], k, γ, A)

    Λ[1] = stableswap_arb_λ(v[1] / v[2], R[1], k, γ, A)
    Λ[2] = stableswap_arb_λ(v[2] / v[1], R[2], k, γ, A)
    return nothing
end