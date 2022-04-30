export CFMM, ProductTwoCoin, GeometricMeanTwoCoin, StableswapTwoCoin
export find_arb!
export stableswap_arb_lambda, stableswap_arb_delta

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
    varphi = d^(1 / 3) + a * b / (3 * d^(1 / 3))

    return varphi
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

function stableswap_arb_lambda(m_p, r, k, gamma, A, tolerance=1e-9, epsilon=1e-5, max_iter=256)
    @inline x(d) = r - d / gamma
    @inline alpha(d) = x(d) - k + k / (4 * A)
    @inline beta(d) = (k^3) / (4 * A * x(d))

    @inline func(d) = (
        m_p
        -
        (-alpha(d) / gamma + beta(d) / (2 * gamma * x(d)))
        /
        (2 * sqrt(alpha(d)^2 + beta(d)))
        -
        1 / (2 * gamma)
    )
    @inline deriv(d) = (
        (2 * alpha(d) * x(d) - beta(d))^2
        -
        4 * (alpha(d)^2 + beta(d)) * (beta(d) + x(d)^2)
    ) / (8 * gamma^2 * x(d)^2 * (alpha(d)^2 + beta(d))^(3 / 2))

    Delta_alpha = 0
    for i = 1:max_iter
        if alpha(Delta_alpha)^2 + beta(Delta_alpha) < 0
            return 0
        end

        val = func(Delta_alpha)
        Delta_alpha = min(Delta_alpha - val / deriv(Delta_alpha), r * gamma * (1 - epsilon))
        println("Lambda ", i, " abs_err:", val, " sln:", Delta_alpha)
        if abs(val) < tolerance
            println()
            return max(Delta_alpha, 0)
        end
        if isnan(Delta_alpha)
            return 0
        end
    end
    # If the solver doesn't converge by this point, return 0
    return 0
end

function stableswap_arb_delta(m_p, r, k, gamma, A, tolerance=1e-9, epsilon=1e-5, max_iter=256)
    @inline y(d) = r + d
    @inline alpha(d) = y(d) - k + k / (4 * A)
    @inline beta(d) = (k^3) / (4 * A * y(d))

    @inline func(d) = (
        gamma
        * m_p
        * (1 / 2 - (alpha(d) - k^3 / (8 * A * y(d)^2)) / (2 * sqrt((alpha(d))^2 + beta(d))))
        -
        1
    )
    @inline deriv(d) = (
        gamma
        * m_p
        * (
            +((alpha(d) - beta(d) / (2 * y(d)))^2 - (1 + beta(d) / y(d)^2) * ((alpha(d))^2 + beta(d)))
            /
            (2 * ((alpha(d))^2 + beta(d))^(3 / 2))
        )
    )

    Delta_beta = 0
    for i = 1:max_iter
        if alpha(Delta_beta)^2 + beta(Delta_beta) < 0
            return 0
        end

        val = func(Delta_beta)
        Delta_beta = max(Delta_beta - val / deriv(Delta_beta), -r * (1 - epsilon))
        println("Delta ", i, " abs_err:", val, " sln:", Delta_beta)
        if abs(val) < tolerance
            println()
            return max(Delta_beta, 0)
        end
        if isnan(Delta_beta)
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

    Δ[1] = stableswap_arb_delta(v[2] / v[1], R[1], k, γ, A)
    Δ[2] = stableswap_arb_delta(v[1] / v[2], R[2], k, γ, A)

    Λ[1] = stableswap_arb_lambda(v[1] / v[2], R[1], k, γ, A)
    Λ[2] = stableswap_arb_lambda(v[2] / v[1], R[2], k, γ, A)
    return nothing
end