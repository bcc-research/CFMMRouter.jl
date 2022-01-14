export CFMM, ProductTwoCoin
export find_arb!

abstract type CFMM{T} end

# Credit: Chris Rackauckas
@def add_generic_fields begin
    R::Vector{T}
    γ::T
    Ai::Vector{Int}                 # idx vector: jth coin in CFMM is Ai[j]
end

@def add_two_coin_fields begin
    R::MVector{2,T}
    γ::T
    Ai::MVector{2,UInt}
end

Base.length(c::CFMM) = length(c.Ai)

# solves min νᵀAᵢ(Λᵢ - Δᵢ) s.t. Λᵢ, Δᵢ ≥ 0 and ϕ(Rᵢ + γΔᵢ - Λᵢ) ≥ ϕ(Rᵢ)
# overwrites Δ and Λ
# returns false if no arb exists (iszero(Δ) && iszero(Λ)), true ow
# THIS DEF IS JUST FOR THE DOCTSTRING
"""
    find_arb!(Δ, Λ, cfmm, ν)

TODO
"""
function find_arb! end


# ------------------------------------------------------------------------------
# |                              CFMM Definitions                              |
# ------------------------------------------------------------------------------
# Each CFMM needs to implement its conjugate function
struct Sum{T} <: CFMM{T}
    @add_generic_fields
end
function find_arb!(Δ::VT, Λ::VT, cfmm::Sum{T}, ν::VT) where {T, VT <: Vector{T}}
    #TODO:
end

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
struct ProductTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
end

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

function ProductTwoCoin(R, γ, idx)
    γ_T, idx_uint, T = two_coin_check_cast(R, γ, idx)

    return ProductTwoCoin{T}(
        MVector{2, T}(R),
        γ_T,
        MVector{2, UInt}(idx_uint)
    )
end

# See App. A of "An Analysis of Uniswap Markets"
@inline prod_arb_δ(m, r, k, γ) = max(sqrt(γ*m*k) - r, 0)/γ
@inline prod_arb_λ(m, r, k, γ) = max(r - sqrt(k/(m*γ)), 0)

# Solves the maximum arbitrage problem for the two-coin constant product case.
# Assumes that v > 0 and γ > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::ProductTwoCoin{T}, v::VT) where {T, VT <: MVector{2, T}}
    R, γ = cfmm.R, cfmm.γ
    k = R[1]*R[2]

    Δ[1] = prod_arb_δ(v[2]/v[1], R[1], k, γ)
    Δ[2] = prod_arb_δ(v[1]/v[2], R[2], k, γ)

    Λ[1] = prod_arb_λ(v[1]/v[2], R[1], k, γ)
    Λ[2] = prod_arb_λ(v[2]/v[1], R[2], k, γ)
    return nothing
end

struct GeometricMeanTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
    w::SVector{2, T}
end

@inline geom_arb_δ(m, r1, r2, η, γ) = max((γ*m*η*r1*r2^η)^(1/(η+1)) - r2, 0)/γ
@inline geom_arb_λ(m, r1, r2, η, γ) = max(r1 - ((r2*r1^(1/η))/(η*γ*m))^(η/(1+η)), 0)

# Solves the maximum arbitrage problem for the two-coin constant product case.
# Assumes that v > 0 and w > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::GeometricMeanTwoCoin{T}, v::VT) where {T, VT <: MVector{2, T}}
    R, γ, w = cfmm.R, cfmm.γ, cfmm.w

    η = w[1]/w[2]

    Δ[1] = geom_arb_δ(v[2]/v[1], R[2], R[1], 1/η, γ)
    Δ[2] = geom_arb_δ(v[1]/v[2], R[1], R[2], η, γ)

    Λ[1] = geom_arb_λ(v[1]/v[2], R[1], R[2], η, γ)
    Λ[2] = geom_arb_λ(v[2]/v[1], R[2], R[1], 1/η, γ)
    return nothing
end