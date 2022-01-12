export CFMM, ProductTwoCoin
export find_arb!

abstract type CFMM{T} end

# Credit: Chris Rackauckas
@def add_generic_fields begin
    R::Vector{T}
    γ::T
    Ai::Vector{Int}
end

@def add_two_coin_fields begin
    R::MVector{2,T}
    γ::T
    Ai::MVector{2,UInt}
end

Base.length(c::CFMM) = length(c.Ai)


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

function ProductTwoCoin(R, γ, idx)
    @assert length(R) == 2
    @assert length(idx) == 2

    T = eltype(R)

    if T <: Integer
        T = Float64
    end

    γ_T = convert(T, γ)
    idx_uint = convert.(UInt, idx)

    return ProductTwoCoin{T}(
        MVector{2, T}(R),
        γ_T,
        MVector{2, UInt}(idx_uint)
    )
end

# Solves the minimum arbitrage problem for the two-coin constant product case.
# Assumes that v > 0.
function find_arb!(Δ::VT, Λ::VT, cfmm::ProductTwoCoin{T}, v::VT) where {T, VT <: MVector{2, T}}
    R, γ = cfmm.R, cfmm.γ
    k = R[1]*R[2]

    Δ[1] = max(sqrt(γ*(v[2]/v[1])*k) - R[1], zero(T))/γ
    Δ[2] = max(sqrt(γ*(v[1]/v[2])*k) - R[2], zero(T))/γ

    Λ[1] = R[1]*max(1 - sqrt((v[2]*R[2])/(γ*v[1])), 0)
    Λ[2] = R[2]*max(1 - sqrt((v[1]*R[1])/(γ*v[2])), 0)
end

struct GeometricMeanTwoCoin{T} <: CFMM{T}
    @add_two_coin_fields
    w::SVector{2, T}
end

# TODO: maybe get rid of this? thought it may be useful for data parallelization 
struct Trade{T}
    cfmm::CFMM{T}
    Λ::Vector{T}
    Δ::Vector{T}
end




# solves min νᵀAᵢ(Λᵢ - Δᵢ) s.t. Λᵢ, Δᵢ ≥ 0 and ϕ(Rᵢ + γΔᵢ - Λᵢ) ≥ ϕ(Rᵢ)
# overwrites Δ and Λ
# returns false if no arb exists (iszero(Δ) && iszero(Λ)), true ow
function find_arb!(trade::Trade{T}, prices::Vector{T}) where {T}
    return find_arb!(trade.Δ, trade.Λ, trade.cfmm, @view(prices[trade.cfmm.Ai]))
end


# Some code to execute trades in CFMMs, maybe useless?
function is_valid(trade::Trade{T}) where {T <: Number}
    #TODO:
    return false
    # return trade.Δ ≥ zero(T) && trade.Λ ≥ zero(T) && trade.cfmm.accept_trade(trade.Δ, trade.Λ) 
    # ϕ(R + γΔ - Λ) == ϕ(R)
end

function execute!(trade::Trade)
    !is_valid(trade) && throw(ArgumentError("Invalid trade"))
    @. trade.cfmm.R += γ*trade.Δ - trade.Λ
end