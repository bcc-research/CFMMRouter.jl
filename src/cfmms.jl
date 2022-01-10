abstract type CFMM{T} end

# Credit: Chris Rackauckas
@def add_generic_fields begin
    R::Vector{T}
    γ::T
    Ai::Vector{Int}
end

Base.length(c::CFMM) = length(c.Ai)


# ------------------------------------------------------------------------------
# |                              CFMM Definitions                              |
# ------------------------------------------------------------------------------
# Each CFMM needs to implement its conjugate function
struct CSMM{T} <: CFMM{T}
    @add_generic_fields
end
function find_arb!(Δ::VT, Λ::VT, cfmm::CSMM{T}, ν::VT) where {T, VT <: Vector{T}}
    #TODO:
end

struct GMMM{T} <: CFMM{T}
    @add_generic_fields
end

struct WGMMM{T} <: CFMM{T}
    @add_generic_fields
    w::Vector{T}
end

struct Curve{T} <: CFMM{T}
    @add_generic_fields
    α::T
    β::T
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