abstract type CFMM end

struct GeometricMean{T <: Real, FeeType <: Real} <: CFMM
    weights::Vector{T}
    reserves::Vector{T}
    fee::FeeType
    delta::Vector{T}
    lambda::Vector{T}
    nu::Vector{T}
    amount::T
end

struct Product{T <: Real, FeeType <: Real} <: CFMM
    reserves::Vector{T}
    fee::FeeType
    delta::Vector{T}
    lambda::Vector{T}
    nu::Vector{T}
    amount::T
end

function Product(reserves::Vector{T}, fee::U) where {T, U}
    @assert length(reserves) == 2

    return GeometricMean{T, U}(
        ones(T, 2)/2,
        reserves,
        fee,
        zero(reserves),
        zero(reserves),
        zero(reserves),
        zero(T),
    )
end

# I should generalize to arbitrary two-coins, but this will do
function fenchel_conj!(p::Product{T, U}) where {T, U}

end