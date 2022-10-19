export Objective, LinearNonnegative, BasketLiquidation, Swap

abstract type Objective end


@doc raw"""
    f(obj::Objective, v)

Evaluates the conjugate of the utility function of `objective` at `v`.
Specifically,
```math
    f(\nu) = \sup_\Psi \left(U(\Psi) - \nu^T \Psi \right).
```
"""
function f end

@doc raw"""
    grad!(g, obj::Objective, v)

Computes the gradient of [`f(obj, v)`](@ref) at v.
"""
function grad! end

@doc raw"""
    lower_limit(obj)

Componentwise lower bound on argument `v` for objective [`f`](@ref).  
Returns a vector with length `length(v)` (number of tokens).
"""
function lower_limit end

@doc raw"""
    upper_limit(obj)

Componentwise upper bound on argument `v` for objective [`f`](@ref).  
Returns a vector with length `length(v)` (number of tokens).
"""
function upper_limit end

# ----- Objective definitions below

@doc raw"""
    LinearNonnegative(c)

Linear objective for the routing problem,
```math
    U(\Psi) = c^T\Psi - \mathbf{I}(\Psi \geq 0),
```
where `c` is a positive price vector.
"""
struct LinearNonnegative{T} <: Objective
    c::AbstractVector{T}
    function LinearNonnegative(c::Vector{T}) where {T<:AbstractFloat}
        all(c .> 0) || throw(ArgumentError("all elements must be strictly positive"))
        return new{T}(
            c,
        )
    end
end
LinearNonnegative(c::Vector{T}) where {T<:Real} = LinearNonnegative(Float64.(c))

function f(obj::LinearNonnegative{T}, v) where {T}
    if all(obj.c .<= v)
        return zero(T)
    end
    return convert(T, Inf)
end

function grad!(g, obj::LinearNonnegative{T}, v) where {T}
    if all(obj.c .<= v)
        g .= zero(T)
    else
        g .= convert(T, Inf)
    end
    return nothing
end

@inline lower_limit(o::LinearNonnegative{T}) where {T} = o.c .+ 1e-8
@inline upper_limit(o::LinearNonnegative{T}) where {T} = convert(T, Inf) .+ zero(o.c)



@doc raw"""
    BasketLiquidation(i, Δin)

Liquidation objective for the routing problem,
```math
    \Psi_i - \mathbf{I}(\Psi_{-i} + Δ^\mathrm{in}_{-i} = 0, ~ \Psi_i \geq 0),
```
where `i` is the desired output token and `Δin` is the basket of tokens to be liquidated.
"""
struct BasketLiquidation{T} <: Objective
    i::Int
    Δin::Vector{T}
    
    function BasketLiquidation(i::Integer, Δin::Vector{T}) where {T<:AbstractFloat}
        !(i > 0 && i <= length(Δin)) && throw(ArgumentError("Invalid index i"))
        return new{T}(
            i,
            Δin,
        )
    end
end
BasketLiquidation(i::Integer, Δin::Vector{T}) where {T<:Real} = BasketLiquidation(i, Float64.(Δin))

function f(obj::BasketLiquidation{T}, v) where {T}
    if v[obj.i] >= 1.0
        return sum(i->(i == obj.i ? 0.0 : obj.Δin[i]*v[i]), 1:length(v))
    end
    return convert(T, Inf)
end

function grad!(g, obj::BasketLiquidation{T}, v) where {T}
    if v[obj.i] >= 1.0
        g .= obj.Δin
        g[obj.i] = zero(T)
    else
        g .= convert(T, Inf)
    end
    return nothing
end

@inline function lower_limit(o::BasketLiquidation{T}) where {T}
    ret = Vector{T}(undef, length(o.Δin))
    fill!(ret, eps())
    ret[o.i] = one(T) + eps()
    return ret
end
@inline upper_limit(o::BasketLiquidation{T}) where {T} = convert(T, Inf) .+ zero(o.Δin)


@doc raw"""
    Swap(i, j, δ, n)

Swap objective for the routing problem with `n` tokens:
```math
    \Psi_i - \mathbf{I}(\Psi_{[n]\setminus\{i,j\}} = 0,\; {\Psi_j = -\delta})
```
where `i` is the desired output token, `j` is the input token, and `δ` the amount input.
Note that this is shorthand for a BasketLiquidation objective where `Δin` is a one-hot vector. 
"""
function Swap(i::Int, j::Int, δ::T, n::Int) where {T<:AbstractFloat}
    Δin = zeros(T, n)
    Δin[j] = δ
    return BasketLiquidation(i, Δin)
end