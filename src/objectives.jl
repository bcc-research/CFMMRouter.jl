export Objective, LinearNonnegative

abstract type Objective end

f(::Objective, v) = error("unimplemented")
grad!(g, ::Objective, v) = error("unimplemented")

# ----- Objective definitions below

struct LinearNonnegative{T} <: Objective where {T}
    c::Vector{T}
end

function LinearNonnegative(c)
    all(c .> 0) || throw(ArgumentError("all elements must be strictly positive"))
    T = eltype(c)

    return LinearNonnegative{T}(
        c,
    )
end

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