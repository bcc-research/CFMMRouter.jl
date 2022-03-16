export Objective, LinearNonnegative

abstract type Objective end

f(::Objective, v) = error("unimplemented")
grad!(g, ::Objective, v) = error("unimplemented")

# ----- Objective definitions below

struct LinearNonnegative{T} <: Objective where {T}
    c::Vector{T}

    lower::Vector{T}
end

function LinearNonnegative(c)
    all(c .> 0) || throw(ArgumentError("all elements must be strictly positive"))
    T = eltype(c)

    return LinearNonnegative{T}(
        c,
        zero(c)
    )
end

function f(obj::LinearNonnegative{T}, v) where {T}
    if all(0 .<= v) && all(v .<= obj.c)
        return zero(T)
    end
    return convert(T, Inf)
end

function grad!(g, obj::LinearNonnegative{T}, v) where {T}
    if all(0 .<= v) && all(v .<= obj.c)
        g .= zero(T)
    end
    return g .= convert(T, Inf)
end

@inline lower_limit(o::LinearNonnegative{T}) where {T} = o.c
@inline upper_limit(o::LinearNonnegative{T}) where {T} = convert(T, Inf) .+ zero(o.c)