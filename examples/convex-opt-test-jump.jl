using JuMP, LinearAlgebra, Ipopt
# using ForwardDiff

R = [100.0, 80.0, 50.0]
π = [0.9, 1, 1.3]
A = 10.0
γ = 0.996

# ϕ(R, α) = sum(R) + α / prod(R)
ϕ(α, R...) = sum(R) + α / prod(R)

# function ∇ϕ(R, α)
#     return ForwardDiff.gradient(x -> ϕ(x, α), R)
# end

# Solve D for given reserves. Taken from the Python mockup:
# https://github.com/curvefi/curve-contract/blob/b0bbf77f8f93c9c5f4e415bce9cd71f0cdee960e/tests/simulation.py#L30-L52
function solve_D(R, A, tol=1e-12)
    Dprev = 0
    n = length(R)
    S = sum(R)
    D = S
    Ann = A * n^n
    while abs(D - Dprev) > tol
        D_P = D
        for x in R
            D_P = D_P * D / (n * x)
        end
        Dprev = D
        D = (Ann * S + D_P * n) * D / ((Ann - 1) * D + (n + 1) * D_P)
    end
    return D
end

# Solves the maximum arbitrage problem for the generic Stableswap case.
function find_arb(R, π, A, γ)
    n = length(R)
    D = solve_D(R, A)
    α = D^(n + 1) / (A * n^(2 * n))

    model = Model(Ipopt.Optimizer)
    register(model, :ϕ, n + 1, ϕ; autodiff=true)

    @variable(model, Δ[1:n] >= 0)
    @variable(model, Λ[1:n] >= 0)

    ex1 = @NLexpression(model, [i = 1:n], R[i] + Δ[i] - Λ[i] / γ)
    # @NLconstraint(model, ϕ(R + Δ - Λ / γ, α) >= ϕ(R, α))
    @NLconstraint(model, ϕ(α, ex1...) >= ϕ(α, R...))
    @NLconstraint(model, [i=1:n], ex1[i] >= 0)
    @NLconstraint(model, [i=1:n], Δ[i] * Λ[i] == 0)

    @objective(model, Max, sum(π[i] * (Λ[i] - Δ[i]) for i = 1:n))
    optimize!(model)

    print(model)
    println("Objective value: ", objective_value(model))
    println("Δ = ", value.(Δ))
    println("Λ = ", value.(Λ))

    R_new = R + value.(Δ) - value.(Λ) / γ
    println("ϕ(R)             = ", ϕ(α, R...))
    println("ϕ(R + Δ - Λ / γ) = ", ϕ(α, R_new...))
end

find_arb(R, π, A, γ)