using Pkg
cd(@__DIR__)
Pkg.activate("..")
using CFMMRouter
using Random
using StatsBase

cosangle(u, v) = dot(u, v)/(norm(u)*norm(v))

TOL = 1e-4

function check_primal_feasibility(r::Router)
    all_flows = zero(r.v)

    for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
        @assert all(Δ .>= -TOL)
        @assert all(Λ .>= -TOL)
        @assert CFMMRouter.ϕ(c, R=c.R + c.γ*Δ - Λ) >= CFMMRouter.ϕ(c) - sqrt(eps())

        all_flows[c.Ai] .+= Λ - Δ
    end

    @show all_flows
    @assert all(all_flows .>= -TOL)
end

function check_dual_feasibility(r::Router)
    @assert all(r.v .>= CFMMRouter.lower_limit(r.objective) .- TOL)
    @assert all(r.v .<= CFMMRouter.upper_limit(r.objective) .+ TOL)
end

function check_opt_conditions_no_fee!(r::Router)
    update_reserves!(r)

    for c in r.cfmms
        p = zero(c.R)
        CFMMRouter.∇ϕ!(p, c)
        @assert cosangle(p, r.v[c.Ai]) ≈ 1.0
    end
end

pools = Vector{CFMM{Float64}}(undef, 100)
Random.seed!(1234)
coins = collect(1:10)

for i in 1:100
    Ai = sample(coins, 2, replace=false)
    pools[i] = ProductTwoCoin(
        1000*rand(2),
        1.0,
        Ai
    )
end

router = Router(
    LinearNonnegative(rand(10)),
    pools,
    10
)

route!(router)

check_primal_feasibility(router)
check_dual_feasibility(router)
check_opt_conditions_no_fee!(router)