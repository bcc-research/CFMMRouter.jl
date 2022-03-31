cosangle(u, v) = dot(u, v)/(norm(u)*norm(v))

TOL = 1e-4

function check_primal_feasibility(r::Router; arb=true)
    all_flows = zero(r.v)

    for (Δ, Λ, c) in zip(r.Δs, r.Λs, r.cfmms)
        @test all(Δ .>= -TOL)
        @test all(Λ .>= -TOL)
        @test CFMMRouter.ϕ(c, R=c.R + c.γ*Δ - Λ) >= CFMMRouter.ϕ(c) - sqrt(eps())

        all_flows[c.Ai] .+= Λ - Δ
    end

    @test all(all_flows .== netflows(r))

    if arb
        @test all(all_flows .>= -TOL)
    else
        @test sum(all_flows .>= -TOL) == 1
    end
end

function check_dual_feasibility(r::Router)
    @test all(r.v .>= CFMMRouter.lower_limit(r.objective) .- TOL)
    @test all(r.v .<= CFMMRouter.upper_limit(r.objective) .+ TOL)
end

function check_opt_conditions_no_fee!(r::Router)
    update_reserves!(r)

    for c in r.cfmms
        p = zero(c.R)
        CFMMRouter.∇ϕ!(p, c)
        @test cosangle(p, r.v[c.Ai]) ≈ 1.0
    end
end

@testset "arbitrage markets" begin
    @testset "simple case" begin
        equal_pool = ProductTwoCoin([100, 100], 1, [1, 2])
        unequal_small_pool = ProductTwoCoin([1, 2], 1, [1, 2])

        router = Router(
            LinearNonnegative(ones(2)),
            [equal_pool, unequal_small_pool],
            2,
        )
        append!(router.cfmms, [equal_pool, unequal_small_pool])

        route!(router)

        check_primal_feasibility(router)
        check_dual_feasibility(router)
        check_opt_conditions_no_fee!(router)
    end

    @testset "random markets, no fee" begin
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
    end
end