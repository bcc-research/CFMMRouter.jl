@testset "swap markets" begin
    @testset "simple case" begin
        equal_pool = ProductTwoCoin([100, 100], 1, [1, 2])
        unequal_small_pool = ProductTwoCoin([1, 2], 1, [1, 2])
        Δin = [5.0, 0.0]

        router = Router(
            BasketLiquidation(1, Δin),
            [equal_pool, unequal_small_pool],
            2,
        )

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

        Δin = vcat([0.0], 100*rand(9))
        router = Router(
            BasketLiquidation(1, Δin),
            pools,
            10
        )

        route!(router)

        check_primal_feasibility(router; arb=false)
        check_dual_feasibility(router)
        check_opt_conditions_no_fee!(router)
    end
end