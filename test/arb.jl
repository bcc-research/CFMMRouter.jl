@testset "arbitrage checks" begin
    @testset "product" begin
        equal_pool = ProductTwoCoin([1, 1], 1, [1, 2])

        Δ = MVector{2, Float64}(undef)
        Λ = MVector{2, Float64}(undef)

        # No arb in the fee-less case
        v = MVector(1.0, 1.0)
        find_arb!(Δ, Λ, equal_pool, v)
        @test iszero(Δ) && iszero(Λ)

        v .*= 2
        find_arb!(Δ, Λ, equal_pool, v)
        @test iszero(Δ) && iszero(Λ)

        # Easy arb in the fee-less case
        v[2] = 1
        find_arb!(Δ, Λ, equal_pool, v)
        @test Δ[1] ≈ 0 && Δ[2] ≈ sqrt(2) - 1
        @test Λ[1] ≈ 1 - sqrt(1/2) && Λ[2] ≈ 0

        # No-arb in the fee case
        equal_pool_fee = ProductTwoCoin([1, 1], .9, [1, 2])
        v = MVector(.95, 1.0)
    end
end