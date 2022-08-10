init_opt_cache(R) = (R⁺=similar(R), ∇ϕR = similar(R))

function optimality_conditions_met(c, Δ, Λ, cfmm; cache=nothing)
    R, γ = cfmm.R, cfmm.γ
    if isnothing(cache)
        cache = init_opt_cache(R)
    end
    R⁺ = cache.R⁺
    ∇ϕR = cache.∇ϕR
    n = length(c)

    # Primal feasibility
    R⁺ = R + γ * Δ - Λ
    pfeas = all(Δ .≥ 0) && all(Λ .≥ 0)
    ϕR = CR.ϕ(cfmm)
    ϕR⁺ = CR.ϕ(cfmm; R=R⁺)
    CR.∇ϕ!(∇ϕR, cfmm; R=R⁺)
    cfmm_sat = ≈(ϕR, ϕR⁺) && ϕR⁺ ≥ ϕR - sqrt(eps())

    opt = maximum(i -> γ * ∇ϕR[i] / c[i], 1:n) ≤ minimum(i -> ∇ϕR[i] / c[i], 1:n) + sqrt(eps())
    return pfeas && cfmm_sat && opt
end

@testset "CFMMs" begin
@testset "arbitrage checks: two coins" begin
    Δ = MVector{2, Float64}(undef)
    Λ = MVector{2, Float64}(undef)
    
    n = 3
    Random.seed!(1234)
    γs = [rand() for i in 1:n]
    Rs = [rand(2)*10 for i in 1:n]
    νs = [@MVector rand(2) for i in 1:n]
    cache = init_opt_cache(Rs[1])

    @testset "product" begin
        equal_pool = ProductTwoCoin([1, 1], 1, [1, 2])

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
        @test length(ProductTwoCoin([1, 1], .9, [1, 2])) == 2
        @test_throws ArgumentError ProductTwoCoin([1, 1], .9, [1])

        for R in Rs, γ in γs, ν in νs
            cfmm = ProductTwoCoin(R, γ, [1, 2])
            find_arb!(Δ, Λ, cfmm, ν)
            @test optimality_conditions_met(ν, Δ, Λ, cfmm; cache=cache)
        end

    end

    @testset "geo mean" begin
        ws = [(w1=rand(); [w1; 1-w1]) for _ in 1:n]
        for R in Rs, γ in γs, ν in νs, w in ws
            cfmm = GeometricMeanTwoCoin(R, w, γ, [1, 2])
            find_arb!(Δ, Λ, cfmm, ν)
            @test optimality_conditions_met(ν, Δ, Λ, cfmm; cache=cache)
        end
    end


end
end
@testset "univ3" begin


    #this test is going to be one were there is no arbitrage because the current_price = target_price
    #there can be some slight numerical error when calculating virtual_reserves

    # #reserves dont matter for univ3 pools, only tick data so just setting to 1
    R = [1.0,1.0]
    # #no fees for now
    γ = .997
    # # 
    ids = [1.0,2.0]
    current_price = 15.0
    current_tick_index = 3.0
    ticks = [Dict("price" => 1/2.0^64, "liquidity" => 0.0),Dict("price" => 5.0, "liquidity" => 1.0), Dict("price" => 10.0, "liquidity" => 2.0), Dict("price" => 20.0, "liquidity" => 1.0), Dict("price" => 30.0, "liquidity" => 0.0), Dict("price" => 2.0^64, "liquidity" => 0)]
    cfmm = UniV3(R,γ,ids, current_price, current_tick_index, ticks)
    
    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [15.0,1.0]
    find_arb!(Δ,Λ,cfmm,v)

    @test (norm(Δ) <= 1e-12) && (norm(Λ) <= 1e-12) #nothing happens


    #In this test the target price is higher than the current_price but we are still within one tick

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [18.0,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    
    @test (Δ[1] == 0) # no token 0 in
    @test (Δ[2] > 0)  # posiive token 1 in
    @test (Λ[1] > 0)  # positive token 0 out
    @test (Λ[2] == 0) # no token 1 out
    @test (16 <= Δ[2]/Λ[1]) && (Δ[2]/Λ[1] <= 17) #Average Price paid is between 16 and 17 


    # now we cross a tick going up 

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [22.0,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    
    @test (Δ[1] == 0) # no token 0 in
    @test (Δ[2] > 0)  # posiive token 1 in
    @test (Λ[1] > 0)  # positive token 0 out
    @test (Λ[2] == 0) # no token 1 out

    @test (17 <= Δ[2]/Λ[1]) && (Δ[2]/Λ[1] <= 18) #Average Price paid is between 17 and 18 


    #out of bounds test upper
    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [29.99,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    
    @test (Δ[1] == 0) # no token 0 in
    @test (Δ[2] > 0)  # posiive token 1 in
    @test (Λ[1] > 0)  # positive token 0 out
    @test (Λ[2] == 0) # no token 1 out


    @test (19.5 <= Δ[2]/Λ[1]) && (Δ[2]/Λ[1] <= 20.5)

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [35,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    
    @test (Δ[1] == 0) # no token 0 in
    @test (Δ[2] > 0)  # posiive token 1 in
    @test (Λ[1] > 0)  # positive token 0 out
    @test (Λ[2] == 0) # no token 1 out


    @test (19.5 <= Δ[2]/Λ[1]) && (Δ[2]/Λ[1] <= 20.5)


    #Now we test when target price is below current price but within current tick

    #starting with just below current price
    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [14.99999,1.0] #basically current price but less than so we check the other branch
    find_arb!(Δ,Λ,cfmm,v)



    @test (norm(Δ) <= 1e-3) && (norm(Λ) <= 1e-3) #nothing happens


    #now below current price but not below current tick price

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [12.0,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    
    @test (Δ[1] > 0)  # positive token 0 in
    @test (Δ[2] == 0) # no token 1 in
    @test (Λ[1] == 0)  # no 0 out
    @test (Λ[2] > 0) # positive 1 out

    @test (13 <= Λ[2]/Δ[1]) && (Λ[2]/Δ[1] <= 14) #Average Price paid is between 13 and 14

    #now below current price and crossing a tick

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [5.01,1.0]
    find_arb!(Δ,Λ,cfmm,v)


    @test (Δ[1] > 0)  # positive token 0 in
    @test (Δ[2] == 0) # no token 1 in
    @test (Λ[1] == 0)  # no 0 out
    @test (Λ[2] > 0) # positive 1 out

    @test (9 <= Λ[2]/Δ[1]) && (Λ[2]/Δ[1] <= 10) #Average Price paid is between 9 and 10

    #now out of bounds test

    Δ = [0.0,0.0]
    Λ = [0.0,0.0]
    v = [4.0,1.0]
    find_arb!(Δ,Λ,cfmm,v)
    

    @test (Δ[1] > 0)  # positive token 0 in
    @test (Δ[2] == 0) # no token 1 in
    @test (Λ[1] == 0)  # no 0 out
    @test (Λ[2] > 0) # positive 1 out

    @test (9 <= Λ[2]/Δ[1]) && (Λ[2]/Δ[1] <= 10) #Average Price paid is between 13 and 14

    
end