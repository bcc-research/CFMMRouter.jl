#=
# Arbitrage
This example illustrates how to use CFMMRouter.jl to solve the multi-market
arbitrage problem
=#
using CFMMRouter
using LinearAlgebra

## Create two pools of the same tokens, no fees (γ=1)
equal_pool = ProductTwoCoin([1e6, 1e6], 1, [1, 2])
unequal_small_pool = ProductTwoCoin([1e3, 2e3], 1, [1, 2])

## Build a routing problem with price vector = [1.0, 1.0]
prices = ones(2)
router = Router(
    LinearNonnegative(prices),
    [equal_pool, unequal_small_pool],
    2,
)

## Optimize!
route!(router)

## Print results
Ψ = round.(Int, netflows(router))
println("Net trade: $Ψ")
println("Profit: $(dot(prices, Ψ))")

#=
We can also see the list of individual trades with each CFMM:
=#
## Print individual trades
for (i, (Δ, Λ)) in enumerate(zip(router.Δs, router.Λs))
    tokens = router.cfmms[i].Ai
    println("CFMM $i:")
    println("\tTendered basket:")
    for (ind, δ) in enumerate(Δ)
        if δ > eps()
            print("\t  ")
            print("$(tokens[ind]): $(round(Int, δ)), ")
        end
    end
    println("\n\tRecieved basket:")
    for (ind, λ) in enumerate(Λ)
        if λ > eps()
            print("\t  ")
            print("$(tokens[ind]): $(round(λ, digits=0)), ")
        end
    end
    print("\n")
end