#=
# Liquidating a basket of tokens
This example illustrates how to use CFMMRouter.jl to liquidate a basket of tokens.
=#
using CFMMRouter
using LinearAlgebra

## Create CFMMs
cfmms = [
    ProductTwoCoin([1e3, 1e4], 0.997, [1, 2]),
    ProductTwoCoin([1e3, 1e2], 0.997, [2, 3]),
    ProductTwoCoin([1e3, 2e4], 0.997, [1, 3])
]

## We want to liquidate a basket of tokens 2 & 3 into token 1
Δin = [0, 1e1, 1e2]

## Build a routing problem with liquidation objective
router = Router(
    BasketLiquidation(1, Δin),
    cfmms,
    maximum([maximum(cfmm.Ai) for cfmm in cfmms]),
)

## Optimize!
route!(router)

## Print results
Ψ = round.(Int, netflows(router))
println("Input Basket: $(round.(Int, Δin))")
println("Net trade: $Ψ")
println("Amount recieved: $(Ψ[1])")

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
            print("\t  $(tokens[ind]): $(round(Int, δ)), ")
        end
    end
    println("\n\tRecieved basket:")
    for (ind, λ) in enumerate(Λ)
        if λ > eps()
            print("\t  $(tokens[ind]): $(round(Int, λ)), ")
        end
    end
    print("\n")
end


#=
## Special Case: Trade Token 1 -> Token 2
=#
Δin = [1e1, 0.0, 0.0]

## Build a routing problem with liquidation objective (output to token 2)
router = Router(
    BasketLiquidation(2, Δin),
    cfmms,
    maximum([maximum(cfmm.Ai) for cfmm in cfmms]),
)

## Optimize!
route!(router)

## Print results
Ψ = round.(Int, netflows(router))
println("Input Basket: $(round.(Int, Δin))")
println("Net trade: $Ψ")
println("Amount recieved: $(Ψ[1])")

#=
List of individual trades with each CFMM:
=#
## Print individual trades
for (i, (Δ, Λ)) in enumerate(zip(router.Δs, router.Λs))
    tokens = router.cfmms[i].Ai
    println("CFMM $i:")
    println("\tTendered basket:")
    for (ind, δ) in enumerate(Δ)
        if δ > eps()
            print("\t  $(tokens[ind]): $(round(Int, δ)), ")
        end
    end
    println("\n\tRecieved basket:")
    for (ind, λ) in enumerate(Λ)
        if λ > eps()
            print("\t  $(tokens[ind]): $(round(Int, λ)), ")
        end
    end
    print("\n")
end