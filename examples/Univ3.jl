#=
# Uniswap V3 Router
This example illustrates how to setup uniswapv3 pools
=#

using CFMMRouter
using LinearAlgebra


## UniV3 pool parameters
current_price = 15.0
lower_ticks = [30.0, 20, 10, 5]
liquidity = [1.0, 2.0, 1.5, 0.0]
Ai = [1, 2]
γ = 0.997

## Create pool
cfmm = UniV3(current_price, lower_ticks, liquidity, γ, Ai)


#=
We find arbitrage assuming a true price of 25.0 (ν1/ν2 = 25.0) and display
the results.
=#
Δ = zeros(2)
Λ = zeros(2)
find_arb!(Δ, Λ, cfmm, [25.0, 1.0])

## Print individual trades
tokens = Ai
println("\tTendered basket:")
for (ind, δ) in enumerate(Δ)
    if δ > eps()
        print("\t  $(tokens[ind]): $(round(δ, digits=3)), ")
    end
end
println("\n\tRecieved basket:")
for (ind, λ) in enumerate(Λ)
    if λ > eps()
        print("\t  $(tokens[ind]): $(round(λ, digits=3)), ")
    end
end
print("\n")