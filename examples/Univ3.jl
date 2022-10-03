using CFMMRouter
using LinearAlgebra, SparseArrays, StaticArrays


#The price in each tick refers to the lower bound on the interval
#so the intervals are [t_i,t_i+1] with liquidity L_i
#         p2---------
# p1------     |     p3-------
#   |          |         |
#   L1         L2        L3
#   |          |         | 
tick_list = [Dict("price" => .1, "liquidity" => 1.0),
             Dict("price" => .5, "liquidity" => 1.0),
             Dict("price" => 1.0, "liquidity" => 100.0),
             Dict("price" => 2.0, "liquidity" => 1.0),
             Dict("price" => 4.0, "liquidity" => 1.0)]

pool = UniV3([100,100],1,[1,2],1.0,3,tick_list)
Δ = [0.0,0.0]
Λ = [0.0,0.0]

find_arb!(Δ,Λ,pool, [1.0,2.0])

print(Δ,Λ)

update_reserves!(pool, Δ, Λ, [1.0,2.0])

print(pool.current_price, pool.current_tick_index)

