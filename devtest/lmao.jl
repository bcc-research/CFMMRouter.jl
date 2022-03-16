using Pkg
cd(@__DIR__)
Pkg.activate("..")
using CFMMRouter

equal_pool = ProductTwoCoin([100, 100], 1, [1, 2])
unequal_small_pool = ProductTwoCoin([1, 2], 1, [1, 2])

router = Router(
    LinearNonnegative(ones(2)),
    [equal_pool, unequal_small_pool],
    2,
)

v = route!(router)
