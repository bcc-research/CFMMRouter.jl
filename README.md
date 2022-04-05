# CFMMRouter

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bcc-research.github.io/CFMMRouter.jl/dev/)
[![Build Status](https://github.com/bcc-research/CFMMRouter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bcc-research/CFMMRouter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/bcc-research/CFMMRouter.jl/branch/main/graph/badge.svg?token=TYizMgRYNE)](https://codecov.io/gh/bcc-research/CFMMRouter.jl)

## Overview
This package contains a fast solver for the CFMM Routing problem, as defined
by Angeris et al. in [Optimal Routing for Constant Function Market Makers](https://angeris.github.io/papers/cfmm-routing.pdf). 
We partially decompose the problem to enable fast solutions when the number
of CFMMs is large relative to the number of tokens.

For more information, check out the [documentation](https://bcc-research.github.io/CFMMRouter.jl/dev/).

## Quick Start
First, add the package locally.
```julia 
using Pkg; Pkg.add(url="https://github.com/bcc-research/CFMMRouter.jl")
```

Make some swap pools.
```julia
using LinearAlgebra
using CFMMRouter

equal_pool = ProductTwoCoin([1e6, 1e6], 1, [1, 2])
unequal_small_pool = ProductTwoCoin([1e3, 2e3], 1, [1, 2])
prices = ones(2)
```

Build a Router & route.
```julia
router = Router(
    LinearNonnegative(prices),
    [equal_pool, unequal_small_pool],
    2,
)
route!(router)
```

Check out the results.
```julia
Ψ = round.(Int, netflows(router))
println("Profit: $(dot(prices, Ψ))")
```

## Performance
This routing algorithm scales approximately linearly in the number of swap pools
for the [arbitrage problem](https://bcc-research.github.io/CFMMRouter.jl/dev/examples/arbitrage/).
These tests were run on a MacBook Pro with a 2.3GHz 8-Core Intel i9 processor.
Several performance improvements are possible.
![alt text](https://github.com/bcc-research/CFMMRouter.jl/blob/main/benchmark/router_scaling.png)

## References
G Angeris, T Chitra, A Evans, and S Boyd. [Optimal Routing for Constant Function Market Makers](https://angeris.github.io/papers/cfmm-routing.pdf)
