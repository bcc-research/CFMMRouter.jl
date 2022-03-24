# CFMMRouter

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tjdiamandis.github.io/CFMMRouter.jl/dev)
[![Build Status](https://github.com/bcc-research/CFMMRouter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bcc-research/CFMMRouter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/bcc-research/CFMMRouter.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bcc-research/CFMMRouter.jl)

## Overview
This package contains a fast solver for the CFMM Routing problem, as defined
by Angeris et al. in [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf). 
We partially decompose the problem to enable fast solutions when the number
of CFMMs is large relative to the number of tokens.

For more information, check out the [documentation](https://tjdiamandis.github.io/CFMMRouter.jl/dev).



## References
G Angeris, T Chitra, A Evans, and S Boyd. [Optimal Routing for Constant Function Market Makers](https://web.stanford.edu/~guillean/papers/cfmm-routing.pdf)