using Pkg
cd(@__DIR__)
Pkg.activate("..")
using CFMMRouter

cfmms = Vector{CFMM}([
    ProductTwoCoin()
])