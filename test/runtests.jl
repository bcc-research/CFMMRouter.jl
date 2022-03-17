using CFMMRouter
using Test

using LinearAlgebra, Random
using StatsBase
# using Convex, SCS
using StaticArrays
const CR = CFMMRouter

include("newton.jl")
include("arb.jl")