using CFMMRouter
using Test

using LinearAlgebra, Random
using StatsBase
using StaticArrays
using ForwardDiff
const CR = CFMMRouter

include("objectives.jl")
include("cfmms.jl")
include("arb.jl")
include("swap.jl")