module CFMMRouter

using LinearAlgebra, SparseArrays, StaticArrays
using Optim
using Printf

include("utils.jl")
include("cfmms.jl")
include("objectives.jl")
include("router.jl")

end
