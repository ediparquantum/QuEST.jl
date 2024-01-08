using Pkg
Pkg.activate(".")

include("src/QuEST.jl")
using .QuEST

env = createQuESTEnv()