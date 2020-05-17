using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

using Test


include("MOI/solvers.jl")

include("MOI/minlptests.jl")
GC.gc()
include("MOI/wrapper.jl")
