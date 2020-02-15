__precompile__()
module ReSHOP

import MathOptInterface
const MOI = MathOptInterface

import MathProgBase
const MPB = MathProgBase

using Compat
using Compat.LinearAlgebra, Compat.SparseArrays

export ReSHOPSolver, getsolvername, getsolveresult, getsolveresultnum, getsolvemessage, getsolveexitcode, LinearQuadraticModel

# Load in `deps.jl`, complaining if it does not exist
const depsjl_path = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    error("ReSHOP not installed properly, run Pkg.build(\"ReSHOP\"), restart Julia and try again")
end

include(depsjl_path)
reshop_valid_index_max = Csize_t(0)

# This is will be called immediately after the module is loaded
function __init__()
   global reshop_valid_index_max = ccall((:reshop_get_index_max, libreshop), Csize_t, (),)
end

const CONFIG = Dict(
:debug => false,
:export_gms => false,
:solver_log => false
)

const solverdata_dir = joinpath(dirname(@__DIR__), ".solverdata")

include("reshop_data.jl")
include("reshop_linearity.jl")
include("reshop_expr.jl")
include("reshop_params.jl")
include("reshop_convert.jl")
include("reshop_fun.jl")
include("reshop_nlp.jl")
include("reshop_gams.jl")

mutable struct ReSHOPSolver <: MPB.AbstractMathProgSolver
    solver_name::String
    options::Dict{String,Any}
    emp::Union{Nothing,Function}
end

# TODO(xhub) write a better struct/enum here
@enum MODEL_TYPE qcp=7 nlp=2 miqcp=6 minlp=5 emp=10

"change the debug state"
setdebug(b::Bool) = global CONFIG[:debug] = b
"export the problem to a GAMS Model file (.gms)"
setexport(b::Bool) = global CONFIG[:export_gms] = b
"printout the log from the solver"
setsolverlog(b::Bool) = global CONFIG[:solver_log] = b


"""
Create a ReSHOPSolver Solver for MPB. The optional arguments are:

# Optional Arguments
- `solver_name::String=""`: solver used for this problem
- `options::Dict{String,Any}=Dict{String,Any}()`: the ReSHOP options

"""
function ReSHOPSolver(solver_name::String="",
                     options::Dict{String,Any}=Dict{String,Any}())
    ReSHOPSolver(solver_name, options, nothing)
end


include("reshop_utils.jl")

include("MBP_wrapper.jl")
include("MOI_wrapper.jl")


include("reshop_mathprgm.jl")
include("reshop_ovf.jl")
include("reshop_solve.jl")

function reshop_options_set(opt::Dict{String,Any})
    jopt = reshop_options_alloc()

    for (k,v) in opt
        reshop_option_set(jopt, k, v)
    end

    return jopt
end

end


