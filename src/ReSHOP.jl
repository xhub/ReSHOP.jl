__precompile__()
module ReSHOP

import MathOptInterface
const MOI = MathOptInterface

import MathProgBase
const MPB = MathProgBase

using Compat
using LinearAlgebra, SparseArrays

using DataDeps

export ReSHOPSolver, getsolvername, getsolveresult, getsolveresultnum, getsolvemessage, getsolveexitcode, LinearQuadraticModel,
       ovf_setreformulation, set_solverstack, get_solverstacks

if VERSION <= v"1.3.0"
   # Load in `deps.jl`, complaining if it does not exist
   const depsjl_path = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
   if !isfile(depsjl_path)
      error("ReSHOP not installed properly, run Pkg.build(\"ReSHOP\"), restart Julia and try again")
   end
   include(depsjl_path)
else
   using ReSHOP_jll
end

include("reshop_utils.jl")

global_solverstack = "RESHOP"
solverstacks = ["RESHOP"]
reshop_valid_index_max = Csize_t(0)

# This is will be called immediately after the module is loaded
function __init__()
    reshop_set_printops(stdout)
    global reshop_valid_index_max = ccall((:rhp_getidxmax, libreshop), Csize_t, (),)

    gamscntr_template_file = joinpath(solverdata_dir, "gamscntr.dat")

    if isfile(gamscntr_template_file)
        input = read(gamscntr_template_file, String)
        gamsdir = gamscntr_getgamsdir(input)
        res = ccall((:rhp_gms_setgamsdir, libreshop), Cint, (Cstring,), gamsdir)
        res != 0 && error("return code $res from ReSHOP")

        ENV["PATH"] *= ":" * gamsdir
    
        # This is needed to prevent the listing of the Process directory
        ENV["DEBUG_PGAMS"] = '0'

        global global_solverstack = "GAMS"
        push!(solverstacks, "GAMS")
    end

    # This comes from https://github.com/chkwon/PATHSolver.jl
    platform = if Sys.iswindows()
        "windows"
    elseif Sys.isapple()
        "macos"
    elseif Sys.islinux()
        "linux"
    else
        error("Unsupported platform.")
    end
    libpath = Dict(
        "windows" => (
            "path50.dll",
            "e227d19109f56628fccfdfedd7ecbdfd1667a3c975dd1d16a160a69d374d5474",
        ),
        "macos" => (
            "libpath50.dylib",
            "8787de93d21f49a46146ebe2ef5844d1c20a80f934a85f60164f9ddc670412f8",
        ),
        "linux" => (
            "libpath50.so",
            "8c36baaea0952729788ec8d964253305b04b0289a1d74ca5606862c9ddb8f2fd",
        ),
    )
    libpath_filename, libpath_sha256 = libpath[platform]
    DataDeps.register(
        DataDeps.DataDep(
            "libpath50",
            "The libpath50 binary from http://pages.cs.wisc.edu/~ferris",
            "http://pages.cs.wisc.edu/~ferris/path/julia/$(libpath_filename)",
            libpath_sha256,
        ),
    )
    lusol = Dict(
        "windows" => (
            "lusol.dll",
            "2e1f0ed17914ddcf1b833898731ff4b85afab0cf914e0707dcff9e4e995cebd8",
        ),
        "macos" => (
            "liblusol.dylib",
            "52d631fd3d753581c62d5b4b636e9cb3f8cc822738fe34c6879443d5b5092f12",
        ),
        "linux" => (
            "liblusol.so",
            "ca87167853cdac9d4697a51a588d13ed9a7c093219743efa1d250cb62ac3dcb7",
        ),
    )
    liblusol_filename, liblusol_sha256 = lusol[platform]
    DataDeps.register(
        DataDeps.DataDep(
            "liblusol",
            "The lusol binary for use with PATH from http://pages.cs.wisc.edu/~ferris",
            "http://pages.cs.wisc.edu/~ferris/path/julia/$(liblusol_filename)",
            liblusol_sha256,
        ),
    )
    if haskey(ENV, "PATH_JL_LOCATION")
        global PATH_FNAME = ENV["PATH_JL_LOCATION"]
        global LUSOL_FNAME = ""
    else
        current = get(ENV, "DATADEPS_ALWAYS_ACCEPT", "false")
        ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
        global PATH_FNAME =
            joinpath(DataDeps.datadep"libpath50", libpath_filename)
        global LUSOL_FNAME =
            joinpath(DataDeps.datadep"liblusol", liblusol_filename)
        ENV["DATADEPS_ALWAYS_ACCEPT"] = current
    end
    return
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
include("reshop_ownsolver.jl")

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


include("MBP_wrapper.jl")
include("MOI_wrapper.jl")


include("reshop_mathprgm.jl")
include("reshop_ovf.jl")
include("reshop_solve.jl")

function set_solverstack(solverstack::String)
  global global_solverstack = solverstack
  return 
end

function get_solverstacks()
  return solverstacks
end

end


