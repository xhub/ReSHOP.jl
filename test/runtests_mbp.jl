using ReSHOP, Compat, JuMP, MathProgBase
using Compat.Test, Compat.LinearAlgebra, Compat.SparseArrays

# This doesn't appears to be in Compat?
@static if VERSION >= v"0.7.0"
    function isdef(symbol)
        return isdefined(@__MODULE__, symbol)
    end
else
    function isdef(symbol)
        return isdefined(symbol)
    end
end


@testset "Basic Functionality" begin
    include("MBP/nl_write.jl")
end

#solver = JuMP.UnsetSolver()
solvers = Any[]
push!(solvers, ReSHOPSolver())

@testset "Examples" begin
    examples_path = joinpath(dirname(dirname(@__FILE__)), "examples")
    for example in ["jump_nltrig.jl",
                    "jump_nlexpr.jl",
                    "jump_minlp.jl",
                    "jump_nonlinearbinary.jl",
                    "jump_no_obj.jl",
                    "jump_const_obj.jl",
                    "jump_pruning.jl"
                    ]
        include(joinpath(examples_path, example))
    end
end
# need to support n-ary max and min operators
#for example in ["jump_maxmin.jl"]
#    include(joinpath(examples_path, example))
#end

# The following list is obtained though
#    grep "solvers = Any" JuMP/test/solvers.jl
lp_solvers = Any[]
ip_solvers = Any[]
ip_dual_solvers = Any[]
semi_solvers = Any[]
sos_solvers = Any[]
lazy_solvers, lazy_soc_solvers, lazylocal_solvers, cut_solvers, cutlocal_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[], Any[], Any[], Any[]
quad_solvers = Any[]
quad_soc_solvers = Any[]
rsoc_solvers = Any[]
nlp_solvers = Any[]
convex_nlp_solvers = Any[]
minlp_solvers = Any[]
sdp_solvers = Any[]
global_solvers = Any[]

push!(lp_solvers, ReSHOPSolver())
push!(lp_solvers, ReSHOPSolver("cplex"))
#push!(lp_solvers, ReSHOPSolver("gurobi"))
push!(lp_solvers, ReSHOPSolver("mosek"))
push!(lp_solvers, ReSHOPSolver("xpress"))

push!(quad_solvers, ReSHOPSolver())
push!(quad_solvers, ReSHOPSolver("cplex"))
#push!(quad_solvers, ReSHOPSolver("gurobi"))
push!(quad_solvers, ReSHOPSolver("mosek"))
push!(quad_solvers, ReSHOPSolver("xpress"))

quad_mip_solvers = copy(quad_solvers)

# Investigate if we can solve those
#soc_solvers = copy(quad_solvers)
soc_solvers = Any[]

push!(convex_nlp_solvers, ReSHOPSolver())
push!(convex_nlp_solvers, ReSHOPSolver("conopt"))
# buggy on my system
#push!(convex_nlp_solvers, ReSHOPSolver("ipopt"))
push!(convex_nlp_solvers, ReSHOPSolver("minos"))
push!(convex_nlp_solvers, ReSHOPSolver("snopt"))
push!(convex_nlp_solvers, ReSHOPSolver("mosek"))

# push!(nlp_solvers, ReSHOPSolver("baron"))
push!(nlp_solvers, ReSHOPSolver("knitro"))

push!(global_solvers, ReSHOPSolver("knitro"))

append!(convex_nlp_solvers, nlp_solvers)

push!(minlp_solvers, ReSHOPSolver())

if VERSION < v"0.7"
    include(Pkg.dir("JuMP","test","qcqpmodel.jl"))
else
    include(joinpath(dirname(pathof(JuMP)), "..", "test", "qcqpmodel.jl"))
end

#
ipt = false

include("jump_test_nonlinear.jl")
