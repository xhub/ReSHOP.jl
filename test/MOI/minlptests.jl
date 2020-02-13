## Test Knitro with MINLPTests
#
using MINLPTests, JuMP, Test

using ReSHOP

const SOLVER = JuMP.with_optimizer(ReSHOP.Optimizer; solver="knitro")

const NLP_SOLVERS = [SOLVER]
const MINLP_SOLVERS = [SOLVER]
const POLY_SOLVERS = []
const MIPOLY_SOLVERS = []

const TERMINATION_TARGET_LOCAL = Dict(
    MINLPTests.FEASIBLE_PROBLEM => JuMP.MOI.ALMOST_LOCALLY_SOLVED,
    MINLPTests.INFEASIBLE_PROBLEM => JuMP.MOI.INFEASIBLE,
)

const TERMINATION_TARGET_NLP_CVX_LOCAL = Dict(
    MINLPTests.FEASIBLE_PROBLEM => JuMP.MOI.LOCALLY_SOLVED,
    MINLPTests.INFEASIBLE_PROBLEM => JuMP.MOI.OTHER_ERROR,
)

const PRIMAL_TARGET_LOCAL = Dict(
    MINLPTests.FEASIBLE_PROBLEM => JuMP.MOI.FEASIBLE_POINT,
    MINLPTests.INFEASIBLE_PROBLEM => JuMP.MOI.INFEASIBLE_POINT,
)


#exclude_solver = Dict
#003_014
#003_015
#206_010 # maybe division by 0 from some solvers


exclude_nlp = [
               "005_010",
               "005_011",  # Uses the function `\`
               "006_010", # handling of user-defined functions
               "008_011", # 
              ]

if Sys.isapple()
    push!(exclude_nlp, "008_010")
end

#nlp_mi_002_010 # no obj function
@testset "JuMP Model Tests" begin
    @testset "$(solver.constructor): nlp" for solver in NLP_SOLVERS
        MINLPTests.test_nlp(solver, exclude = exclude_nlp)
        # For 005_010, Knitro founds a different solution, close
        # to those of MINLPTests.
        MINLPTests.nlp_005_010(solver, 1e-5, 1e-5, 1e-5)
        MINLPTests.test_nlp_cvx(solver, exclude = [
            "001_010", # solution is locally optimal
            "002_010", # ditto
        ])
#        MINLPTests.nlp_cvx_001_010(solver, 
    end
    @testset "$(solver.constructor): nlp_mi" for solver in MINLP_SOLVERS
        MINLPTests.test_nlp_mi(solver, exclude = [
            "005_010",
            "005_011",  # Uses the function `\`
            "006_010",  # handling of user-defined functions.
            "007_010",  # Not consistent across solver version
            "007_020",  # Not consistent across solver version
        ], termination_target=TERMINATION_TARGET_LOCAL, primal_target=PRIMAL_TARGET_LOCAL)
        MINLPTests.test_nlp_mi_cvx(solver)
    end
end
