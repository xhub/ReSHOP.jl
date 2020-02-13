################################################################################################
#
# This comes from Knitro


# Default configuration.
const config = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                               optimal_status=MOI.OPTIMAL,
                               query=false,
                               infeas_certificates=false, # Do not ask for infeasibility certificates.
                               modify_lhs=false)

const config_local = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                                     optimal_status=MOI.LOCALLY_SOLVED,
                                     query=false,
                                     duals=false, # well presolve give the solution, hence no duals ...
                                     infeas_certificates=false, # Do not ask for infeasibility certificates.
                                     modify_lhs=false)

const config_noduals = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                                       optimal_status=MOI.OPTIMAL,
                                       query=false,
                                       duals=false,
                                       infeas_certificates=false,
                                       modify_lhs=false)

configs = Dict(
               "local" => config_local,
               "nodual" => config_noduals
              )

optimizer_kw = Dict(:rtol => 1e-9)

const OPTIMIZER = ReSHOP.Optimizer(; optimizer_kw...)

# MOIT.nlptest does not support :ExprGraph
#@testset "MOI NLP tests" begin
#    @testset "with $nlp_solver" for nlp_solver in nlp_solvers
#        MOIT.nlptest(ReSHOP.Optimizer(solver=nlp_solver; optimizer_kw...), get(config_solver, nlp_solver, config))
#    end
#end

@testset "MOI utils" begin
    @testset "SolverName" begin
        optimizer = ReSHOP.Optimizer()
        @test MOI.get(optimizer, MOI.SolverName()) == "ReSHOP"
    end
    @testset "supports_default_copy_to" begin
        optimizer = ReSHOP.Optimizer()
        @test MOIU.supports_default_copy_to(optimizer, false)
        # Use `@test !...` if names are not supported
        @test MOIU.supports_default_copy_to(optimizer, true)
    end
    @testset "MOI.Silent" begin
        optimizer = ReSHOP.Optimizer()
        @test MOI.supports(optimizer, MOI.Silent())
        MOI.set(optimizer, MOI.Silent(), true)
        @test MOI.get(optimizer, MOI.Silent()) == true
        MOI.set(optimizer, MOI.Silent(), false)
        @test MOI.get(optimizer, MOI.Silent()) == false
    end
#    @testset "MOI.TimeLimitSec" begin
#        optimizer = ReSHOP.Optimizer()
#        @test MOI.supports(optimizer, MOI.TimeLimitSec())
#        # TimeLimitSec is set to 1e8 by default in Knitro.
#        @test MOI.get(optimizer, MOI.TimeLimitSec()) == 1e8
#        my_time_limit = 10.
#        MOI.set(optimizer, MOI.TimeLimitSec(), my_time_limit)
#        @test MOI.get(optimizer, MOI.TimeLimitSec()) == my_time_limit
#    end
end

@testset "MOI Linear tests" begin
    @testset "using $lp_solver" for lp_solver in lp_solvers
    optimizer = ReSHOP.Optimizer(solver=lp_solver; optimizer_kw...)
    exclude = [
               "linear1",
               "linear2", # DualObjectivevalue not supported
               "linear4", # ConstraintSet not supported
               "linear6", # ConstraintSet not supported
               "linear10", # No support of SAF in  Interval for now
               "linear10b", # No support of SAF in  Interval for now
               "linear11", # problem accessing constraint function
               "linear12", # Same as above.
               "linear14", # Delete not allowed
               "linear15", # DualObjectivevalue not supported
               ]
    append!(exclude, get(lp_exclude, lp_solver, []))
    MOIT.contlineartest(optimizer, config, exclude)
    # Test linear2 and linear15 without querying the dual solution.
    MOIT.linear15test(optimizer, config_noduals)
    MOIT.linear2test(optimizer, config_noduals)
    end
end

@testset "MOI QP/QCQP tests" begin
    @testset "with $quad_solver" for (quad_solver, config_str) in quad_solvers
    # Exclude NCQP
    # Exclude qp2 and qp3, no solver supports it
    exclude = String["qp2", "qp3", "ncqcp1", "ncqcp2"]
    append!(exclude, get(quad_exclude, quad_solver, []))
    solver_config = get(configs, config_str, config)
    MOIT.contquadratictest(ReSHOP.Optimizer(solver=quad_solver; optimizer_kw...), get(config_solver, quad_solver, solver_config), exclude)
    end
end

@testset "MOI non-convex QP/QCQP tests" begin
    @testset "with $global_solver" for (global_solver, config_str) in global_solvers
    solver_config = get(configs, config_str, config)
    moi_solver = ReSHOP.Optimizer(solver=global_solver; optimizer_kw...)
    if global_solver == "baron"
        MOIT.ncqcp1test(moi_solver, config_local)
    else
        MOIT.ncqcp1test(moi_solver, get(config_solver, global_solver, solver_config))
    end
    MOIT.ncqcp2test(moi_solver, get(config_solver, global_solver, solver_config))
    end
end

#= @testset "MOI SOCP tests" begin =#
#=     # TODO: DualObjectivevalue not supported =#
#=     # Presolve must be switch off to get proper dual variables. =#
#=     config2 = MOIT.TestConfig(atol=1e-4, rtol=1e-4, infeas_certificates=false, =#
#=                               optimal_status=MOI.LOCALLY_SOLVED, query=false) =#
#=     # Behavior in infeasible case doesn't match test. =#
#=     exclude = ["lin4"] =#
#=     BRIDGED2 = MOIB.full_bridge_optimizer(ReSHOP.Optimizer(outlev=0, presolve=0), Float64) =#
#=     MOIT.lintest(BRIDGED2, config2) =#
#= end =#

# We don't support ConstraintFunction for now
#@testset "MOI MILP test with $mip_solver" for mip_solver in mip_solvers
#    exclude = [
#               "int1",        # ObjectiveBound()
#               "int3",        # Need SAF in Interval
#               "indicator1",  # Sadly not support yet
#               "indicator2",  # Sadly not support yet
#               "indicator3",  # Need SAF in Interval
#              ]
#    MOIT.intlineartest(ReSHOP.Optimizer(solver=mip_solver; optimizer_kw...), config, exclude)
#end
