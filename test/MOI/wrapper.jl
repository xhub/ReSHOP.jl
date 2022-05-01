################################################################################################
#
# This comes from Knitro


# MOIT.nlptest does not support :ExprGraph
#@testset "MOI NLP tests" begin
#    @testset "with $nlp_solver" for nlp_solver in nlp_solvers
#        MOIT.nlptest(ReSHOP.Optimizer(solver=nlp_solver; optimizer_kw...), get(config_solver, nlp_solver, config))
#    end
#end

#@testset "MOI Basic constraint test" begin
#    MOIT.basic_constraint_tests(ReSHOP.Optimizer(), config, get_constraint_function=false, get_constraint_set=false)
#end

# TODO: uncomment
#@testset "MOI unit test" begin
#
#    model = ReSHOP.Optimizer(solver="cplex")
#    exclude = ["number_threads",
#			  "solve_with_lowerbound", # c1 is not registered
#				"solve_singlevariable_obj", # problem with dual constraint value
#			  ]
#
#    MOI.Test.unittest(model, config, exclude)
#end

@testset "MOI Linear tests" begin
    @testset "using $lp_solver" for lp_solver in lp_solvers
    optimizer = ReSHOP.Optimizer(solver=lp_solver; optimizer_kw...)
    exclude = [
               "linear1",
               "linear2", # DualObjectivevalue not supported
               "linear4", # ConstraintSet not supported
               "linear5", # myo_exportmodel_gams :: rosetta arrays are already present. This is not possible!
               "linear6", # ConstraintSet not supported
               "linear7", # myo_exportmodel_gams :: rosetta arrays are already present. This is not possible!
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
    MOIT.linear2test(optimizer, config_noduals)
    # TODO URG this test now fails a lot
#    MOIT.linear15test(optimizer, config_noduals)
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
    MOIT.ncqcp1test(moi_solver, get(config_solver, global_solver, solver_config))
    MOIT.ncqcp2test(moi_solver, get(config_solver, global_solver, solver_config))
    end
end

#TODO: fix problem type issue
#@testset "MOI MPCC" begin
#   @testset "with knitro" begin
#   solver_config = config
#   moi_solver = ReSHOP.Optimizer(solver="knitro"; optimizer_kw...)
#   MOIT.test_qp_complementarity_constraint(moi_solver, get(config_solver, moi_solver, solver_config))
#   end
#end

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
