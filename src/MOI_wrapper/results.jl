# MathOptInterface results
MOI.get(model::Optimizer, ::MOI.RawStatusString) = string(get_status(model.ctx))

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.number_solved == 0
        return MOI.OPTIMIZE_NOT_CALLED
    end
    (solver_code, model_code) = get_status(model.ctx)
    if solver_code == :Optimal
        if model_code == :OptimalGlobal
            return MOI.OPTIMAL
        elseif model_code == :OptimalLocal || model_code == :Feasible
            return MOI.LOCALLY_SOLVED
        elseif model_code == :Unbounded
            return MOI.INFEASIBLE_OR_UNBOUNDED
        elseif model_code == :InfeasibleGlobal
            return MOI.INFEASIBLE
        elseif model_code == :InfeasibleLocal
            return MOI.LOCALLY_INFEASIBLE
        # :InfeasibleIntermed should not match here
        elseif model_code == :Integer
            # This actually means that a feasible (integer) solution has been found
            return MOI.ALMOST_LOCALLY_SOLVED
        elseif model_code == :IntegerInfeasible
            return MOI.INFEASIBLE
        elseif model_code == :InfeasibleNoSolution
            return MOI.INFEASIBLE_OR_UNBOUNDED
        elseif model_code == :UnboundedNoSolution
            return MOI.INFEASIBLE_OR_UNBOUNDED
# DUAL_INFEASIBLE
        else
            MOI.OTHER_ERROR
        end
    # Solved to relaxed tolerances
# ALMOST_OPTIMAL
# ALMOST_INFEASIBLE
# ALMOST_DUAL_INFEASIBLE
# ALMOST_LOCALLY_SOLVED
    # Limits
    elseif solver_code == :IterationLimit
        return MOI.ITERATION_LIMIT
    elseif solver_code == :TimeLimit
        return MOI.TIME_LIMIT
# NODE_LIMIT
# SOLUTION_LIMIT
# MEMORY_LIMIT
# OBJECTIVE_LIMIT
# NORM_LIMIT
# OTHER_LIMIT
    # Problematic
    elseif solver_code == :SolverInterrupt
        # TODO Check that it is correct
        return MOI.SLOW_PROGRESS
    elseif solver_code == :EvalError
        return MOI.NUMERICAL_ERROR
# INVALID_MODEL
# INVALID_OPTION
    elseif solver_code == :UserInterrupt
        return MOI.INTERRUPTED
    elseif solver_code == :Capability
        gams_solver = ctx_get_solvername(m.reshop_ctx_dest)
        println("ReSHOP: GAMS solver $(gams_solver) cannot solve the specified problem")
        return MOI.OTHER_ERROR
    elseif solver_code in (:SystemError, :LicenseIssue, :SetupError)
        return MOI.OTHER_ERROR
    elseif solver_code == :SolverError && model_code == :ErrorNoSolution
        return MOI.OTHER_ERROR
    else
        println("ReSHOP: unhandle case: solver code is $(solver_code) and model code is $(model_code)")
        return MOI.OTHER_ERROR
    end
end

# TODO
function MOI.get(model::Optimizer, ::MOI.ResultCount)
    has_result = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED, MOI.ITERATION_LIMIT, MOI.TIME_LIMIT]
    if (MOI.get(model, MOI.TerminationStatus()) in has_result)
        return 1
    else
        return 0
    end
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    if model.number_solved == 0
        return MOI.NO_SOLUTION
    end

    (solver_code, model_code) = get_status(model.ctx)
    if solver_code == :Optimal
        if model_code == :OptimalGlobal
            return MOI.FEASIBLE_POINT
        elseif model_code == :OptimalLocal || model_code == :Feasible
            return MOI.FEASIBLE_POINT
        elseif model_code == :Unbounded
            return MOI.NO_SOLUTION
        elseif model_code == :InfeasibleGlobal
            return MOI.INFEASIBLE_POINT
        elseif model_code == :InfeasibleLocal
            return MOI.INFEASIBLE_POINT
        # :InfeasibleIntermed should not match here
        elseif model_code == :Integer
            # This actually means that a feasible (integer) solution has been found
            return MOI.FEASIBLE_POINT
        elseif model_code == :IntegerInfeasible
            return MOI.INFEASIBLE_POINT
        elseif model_code == :InfeasibleNoSolution
            return MOI.INFEASIBLE_POINT
        elseif model_code == :UnboundedNoSolution
            return MOI.NO_SOLUTION
# DUAL_INFEASIBLE
        else
            MOI.UNKNOWN_RESULT_STATUS
        end
    # Solved to relaxed tolerances
# ALMOST_OPTIMAL
# ALMOST_INFEASIBLE
# ALMOST_DUAL_INFEASIBLE
# ALMOST_LOCALLY_SOLVED
    # Limits
    elseif solver_code == :IterationLimit
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :TimeLimit
        return MOI.UNKNOWN_RESULT_STATUS
# NODE_LIMIT
# SOLUTION_LIMIT
# MEMORY_LIMIT
# OBJECTIVE_LIMIT
# NORM_LIMIT
# OTHER_LIMIT
    # Problematic
    elseif solver_code == :SolverInterrupt
        # TODO Check that it is correct
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :EvalError
        return MOI.UNKNOWN_RESULT_STATUS
# INVALID_MODEL
# INVALID_OPTION
    elseif solver_code == :UserInterrupt
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :Capability
        gams_solver = ctx_get_solvername(m.reshop_ctx_dest)
        println("ReSHOP: GAMS solver $(gams_solver) cannot solve the specified problem")
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code in (:SystemError, :LicenseIssue, :SetupError)
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :SolverError && model_code == :ErrorNoSolution
        return MOI.NO_SOLUTION
    else
        println("ReSHOP: unhandle case: solver code is $(solver_code) and model code is $(model_code)")
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(model::Optimizer, ::MOI.DualStatus)
    if model.number_solved == 0
        return MOI.NO_SOLUTION
    end

    (solver_code, model_code) = get_status(model.ctx)
    if solver_code == :Optimal
        if model_code == :OptimalGlobal
            return MOI.FEASIBLE_POINT
        elseif model_code == :OptimalLocal || model_code == :Feasible
            return MOI.FEASIBLE_POINT
        elseif model_code == :Unbounded
            return MOI.FEASIBLE_POINT
        elseif model_code == :InfeasibleGlobal
            return MOI.NO_SOLUTION
        elseif model_code == :InfeasibleLocal
            return MOI.NO_SOLUTION
        # :InfeasibleIntermed should not match here
        elseif model_code == :Integer
            # This actually means that a feasible (integer) solution has been found
            return MOI.FEASIBLE_POINT
        elseif model_code == :IntegerInfeasible
            return MOI.NO_SOLUTION
        elseif model_code == :InfeasibleNoSolution
            return MOI.NO_SOLUTION
        elseif model_code == :UnboundedNoSolution
            return MOI.FEASIBLE_POINT
# DUAL_INFEASIBLE
        else
            MOI.UNKNOWN_RESULT_STATUS
        end
    # Solved to relaxed tolerances
# ALMOST_OPTIMAL
# ALMOST_INFEASIBLE
# ALMOST_DUAL_INFEASIBLE
# ALMOST_LOCALLY_SOLVED
    # Limits
    elseif solver_code == :IterationLimit
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :TimeLimit
        return MOI.UNKNOWN_RESULT_STATUS
# NODE_LIMIT
# SOLUTION_LIMIT
# MEMORY_LIMIT
# OBJECTIVE_LIMIT
# NORM_LIMIT
# OTHER_LIMIT
    # Problematic
    elseif solver_code == :SolverInterrupt
        # TODO Check that it is correct
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :EvalError
        return MOI.UNKNOWN_RESULT_STATUS
# INVALID_MODEL
# INVALID_OPTION
    elseif solver_code == :UserInterrupt
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :Capability
        gams_solver = ctx_get_solvername(m.reshop_ctx_dest)
        println("ReSHOP: GAMS solver $(gams_solver) cannot solve the specified problem")
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code in (:SystemError, :LicenseIssue, :SetupError)
        return MOI.UNKNOWN_RESULT_STATUS
    elseif solver_code == :SolverError && model_code == :ErrorNoSolution
        return MOI.NO_SOLUTION
    else
        println("ReSHOP: unhandle case: solver code is $(solver_code) and model code is $(model_code)")
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(model::Optimizer, ::S) where S <: MOI.ObjectiveValue
    if model.number_solved == 0
        error("ObjectiveValue not available.")
    end
    obj_equ = ctx_getobjequ(model.ctx)
    if is_valid_index(obj_equ)
        return ctx_getequval(model.ctx, obj_equ)
    else
        objvar = ctx_getobjvar(model.ctx)
        if is_valid_index(objvar)
            return ctx_getvarval(model.ctx, objvar)
        else
            error("No valid objective function of variable found")
        end
    end
end

function MOI.get(model::Optimizer,
                 ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.number_solved == 0
        error("VariablePrimal not available.")
    end
    check_inbounds(model, vi)
    return ctx_getvarval(model.ctx, vi.value-1)
end
function MOI.get(model::Optimizer,
                 ::MOI.VariablePrimal, vi::Vector{MOI.VariableIndex})
    if model.number_solved == 0
        error("VariablePrimal not available.")
    end
    return [ctx_getvarval(model.ctx, v.value-1) for v in vi]
end

macro checkcons(model, ci)
    quote
        if $(esc(model)).number_solved == 0
            error("Solve problem before accessing solution.")
        end
        if !(1 <= $(esc(ci)).value <= number_constraints($(esc(model))))
            error("Invalid constraint index ", $(esc(ci)).value)
        end
    end
end

##################################################
## ConstraintPrimal
## This is the value of the function at the current iterate
function _get_equval(ctx, eidx)
    return ctx_getequval(ctx, eidx)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: SF, T <: SS}
    @checkcons(model, ci)
    return _get_equval(model.ctx, ci.value-1)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VAF, T <: Union{MOI.Nonnegatives, MOI.Nonpositives}}
    @checkcons(model, ci)
    return [_get_equval(model.ctx, eidx) + ctx_getcst(model.ctx, eidx) for eidx in model.vaf_mapping[ci]]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VOV, T <: Union{MOI.Nonnegatives, MOI.Nonpositives}}
    @checkcons(model, ci)
    return [ctx_getvarval(model.ctx, vidx) for vidx in model.vaf_mapping[ci]]
end

# delete
function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VAF, T <: MOI.Zeros}
    @checkcons(model, ci)
    ncons = length(model.vaf_mapping[ci])
    return zeros(ncons)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VOV, T <: MOI.Zeros}
    @checkcons(model, ci)
    ncons = length(model.vov_mapping[ci])
    return zeros(ncons)
end

#function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
#                 ci::MOI.ConstraintIndex{S, T}) where {S <: Union{VAF, VOV}, T <: MOI.SecondOrderCone}
#    @checkcons(model, ci)
#    x = get_solution(model.inner)
#    index = model.constraint_mapping[ci] .+ 1
#    return x[index]
#end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S <: SS
    if model.number_solved == 0
        error("ConstraintPrimal not available.")
    end
    return ctx_getvarval(model.ctx, ci.value-1)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal,
                 ci::Vector{MOI.ConstraintIndex{MOI.SingleVariable, S}}) where S <: SS
    if model.number_solved == 0
        error("ConstraintPrimal not available.")
    end
    return [ctx_getvarval(model.ctx, c.value-1) for c in ci]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintPrimal, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return ctx_getvarval(model.ctx, vi.value-1)
end

##################################################
## ConstraintDual
#

sense_to_sign(model::Optimizer) = model.sense == MOI.MAX_SENSE ? -1. : 1.

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: SF, T <: LS}
    @checkcons(model, ci)
    return sense_to_sign(model) * ctx_getequmult(model.ctx, ci.value-1)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VAF, T <: VLS}
    @checkcons(model, ci)
    return sense_to_sign(model) * [ctx_getequmult(model.ctx, eidx) for eidx in model.vaf_mapping[ci]]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{S, T}) where {S <: VOV, T <: VLS}
    @checkcons(model, ci)
    return sense_to_sign(model) * [ctx_getvarmult(model.ctx, vidx) for vidx in model.vov_mapping[ci]]
end

function MOI.get(model::Optimizer, ::MOI.ConstraintDual, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return sense_to_sign(model) * ctx_getvarmult(model.ctx, vi.value-1)
end

###
# Get constraint of a SOC constraint.
#
# Use the following mathematical property.  Let
#
#   ||u_i || <= t_i      with dual constraint    || z_i || <= w_i
#
# At optimality, we have
#
#   w_i * u_i  = - t_i z_i
#
###
## TODO investigate
#function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
#                 ci::MOI.ConstraintIndex{S, T}) where {S <: Union{VAF, VOV}, T <: MOI.SecondOrderCone}
#    @checkcons(model, ci)
#    index_var = model.constraint_mapping[ci] .+ 1
#    index_con = ci.value
#    x =  get_solution(model.inner)[index_var]
#    # By construction.
#    t_i = x[1]
#    u_i = x[2:end]
#    w_i = get_dual(model.inner)[index_con]
#
#    dual = [-w_i; 1/t_i * w_i * u_i]
#
#    return sense_to_sign(model) * dual
#end

## Reduced costs.
function MOI.get(model::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, <: SS})
    if model.number_solved == 0
        error("ConstraintDual not available.")
    end
    if is_constraint_active(model.ctx, ci)
        return sense_to_sign(model) * ctx_getvarmult(model.ctx, ci.value-1)
    else
        return 0.
    end
end

function MOI.get(model::Optimizer, ::MOI.NLPBlockDual)
    if model.number_solved == 0
        error("NLPBlockDual not available.")
    end
    return [sense_to_sign(model) * ctx_getequmult(model.ctx, eidx) for eidx in model.start_nl_cons:model.start_nl_cons+model.len_nl_cons]
end

###
# TODO
# MOI.get(model::Optimizer, ::MOI.SolveTime) = TODO
# Additional getters
#MOI.get(model::Optimizer, ::MOI.NodeCount)  = TODO
#MOI.get(model::Optimizer, ::MOI.BarrierIterations)  = TODO
#MOI.get(model::Optimizer, ::MOI.RelativeGap)  = TODO
#MOI.get(model::Optimizer, ::MOI.ObjectiveBound)  = TODO
