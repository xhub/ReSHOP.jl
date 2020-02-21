function MOI.add_variable(model::Optimizer)
    avar = _ensure_avar(model)
    rhp_add_var(model.ctx, avar)
    return MOI.VariableIndex(1+reshop_avar_get(avar, UInt32(0)))
end

function MOI.add_variables(model::Optimizer, n::Integer)
    avar = _ensure_avar(model)
    rhp_add_var(model.ctx, n, avar)
    # We have faith that we have a compact avar
    vidx = 1+reshop_avar_get(avar, UInt32(0))
    return [MOI.VariableIndex(i) for i=vidx:vidx+n-1]
end


# MathOptInterface Variables

##################################################
## Getters
MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = ctx_numvar(model.ctx)

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:ctx_numvar(model.ctx)]
end

function check_inbounds(model::Optimizer, vi::MOI.VariableIndex)
    num_variables = ctx_numvar(model.ctx)
    if !(1 <= vi.value <= 1+num_variables)
        error("Invalid variable index $vi. ($num_variables variables in the model.)")
    end
    return
end

##################################################
## Check inbounds for safety.
function check_inbounds(model::Optimizer, var::MOI.SingleVariable)
    return check_inbounds(model, var.variable)
end

function check_inbounds(model::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        check_inbounds(model, term.variable_index)
    end
    return
end

function check_inbounds(model::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        check_inbounds(model, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds(model, term.variable_index_1)
        check_inbounds(model, term.variable_index_2)
    end
    return
end

function has_upper_bound(model::Optimizer, vi::MOI.VariableIndex)
    (lb, ub) = ctx_getvarbounds(model.ctx, vi.value-1)
    return isfinite(ub)
end

function has_lower_bound(model::Optimizer, vi::MOI.VariableIndex)
    (lb, ub) = ctx_getvarbounds(model.ctx, vi.value-1)
    return isfinite(lb)
end

function is_fixed(model::Optimizer, vi::MOI.VariableIndex)
    (lb, ub) = ctx_getvarbounds(model.ctx, vi.value-1)
    return isfinite(lb) && isfinite(ub) && abs(ub-lb) < eps()
end

function has_box_bounds(model::Optimizer, vi::MOI.VariableIndex)
    (lb, ub) = ctx_getvarbounds(model.ctx, vi.value-1)
    return isfinite(lb) && isfinite(ub)
end

##################################################
## PrimalStart
function MOI.set(model::Optimizer, ::MOI.VariablePrimalStart,
                 vi::MOI.VariableIndex, value::Union{Real, Nothing})
    check_inbounds(model, vi)
    val = isa(value, Real) ? Cdouble(value) : 0.
    ctx_setvarval(model.ctx, vi.value - 1, val)
end

##################################################
## Naming
function MOI.set(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex, name::String)
    check_inbounds(model, vi)
    ctx_setvarname(model.ctx, vi.value-1, name)
end

##################################################
## getter
function MOI.get(model::Optimizer, ::MOI.VariablePrimalStart, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return ctx_getvarval(model.ctx, vi.value-1)
end

function MOI.get(model::Optimizer, ::MathOptInterface.VariableIndex, name::String)
    return ctx_getvarname(model.ctx, vi.value-1)
end
