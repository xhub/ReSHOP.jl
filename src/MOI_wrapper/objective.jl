# MathOptInterface objective

MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

# Objective definition.
function add_objective!(model::Optimizer, objective::MOI.ScalarQuadraticFunction)
    eidx = rhp_add_equ(model.ctx)
    # We parse the expression passed in arguments.
    qvidx1, qvidx2, qcoefs = canonical_quadratic_reduction(objective)
    lvidx, lcoefs = canonical_linear_reduction(objective)
    rhp_equ_add_quadratic(model.ctx, eidx, qvidx1, qvidx2, qcoefs)
    avar = _ensure_avar(model)
    rhp_avar_set(avar, lvidx)
    rhp_equ_add_linear_chk(model.ctx, eidx, avar, lcoefs)
    reshop_set_rhs(model.ctx, eidx, -objective.constant)
    rhp_set_objeqn(model.ctx, eidx)
    return
end

function add_objective!(model::Optimizer, objective::MOI.ScalarAffineFunction)
    eidx = rhp_add_equ(model.ctx)
    lvidx, lcoefs = canonical_linear_reduction(objective)
    avar = _ensure_avar(model)
    rhp_avar_set(avar, lvidx)
    rhp_equ_add_linear(model.ctx, eidx, avar, lcoefs)
    reshop_set_rhs(model.ctx, eidx, -objective.constant)
    rhp_set_objeqn(model.ctx, eidx)
    return
end

function add_objective!(model::Optimizer, var::MOI.SingleVariable)
    check_inbounds(model, var)
    # TODO(Xhub) because of ovf, we need to always add an equation
#    rhp_set_objvar(model.ctx, var.variable.value - 1)
    eidx = rhp_add_equ(model.ctx)
    avar = _ensure_avar(model)
    rhp_avar_set(avar, var.variable.value - 1)
    rhp_equ_add_linear(model.ctx, eidx, avar, [1.,])
    rhp_set_objeqn(model.ctx, eidx)
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction,
                 func::Union{MOI.SingleVariable, MOI.ScalarAffineFunction,
                             MOI.ScalarQuadraticFunction})
    add_objective!(model, func)
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    model.sense = sense
    reshop_set_objsense(model.ctx, sense)
end
