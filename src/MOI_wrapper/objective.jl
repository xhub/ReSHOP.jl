# MathOptInterface objective

MOI.get(model::Optimizer, ::MOI.ObjectiveSense) = model.sense

# Objective definition.
function add_objective!(model::Optimizer, objective::MOI.ScalarQuadraticFunction)
  avar = _ensure_avar(model)
  eidx = rhp_addequ_nocst(model.ctx, avar, objective)
  rhp_set_objeqn(model.ctx, eidx)
  reshop_set_cst(model.ctx, eidx, objective.constant)
  return
end

function add_objective!(model::Optimizer, objective::MOI.ScalarAffineFunction)
  avar = _ensure_avar(model)
  eidx = rhp_addequ_nocst(model.ctx, avar, objective)
  rhp_set_objeqn(model.ctx, eidx)
  reshop_set_cst(model.ctx, eidx, objective.constant)
  return
end

function add_objective!(model::Optimizer, var::MOI.SingleVariable)
    check_inbounds(model, var)
    # TODO(Xhub) because of ovf, we need to always add an equation
#    rhp_set_objvar(model.ctx, var.variable.value - 1)
    avar = _ensure_avar(model)
    eidx = rhp_addequ_nocst(model.ctx, avar, var)
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
