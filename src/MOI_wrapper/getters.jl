const EMPTYSTRING = ""

##################################################
#
##################################################
## Getters
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints) = ctx_m(model.ctx)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.ZeroOne}) =
    rhp_get_nb_bin_var(model.ctx)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.Integer}) =
    rhp_get_nb_int_var(model.ctx)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.SingleVariable, S}) where S <: SS =
    rhp_get_nb_var(model.ctx, S)
# TODO: a bit hacky, but that should work for MOI Test.
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VAF, MOI.Nonnegatives}) =
    sum(typeof.(collect(keys(model.vaf_mapping))) .== MOI.ConstraintIndex{VAF, MOI.Nonnegatives})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VAF, MOI.Nonpositives}) =
    sum(typeof.(collect(keys(model.vaf_mapping))) .== MOI.ConstraintIndex{VAF, MOI.Nonpositives})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VAF, MOI.Zeros}) =
    sum(typeof.(collect(keys(model.vaf_mapping))) .== MOI.ConstraintIndex{VAF, MOI.Zeros})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VAF, MOI.SecondOrderCone}) =
    sum(typeof.(collect(keys(model.vov_mapping))) .== MOI.ConstraintIndex{VAF, MOI.SecondOrderCone})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VOV, T}) where T <: VLS =
    sum(typeof.(collect(keys(model.vov_mapping))) .== MOI.ConstraintIndex{VOV, T})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{VOV, S}) where S <: Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}} =
    rhp_get_nb_var(model.ctx, S)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, S}) where S <: LS  =
    rhp_get_nb_lequ(model.ctx, S)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarQuadraticFunction{Float64}, S}) where S <: LS  =
    sum(typeof.(collect(keys(model.quadfn_mapping))) .== MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S})

function MOI.get(model::Optimizer, ::MOI.ConstraintSet, ci::MOI.ConstraintIndex{VOV, S}) where S <: Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}
    return model.sos_sets[ci]
end

function MOI.get(model::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    return ctx_getvarname(model.ctx, vi, name)
end

function MOI.get(model::Optimizer, ::MOI.VariableName, vi::MOI.ConstraintIndex)
    check_inbounds(model, vi)
    return ctx_getvarname(model.ctx, vi, name)
end
