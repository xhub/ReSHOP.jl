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
    sum(typeof.(collect(keys(model.vov_mapping))) .== MOI.ConstraintIndex{VOV, S})
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, S}) where S <: LS  =
    rhp_get_nb_lequ(model.ctx, S)
MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarQuadraticFunction{Float64}, S}) where S <: LS  =
    sum(typeof.(collect(keys(model.quadfn_mapping))) .== MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S})

function MOI.get(model::Optimizer, ::MOI.ConstraintSet, ci::MOI.ConstraintIndex{VOV, S}) where S <: Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}
    return model.sos_sets[ci]
end

function MOI.get(model::Optimizer, ::MOI.VariableName, ci::MOI.ConstraintIndex)
    check_inbounds(model, ci)
    return ctx_getequname(model.ctx, ci.value-1)
end

_chk_citype = Dict(
     MOI.LessThan{Float64} => (false, true),
     MOI.GreaterThan{Float64} => (true, false),
     MOI.EqualTo{Float64} => (true, true),
     MOI.Interval{Float64} => (true, true),
    )

_cone_moi_to_rhp = Dict(
                        MOI.LessThan{Float64} => RHP_CONE_R_MINUS,
                        MOI.GreaterThan{Float64} => RHP_CONE_R_PLUS,
                        MOI.EqualTo{Float64} => RHP_CONE_0
                       )

function chk_equ_citype(ctx, eidx, ::Type{MOI.ConstraintIndex{MOI.SingleVariable, S}}) where S<:SS
    (lb, ub) = ctx_getvarbounds(ctx, eidx)
    has_lb = isfinite(lb)
    has_ub = isfinite(ub)
    return (has_lb, has_ub) == _chk_citype[S]
end

function chk_equ_citype(ctx, eidx, ::Type{MOI.ConstraintIndex{F, S}}) where  {F<:SF, S<:LS}
    type, cone = reshop_get_equtype(ctx, eidx)
    res = (type == 2) && (cone == get(_cone_moi_to_rhp, S, RHP_CONE_NONE))
    return res
end

function MOI.get(model::Optimizer, ci_type::Type{MOI.ConstraintIndex{MOI.SingleVariable, S}}, cstr::String) where S<:SS
    eidx = ctx_getequbyname(model.ctx, cstr)
    if eidx >= 0
        if (!chk_equ_citype(model.ctx, eidx, ci_type))
            return nothing
        end
        return MOI.ConstraintIndex{MOI.SingleVariable, S}(eidx+1)
    elseif eidx == -2
        error("Multiple variables have the name $(name).")
    else
        return nothing
    end
end

function MOI.get(model::Optimizer, ci_type::Type{MOI.ConstraintIndex{F, S}}, cstr::String) where {F<:SF, S<:LS}
    eidx = ctx_getequbyname(model.ctx, cstr)
    if eidx >= 0
        chk_equ_citype(model.ctx, eidx, ci_type)
        return MOI.ConstraintIndex{F, S}(eidx+1)
    elseif eidx == -2
        error("Multiple variables have the name $(name).")
    else
        return nothing
    end
end


