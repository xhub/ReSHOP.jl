function MOI.delete(model::Optimizer, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
    rhp_delete_var(model.ctx, vi.value-1)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{<: Union{SF, SF}, <: LS})
    rhp_delete_eqn(model.ctx, ci.value-1)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{VAF, <: VLS})
    if haskey(model.vaf_mapping, ci)
        for eidx in model.vaf_mapping[ci]
            rhp_delete_eqn(model.ctx, eidx)
        end
    else
        error("Unknown contraint $ci")
    end
    return
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    ctx_setvarub(ctx, ci.value-1, Inf64)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}})
    ctx_setvarlb(ctx, ci.value-1, -Inf64)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}})
    ctx_setvarub(ctx, ci.value-1, Inf64)
    ctx_setvarlb(ctx, ci.value-1, -Inf64)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
    ctx_setvarub(ctx, ci.value-1, Inf64)
    ctx_setvarlb(ctx, ci.value-1, -Inf64)
end