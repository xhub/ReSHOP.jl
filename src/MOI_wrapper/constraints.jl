# MathOptInterface constraints

##################################################
# Generic constraint definition
#--------------------------------------------------
# Bound constraint on variables.

# TODO(Xhub) is this sufficient?
function chk_var_cons(model::Optimizer, vi::MOI.VariableIndex)
    check_inbounds(model, vi)
end


function chk_var_cons(model::Optimizer, vi::MOI.VariableIndex, rel::SS)
    check_inbounds(model, vi)

    if isa(rel, MOI.LessThan{Float64})
        if isnan(rel.upper)
           error("Invalid upper bound value $(rel.upper).")
       end
        if has_upper_bound(model, vi)
            error("Upper bound on variable $vi already exists.")
        end
        if is_fixed(model, vi)
            error("Variable $vi is fixed. Cannot also set upper bound.")
        end
    elseif isa(rel, MOI.GreaterThan{Float64})
        if isnan(rel.lower)
            error("Invalid lower bound value $(rel.lower).")
        end
        if has_lower_bound(model, vi)
            error("Lower bound on variable $vi already exists.")
        end
        if is_fixed(model, vi)
            error("Variable $vi is fixed. Cannot also set lower bound.")
        end
    elseif isa(rel, MOI.EqualTo{Float64})
        if isnan(rel.value)
            error("Invalid fixed value $(eq.value).")
        end
        if has_lower_bound(model, vi)
            error("Variable $vi has a lower bound. Cannot be fixed.")
        end
        if has_upper_bound(model, vi)
            error("Variable $vi has an upper bound. Cannot be fixed.")
        end
        if is_fixed(model, vi)
            error("Variable $vi is already fixed.")
        end
    elseif isa(rel, MOI.Interval{Float64})
        if isnan(rel.lower)
            error("Invalid lower value $(eq.lower).")
        end
        if isnan(rel.upper)
            error("Invalid upper value $(eq.upper).")
        end

        if has_lower_bound(model, vi)
            error("Variable $vi already has a lower bound.")
        end
        if has_upper_bound(model, vi)
            error("Variable $vi already has an upper bound.")
        end
        if is_fixed(model, vi)
            error("Variable $vi is already fixed.")
        end
    end


end

function reshop_set_bound(ctx::Ptr{context}, vidx::Integer, rel::MOI.LessThan{Float64})
    ctx_setvarub(ctx, vidx, rel.upper)
end

function reshop_set_bound(ctx::Ptr{context}, vidx::Integer, rel::MOI.GreaterThan{Float64})
    ctx_setvarlb(ctx, vidx, rel.lower)
end

function reshop_set_bound(ctx::Ptr{context}, vidx::Integer, rel::MOI.EqualTo{Float64})
    ctx_setvarfx(ctx, vidx, rel.value)
end

function reshop_set_bound(ctx::Ptr{context}, vidx::Integer, rel::MOI.Interval{Float64})
    ctx_setvarub(ctx, vidx, rel.upper)
    ctx_setvarlb(ctx, vidx, rel.lower)
end

function MOI.add_constraint(model::Optimizer, v::MOI.SingleVariable, rel::SS)
    vi = v.variable
    chk_var_cons(model, vi, rel)
    reshop_set_bound(model.ctx, vi.value-1, rel)
    return MOI.ConstraintIndex{MOI.SingleVariable, typeof(rel)}(vi.value)
end

function MOI.add_constraint(model::Optimizer,
                            func::MOI.ScalarAffineFunction{Float64}, set::LS)
    check_inbounds(model, func)
    eidx = rhp_add_equ(model.ctx)
    # Parse structure of constraint.
    vidx, coefs = canonical_linear_reduction(func)
    avar = _ensure_avar(model)
    rhp_avar_set(avar, vidx)
    rhp_equ_add_linear(model.ctx, eidx, avar, coefs)
    _set_rhs(model, eidx, func.constant, set)
    # Add constraint to index.
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(eidx+1)
end

function MOI.add_constraint(model::Optimizer,
                            func::MOI.ScalarQuadraticFunction{Float64}, set::LS)
    check_inbounds(model, func)
    eidx = rhp_add_equ(model.ctx)

    # Parse linear structure of constraint.
    # Parse quadratic term.
    qvar1, qvar2, qcoefs = canonical_quadratic_reduction(func)

    rhp_equ_add_quadratic(model.ctx, eidx, qvar1, qvar2, qcoefs)
    vidx, coefs = canonical_linear_reduction(func)

    avar = _ensure_avar(model)
    rhp_avar_set(avar, vidx)
    rhp_equ_add_lin_tree(model.ctx, eidx, coefs, avar, 1.)
    _set_rhs(model, eidx, func.constant, set)
    # Add constraints to index.
    ci = MOI.ConstraintIndex{typeof(func), typeof(set)}(eidx+1)
    model.quadfn_mapping[ci] = eidx
    return ci
end

function MOI.add_constraint(model::Optimizer,
                            func::MOI.VectorAffineFunction, set::MOI.AbstractVectorSet)
    # TODO: add check inbounds for VectorAffineFunction.
    current_m = ctx_numequ(model.ctx)
    num_cons = length(func.constants)
    # Add constraints inside KNITRO.
    rhp_add_equs(model.ctx, num_cons)

    # Parse vector affine expression.
    # The indices are already 0-based
    cidx, vidx, coeffs = canonical_vector_affine_reduction(func)
    # Get true indices
    cidx .+= current_m

    for (idx, eidx) in enumerate(cidx)
        ctx_add_lin_var(model.ctx, eidx, vidx[idx], coeffs[idx])
    end

    for (idx, cst) in enumerate(func.constants)
        _set_rhs(model, current_m-1+idx, cst, set)
    end

    # Add constraints to index.
    ci = MOI.ConstraintIndex{typeof(func), typeof(set)}(cidx[1]+1)
    model.vaf_mapping[ci] = unique(cidx)
    return ci
end

# Add second order cone constraint.
#function MOI.add_constraint(model::Optimizer,
#                            func::MOI.VectorAffineFunction, set::MOI.SecondOrderCone)
#    # Add constraints inside KNITRO.
#    rhp_add_equ(model.ctx)
#    eidx = ctx_numequ(model.ctx)
#
#    # Parse vector affine expression.
#    indexcoords, vidx, coefs = canonical_vector_affine_reduction(func)
#    constants = func.constants
#    # Distinct two parts of secondordercone.
#    # First row corresponds to linear part of SOC.
#    indlinear = indexcoords .== 0
#    indcone = indexcoords .!= 0
#    ncoords = length(constants) - 1
#    @assert ncoords == set.dimension - 1
#
#    # Load Second Order Conic constraint.
#    ## i) linear part
#    KN_set_con_upbnd(model.ctx, eidx, constants[1])
#    rhp_equ_add_linear_struct(model.ctx, eidx,
#                             vidx[indlinear], -coefs[indlinear])
#
#    ## ii) soc part
#    index_var_cone = vidx[indcone]
#    nnz = length(index_var_cone)
#    index_coord_cone = convert.(Cint, indexcoords[indcone] .- 1)
#    coefs_cone = coefs[indcone]
#    const_cone = constants[2:end]
#
#    rhp_equ_add_L2norm(model.ctx,
#                      eidx, ncoords, nnz,
#                      index_coord_cone,
#                      index_var_cone,
#                      coefs_cone,
#                      const_cone)
#
#    # Add constraints to index.
#    ci = MOI.ConstraintIndex{typeof(func), typeof(set)}(eidx)
#    model.constraint_mapping[ci] = vidx
#    return ci
#end

function MOI.add_constraint(model::Optimizer, vov::VOV, set::T) where T <: VLS
    chk_var_cons(model, vv)
    vidx = convert.(Cint, [v.value for v in vov.variables] .- 1)
    bnd = zeros(Float64, length(vidx))

    if isa(set, MOI.Zeros)
        ctx_setvarfx(model.ctx, vidx, bnd)
    elseif isa(set, MOI.Nonnegatives)
        ctx_setvarlb(model.ctx, vidx, bnd)
    elseif isa(set, MOI.Nonpositives)
        ctx_setvarub(model.ctx, vidx, bnd)
    end

    # TODO is this necessary?
#    ncons = MOI.get(model, MOI.NumberOfConstraints{VOV, T}())
    ci = MOI.ConstraintIndex{VOV, T}(vidx[1]+1)
    model.vov_mapping[ci] = vidx
    return ci
end

function MOI.add_constraint(model::Optimizer, vov::VOV, set::MOI.SecondOrderCone)

    # the variables have the following structure
    # (t, x) where t ≥ ||x||_2
    #

    vidx = [v.value - 1 for v in func.variables]

    soc_data = rhp_cone_soc(vidx[2:end])
    reshop_set_vartype(model.ctx, vidx[1], RHP_VARTYPE_SOC, soc_data)

    # Add constraints to index.
    # The index of the t variable is used, since the multiplier
    # would only be defined for it.
    ci = MOI.ConstraintIndex{typeof(func), typeof(set)}(vidx[1]+1)
    model.vov_mapping[ci] = vidx
    return ci
end

function MOI.add_constraint(model::Optimizer, vov::VOV,
                            set::Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}})
    @assert length(vov.variables) == length(set.weights)
    vidx = convert.(Cint, [v.value for v in vov.variables] .- 1)
    avar = _ensure_avar(model)
    rhp_avar_set(avar, vidx)
    reshop_set_vartype(model.ctx, avar, set, set.weights)
    ci = MOI.ConstraintIndex{VOV, typeof(set)}(vidx[1]+1)
    model.vov_mapping[ci] = vidx
    model.sos_sets[ci] = set
    return ci
end

##################################################
## Binary & Integer constraints.

# Define integer and boolean constraints.
# Those don't need to be stored in the constraint index
function MOI.add_constraint(model::Optimizer,
                            v::MOI.SingleVariable, S::Union{MOI.ZeroOne, MOI.Integer})
    vi = v.variable
    check_inbounds(model, vi)
    reshop_set_vartype(model.ctx, vi.value - 1, set_to_reshop[typeof(S)])
    return MOI.ConstraintIndex{MOI.SingleVariable, typeof(S)}(vi.value)
end

##################################################
## Constraint DualStart
function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex, value::Union{Real, Nothing})
    check_inbounds(model, ci)
    val = isa(value, Real) ? Cdouble(value) : 0.
    if ci in model.vaf_mapping
        ctx_setequmult(model.ctx, ci.value - 1, val)
    elseif check_inbounds(model, ci.value)
        ctx_setvarmult(model.ctx, ci.value - 1, val)
    else
        error("Unknown contraint $ci")
    end
    return
end

function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex{VOV, T}, value::Union{Real, Nothing}) where T
    check_inbounds(model, ci)
    val = isa(value, Real) ? Cdouble(value) : 0.

    if ci in model.vov_mapping
        for vidx in model.vov_mapping[ci]
            ctx_setvarmult(model.ctx, vidx, val)
        end
    else
        error("Unknown contraint $ci")
    end
    return
end

function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex{VOV, T}, value::AbstractVector{Real}) where T
    check_inbounds(model, ci)

    if ci in model.vov_mapping
        @assert length(value) == length(model.vov_mapping[ci])
        for (idx, vidx) in enumerate(model.vov_mapping[ci])
            ctx_setvarmult(model.ctx, vidx, Cdouble(value[idx]))
        end
    else
        error("Unknown contraint $ci")
    end
    return
end

function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex{VAF, T}, value::Union{Real, Nothing}) where T
    check_inbounds(model, ci)
    val = isa(value, Real) ? Cdouble(value) : 0.

    if ci in model.vaf_mapping
        for eidx in model.vaf_mapping[ci]
            ctx_setequmult(model.ctx, eidx, val)
        end
    else
        error("Unknown contraint $ci")
    end
    return
end

function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex{VAF, T}, value::AbstractVector{Real}) where T
    check_inbounds(model, ci)

    if ci in model.vaf_mapping
        @assert length(value) == length(model.vaf_mapping[ci])
        for (idx, eidx) in enumerate(model.vaf_mapping[ci])
            ctx_setequmult(model.ctx, eidx, Cdouble(value[idx]))
        end
    else
        error("Unknown contraint $ci")
    end
    return
end

##################################################
## Constraint naming
# TODO
#function MOI.set(model::Optimizer, ::MOI.ConstraintName, ci::MOI.VariableIndex{:< ALLV, <: Union{SS,LS,VLS}}, name::String)
#    ctx_set_equname(model.ctx, ci.value-1, name)
#end
#
#function MOI.set(model::Optimizer, ::MOI.ConstraintName, ci::MOI.VariableIndex{:< ALLV, <: Union{}}, name::String)
#    get(model.fake_cons_name, ci, EMPTYSTRING)
#end
