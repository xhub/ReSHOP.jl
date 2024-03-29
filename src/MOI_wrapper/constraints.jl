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
    _set_cst(model, eidx, func.constant, set)
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
    # TODO(Xhub) change this when quadratic constrains are supported
    rhp_equ_add_lin_tree(model.ctx, eidx, coefs, avar, 1.)
    _set_cst(model, eidx, func.constant, set)
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
    # Add constraints inside ReSHOP.
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
        _set_cst(model, current_m-1+idx, cst, set)
    end

    # Add constraints to index.
    ci = MOI.ConstraintIndex{typeof(func), typeof(set)}(cidx[1]+1)
    model.vaf_mapping[ci] = unique(cidx)
    return ci
end

# Add second order cone constraint.
#function MOI.add_constraint(model::Optimizer,
#                            func::MOI.VectorAffineFunction, set::MOI.SecondOrderCone)
#    # Add constraints inside ReSHOP.
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
    vidx = convert.(Cint, [v.value for v in vov.variables] .- 1)
    chk_var_cons.(model, vov.variables)
    bnd = zeros(Float64, length(vidx))

    if isa(set, MOI.Zeros)
        ctx_setvarfx.(model.ctx, vidx, bnd)
    elseif isa(set, MOI.Nonnegatives)
        ctx_setvarlb.(model.ctx, vidx, bnd)
    elseif isa(set, MOI.Nonpositives)
        ctx_setvarub.(model.ctx, vidx, bnd)
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

function MOI.add_constraint(model::Optimizer, vaf::VAF, set::MOI.Complements)

    # the vaf needs to be parsed
    idxs, vars, coeffs = canonical_vector_affine_reduction(vaf)
    # TODO: for 0.22, we must set dim to set.dimension/2
    dim = set.dimension

    rhp_vars = Vector{RHP_IDXT}(undef, dim)
    equs = Vector{Tuple{Vector{RHP_IDXT},Vector{Float64}}}(undef, dim)

    for i in 1:dim
        equs[i] = (Vector{RHP_IDXT}(undef, 0), Vector{Float64}(undef, 0))
    end

    for (idx, var, coeff) in zip(idxs, vars, coeffs)
        # past dim, we have the complementary variable
        if idx >= dim
            rhp_vars[idx-dim+1] = var
         else
            push!(equs[idx+1][1], var)
            push!(equs[idx+1][2], coeff)
        end
    end

    # Now we define all the equations
    avar = _ensure_avar(model)
    current_m = ctx_numequ(model.ctx)
    # Add constraints inside ReSHOP.
    rhp_add_equs(model.ctx, dim)

    # Add all relations
    for i in 1:dim
        ei = current_m+i-1
        vi = rhp_vars[i]
        vis_fn, coeffs = equs[i]
        rhp_avar_set(avar, vis_fn)
        rhp_equ_add_linear(model.ctx, ei, avar, coeffs)
        rhp_set_equasmapping(model.ctx, ei)
        rhp_set_perp(model.ctx, ei, vi)
        reshop_set_cst(model.ctx, ei, vaf.constants[i])
    end


   # Add constraints to index.
    # The index of the t variable is used, since the multiplier
    # would only be defined for it.
    ci = MOI.ConstraintIndex{typeof(vaf), typeof(set)}(current_m+1)

    if dim > 1
      model.vaf_mapping[ci] = collect(current_m:current_m+dim-1)
    end

    return ci
end

function MOI.add_constraint(model::Optimizer, vov::VOV, set::MOI.Complements)

    # TODO: for 0.22, we must set dim to set.dimension/2
    dim = set.dimension

    rhp_vars = RHP_IDXT[vi.value-1 for vi in vov.variables]

    # Now we define all the equations
    current_m = ctx_numequ(model.ctx)

    # Add constraints inside ReSHOP
    rhp_add_equs(model.ctx, dim)

    # Add all relations
    for i in 1:dim
        ei = current_m+i-1
        vi_fn = rhp_vars[i]
        vi = rhp_vars[i+dim]
        ctx_add_lin_var(model.ctx, ei, vi_fn, 1.)
        rhp_set_equasmapping(model.ctx, ei)
        rhp_set_perp(model.ctx, ei, vi)
    end


    # Add constraints to index.
    # The index of the t variable is used, since the multiplier
    # would only be defined for it.
    ci = MOI.ConstraintIndex{VOV, typeof(set)}(current_m+1)

    if dim > 1
      model.vov_mapping[ci] = collect(current_m:current_m+dim-1)
    end

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
## Change some constraints
function MOI.set(model::Optimizer, ::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.SingleVariable, S}, new_set::S) where S <: SS
    reshop_set_bound(model.ctx, ci.value-1, new_set);
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S <: SS
    vidx = ci.value-1
    ctx_setvarlb(model.ctx, vidx, -Inf64)
    ctx_setvarub(model.ctx, vidx, Inf64)
end

function MOI.delete(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, <:Union{MOI.ZeroOne, MOI.Integer}})
    reshop_set_vartype(model.ctx, ci.value - 1, RHP_VARTYPE_X)
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
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, S}, value::Union{Real, Nothing}) where S <: SS

    val = isa(value, Real) ? Cdouble(value) : 0.
    ctx_setvarmult(model.ctx, ci.value - 1, val)
end

function MOI.set(model::Optimizer, ::MOI.ConstraintDualStart,
                 ci::MOI.ConstraintIndex, value::Union{Real, Nothing})
    check_inbounds(model, ci)
    val = isa(value, Real) ? Cdouble(value) : 0.
    if ci in model.vaf_mapping
        ctx_setequmult(model.ctx, ci.value - 1, val)
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
# getters
#

function MOI.get(model::Optimizer, loc::MOI.ListOfConstraints)
    @warn "MOI.get(model::Optimizer, loc::MOI.ListOfConstraints) is not yet implemented for ReSHOP"
    list = Vector{Tuple{DataType,DataType}}()
    return list
    for S in (
        MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64},
        MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}, MOI.Integer, MOI.ZeroOne
       )
        nb = MOI.get(model, MOI.NumberOfConstraints{MOI.SingleVariable, S}())
        if nb > 0
            push!(list, (MOI.SingleVariable, S))
        end
    end
    for S in (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros, MOI.SecondOrderCone)
        nb = MOI.get(model, MOI.NumberOfConstraints{VAF, S}())
        if nb > 0
            push!(list, (VAF, S))
        end
    end
    for S in (MOI.Nonnegatives,  MOI.Nonpositives, MOI.Zeros, MOI.SOS1, MOI.SOS2)
         nb = MOI.get(model, MOI.NumberOfConstraints{VOV, S}())
         if nb > 0
             push!(list, (VOV, S))
         end
    end
    for F in (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction)
        for S in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64})
            nb = MOI.get(model, MOI.NumberOfConstraints{F, S}())
            if nb > 0
                push!(list, (F, S))
            end
        end
    end

    return list
#    for S in (MOI.SecondOrderCone, MOI.RotatedSecondOrderCone, MOI.ExponentialCone,
#              MOI.DualExponentialCone, MOI.PowerCone, MOI.DualPowerCone)
#        nb = MOI.get(model, MOI.ScalarAffineFunction
end

function _get_var_set(model::Optimizer, S::Type{<:Union{MOI.EqualTo, MOI.GreaterThan}}, lb::Float64, ub::Float64)
    return S(lb)
end

function _get_var_set(model::Optimizer, S::Type{<:MOI.LessThan}, lb::Float64, ub::Float64)
    return S(ub)
end

# TODO semiinteger should have a different codepath
function _get_var_set(model::Optimizer, S::Type{<:Union{MOI.Interval, MOI.Semicontinuous, MOI.Semiinteger}}, lb::Float64, ub::Float64)
    return S(lb, ub)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S
    lb, ub = ctx_getvarbounds(model.ctx, ci.value-1)
    return _get_var_set(model, S, lb, ub)
end

##################################################
#Constraint validity check
#

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex)
    error("MOI.is_valid not implemented for type $(ci)")
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}})
    return has_upper_bound(model, MOI.VariableIndex(ci.value)) && !is_fixed(model, MOI.VariableIndex(ci.value))
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    return has_lower_bound(model, MOI.VariableIndex(ci.value)) && !is_fixed(model, MOI.VariableIndex(ci.value))
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}})
    return is_fixed(model, MOI.VariableIndex(ci.value))
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
    return has_box_bounds(model, MOI.VariableIndex(ci.value))
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer})
    return reshop_getvartype(model.ctx, ci.value-1) == RHP_VARTYPE_I
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne})
    return reshop_getvartype(model.ctx, ci.value-1) == RHP_VARTYPE_B
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semicontinuous{Float64}})
    return reshop_getvartype(model.ctx, ci.value-1) == RHP_VARTYPE_SC
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semiinteger{Float64}})
    return reshop_getvartype(model.ctx, ci.value-1) == RHP_VARTYPE_SI
end

function MOI.is_valid(model::Optimizer, ci::MOI.ConstraintIndex{VOV, <:VLS})
    allvars = get(model.vov_mapping, ci, Cuint[])
    if length(allvars) == 0
        return false
    end

    return !any([MOI.VariableIndex(1+vi) in model.delvars for vi in allvars])
end

##################################################
## Constraint naming
# TODO
function MOI.set(model::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex{<:SF,<:LS}, name::String)
    ctx_setequname(model.ctx, ci.value-1, name)
end

function MOI.set(model::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex{MOI.SingleVariable,<:SS}, name::String)
    ctx_setequname(model.ctx, ci.value-1, name)
end

function MOI.set(model::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex, name::String)
    println("unsupported constraint named $name at index $(ci.value) of type $ci")
    model.fake_cons_name[ci] = name
end

function MOI.get(model::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex{<:SF, <:LS})
    return ctx_getequname(model.ctx, ci.value-1)
end

function MOI.get(model::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex)
    return get(model.fake_cons_name, ci, "")
end

function MOI.get(model::Optimizer, ::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{MOI.SingleVariable,<: SS})
    MOI.throw_if_not_valid(model, ci)
    return MOI.SingleVariable(MOI.VariableIndex(ci.value))
end
