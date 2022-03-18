function _ensure_avar(model)
    if (model.avar_cache != nothing)
        return model.avar_cache
    end
    model.avar_cache = reshop_avar()
    return model.avar_cache
end


function _set_cst(model::Optimizer, eidx, fn_cst, set::LS)
    # Add bound to constraint.
    if isa(set, MOI.LessThan{Float64})
        val = set.upper
    elseif isa(set, MOI.GreaterThan{Float64})
        val = set.lower
    elseif isa(set, MOI.EqualTo{Float64})
        val = set.value
    end

    # Set the rhs
    reshop_set_equtype(model.ctx, eidx, set_to_reshop[typeof(set)])
    reshop_set_cst(model.ctx, eidx, fn_cst - val)

end

function _set_cst(model::Optimizer, eidx, fn_cst, set::VLS)
    # Add bound to constraint.

    # Set the rhs
    reshop_set_equtype(model.ctx, eidx, set_to_reshop[typeof(set)])
    reshop_set_cst(model.ctx, eidx, fn_cst)
end

function get_status(ctx::Ptr{context})
    # beaware of the 1-index!
    solvestat = rhp_get_solvestat(ctx)
    modelstat = rhp_get_modelstat(ctx)

    if solvestat < 1 || solvestat > length(solver_stat)
      println("solvestat = $(solvestat) is outside of range")
      solvestatsymb = :InternalError
    else
      solvestatsymb = solver_stat[rhp_get_solvestat(ctx)]
    end

    if modelstat < 1 || modelstat > length(model_stat)
      println("modelstat = $(modelstat) is outside of range")
      modelstatsymb = :ErrorUnknown
    else
      modelstatsymb = model_stat[rhp_get_modelstat(ctx)]
    end

    return (solvestatsymb, modelstatsymb)
end

function is_valid_index(idx::RHP_IDXT)
    return idx >= 0
end

# Adding the constant part requires setting the equation type
function rhp_addequ_nocst(ctx, avar, func::MOI.ScalarAffineFunction)
  eidx = rhp_add_equ(ctx)
  lvidx, lcoefs = canonical_linear_reduction(func)
  rhp_avar_set(avar, lvidx)
  rhp_equ_add_linear(ctx, eidx, avar, lcoefs)
  return eidx
end

function rhp_addequ_nocst(ctx, avar, func::MOI.ScalarQuadraticFunction)
    eidx = rhp_add_equ(ctx)
    # We parse the expression passed in arguments.
    qvidx1, qvidx2, qcoefs = canonical_quadratic_reduction(func)
    lvidx, lcoefs = canonical_linear_reduction(func)
    rhp_equ_add_quadratic(ctx, eidx, qvidx1, qvidx2, qcoefs)
    rhp_avar_set(avar, lvidx)
    rhp_equ_add_linear_chk(ctx, eidx, avar, lcoefs)
    return eidx
end

function rhp_addequ_nocst(ctx, avar, var::MOI.SingleVariable)
    eidx = rhp_add_equ(ctx)
    rhp_avar_set(avar, var.variable.value - 1)
    rhp_equ_add_linear(ctx, eidx, avar, [1.,])
    return eidx
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}})
     return true
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}})
     return ctx_getvarbstat(ctx, ci.value-1) == RHP_BASIS_STATUS_UPPER
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}})
     return ctx_getvarbstat(ctx, ci.value-1) == RHP_BASIS_STATUS_LOWER
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
     bstat = ctx_getvarbstat(ctx, ci.value-1)
     return bstat == RHP_BASIS_STATUS_LOWER || bstat == RHP_BASIS_STATUS_UPPER
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{<:SF, MOI.EqualTo{Float64}})
     return true
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{<:SF, MOI.LessThan{Float64}})
     return ctx_getequbstat(ctx, ci.value-1) == RHP_BASIS_STATUS_UPPER
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{<:SF, MOI.GreaterThan{Float64}})
     return ctx_getequbstat(ctx, ci.value-1) == RHP_BASIS_STATUS_LOWER
end

function is_constraint_active(ctx::Ptr{context}, ci::MOI.ConstraintIndex{<:SF, MOI.Interval{Float64}})
     bstat = ctx_getequbstat(ctx, ci.value-1)
     return bstat == RHP_BASIS_STATUS_LOWER || bstat == RHP_BASIS_STATUS_UPPER
end

# TODO: fix this
#function reset_var_bnd(ctx::Ptr{context}, ci::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}})
#     ctx
#end

# Utils for MathOptInterface
#
##################################################
# Import legacy from LinQuadOptInterface to ease the integration
# of ReSHOP quadratic and linear facilities.
##################################################
# URL: https://github.com/JuliaOpt/LinQuadOptInterface.jl
#
# LICENSE:
# MIT License
# Copyright (c) 2017 Oscar Dowson, Joaquim Dias Garcia and contributors
##################################################

function reduce_duplicates(rows::Vector{T}, cols::Vector{T}, vals::Vector{S}) where T where S
    @assert length(rows) == length(cols) == length(vals)
    for i in 1:length(rows)
        if rows[i] > cols[i]
            tmp = rows[i]
            rows[i] = cols[i]
            cols[i] = tmp
        end
    end
    findnz(dropzeros!(sparse(rows, cols, vals)))
end

"""
    canonical_quadratic_reduction(func::ScalarQuadraticFunction)

Reduce a ScalarQuadraticFunction into three arrays, returned in the following
order:
 1. a vector of quadratic row indices
 2. a vector of quadratic column indices
 3. a vector of quadratic coefficients

Warning: we assume in this function that all variables are correctly
ordered, that is no deletion or swap has occured.
"""
function canonical_quadratic_reduction(func::MOI.ScalarQuadraticFunction)
    quad_columns_1, quad_columns_2, quad_coefficients = (
        RHP_IDXT[term.variable_index_1.value for term in func.quadratic_terms],
        RHP_IDXT[term.variable_index_2.value for term in func.quadratic_terms],
        [term.coefficient for term in func.quadratic_terms]
    )
    # MOI stores the diagonal element with a coefficient of 2.
    # This is most likely very smart for derivative, but it is not a panacee
    for i in 1:length(quad_coefficients)
        if quad_columns_1[i] == quad_columns_2[i]
            quad_coefficients[i] *= .5
        end
    end
    r, c, v = reduce_duplicates(quad_columns_1, quad_columns_2, quad_coefficients)
    return (r .- 1., c .- 1., v)
end

"""
    canonical_linear_reduction(func::Quad)

Reduce a ScalarQuadraticFunction into two arrays, returned in the following
order:
 1. a vector of linear column indices
 2. a vector of linear coefficients

Warning: we assume in this function that all variables are correctly
ordered, that is no deletion or swap has occured.
"""
function canonical_linear_reduction(func::MOI.ScalarQuadraticFunction)
    f = MOIU.canonical(func)
    affine_columns = RHP_IDXT[term.variable_index.value - 1 for term in f.affine_terms]
    affine_coefficients = [term.coefficient for term in f.affine_terms]
    return affine_columns, affine_coefficients
end

function canonical_linear_reduction(func::MOI.ScalarAffineFunction)
    f = MOIU.canonical(func)
    affine_columns = RHP_IDXT[term.variable_index.value - 1 for term in f.terms]
    affine_coefficients = [term.coefficient for term in f.terms]
    return affine_columns, affine_coefficients
end

function canonical_vector_affine_reduction(func::MOI.VectorAffineFunction)
    index_cols = RHP_IDXT[]
    index_vars = RHP_IDXT[]
    coeffs = Float64[]

    for t in func.terms
        push!(index_cols, t.output_index - 1)
        push!(index_vars, t.scalar_term.variable_index.value - 1)
        push!(coeffs, t.scalar_term.coefficient)
    end
    return index_cols, index_vars, coeffs
end

function get_solverstack(model::Optimizer)
  solver_stack = model.solver_stack
  if solver_stack == ""
    solver_stack = global_solverstack
  end
end
