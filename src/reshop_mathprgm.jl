# TODO: delete?
function reshop_set_mathprgm_modeltype(m::ReSHOPMathProgBaseModel, idx)
    error("to implement")
#    discrete = any((m.vartypes .== :Int) + (m.vartypes .== :Bin) .> 0)
#    if discrete
#        if m.model_type == qcp
#            m.model_type = miqcp
#        elseif m.model_type == nlp
#             m.model_type = minlp
#        end
#    end
end

function _get_true_idx(ci::MOI.ConstraintIndex{T, <: MOI.AbstractSet}, NLoffset::Int) where T <: Union{MOI.ScalarAffineFunction{Float64},MOI.VectorAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64},MOI.VectorQuadraticFunction{Float64}}
  return ci.value-1
end

# TODO: put type info here
function _get_rhp_moi_mdl(mp)
  return mp.emp.backend.moi_backend
end

function _mp_addcons(mp, reshop_mp, ci::MOI.ConstraintIndex, NLoffset::Int)
  rhp_moi_mdl = _get_rhp_moi_mdl(mp)

  if ci in keys(rhp_moi_mdl.vaf_mapping)
    for ei in rhp_moi_mdl.vaf_mapping[ci]
      emp_mp_constraint(reshop_mp, ei)
    end
  else
    ei = _get_true_idx(ci, NLoffset)
    emp_mp_constraint(reshop_mp, ei)
  end
end

function _mp_addequs(mp, reshop_mp, ci::MOI.ConstraintIndex, NLoffset::Int)
  rhp_moi_mdl = _get_rhp_moi_mdl(mp)

  if ci in keys(rhp_moi_mdl.vaf_mapping)
    for ei in rhp_moi_mdl.vaf_mapping[ci]
      emp_mp_addequ(reshop_mp, ei)
    end
  else
    ei = _get_true_idx(ci, NLoffset)
    emp_mp_addequ(reshop_mp, ei)
  end
end

function _mp_addequvar(mp, NLoffset::Int, reshop_mp)
  # Iterate over the constraint and add them.
  for ci in mp.cons
      _mp_addcons(mp, reshop_mp, ci, NLoffset)
  end

  # These equations already have some metadata set
  for ci in mp.equs
      _mp_addequs(mp, reshop_mp, ci, NLoffset)
  end

  # Now the nonlinear equations
  if !isnothing(mp.nlp_data)
    for nlidx in mp.nlp_data.indices
      emp_mp_constraint(reshop_mp, nlidx+NLoffset)
    end
  end

  for idx in mp.vars
    emp_mp_var(reshop_mp, idx-1)
  end
end

function _mp_opt(mp, mdl::Ptr{reshop_model}, NLoffset::Int, reshop_mp,
                 objective::Union{MOI.ScalarQuadraticFunction,MOI.ScalarAffineFunction})
  typ = 0

  emp_mp_start(reshop_mp, typ)
  emp_mp_setsense(reshop_mp, mp.objectivesense)

  moi_mdl = _get_rhp_moi_mdl(mp)
  avar = _ensure_avar(moi_mdl)
  objequ = rhp_addequ_nocst(moi_mdl.ctx, avar, objective)
  # objequ is already a reshop index
  rhp_set_equasmapping(moi_mdl.ctx, objequ)
  reshop_set_cst(moi_mdl.ctx, objequ, objective.constant)
  emp_mp_objequ(reshop_mp, objequ)

  # Add all the variable and equations of this MP
  _mp_addequvar(mp, NLoffset, reshop_mp)
end

function _mp_opt(mp, mdl::Ptr{reshop_model}, NLoffset::Int, reshop_mp,
                 objvar::MOI.SingleVariable)
  typ = 0

  emp_mp_start(reshop_mp, typ)
  emp_mp_setsense(reshop_mp, mp.objectivesense)

  moi_mdl = _get_rhp_moi_mdl(mp)
  avar = _ensure_avar(moi_mdl)
  objequ = rhp_addequ_nocst(moi_mdl.ctx, avar, objvar)
  # objequ is already a reshop index
  rhp_set_equasmapping(moi_mdl.ctx, objequ)
  emp_mp_objequ(reshop_mp, objequ)

  # Add all the variable and equations of this MP
  _mp_addequvar(mp, NLoffset, reshop_mp)
end

function _mp_opt(mp, mdl::Ptr{reshop_model}, NLoffset::Int, reshop_mp, objequ::Int)
  typ = 0

  emp_mp_start(reshop_mp, typ)
  emp_mp_setsense(reshop_mp, mp.objectivesense)
  if (objequ <= 0)
    error("wrong objective equation value $(mp.objequ)")
  end

  # objequ is already a reshop index
  emp_mp_objequ(reshop_mp, objequ + NLoffset - 1)

  _mp_addequvar(mp, NLoffset, reshop_mp)
end

# JuMP.AbstractJuMPScalar

function reshop_declare_mathprgm(mp, mdl::Ptr{reshop_model}, NLoffset::Int)

    reshop_mp = emp_mp_alloc(mdl)

    CONFIG[:debug] && println("DEBUG: reshop_declare_mathprgm: tackling MP sense = $(mp.sense) objequ = $(mp.objequ) equs = $(mp.equs) vars = $(mp.vars)")

    objequ = mp.objective_function


    if mp.objectivesense != MOI.FEASIBILITY_SENSE #|| mp.objvar > 0
      _mp_opt(mp, mdl, NLoffset, reshop_mp, objequ)
    else
        typ = 2

        emp_mp_start(reshop_mp, typ)

        _mp_addequvar(mp, NLoffset, reshop_mp)

        # TODO: delete
#        VIvars = keys(mp.matching)
#        # TODO this should be an array of booleans ... --xhub
#        equ_seen = fill((-1,false), length(mp.equs))
#        sidx = 1
#        for vidx in mp.vars
#            if vidx in VIvars
#                (eidx, isnl) = mp.matching[vidx]
#                if isnl
#                    emp_mp_vipair(reshop_mp, eidx-1+NLoffset, vidx-1)
#                else
#                    emp_mp_vipair(reshop_mp, eidx-1, vidx-1)
#                end
#                equ_seen[sidx] = (eidx, isnl)
#                sidx += 1
#                CONFIG[:debug] && println("DEBUG: eqn $eidx perp x[$vidx]")
#            else
#                # QVI case
#                emp_mp_var(reshop_mp, vidx-1)
#            end
#        end
#
#        for (eidx, isnl) in mp.equs
#            if !((eidx, isnl) in equ_seen)
#                CONFIG[:debug] && println("DEBUG: eqn $eidx is a constraint")
#                if isnl
#                    emp_mp_constraint(reshop_mp, eidx-1+NLoffset)
#                else
#                    emp_mp_constraint(reshop_mp, eidx-1)
#                end
#            end
#        end
    end

    return reshop_mp
end
