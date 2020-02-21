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

function reshop_declare_mathprgm(mp, mdl::Ptr{reshop_model}, NLoffset::Int)

    reshop_mp = emp_mp_alloc(mdl)

    CONFIG[:debug] && println("DEBUG: reshop_declare_mathprgm: tackling MP sense = $(mp.sense) objequ = $(mp.objequ) equs = $(mp.equs) vars = $(mp.vars)")

    (objequ, isNL) = mp.objequ
    if objequ > 0 || mp.objvar > 0
        typ = 0

        emp_mp_start(reshop_mp, typ)
        emp_mp_objdir(reshop_mp, mp.sense)
        if (objequ <= 0)
            error("wrong objective equation value $(mp.objequ)")
        end

        if isNL
            emp_mp_objequ(reshop_mp, objequ-1 + NLoffset)
        else
            emp_mp_objequ(reshop_mp, objequ-1)
        end

        emp_mp_objvar(reshop_mp, mp.objvar-1)

        for (eidx, isnl) in mp.equs
            if (eidx, isnl) != mp.objequ
                if isnl
                    emp_mp_constraint(reshop_mp, eidx-1+NLoffset)
                else
                    emp_mp_constraint(reshop_mp, eidx-1)
                end
            end
        end

        for vidx in mp.vars
            emp_mp_var(reshop_mp, vidx-1)
        end
    else
        typ = 2

        emp_mp_start(reshop_mp, typ)
        VIvars = keys(mp.matching)
        # TODO this should be an array of booleans ... --xhub
        equ_seen = fill((-1,false), length(mp.equs))
        sidx = 1
        for vidx in mp.vars
            if vidx in VIvars
                (eidx, isnl) = mp.matching[vidx]
                if isnl
                    emp_mp_vipair(reshop_mp, eidx-1+NLoffset, vidx-1)
                else
                    emp_mp_vipair(reshop_mp, eidx-1, vidx-1)
                end
                equ_seen[sidx] = (eidx, isnl)
                sidx += 1
                CONFIG[:debug] && println("DEBUG: eqn $eidx perp x[$vidx]")
            else
                # QVI case
                emp_mp_var(reshop_mp, vidx-1)
            end
        end

        for (eidx, isnl) in mp.equs
            if !((eidx, isnl) in equ_seen)
                CONFIG[:debug] && println("DEBUG: eqn $eidx is a constraint")
                if isnL
                    emp_mp_constraint(reshop_mp, eidx-1+NLoffset)
                else
                    emp_mp_constraint(reshop_mp, eidx-1)
                end
            end
        end
    end

    return reshop_mp
end
