if VERSION < v"0.7"
    Cvoid = Void
end

function reshop_set_modeltype(m::ReSHOPMathProgModel)
    discrete = any((m.vartypes .== :Int) + (m.vartypes .== :Bin) .> 0)
    if discrete
        if m.model_type == qcp
            m.model_type = miqcp
        elseif m.model_type == nlp
             m.model_type = minlp
        end
    end
    reshop_set_modeltype(m.reshop_ctx, m.model_type)
    res = ccall((:ctx_setobjsense, libreshop), Cint, (Ptr{context}, Cint), m.reshop_ctx, sense_to_reshop[m.sense])
    res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_modeltype(ctx::Ptr{context}, model_type)
    res = ccall((:ctx_setmodeltype, libreshop), Cint, (Ptr{context}, Cint), ctx, model_type)
    res != 0 && error("return code $res from ReSHOP")
end

# TODO(xhub) reimplement writemodel using CONVERTD
function reshop_solve(ctx::Ptr{context}, ctx_dest::Ptr{context}, solver_name::String, emp::Ptr{empinfo}=Ptr{empinfo}(C_NULL))

    CONFIG[:solver_log] && hack_solver_log()

    if emp != C_NULL
        res = ccall((:reshop_transform, libreshop), Cint, (Ptr{empinfo}, Ptr{context}), emp, ctx_dest)
        res != 0 && error("return code $res from ReSHOP")

        if CONFIG[:export_gms]
#            ccall((:ctx_writemodel, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "validation.gms")
            ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
            ccall((:ctx_callsolver, libreshop), Cint, (Ptr{context},), ctx_dest)
        end

        ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{ReSHOP.context}, Cstring), ctx_dest, solver_name)

        res = ccall((:reshop_solve, libreshop), Cint, (Ptr{empinfo},), emp)
        return res
    else
        res = ccall((:model_compress, libreshop), Cint, (Ptr{context}, Ptr{context}, Ptr{empinfo}, Ptr{Cvoid}), ctx, ctx_dest, emp, C_NULL)
        res != 0 && error("return code $res from ReSHOP")
        res = ccall((:ctx_exportmodel, libreshop), Cint, (Ptr{context}, Ptr{context}), ctx, ctx_dest)
        res != 0 && error("return code $res from ReSHOP")

        if CONFIG[:export_gms]
#            ccall((:ctx_writemodel, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "validation.gms")
            ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
            ccall((:ctx_callsolver, libreshop), Cint, (Ptr{context},), ctx_dest)
        end

        # switch back to the default solver
        ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{ReSHOP.context}, Cstring), ctx_dest, solver_name)

        return ccall((:ctx_callsolver, libreshop), Cint, (Ptr{context},), ctx_dest)
    end
end

function reshop_setup_gams()
    ctx = ccall((:ctx_alloc, libreshop), Ptr{context}, (Cuint,), 0)

    gamscntr_template_file = joinpath(solverdata_dir, "gamscntr.dat")

    if !isfile(gamscntr_template_file)
        reshop_init_gams_solverdata()
    end
    if !isfile(gamscntr_template_file)
        error("Could not create template GAMS control file! Make sure that GAMS is properly installed and available via the system path")
    end

    gamscntr_template = read(gamscntr_template_file, String)

    gams_dir = mktempdir(pwd())
    cur_dir = gams_dir

    open(joinpath(gams_dir, "gamscntr.dat"), "w") do gamscntr_file
        println(gamscntr_file, replace(gamscntr_template, r"@@SUB@@" => cur_dir))
    end

    # we need an empty Matrixfile
    touch(joinpath(cur_dir, "gamsmatr.dat"))

    res = ccall((:gams_set_gamscntr, libreshop), Cint, (Ptr{context}, Cstring), ctx, joinpath(cur_dir, "gamscntr.dat"))
    res != 0 && error("return code $res from ReSHOP")

    # hm bad hack
    gamsdir = split(gamscntr_template, term_str)[29]

    CONFIG[:debug] && println("DEBUG: gamsdir is ``$gamsdir''")

    res = ccall((:gams_set_gamsdir, libreshop), Cint, (Ptr{context}, Cstring), ctx, gamsdir)
    res != 0 && error("return code $res from ReSHOP")

    ENV["PATH"] *= ":" * gamsdir

    # This is needed to prevent the listing of the Process directory
    ENV["DEBUG_PGAMS"] = '0'

    return (ctx, gams_dir)
end

function reshop_init_gams_solverdata()
    if !isdir(solverdata_dir)
        error("No directory named $solverdata_dir. ")
    end

    gamscntr = joinpath(solverdata_dir, "gamscntr.dat")
    if isfile(gamscntr)
        return
    end

    gms_file = joinpath(solverdata_dir, "dummy.gms")

    mktempdir(solverdata_dir) do substr
        # use run(`gams $gms_file scrdir=$substr lo=0 curdir=$substr`) ?
        run(`gams $gms_file scrdir=$substr lo=0`)

        open(gamscntr, "w") do out_gamscntr
            input = read(joinpath(substr, "gamscntr.dat"), String)
            input = replace(input, pwd() => "@@SUB@@")
            println(out_gamscntr, replace(input, substr => "@@SUB@@"))
        end
    end
end
