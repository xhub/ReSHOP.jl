function reshop_solve(mdl::Ptr{reshop_model}, mdl_solver::Ptr{reshop_model}, ctx_dest::Ptr{context}, solver_name::String)

    CONFIG[:solver_log] && hack_solver_log()

    # Set the modeltype
    res = ccall((:myo_analyze_modeltype, libreshop), Cint, (Ptr{reshop_model}, Ptr{filter_ops}), mdl, C_NULL)
    res != 0 && error("return code $res from ReSHOP")

    # TODO(xhub) HIGH why is true here?
    if true
        res = ccall((:reshop_transform, libreshop), Cint, (Ptr{reshop_model}, Ptr{reshop_model}), mdl, mdl_solver)
        res != 0 && error("return code $res from ReSHOP")

        if CONFIG[:export_gms]
#            ccall((:ctx_writemodel, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "validation.gms")
            ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
            ccall((:ctx_callsolver, libreshop), Cint, (Ptr{context},), ctx_dest)
        end

        ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{ReSHOP.context}, Cstring), ctx_dest, solver_name)

        res = ccall((:reshop_solve, libreshop), Cint, (Ptr{reshop_model},), mdl_solver)
    else
        res = ccall((:model_compress, libreshop), Cint, (Ptr{reshop_model}, Ptr{reshop_model}, Ptr{Cvoid}), mdl, mdl_solver, C_NULL)
        res != 0 && error("return code $res from ReSHOP")
        res = ccall((:ctx_exportmodel, libreshop), Cint, (Ptr{context}, Ptr{context}), ctx, ctx_dest)
        res != 0 && error("return code $res from ReSHOP")

        res = ccall((:ctx_callsolver, libreshop), Cint, (Ptr{context},), ctx_dest)
    end

    return res
end
