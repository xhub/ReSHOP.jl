function reshop_solve(mdl::Ptr{reshop_model}, mdl_solver::Ptr{reshop_model}, ctx_dest::Ptr{context}, solver_name::String)

    CONFIG[:solver_log] && hack_solver_log()

    # Set the modeltype
    res = ccall((:myo_analyze_modeltype, libreshop), Cint, (Ptr{reshop_model}, Ptr{filter_ops}), mdl, C_NULL)
    res != 0 && error("return code $res from ReSHOP")

    # TODO(xhub) HIGH why is true here?
    res = ccall((:reshop_transform, libreshop), Cint, (Ptr{reshop_model}, Ptr{reshop_model}), mdl, mdl_solver)
    res != 0 && error("return code $res from ReSHOP")

    if CONFIG[:export_gms]
        ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
        ccall((:reshop_solve, libreshop), Cint, (Ptr{reshop_model},), mdl_solver)
    end

    ccall((:gams_set_solverstr, libreshop), Cint, (Ptr{ReSHOP.context}, Cstring), ctx_dest, solver_name)
    res = ccall((:reshop_solve, libreshop), Cint, (Ptr{reshop_model},), mdl_solver)

    return res
end
