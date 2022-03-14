# Use solver internal to ReSHOP

function reshop_setup_ownsolver()
    ctx = ctx_alloc(RHP_MDL_RHP)

    rhp_PATH_setfname(PATH_FNAME)

    return ctx
end

