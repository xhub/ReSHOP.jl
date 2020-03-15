##############################################################################
# GAMS specific functions
##############################################################################

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

    res = ccall((:gams_setgamscntr, libreshop), Cint, (Ptr{context}, Cstring), ctx, joinpath(cur_dir, "gamscntr.dat"))
    res != 0 && error("return code $res from ReSHOP")

    # hm bad hack
    gamsdir = split(gamscntr_template, term_str)[29]

    CONFIG[:debug] && println("DEBUG: gamsdir is ``$gamsdir''")

    res = ccall((:gams_setgamsdir, libreshop), Cint, (Ptr{context}, Cstring), ctx, gamsdir)
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
