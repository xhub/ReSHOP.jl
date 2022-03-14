##############################################################################
# GAMS specific functions
##############################################################################

# Remove the stored gamscnrt.dat file
function reset_gamscntr()
  gamscntr_template_file = joinpath(solverdata_dir, "gamscntr.dat")
  if isfile(gamscntr_template_file)
    rm(gamscntr_template_file, force=true)
  end
end

function reshop_setup_gams()
    ctx = ctx_alloc(RHP_MDL_GAMS)

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

function reshop_init_gams_solverdata(force=false)
    curdir = pwd()

    if !isdir(solverdata_dir)
        error("No directory named $solverdata_dir. ")
    end

    gamscntr = joinpath(solverdata_dir, "gamscntr.dat")
    if isfile(gamscntr)
        !force  && return
    end

    gms_file = joinpath(solverdata_dir, "dummy.gms")

#    println("Running the gams command to initial bootstrap file\n")

    mktempdir(solverdata_dir) do substr
        cd(solverdata_dir)
        # use run(`gams $gms_file scrdir=$substr lo=0 curdir=$substr`) ?
        run(`gams $gms_file scrdir=$substr`)

        open(gamscntr, "w") do out_gamscntr
            input = read(joinpath(substr, "gamscntr.dat"), String)
            # the order of subsitution matters here
            input = replace(input, substr => "@@SUB@@")
            input = replace(input, pwd() => "@@SUB@@")
            println(out_gamscntr, input)
        end
    end

#    println("\nGAMS boostrap file successfully generated")
    cd(curdir)

end
