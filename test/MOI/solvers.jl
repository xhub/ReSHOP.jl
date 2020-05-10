lp_solvers = Any[]
ip_solvers = Any[]
ip_dual_solvers = Any[]
semi_solvers = Any[]
sos_solvers = Any[]
lazy_solvers, lazy_soc_solvers, lazylocal_solvers, cut_solvers, cutlocal_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[], Any[], Any[], Any[]
quad_solvers = Any[]
quad_soc_solvers = Any[]
rsoc_solvers = Any[]
nlp_solvers = Any[]
convex_nlp_solvers = Any[]
minlp_solvers = Any[]
sdp_solvers = Any[]

#lp_solvers = ["baron", "cbc", "conopt", "cplex", "knitro", "minos", "pathnlp", "mosek", "xpress"]
lp_solvers = ["baron", "conopt", "knitro", "minos", "pathnlp", "mosek", "xpress"]

if Sys.islinux() || Sys.iswindows()
    push!(lp_solvers, "xa")
end

lp_exclude = Dict("knitro" => ["linear8c", "linear8b"])

quad_solvers = [("baron", "local"),
                ("conopt", "local"),
#                ("cplex", "nodual"),
                ("knitro", ""),
                ("minos", "local"),
                ("pathnlp", "local"),
                ("mosek", ""),
                ("xpress", "")]
quad_exclude = Dict("baron" => ["qcp1"], # SLOW_PROGRESS
                    "conopt" => ["socp1", "qcp1"],
#                    "cplex" => ["socp1"],
                    "minos" => ["socp1"],
                    "pathnlp" => ["socp1"],
                    "xpress" => ["qcp1"])
# mosek fails?
mip_solvers = ["xpress", "scip"]#, "cplex", "cbc"]

global_solvers = [("knitro", "")]

# fails with "'gmsbaxnx.exe' is not recognized as an internal or external command, operable program or batch file"
if !Sys.iswindows()
    push!(global_solvers, ("baron", ""))
end

# Investigate if we can solve those
#soc_solvers = copy(quad_solvers)
soc_solvers = Any[]

push!(convex_nlp_solvers, "")
push!(convex_nlp_solvers, "conopt")
# buggy on my system
#push!(convex_nlp_solvers, "ipopt")
push!(convex_nlp_solvers, "minos")
push!(convex_nlp_solvers, "snopt")
push!(convex_nlp_solvers, "mosek")

# push!(nlp_solvers, "baron")
push!(nlp_solvers, "knitro")

append!(convex_nlp_solvers, nlp_solvers)

push!(minlp_solvers, "")

# Default configuration.
const config_mosek = MOIT.TestConfig(atol=1e-4, rtol=1e-4,
                               optimal_status=MOI.OPTIMAL,
                               query=false,
                               infeas_certificates=false, # Do not ask for infeasibility certificates.
                               modify_lhs=false)

const config_mosek_noduals = MOIT.TestConfig(atol=1e-4, rtol=1e-4,
                                       optimal_status=MOI.OPTIMAL,
                                       query=false,
                                       duals=false,
                                       infeas_certificates=false,
                                       modify_lhs=false)

# Default configuration.
const config_knitro = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                               optimal_status=MOI.LOCALLY_SOLVED,
                               query=false,
                               infeas_certificates=false, # Do not ask for infeasibility certificates.
                               modify_lhs=false)

const config_knitro_noduals = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                                       optimal_status=MOI.LOCALLY_SOLVED,
                                       query=false,
                                       duals=false,
                                       infeas_certificates=false,
                                       modify_lhs=false)

config_solver = Dict(
"mosek"   => config_mosek,
"knitro" => config_knitro
)

config_solver_noduals = Dict(
"mosek"   => config_mosek_noduals,
"knitro" => config_knitro_noduals
)
