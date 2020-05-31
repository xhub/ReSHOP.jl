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

lp_solvers = ["baron" "knitro"]

lp_exclude = Dict("knitro" => ["linear8c", "linear8b"])

quad_solvers = [("baron", "nodual"),
                ("knitro", ""),
               ]
quad_exclude = Dict("baron" => ["qcp1"], # SLOW_PROGRESS
                    )
# mosek fails?
global_solvers = [("knitro", "")]

# fails with "'gmsbaxnx.exe' is not recognized as an internal or external command, operable program or batch file"
if !Sys.iswindows()
    push!(global_solvers, ("baron", ""))
end


push!(nlp_solvers, "baron")
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
"knitro" => config_knitro
)

config_solver_noduals = Dict(
"knitro" => config_knitro_noduals
)
