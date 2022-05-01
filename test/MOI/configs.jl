# Default configuration.
const config = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                               optimal_status=MOI.OPTIMAL,
                               query=false,
                               infeas_certificates=false, # Do not ask for infeasibility certificates.
                               modify_lhs=false)

const config_local = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                                     optimal_status=MOI.LOCALLY_SOLVED,
                                     query=false,
                                     duals=false, # well presolve give the solution, hence no duals ...
                                     infeas_certificates=false, # Do not ask for infeasibility certificates.
                                     modify_lhs=false)

const config_noduals = MOIT.TestConfig(atol=1e-5, rtol=1e-8,
                                       optimal_status=MOI.OPTIMAL,
                                       query=false,
                                       duals=false,
                                       infeas_certificates=false,
                                       modify_lhs=false)

const configs = Dict(
                     "local" => config_local,
                     "nodual" => config_noduals
                    )

optimizer_kw = Dict(:rtol => 1e-9)

const OPTIMIZER = ReSHOP.Optimizer(; optimizer_kw...)
