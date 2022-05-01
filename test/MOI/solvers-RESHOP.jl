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
global_solvers = Any[]

lp_solvers = ["path"]
lp_exclude = Dict()

quad_solvers = [("path", "local")]
quad_exclude = Dict( "path" => [])
mip_solvers = []

# Investigate if we can solve those
#soc_solvers = copy(quad_solvers)
soc_solvers = Any[]

push!(convex_nlp_solvers, "path")

append!(convex_nlp_solvers, nlp_solvers)

#push!(minlp_solvers, "")
