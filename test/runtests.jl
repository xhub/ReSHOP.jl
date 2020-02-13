using ReSHOP
using Pkg
using Compat.Test

pkgs = Pkg.installed()

@testset "Common basic tests" begin
	include("nl_convert.jl")
	include("nl_linearity.jl")
end

if pkgs["JuMP"] >= v"0.19"
	include("runtests_moi.jl")
	include("export_gms.jl")
	export_gms(mktempdir())
else
	include("runtests_mbp.jl")
end
