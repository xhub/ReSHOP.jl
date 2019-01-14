@testset "nl_write" begin
    # Turn on debug mode so files persist
    old_debug = ReSHOP.CONFIG[:debug]
    ReSHOP.setdebug(true)

#    filename = "test"
#    filepath = joinpath(ReSHOP.solverdata_dir, filename)
#    ReSHOP.clean_solverdata()

#    context("all temp files deleted successfully") do
#        @fact length(readdir(ReSHOP.solverdata_dir)) --> 1
#    end

    m = Model(solver=ReSHOP.ReSHOPSolver())
    @variable(m, x >= 0)
    @objective(m, Min, x)
    solve(m)

    # Reset debug mode and clean up
    ReSHOP.setdebug(old_debug)
    ReSHOP.clean_solverdata()
end
