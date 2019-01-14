using JuMP, Compat.Test, ReSHOP

## Solve test problem with non-linear binary variables
 #
 #  min   100 * (x2 - (0.5 + x1) ^ 2) ^ 2 + (1 - x1) ^ 2
 #  s.t.  x1, x2 binary
 #
 #  The solution is (0, 0).
 ##

if !isdef(:solver_sbb); solver_sbb = ReSHOP.ReSHOPSolver("sbb"); end

@testset "example: jump_nonlinearbinary" begin
    m = Model(solver=solver_sbb)
    @variable(m, x[1:2], Bin)

    # Set some non-binary bounds on x1 and x2. These should be ignored.
    # The optimal solution if x is Int is (1, 2) which is allowed by these bounds
    setupperbound(x[1], 2)
    setupperbound(x[2], 2)

    @NLobjective(m, Min, 100*(x[2] - (0.5 + x[1])^2)^2 + (1 - x[1])^2)

    @test solve(m) == :Optimal
    @test getvalue(x)[:] == [0.0, 0.0]
    @test getobjectivevalue(m) == 7.25
end
