using JuMP, Compat.Test, ReSHOP

## Solve test problem 1 (Synthesis of processing system) in
 #  M. Duran & I.E. Grossmann, "An outer approximation algorithm for
 #  a class of mixed integer nonlinear programs", Mathematical
 #  Programming 36, pp. 307-339, 1986.  The problem also appears as
 #  problem synthes1 in the MacMINLP test set.
 #
 #  min   5 x4 + 6 x5 + 8 x6 + 10 x1 - 7 x3 -18 math.log(x2 + 1)
 #       - 19.2 math.log(x1 - x2 + 1) + 10
 #  s.t.  0.8 math.log(x2 + 1) + 0.96 math.log(x1 - x2 + 1) - 0.8 x3 >= 0
 #        math.log(x2 + 1) + 1.2 math.log(x1 - x2 + 1) - x3 - 2 x6 >= -2
 #        x2 - x1 <= 0
 #        x2 - 2 x4 <= 0
 #        x1 - x2 - 2 x5 <= 0
 #        x4 + x5 <= 1
 #        0 <= x1 <= 2
 #        0 <= x2 <= 2
 #        0 <= x3 <= 1
 #        x1, x2, x3 continuous
 #        x4, x5, x6 binary
 #
 #
 #  The solution is (1.30098, 0, 1, 0, 1, 0).
 ##

if !isdef(:solver); solver = ReSHOP.ReSHOPSolver(); end

@testset "example: jump_minlp" begin
    m = Model(solver=solver)
    x_U = [2,2,1]
    @variable(m, x_U[i] >= x[i=1:3] >= 0)
    @variable(m, y[4:6], Bin)

    @NLobjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(x[1]-x[2]+1))
    @NLconstraint(m, 0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0)
    @NLconstraint(m, log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[6] >= -2)
    @NLconstraint(m, x[2] - x[1] <= 0)
    @NLconstraint(m, x[2] - 2*y[4] <= 0)
    @NLconstraint(m, x[1] - x[2] - 2*y[5] <= 0)
    @NLconstraint(m, y[4] + 2*y[5]*0.5 <= 1)

    @test solve(m) == :Optimal
    @test isapprox(getvalue(x)[:], [1.30098, 0.0, 1.0], atol=1e-5)
    @test isapprox(getvalue(y)[:], [0.0, 1.0, 0.0], atol=1e-5)
    @test isapprox(getobjectivevalue(m), 6.00975, atol=1e-5)
end

# this test fails in the Jump NL series with Julia 0.6 on Mac.
# Put it here to check that it can solve and that it is likely a memory corruption issue
@testset "Mixed integer nonlinear problems v2" begin
    # Solve test problem 1 (Synthesis of processing system) in
    # M. Duran & I.E. Grossmann, "An outer approximation algorithm for
    # a class of mixed integer nonlinear programs", Mathematical
    # Programming 36, pp. 307-339, 1986.  The problem also appears as
    # problem synthes1 in the MacMINLP test set.
    m = Model(solver=solver)
    x_U = [2,2,1]
    @variable(m, x_U[i] >= x[i=1:3] >= 0)
    @variable(m, y[4:6], Bin)
    @NLobjective(m, Min, 10 + 10*x[1] - 7*x[3] + 5*y[4] + 6*y[5] + 8*y[6] - 18*log(x[2]+1) - 19.2*log(x[1]-x[2]+1))
    @NLconstraints(m, begin
        0.8*log(x[2] + 1) + 0.96*log(x[1] - x[2] + 1) - 0.8*x[3] >= 0
        log(x[2] + 1) + 1.2*log(x[1] - x[2] + 1) - x[3] - 2*y[6] >= -2
        x[2] - x[1] <= 0
        x[2] - 2*y[4] <= 0
        x[1] - x[2] - 2*y[5] <= 0
        y[4] + y[5] <= 1
    end)
    status = solve(m)

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), 6.00976, atol=1e-5)
    @test isapprox(getvalue(x), [1.30098, 0.0, 1.0], atol=1e-5)
    @test isapprox(getvalue(y)[:], [0.0, 1.0, 0.0], atol=1e-5)
end


