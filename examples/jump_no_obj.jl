using JuMP, Compat.Test, ReSHOP

# Example with no objective (#50)

if !isdef(:solver);global solver; solver = ReSHOP.ReSHOPSolver(); end

@testset "example: jump_no_obj" begin
    m = Model(solver=solver)
    @variable(m, 0 <= yp <= 1, Int)
    @variable(m, 0 <= l <= 1000.0)
    @variable(m, 0 <= f <= 1000.0)
    @NLconstraint(m, .087 * l >= f ^ 2)
    @constraint(m, l <= yp * 1000.0)

    @test solve(m) == :Optimal
    @test getobjectivevalue(m) == 0
end
