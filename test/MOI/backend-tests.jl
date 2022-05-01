@testset "MOI utils" begin
    @testset "SolverName" begin
        optimizer = ReSHOP.Optimizer()
        @test MOI.get(optimizer, MOI.SolverName()) == "ReSHOP"
    end
    @testset "supports_default_copy_to" begin
        optimizer = ReSHOP.Optimizer()
        @test MOIU.supports_default_copy_to(optimizer, false)
        # Use `@test !...` if names are not supported
        @test MOIU.supports_default_copy_to(optimizer, true)
    end
    @testset "MOI.Silent" begin
        optimizer = ReSHOP.Optimizer()
        @test MOI.supports(optimizer, MOI.Silent())
        MOI.set(optimizer, MOI.Silent(), true)
        @test MOI.get(optimizer, MOI.Silent()) == true
        MOI.set(optimizer, MOI.Silent(), false)
        @test MOI.get(optimizer, MOI.Silent()) == false
    end
#    @testset "MOI.TimeLimitSec" begin
#        optimizer = ReSHOP.Optimizer()
#        @test MOI.supports(optimizer, MOI.TimeLimitSec())
#        # TimeLimitSec is set to 1e8 by default in Knitro.
#        @test MOI.get(optimizer, MOI.TimeLimitSec()) == 1e8
#        my_time_limit = 10.
#        MOI.set(optimizer, MOI.TimeLimitSec(), my_time_limit)
#        @test MOI.get(optimizer, MOI.TimeLimitSec()) == my_time_limit
#    end
end


