function export_gms(path)
    cd(path)
    model = ReSHOP.Optimizer()
    ReSHOP.setexport(true)

    MOI.empty!(model)

    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    t = MOI.add_variable(model)

    c1f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x,y]), 0.0)
    c1 = MOI.add_constraint(model, c1f, MOI.GreaterThan(1.0))

    c2f = MOI.ScalarQuadraticFunction(
        MOI.ScalarAffineTerm{Float64}[],
        MOI.ScalarQuadraticTerm.([2.0, 2.0, -2.0], [x, y, t], [x, y, t]),
        0.0
    )
    c2 = MOI.add_constraint(model, c2f, MOI.LessThan(0.0))

    bound = MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(0.0))

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, t)], 0.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(model)
end
