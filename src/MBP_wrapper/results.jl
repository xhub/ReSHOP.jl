function report_results(m::ReSHOPMathProgBaseModel)
    # TODO(Xhub) fix this hack
    reshop_postprocess(m.reshop_mdl_solver)

    # Next, read for the variable values
    report_results_common(m)

    # TODO(xhub) this should not be necessary
    if m.solve_exitcode == 0
        if m.objlinearity == :Nonlin
            # Try to use NLPEvaluator if we can.
            # Can fail due to unsupported functions so fallback to eval
            try
                m.objval = eval_f(m.d, m.solution)
            catch
                CONFIG[:debug] && println("Error: could not evaluate the objective function")
            end
        end

        # Calculate objective value from nonlinear and linear parts
        obj_nonlin = eval(substitute_vars!(deepcopy(m.obj), m.solution))
        obj_lin = evaluate_linear(m.lin_obj, m.solution)
        if (length(m.quad_obj) == 3)
           ridx, cidx, vals = m.quad_obj
           obj_quad = evaluate_quad(ridx, cidx, vals, m.solution)
        else
           obj_quad = 0.
        end
        m.objval = obj_nonlin + obj_lin + obj_quad
    end
end


