# MathOptInterface NLP

include("reshop_nlp_moi.jl")

function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    MOI.initialize(nlp_data.evaluator, [:ExprGraph])
    # Process constraints
    load_nlp_constraints(model, nlp_data)

    if nlp_data.has_objective
        load_nlp_objective(model, nlp_data)
    end


end

# Keep loading of NLP constraints apart to load all NLP model all in once
# inside Knitro.
function load_nlp_constraints(model::Optimizer, nlp_data::MOI.NLPBlockData)
    # This heavily relies on some convention
    model.start_nl_cons = Int(ctx_numequ(model.ctx))
    model.len_nl_cons = length(nlp_data.constraint_bounds)
    for i=1:length(nlp_data.constraint_bounds)
        c = MOI.constraint_expr(nlp_data.evaluator, i)

        ub = +Inf
        lb = -Inf

        # Remove relations and bounds from constraint expressions
        if length(c.args) == 3
            expected_head = :call
            expr_index = 2

            @assert c.head == expected_head
            # Single relation constraint: expr rel bound
            rel = c.args[1]
            if rel in [:<=, :(==)]
                ub = c.args[3]
            end
            if rel in [:>=, :(==)]
                lb = c.args[3]
            end
            c = c.args[expr_index]
            @assert isfinite(lb) || isfinite(ub)
        else
            # Double relation constraint: bound <= expr <= bound
            @assert c.head == :comparison
            error("doubly bounded NLP constraint is not yet supported")
        end

        # Convert non-linear expression to non-linear, linear and constant
        lin_part = Dict{Int32, Float64}()
        c, constant, conlinearities = process_expression!(c, lin_part)

        # Update bounds on constraint
        lb -= constant
        ub -= constant

        eidx = rhp_add_equ(model.ctx)

        tree, node = reshop_get_treedata(model.ctx, eidx)
        reshop_add_nlexpr(model.ctx, tree, node, c)
        lin_part = filter( x -> (abs(last(x)) > eps(0.)), lin_part)
        if length(lin_part) > 0
            avar = _ensure_avar(model)
            rhp_avar_set(avar, Vector{Int32}(keys(lin_part) .- 1))
            rhp_equ_add_linear(model.ctx, eidx, avar, collect(values(lin_part)))
        end
        reshop_set_equtype(model.ctx, eidx, relation_to_reshop[rel])
        isfinite(lb) && reshop_set_rhs(model.ctx, eidx, lb)
        isfinite(ub) && reshop_set_rhs(model.ctx, eidx, ub)
    end

 end

function load_nlp_objective(model::Optimizer, nlp_data::MOI.NLPBlockData)
# Process objective
    obj = MOI.objective_expr(nlp_data.evaluator)

    if length(obj.args) < 2
        obj = 0
    else
        # Convert non-linear expression to non-linear, linear and constant
        lin_part = Dict{Int32, Float64}()
        obj, constant, objlinearity = process_expression!(obj, lin_part)

        # Add constant back into non-linear expression
        if constant != 0
            obj = add_constant(obj, constant)
        end
    end

    eidx = rhp_add_equ(model.ctx)
    tree, node = reshop_get_treedata(model.ctx, eidx)
    reshop_add_nlexpr(model.ctx, tree, node, obj)
    lin_part = filter( x -> (abs(last(x)) > eps(0.)), lin_part)
    if length(lin_part) > 0
       avar = _ensure_avar(model)
       rhp_avar_set(avar, Vector{Int32}(keys(lin_part) .- 1))
       rhp_equ_add_linear(model.ctx, eidx, avar, collect(values(lin_part)))
    end
    rhp_set_objeqn(model.ctx, eidx)
end
