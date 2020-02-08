##############################################################################
# This files constains utilities to manipulate Julia Expr
##############################################################################



function process_expression!(nonlin_expr::Expr, lin_expr::Dict{Int32, Float64},
                             varlinearities::Vector{Symbol})
    # Get list of all variables in the expression
    extract_variables!(lin_expr, nonlin_expr)
    # Extract linear and constant terms from non-linear expression
    tree = LinearityExpr(nonlin_expr)
    tree = pull_up_constants(tree)
    _, tree, constant = prune_linear_terms!(tree, lin_expr)
    # Make sure all terms remaining in the tree are .nl-compatible
    nonlin_expr = convert_formula(tree)

    # Track which variables appear nonlinearly
    nonlin_vars = Dict{Int32, Float64}()
    extract_variables!(nonlin_vars, nonlin_expr)
    for j in keys(nonlin_vars)
        varlinearities[j] = :Nonlin
    end

    # Remove variables at coeff 0 that aren't also in the nonlinear tree
    for (j, coeff) in lin_expr
        if coeff == 0 && !(j in keys(nonlin_vars))
            delete!(lin_expr, j)
        end
    end

    # Mark constraint as nonlinear if anything is left in the tree
    linearity = nonlin_expr != 0 ? :Nonlin : :Lin

    return nonlin_expr, constant, linearity
end
function process_expression!(nonlin_expr::Real, lin_expr, varlinearities)
    # Special case where body of constraint is constant
    # Return empty nonlinear and linear parts, and use the body as the constant
    0, nonlin_expr, :Lin
end

# For MOI
function process_expression!(nonlin_expr::Expr, lin_expr::Dict{Int32, Float64})
    # Get list of all variables in the expression
    extract_variables!(lin_expr, nonlin_expr)
    # Extract linear and constant terms from non-linear expression
    tree = LinearityExpr(nonlin_expr)
    tree = pull_up_constants(tree)
    _, tree, constant = prune_linear_terms!(tree, lin_expr)
    # Make sure all terms remaining in the tree are .nl-compatible
    nonlin_expr = convert_formula(tree)

    # Track which variables appear nonlinearly
    nonlin_vars = Dict{Int32, Float64}()
    extract_variables!(nonlin_vars, nonlin_expr)

    # Remove variables at coeff 0 that aren't also in the nonlinear tree
    for (j, coeff) in lin_expr
        if abs(coeff) < eps(0.) && !(j in keys(nonlin_vars))
            delete!(lin_expr, j)
        end
    end

    # Mark constraint as nonlinear if anything is left in the tree
    linearity = nonlin_expr != 0 ? :Nonlin : :Lin

    return nonlin_expr, constant, linearity
end

function process_expression!(nonlin_expr::Real, lin_expr)
    # Special case where body of constraint is constant
    # Return empty nonlinear and linear parts, and use the body as the constant
    0, nonlin_expr, :Lin
end

# We need to track linear coeffs of all variables present in the expression tree
extract_variables!(lin_constr::Dict{Int32, Float64}, c) = c
extract_variables!(lin_constr::Dict{Int32, Float64}, c::LinearityExpr) =
    extract_variables!(lin_constr, c.c)
function extract_variables!(lin_constr::Dict{Int32, Float64}, c::Expr)
    if c.head == :ref
        CONFIG[:debug] && println("DEBUG: extract_variables :: variable case: $c")
        if c.args[1] == :x
            varidx = getvidx(c.args[2])
            lin_constr[varidx] = 0
        else
            error("Unrecognized reference expression $c")
        end
    else
        map(arg -> extract_variables!(lin_constr, arg), c.args)
    end
end

add_constant(c, constant::Real) = c + constant
add_constant(c::Expr, constant::Real) = Expr(:call, :+, c, constant)

substitute_vars!(c, x::Array{Float64}) = c
function substitute_vars!(c::Expr, x::Array{Float64})
    if c.head == :ref
        if c.args[1] == :x
            varidx = getvidx(c.args[2])
            c = x[varidx]
        else
            error("Unrecognized reference expression $c")
        end
    else
        if c.head == :call
            # Convert .nl unary minus (:neg) back to :-
            if c.args[1] == :neg
                c.args[1] = :-
            # Convert .nl :sum back to :+
            elseif c.args[1] == :sum
                c.args[1] = :+
            end
        end
        map!(arg -> substitute_vars!(arg, x), c.args, c.args)
    end
    c
end

function evaluate_linear(linear_coeffs::Dict{Int32, Float64}, x::Array{Float64})
    total = 0.0
    for (i, coeff) in linear_coeffs
        total += coeff * x[i]
    end
    total
end

function evaluate_quad(rowidx, colidx, qvals, x::Array{Float64})
   n = length(x)
   mat = sparse(rowidx, colidx, qvals, n, n)
   # This is soooo ugly --xhub
   Q = (mat + mat') - Diagonal(diag(mat))
   total = .5*x'*Q*x
   return total
end


