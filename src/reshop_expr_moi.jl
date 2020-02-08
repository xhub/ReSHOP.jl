##############################################################################
# This files constains utilities to manipulate Julia Expr
##############################################################################



# We need to track linear coeffs of all variables present in the expression tree
extract_variables!(lin_constr::Dict{Int, Float64}, c) = c
extract_variables!(lin_constr::Dict{Int, Float64}, c::LinearityExpr) =
    extract_variables!(lin_constr, c.c)
function extract_variables!(lin_constr::Dict{Int, Float64}, c::Expr)
    if c.head == :ref
        CONFIG[:debug] && println("DEBUG: extract_variables :: variable case: $c")
        if c.args[1] == :x
            @assert isa(c.args[2], Int)
            lin_constr[c.args[2]] = 0
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
            index = c.args[2]
            @assert isa(index, Int)
            c = x[index]
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

function evaluate_linear(linear_coeffs::Dict{Int, Float64}, x::Array{Float64})
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


