##############################################################################
# The code is adapted from AMPLNlWriter 
# Many thanks for all contributors there
##############################################################################

##############################################################################
# Low level functions
##############################################################################

function write_arithm_op(ctx, tree, node, m, fn, args)
    # TODO(xhub) implement better management of variable/constant base
    len = length(args)
    equtree_arithm(tree, node, arithm_ops[fn], len)

    for i in 0:(len-1)
        child = equnode_get_child_addr(node, i)
        reshop_add_nlexpr(ctx, tree, child, m, args[len-i])
    end
end


# Convert an expression tree
reshop_add_nlexpr(ctx, tree, node, m, c) = error("Unrecognized expression $c")
# Handle numerical constants e.g. pi
reshop_add_nlexpr(ctx, tree, node, m, c::Symbol) =  reshop_add_nlexpr(ctx, tree, node, m, float(eval(c)))

# write down constant
function reshop_add_nlexpr(ctx, tree, node, m, c::Real)
#    if abs(c) > 0.
        equtree_cst(ctx, tree, node, c)
#    end
end

reshop_add_nlexpr(ctx, tree, node, m, c::LinearityExpr) = reshop_add_nlexpr(ctx, tree, node, m, c.c)
function reshop_add_nlexpr(ctx, tree, node, m, c::Expr)
    CONFIG[:debug] && println("DEBUG: encoding expr $c")
    if c.head == :ref
        CONFIG[:debug] && println("DEBUG: variable case: $c")
        # This is a variable
        if c.args[1] == :x
            @assert isa(c.args[2], Int)
            equtree_var(ctx, tree, node, m.v_index_map[c.args[2]], 1.)
        else
            error("Unrecognized reference expression $c")
        end
    elseif c.head == :call
        if c.args[1] in keys(arithm_ops)
            write_arithm_op(ctx, tree, node, m, c.args[1], c.args[2:end])
        elseif c.args[1] == :neg
            equtree_umin(ctx, tree, node)
            reshop_add_nlexpr(ctx, tree, node, m, c.args[2])
        # TODO(xhub) support nary_functions
        elseif c.args[1] in nary_functions
            error("Unsupported n-ary function $(c.args[1])")
        else
            equtree_call(ctx, tree, node, func_to_reshop[c.args[1]])
            len = length(c.args)

            for i in 0:(len - 2)
                child = equnode_get_child_addr(node, i)
                reshop_add_nlexpr(ctx, tree, child, m, c.args[i+2])
            end
        end

    elseif c.head == :comparison
        # .nl only handles binary comparison
        @assert length(c.args) == 3
        # Output comparison type first, followed by args
        CONFIG[:debug] && println(f, nl_operator(c.args[2]))
        map(arg -> reshop_add_nlexpr(ctx, tree, node, m, arg), c.args[1:2:end])

    elseif c.head in [:&&, :||]
        error("unsupported binary condition")
    else
        error("Unrecognized expression $c")
    end
end
