##############################################################################
# The code is adapted from AMPLNlWriter 
# Many thanks for all contributors there
##############################################################################

##############################################################################
# Low level functions
##############################################################################

function write_arithm_op(ctx, tree, node, fn, args)
    # TODO(xhub) implement better management of variable/constant base
    len = length(args)
    equtree_arithm(tree, node, arithm_ops[fn], len)
    pnode = equnode_deref(node)
    for i in 0:(len-1)
        child = equnode_get_child_addr(pnode, i)
        reshop_add_nlexpr(ctx, tree, child, args[len-i])
    end
end

function write_power_expr(c, ctx, tree, node)
    len = length(c.args)
    @assert len == 3
    if isa(c.args[2], AbstractFloat) && abs(round(c.args[2]) - c.args[2]) > eps(c.args[2])
        opcode = 74 # fncvpower
    elseif isa(c.args[3], AbstractFloat) && abs(round(c.args[3]) - c.args[3]) > eps(c.args[3])
        opcode = 75 # fnvcpower
    else
        opcode = 21
    end

    equtree_call(ctx, tree, node, (opcode, 2))
    pnode = equnode_deref(node)

    for i in 0:(len - 2)
        child = equnode_get_child_addr(pnode, i)
        reshop_add_nlexpr(ctx, tree, child, c.args[i+2])
    end
end

    # Convert an expression tree
reshop_add_nlexpr(ctx, tree, node, c) = error("Unrecognized expression $c")
# Handle numerical constants e.g. pi
reshop_add_nlexpr(ctx, tree, node, c::Symbol) =  reshop_add_nlexpr(ctx, tree, node, float(eval(c)))

# write down constant
function reshop_add_nlexpr(ctx, tree, node, c::Real)
#    if abs(c) > 0.
        equtree_cst(ctx, tree, node, c)
#    end
end

reshop_add_nlexpr(ctx, tree, node, c::LinearityExpr) = reshop_add_nlexpr(ctx, tree, node, c.c)
function reshop_add_nlexpr(ctx, tree, node, c::Expr)
    CONFIG[:debug] && println("DEBUG: encoding expr $c")
    if c.head == :ref
        CONFIG[:debug] && println("DEBUG: variable case: $c")
        # This is a variable
        if c.args[1] == :x
            @assert isa(c.args[2], MOI.VariableIndex)
            equtree_var(ctx, tree, node, c.args[2].value-1, 1.)
        else
            error("Unrecognized reference expression $c")
        end
    elseif c.head == :call
        if c.args[1] in keys(arithm_ops)
            write_arithm_op(ctx, tree, node, c.args[1], c.args[2:end])
        elseif c.args[1] == :neg
            equtree_umin(ctx, tree, node)
            reshop_add_nlexpr(ctx, tree, node, c.args[2])
        # TODO(xhub) support nary_functions
        elseif c.args[1] in nary_functions
            error("Unsupported n-ary function $(c.args[1])")
        else
            if (c.args[1] == :^) return write_power_expr(c, ctx, tree, node) end

            equtree_call(ctx, tree, node, func_to_reshop[c.args[1]])
            len = length(c.args)
            pnode = equnode_deref(node)

            for i in 0:(len - 2)
                child = equnode_get_child_addr(pnode, i)
                reshop_add_nlexpr(ctx, tree, child, c.args[i+2])
            end
        end

    elseif c.head == :comparison
        # .nl only handles binary comparison
        @assert length(c.args) == 3
        # Output comparison type first, followed by args
        CONFIG[:debug] && println(f, nl_operator(c.args[2]))
        map(arg -> reshop_add_nlexpr(ctx, tree, node, arg), c.args[1:2:end])

    elseif c.head in [:&&, :||]
        error("unsupported binary condition")
    else
        error("Unrecognized expression $c")
    end
end
