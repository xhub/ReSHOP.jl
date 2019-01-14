function reshop_reg_eqns(ctx, m::ReSHOPMathProgModel, equil=false)
    if has_objective(m)
        obj_offset = 1
    else
        obj_offset = 0
    end

    reshop_declare_eqns(ctx, m, equil)

    if has_objective(m)
        reshop_add_obj_nl(ctx, m)
        reshop_add_obj_lin(ctx, m)
        reshop_add_obj_quad(ctx, m)
    end

    if m.ncon > 0 || length(m.quad_equs) > 0
#        write_nl_k_block(f, m)
        if length(m.nonquad_idx) != m.ncon
            error("The number of constraint and the equation index do not match")
        end
        reshop_add_cons_lin(ctx, m, obj_offset)
        reshop_add_cons_nl(ctx, m, obj_offset)
        reshop_add_cons_quad(ctx, m, obj_offset)
        reshop_add_contraint_sense(ctx, m, obj_offset)
    end

    # This is most likely broken, we need another way of getting this data if
    # it is important
#    return m.ncon+obj_offset+length(m.quad_equs)
end


function create_reshop_ctx(m::ReSHOPMathProgModel)
    # Not all model have an objective function
    if has_objective(m)
        obj_offset = 1
    else
        obj_offset = 0
    end

    ctx = ctx_create(m.nvar, length(m.quad_equs)+m.ncon+obj_offset)

    reshop_declare_vars(ctx, m)
    reshop_add_var_guess(ctx, m)

    reshop_reg_eqns(ctx, m)

    return ctx
end

function has_objective(m::ReSHOPMathProgModel)
    CONFIG[:debug] && println("DEBUG: has_objective $(m.lin_obj) $(m.obj) $(m.quad_obj)")
    res = !isempty(m.lin_obj) || (isa(m.obj, Expr)) || length(m.quad_obj) == 3
    if !res return res end
    res = any(v != 0. for v in values(m.lin_obj)) || isa(m.obj, Expr) && (m.obj != 0.) || length(m.quad_obj) == 3
    return res
end

# Nonlinear constraint trees
function reshop_add_cons_nl(ctx, m::ReSHOPMathProgModel, offset)
    for idx in 1:m.ncon
        if isa(m.constrs[idx], Expr) && (m.constrs[idx] != 0.)
            eidx = m.nonquad_idx[idx] + offset - 1 + m.offset
            tree, node = reshop_get_treedata(ctx, eidx)
            CONFIG[:debug] && println("reshop_add_cons_nl: storing expression $(eidx)")
            reshop_add_nlexpr(ctx, tree, node, m, m.constrs[idx])
        end
    end
end

function reshop_quad(ctx, m::ReSHOPMathProgModel, idx, equ, offset, isObj::Bool=false)
    ##########################################################################
    #
    ##########################################################################
    if isObj
        eidx = m.offset
    else
        eidx = m.quad_idx[idx] - 1 + m.offset + offset
    end

    CONFIG[:debug] && println("reshop_quad: $equ")
    lidx, lval, rowidx, colidx, qval = equ[1:5]
    CONFIG[:debug] && println("reshop_quad: $lidx $lval $rowidx $colidx $qval")

    ##########################################################################
    # lin_dict collects all the linear terms
    # 
    # we don't construct the dict directly, since it may have duplicate values
    ##########################################################################

    lin_dict = Dict(zip(lidx, Iterators.repeated(0.)))
    for i in 1:length(lidx)
        lin_dict[lidx[i]] += lval[i]
    end

    ###########################################################################
    # MPB is a bit insane here: via setquadobj! it provides a triangular matrix
    # representing .5 x^TQx ... But this makes no sense since we assume that
    # the multiplication is commutative here. Therefore the values of the off-
    # diagonal element have to be multiplied by 2
    ###########################################################################

    quad_dict = Dict(zip(zip(rowidx, colidx), Iterators.repeated(0.)))

    if isObj
        for i in 1:length(rowidx)
            if rowidx[i] != colidx[i]
                quad_dict[(rowidx[i], colidx[i])] += 2*qval[i]
            else
                quad_dict[(rowidx[i], colidx[i])] += qval[i]
            end
        end
    else
        for i in 1:length(rowidx)
            quad_dict[(rowidx[i], colidx[i])] += qval[i]
        end
    end

    qidxC = keys(quad_dict)
    rowidxU = [elt[1]-1 for elt in qidxC]
    colidxU = [elt[2]-1 for elt in qidxC]
    qvalU = collect(values(quad_dict))


    ##########################################################################
    # Add the linear part <c,x> for the purely linear variables
    ##########################################################################

    lidxS = Set(lidx)
    lin_vars = setdiff(lidxS, union(Set(rowidx), Set(colidx)))
    qidxS = setdiff(lidxS, lin_vars)
    for lindx in lin_vars
        # add as linear variable
        v = lin_dict[lindx]
        if abs(v) > 0.
            CONFIG[:debug] && println("DEBUG: in reshop_quad, lin var $lindx with value $v")
            ctx_add_lin_var(ctx, eidx, m.v_index_map[lindx], v)
        end
    end

    ##########################################################################
    # Add the term <c,x> for the quadratic variables
    ##########################################################################

    if length(lin_dict) > 0 && length(qidxS) > 0
        # TODO(xhub) filter for qvals[i] = 0.
        qvals = collect(lin_dict[qidx] for qidx in qidxS)
        CONFIG[:debug] && println("DEBUG: in reshop_quad, qidxS = $qidxS; qval= $qvals")
        abs_var = reshop_avar(length(qidxS), collect(qidxS) .- 1)
        CONFIG[:debug] && println("DEBUG: in reshop_quad, quad var $(collect(qidxS)) with value $qvals")
        reshop_equ_add_lin_tree(ctx, eidx, qvals, abs_var, 1.)
        reshop_avar_free(abs_var)
    end

    ##########################################################################
    # Add the quadratic terms <x, Mx>
    ##########################################################################

    if length(qvalU) > 0
        CONFIG[:debug] && println("DEBUG: in reshop_quad, quad term $rowidxU $colidxU $qvalU")
        @assert length(qvalU) == length(rowidxU)
        @assert length(qvalU) == length(colidxU)
        mat = reshop_mat_coo(rowidxU, colidxU, qvalU)
        midxS = union(BitSet(rowidx), BitSet(colidx))
        avar = reshop_avar(length(midxS), collect(i-1 for i in midxS))
        if isObj
            c = 1.
        else
            c = 2.
        end
        # the coeff is 2 since this function adds .5 x^TMx
        reshop_equ_add_quadratic(ctx, eidx, mat, avar, c)
        reshop_avar_free(avar)
        reshop_mat_free(mat)
    end
end

function reshop_add_obj_quad(ctx, m::ReSHOPMathProgModel)
    CONFIG[:debug] && println("DEBUG: quad_obj = $(m.quad_obj)")
    if length(m.quad_obj) == 3
        rowidx = m.quad_obj[1]
        colidx = m.quad_obj[2]
        quadval = m.quad_obj[3]
        reshop_quad(ctx, m, 0, (Int[], Float64[], rowidx, colidx, quadval), 0, true)
    elseif length(m.quad_obj) == 0
        #doing noting here
    else
        error("reshop_add_obj_quad :: invalid quad_obj object $(m.quad_obj)")
    end
end

function reshop_add_cons_quad(ctx, m::ReSHOPMathProgModel, offset)
    for (idx, equ) in enumerate(m.quad_equs)
        reshop_quad(ctx, m, idx, equ, offset)
    end
end

# Nonlinear objective tree
function reshop_add_obj_nl(ctx, m::ReSHOPMathProgModel)
    # Get tree, node, ...
    if m.obj != 0.
        tree, node = reshop_get_treedata(ctx, m.offset)
        CONFIG[:debug] && println("reshop_add_obj_nl: storing expression $(m.offset)")
        reshop_add_nlexpr(ctx, tree, node, m, m.obj)
    end
end

# Initial primal guesses
function reshop_add_var_guess(ctx, m::ReSHOPMathProgModel)
    for idx in 0:(m.nvar - 1)
        i = m.v_index_map_rev[idx]
        CONFIG[:debug] && println("Setting level value variable $i to $(m.x_0[i])")
        ctx_setvarlone(ctx, idx, m.x_0[i])
    end
end

# Constraint bounds
function reshop_add_contraint_sense(ctx, m::ReSHOPMathProgModel, offset)
    for idx in 1:m.ncon
        lower = m.g_l[idx]
        upper = m.g_u[idx]
        rel = m.r_codes[idx]
        # TODO(xhub) document
        if rel == -1
            error("Doubly constraint equation is not yet supported")
        elseif rel == 4         # ==
            value = lower
        elseif rel == 1         # >=
            value = lower
        elseif rel == 2         # <=
            value = upper
        elseif rel == 0         # ???
            error("unsupported rel = $(rel) for equation index $idx")
            value == upper
        else
            error("unsupported rel = $(rel) for equation index $idx")
        end
        eidx = m.nonquad_idx[idx] + offset - 1 + m.offset
        CONFIG[:debug] && println("Setting sense and rhs for equation $eidx: $rel $value")
        reshop_set_rhs(ctx, eidx, value)
        reshop_set_equtype(ctx, eidx, rel)
    end

    for (idx, equ) in enumerate(m.quad_equs)
        rel, value = equ[end-1:end]
        eidx = m.quad_idx[idx] + offset - 1 + m.offset
        reshop_set_rhs(ctx, eidx, value)
        reshop_set_equtype(ctx, eidx, relation_to_reshop[rel])
    end
end

function reshop_declare_var(ctx, vtype, lower, upper)
    if lower == -Inf
        if upper == Inf
            reshop_add_free_var(ctx, 1)
        elseif upper == 0
            reshop_add_neg_var(ctx, 1)
        else
            reshop_add_box_var(ctx, lower, upper)
        end
    else
        if lower == upper
            reshop_add_box_var(ctx, lower, upper)
        elseif upper == Inf
            if lower == 0
                reshop_add_pos_var(ctx, 1)
            else
                reshop_add_box_var(ctx, lower, upper)
            end
        else
            reshop_add_box_var(ctx, lower, upper)
        end
    end
    idx = hack_last_vidx(ctx)
    CONFIG[:debug] && println("DEBUG: var $(idx) has type $(vtype)")
    # We have no consistency check
    # TODO SOS, semicont
    if vtype == :Bin
        reshop_set_vartype(ctx, idx, 1)
    elseif vtype == :Int
        reshop_set_vartype(ctx, idx, 2)
    end

end

# Variable bounds
function reshop_declare_vars(ctx, m::ReSHOPMathProgModel)
    for idx in 0:(m.nvar - 1)
        i = m.v_index_map_rev[idx]
        lower = m.x_l[i]
        upper = m.x_u[i]
        reshop_declare_var(ctx, m.vartypes[i], lower, upper)
    end

    # Set variable names
    if m.d.hasvalue
        ctx_setvarnames(ctx, m.d.value.m.colNames)
    end
end

# Jacobian counts
# TODO(xhub) resuse this to prealloc the model_repr?
#function write_nl_k_block(f, m::ReSHOPMathProgModel)
#    println(f, "k$(m.nvar - 1)")
#    total = 0
#    for index = 0:(m.nvar - 2)
#        i = m.v_index_map_rev[index]
#        total += m.j_counts[i]
#        println(f, total)
#    end
#end

function reshop_declare_eqns(ctx, m::ReSHOPMathProgModel, equil=false)
    if has_objective(m)
        cor = 0
    else
        cor = -1
    end

    for idx in 0:m.ncon+cor+length(m.quad_equs)
        reshop_decl_eqn(ctx, idx+m.offset)
    end
    if equil || !has_objective(m)
        reshop_set_objeqn(ctx, -1)
    else
        reshop_set_objeqn(ctx, 0)
    end
#    return m.ncon+cor+length(m.quad_equs)
end

# Linear constraint expressions
function reshop_add_cons_lin(ctx, m::ReSHOPMathProgModel, offset)
    for idx in 1:m.ncon
        num_vars = length(m.lin_constrs[idx])
        if num_vars > 0
#TODO(Xhub)            reshop_reserve_lequ(ctx, idx, num_vars)
            eidx = m.nonquad_idx[idx] + offset - 1 + m.offset
            for (k,v) in m.lin_constrs[idx]
                if abs(v) > 0.
                    CONFIG[:debug] && println("DEBUG: add var $k with value $v in eqn $(eidx)")
                    ctx_add_lin_var(ctx, eidx, m.v_index_map[k], v)
                end
            end
        end
    end
end

# Linear objective expression
function reshop_add_obj_lin(ctx, m::ReSHOPMathProgModel)
    for (k,v) in m.lin_obj
        if abs(v) > 0.
            ctx_add_lin_var(ctx, m.offset, m.v_index_map[k], v)
        end
    end
end

function reshop_add_quad_lin(ctx, m::ReSHOPMathProgModel)
end


##############################################################################
# Low level functions
##############################################################################

function write_arithm_op(ctx, tree, node, m, fn, args)
    # TODO(xhub) implement better management of variable/constant base
    len = length(args)
    equtree_arithm(tree, node, arithm_ops[fn], len)
    pnode = equnode_deref(node)
    for i in 0:(len-1)
        child = equnode_get_child_addr(pnode, i)
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
            pnode = equnode_deref(node)

            for i in 0:(len - 2)
                child = equnode_get_child_addr(pnode, i)
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
