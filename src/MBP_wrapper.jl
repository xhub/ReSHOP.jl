##############################################################################
# Some utilities
##############################################################################

function add_to_index_maps!(forward_map::Dict{Int, Int},
                            backward_map::Dict{Int, Int},
                            inds::Array{Int},
                            offset::Int)
    for i in inds
        index = length(forward_map) + offset
        forward_map[i] = index
        backward_map[index] = i
    end
end

function add_to_index_maps!(forward_map::Dict{Int, Int},
                            inds::Array{Int},
                            offset::Int)
    for i in inds
        index = length(forward_map) + offset
        forward_map[i] = index
    end
end

##############################################################################
# The main structure definition
##############################################################################

mutable struct ReSHOPMathProgBaseModel <: MPB.AbstractMathProgModel
    options::Dict{String, Any}

    solver_name::String

    x_l::Vector{Float64}
    x_u::Vector{Float64}
    g_l::Vector{Float64}
    g_u::Vector{Float64}

    nvar::Int
    ncon::Int

    obj
    constrs::Vector{Any}

    lin_constrs::Vector{Dict{Int32, Float64}}
    lin_obj::Dict{Int32, Float64}

    r_codes::Vector{Int}
    j_counts::Vector{Int}

    vartypes::Vector{Symbol}
    varlinearities_con::Vector{Symbol}
    varlinearities_obj::Vector{Symbol}
    conlinearities::Vector{Symbol}
    objlinearity::Symbol

    v_index_map::Dict{Int, Int}
    v_index_map_rev::Dict{Int, Int}

    # was:
    # c_index_map::Dict{Int, Int}
    # c_index_map_rev::Dict{Int, Int}

    nonquad_idx::Dict{Int, Int}
    quad_idx::Dict{Int, Int}

    sense::Symbol

    x_0::Vector{Float64}

    objval::Float64
    solution::Vector{Float64}

    status::Symbol
    solve_exitcode::Int
    solve_result_num::Int
    solve_result::String
    model_result_num::Int
    model_result::String
    solve_message::String
    solve_time::Float64

    model_type::MODEL_TYPE

    quad_equs::Vector{Any}
#    quad_obj::Tuple{Vector{Int}, Vector{Int}, Vector{Float64}}
    quad_obj::Tuple
    offset::Int
    emp::Nullable{Function}


    d::Nullable{MPB.AbstractNLPEvaluator}

    reshop_ctx::Ptr{context}
    reshop_ctx_dest::Ptr{context}
    reshop_mdl::Ptr{reshop_model}
    reshop_mdl_solver::Ptr{reshop_model}
    reshop_options::Ptr{reshop_options}
    gams_dir::String

    function ReSHOPMathProgBaseModel(solver_name::String,
                                options::Dict{String,Any},
                                model_type::MODEL_TYPE,
                                emp)
        o = new(options,
            solver_name,
            zeros(0),
            zeros(0),
            zeros(0),
            zeros(0),
            0,
            0,
            :(0),
            [],
            Dict{Int, Float64}[],
            Dict{Int, Float64}(),
            Int[],
            Int[],
            Symbol[],
            Symbol[],
            Symbol[],
            Symbol[],
            :Lin,
            Dict{Int, Int}(),
            Dict{Int, Int}(),
            Dict{Int, Int}(),
            Dict{Int, Int}(),
            :Min,
            zeros(0),
            NaN,
            zeros(0),
            :NotSolved,
            -1,
            -1,
            "?",
            -1,
            "?",
            "",
            NaN,
            model_type,
            Vector{Any}(),
            (),
            0,
            emp,
            Nullable{MPB.AbstractNLPEvaluator}(),
            Ptr{context}(C_NULL),
            Ptr{context}(C_NULL),
            Ptr{reshop_model}(C_NULL),
            Ptr{reshop_model}(C_NULL),
            Ptr{reshop_options}(C_NULL),
            "")
        @compat finalizer(reshop_cleanup, o)
        o
    end
end

include(joinpath("MBP_wrapper", "utils.jl"))

struct ReSHOPLinearQuadraticModel <: MPB.AbstractLinearQuadraticModel
    inner::ReSHOPMathProgBaseModel
end

struct ReSHOPNonlinearModel <: MPB.AbstractNonlinearModel
    inner::ReSHOPMathProgBaseModel
end

struct ReSHOPConicModel <: MPB.AbstractNonlinearModel
    inner::ReSHOPMathProgBaseModel
end

MPB.NonlinearModel(s::ReSHOPSolver) = ReSHOPNonlinearModel(
    ReSHOPMathProgBaseModel(s.solver_name, s.options, nlp, s.emp)
)

MPB.LinearQuadraticModel(s::ReSHOPSolver) = ReSHOPLinearQuadraticModel(
    ReSHOPMathProgBaseModel(s.solver_name, s.options, qcp, s.emp)
)

function MPB.ConicModel(s::ReSHOPSolver)
    error("ConicModel is not yet supported")
end

function MPB.loadproblem!(outer::ReSHOPNonlinearModel, nvar::Integer, ncon::Integer,
                      x_l, x_u, g_l, g_u, sense::Symbol,
                      d::MPB.AbstractNLPEvaluator)
    m = outer.inner

    m.nvar, m.ncon = nvar, ncon
    loadcommon!(m, x_l, x_u, g_l, g_u, sense)

    m.d = d

    MPB.initialize(m.d.value, [:ExprGraph])

    # Process constraints
    m.constrs = map(1:m.ncon) do i
        c = MPB.constr_expr(m.d.value, i)

        # Remove relations and bounds from constraint expressions
        if length(c.args) == 3
            if VERSION < v"0.5-"
                expected_head = :comparison
                expr_index = 1
                rel_index = 2
            else
                expected_head = :call
                expr_index = 2
                rel_index = 1
            end

            @assert c.head == expected_head
            # Single relation constraint: expr rel bound
            rel = c.args[rel_index]
            m.r_codes[i] = relation_to_reshop[rel]
            if rel == [:<=, :(==)]
                m.g_u[i] = c.args[3]
            end
            if rel in [:>=, :(==)]
                m.g_l[i] = c.args[3]
            end
            c = c.args[expr_index]
        else
            # Double relation constraint: bound <= expr <= bound
            @assert c.head == :comparison
            m.r_codes[i] = relation_to_reshop[:multiple]
            m.g_u[i] = c.args[5]
            m.g_l[i] = c.args[1]
            c = c.args[3]
        end

        # Convert non-linear expression to non-linear, linear and constant
        c, constant, m.conlinearities[i] = process_expression!(
            c, m.lin_constrs[i], m.varlinearities_con)

        # Update bounds on constraint
        m.g_l[i] -= constant
        m.g_u[i] -= constant

        # Update jacobian counts using the linear constraint variables
        for j in keys(m.lin_constrs[i])
            m.j_counts[j] += 1
        end
        c
    end

    # Process objective
    m.obj = MPB.obj_expr(m.d.value)

    if length(m.obj.args) < 2
        m.obj = 0
    else
        # Convert non-linear expression to non-linear, linear and constant
        m.obj, constant, m.objlinearity = process_expression!(
            m.obj, m.lin_obj, m.varlinearities_obj)

        # Add constant back into non-linear expression
        if constant != 0
            m.obj = add_constant(m.obj, constant)
        end
    end
    m
end

function MPB.loadproblem!(outer::ReSHOPLinearQuadraticModel, A::AbstractMatrix,
                      x_l, x_u, c, g_l, g_u, sense)
    m = outer.inner
    m.ncon, m.nvar = size(A)

    loadcommon!(m, x_l, x_u, g_l, g_u, sense)

    # Load A into the linear constraints
    @assert (m.ncon, m.nvar) == size(A)
    load_A!(m, A)
    m.constrs = zeros(m.ncon)  # Dummy constraint expression trees

    # Load c
    for (index, val) in enumerate(c)
        m.lin_obj[index] = val
    end
    # TODO(xhub) see if we can get rid of that
    m.obj = 0  # Dummy objective expression tree

    # Process variables bounds
    for j = 1:m.ncon
        lower = m.g_l[j]
        upper = m.g_u[j]
        if lower == -Inf
            if upper == Inf
                error("Neither lower nor upper bound on constraint $j")
            else # <=
                m.r_codes[j] = 2
            end
        else
            if lower == upper  # ==
                m.r_codes[j] = 4
            elseif upper == Inf # >=
                m.r_codes[j] = 1
            else # lb <= expr <= ub
                m.r_codes[j] = -1
            end
        end
    end
    m
end

function load_A!(m::ReSHOPMathProgBaseModel, A::SparseMatrixCSC{Float64})
    for var = 1:A.n, k = A.colptr[var] : (A.colptr[var + 1] - 1)
        m.lin_constrs[A.rowval[k]][var] = A.nzval[k]
        m.j_counts[var] += 1
    end
end

function load_A!(m::ReSHOPMathProgBaseModel, A::Matrix{Float64})
    for con = 1:m.ncon, var = 1:m.nvar
        val = A[con, var]
        if val != 0
            m.lin_constrs[A.rowval[k]][var] = A.nzval[k]
            m.j_counts[var] += 1
        end
    end
end

function loadcommon!(m::ReSHOPMathProgBaseModel, x_l, x_u, g_l, g_u, sense)
    m.x_l, m.x_u = x_l, x_u
    m.g_l, m.g_u = g_l, g_u
    setsense!(m, sense)

    m.lin_constrs = [Dict{Int, Float64}() for _ in 1:m.ncon]
    m.j_counts = zeros(Int, m.nvar)

    m.r_codes = Vector{Int}(undef, m.ncon)

    m.varlinearities_con = fill(:Lin, m.nvar)
    m.varlinearities_obj = fill(:Lin, m.nvar)
    m.conlinearities = fill(:Lin, m.ncon)
    m.objlinearity = :Lin

    m.vartypes = fill(:Cont, m.nvar)
    m.x_0 = zeros(m.nvar)
end

getvartype(m::ReSHOPMathProgBaseModel) = copy(m.vartypes)
function setvartype!(m::ReSHOPMathProgBaseModel, cat::Vector{Symbol})
    @assert all(x-> (x in [:Cont,:Bin,:Int,:external]), cat)
    m.vartypes = copy(cat)
end

getsense(m::ReSHOPMathProgBaseModel) = m.sense
function setsense!(m::ReSHOPMathProgBaseModel, sense::Symbol)
    @assert sense == :Min || sense == :Max
    m.sense = sense
end

setwarmstart!(m::ReSHOPMathProgBaseModel, v::Vector{Float64}) = m.x_0 = v

function MPB.addquadconstr!(m::ReSHOPLinearQuadraticModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
   # we have to do a little translation of the sense here ...
   push!(m.inner.quad_equs, tuple(linearidx,linearval,quadrowidx,quadcolidx,quadval, quad_relation_sense[sense], rhs))
end

function MPB.setquadobj!(m::ReSHOPLinearQuadraticModel, rowidx, colidx, quadval)
   m.inner.quad_obj = tuple(rowidx, colidx, quadval)
end

function make_var_index!(m::ReSHOPMathProgBaseModel)
    m.v_index_map = Dict(zip(1:m.nvar, 0:(m.nvar-1)))
    m.v_index_map_rev = Dict(zip(0:(m.nvar-1), 1:m.nvar))
end

#function make_var_index!(m::ReSHOPMathProgBaseModel)
#    nonlin_cont = Int[]
#    nonlin_int = Int[]
#    lin_cont = Int[]
#    lin_int = Int[]
#    lin_bin = Int[]
#
#    # TODO(xhub) we do that multiple times in the EMP context
#    for i in 1:length(m.vartypes)
#        if m.varlinearities_obj[i] == :Nonlin ||
#           m.varlinearities_con[i] == :Nonlin
#            if m.vartypes[i] == :Cont || m.vartypes[i] == :external
#                push!(nonlin_cont, i)
#            else
#                push!(nonlin_int, i)
#            end
#        else
#           if m.vartypes[i] == :Cont || m.vartypes[i] == :external
#                push!(lin_cont, i)
#            elseif m.vartypes[i] == :Int
#                push!(lin_int, i)
#            else
#                push!(lin_bin, i)
#            end
#        end
#    end
#
#    # Index variables in required order
#    for var_list in (nonlin_cont, nonlin_int, lin_cont, lin_bin, lin_int)
#        add_to_index_maps!(m.v_index_map, m.v_index_map_rev, var_list, 0)
#    end
#    CONFIG[:debug] && println("DEBUG: $(m.v_index_map)")
#end

function make_con_index!(m::ReSHOPMathProgBaseModel)
    nonlin_cons = Int[]
    lin_cons = Int[]

    for i in 1:m.ncon
        if m.conlinearities[i] == :Nonlin
            push!(nonlin_cons, i)
        else
            push!(lin_cons, i)
        end
    end
    for con_list in (nonlin_cons, lin_cons)
        add_to_index_maps!(m.nonquad_idx, con_list, 1)
    end

    if length(m.quad_idx) == 0
       m.quad_idx = Dict(enumerate(Vector{Int}(Compat.range(1+m.ncon, length=length(m.quad_equs)))))
    end

    CONFIG[:debug] && println("DEBUG: make_con_index: nonquad_idx: $(m.nonquad_idx)\n quad_idx: $(m.quad_idx)")
end

function optimize!(m::ReSHOPMathProgBaseModel)
    m.status = :NotSolved
    m.solve_exitcode = -1
    m.solve_result_num = -1
    m.solve_result = "?"
    m.model_result_num = -1
    m.model_result = "?"
    m.solve_message = ""

    # There is no non-linear binary type, only non-linear discrete, so make
    # sure binary vars have bounds in [0, 1]
    for i in 1:m.nvar
        if m.vartypes[i] == :Bin
            if m.x_l[i] < 0
                m.x_l[i] = 0
            end
            if m.x_u[i] > 1
                m.x_u[i] = 1
            end
        end
    end

    m.reshop_options = reshop_options_set(m.options)

    # TODO(Xhub) this hack has to go.
    if m.emp.hasvalue
       return m.emp.value()
    end

    make_var_index!(m)
    make_con_index!(m)

    # Run solver and save exitcode
    t = time()
    m.reshop_ctx = create_reshop_ctx(m)
    reshop_set_modeltype(m)
    # Solve via gams for now
    m.reshop_ctx_dest, m.gams_dir = reshop_setup_gams()

    m.reshop_mdl = reshop_alloc(m.reshop_ctx)
    m.reshop_mdl_solver = reshop_alloc(m.reshop_ctx_dest)
    m.solve_exitcode = reshop_solve(m.reshop_mdl, m.reshop_mdl_solver, m.reshop_ctx_dest, m.solver_name)

#    ccall((:print_model, libreshop), Cint, (Ptr{context},), m.reshop_ctx)
    m.solve_time = time() - t

    if m.solve_exitcode == 0
        report_results(m)
    else
        println("ReSHOP: solver failed with status $(m.solve_exitcode)")
        m.status = :Error
        m.solution = fill(NaN, m.nvar)
        m.solve_result = "failure"
        m.solve_result_num = 999
    end
end

function getconstrduals(m::ReSHOPMathProgBaseModel)

    ctx = m.reshop_ctx

    x = fill(NaN, numconstr(m))

   if has_objective(m)
        offset = m.offset
    else
        offset = m.offset - 1
    end

    for idx in 1:m.ncon
       eidx = m.nonquad_idx[idx] + offset
       x[idx] = ctx_getequmult(ctx, eidx)
    end

    return x
end

function getquadconstrduals(quadm::ReSHOPLinearQuadraticModel)

    m = quadm.inner
    ctx = m.reshop_ctx

    x = fill(NaN, numconstr(m))

    if has_objective(m)
        offset = m.offset
    else
        offset = m.offset - 1
    end

    for (idx, equ) in enumerate(m.quad_equs)
       eidx = m.quad_idx[idx] + offset
       x[idx] = ctx_getequmult(ctx, eidx)
    end

    return x
end

function getreducedcosts(m::ReSHOPMathProgBaseModel)

    ctx = m.reshop_ctx

    x = fill(NaN, numvar(m))

    for idx in 1:numvar(m)
       x[idx] = ctx_getvarmult(ctx, idx-1)
    end

    return x
end

function report_results_common(m::ReSHOPMathProgBaseModel)
    x = fill(NaN, m.nvar)
    m.objval = NaN

    for index in 0:(m.nvar - 1)
        i = m.v_index_map_rev[index]
        x[i] = ctx_getvarval(m.reshop_ctx, index)
    end

    m.solution = x

    ###########################################################################
    # Convert solve_result
    #
    # GAMS return two information:
    #  - the solve status
    #  - the model status
    #
    # - :Optimal
    # - :Infeasible
    # - :Unbounded
    # - :UserLimit (iteration limit or timeout)
    # - :Error (and maybe others)
    ###########################################################################

    m.model_result_num = rhp_get_modelstat(m.reshop_ctx)
    m.solve_result_num = rhp_get_solvestat(m.reshop_ctx)

    m.model_result = unsafe_string(ccall((:ctx_getmodelstattxt, libreshop), Cstring, (Ptr{context}, Cint), m.reshop_ctx_dest, m.model_result_num))
    m.solve_result = unsafe_string(ccall((:ctx_getsolvestattxt, libreshop), Cstring, (Ptr{context}, Cint), m.reshop_ctx, m.solve_result_num))

    # GAMS already uses an 1-indices
    solver_code = solver_stat[m.solve_result_num]
    model_code = model_stat[m.model_result_num]

    CONFIG[:debug] && println("solver stat $(m.solve_result) ($(m.solve_result_num)); model stat $(m.model_result) ($(m.model_result_num))")

    if solver_code == :Optimal
       if model_code == :OptimalGlobal || model_code == :OptimalLocal || model_code == :Integer
          m.status = :Optimal
       elseif model_code == :Unbounded || model_code == :UnboundedNoSolution
          m.status = :Unbounded
       elseif model_code == :InfeasibleGlobal || model_code == :InfeasibleLocal || model_code == :InfeasibleIntermed || model_code == :InfeasibleNoSolution
          m.status = :Infeasible
       elseif model_code == :Feasible
          # TODO investigate that. Baron is weird
          m.status = :Optimal
       elseif model_code == :NoSolutionReturned
          gams_solver = ctx_get_solvername(m.reshop_ctx_dest)
          m.status = :Optimal
#          if gams_solver == "jams" || gams_solver == "JAMS"
             # This is fine, we have a kludge in the code
#          else
#             println("ReSHOP: Solve successed, but no solution was returned by solver $(gams_solver)!")
#             m.status = :Error
#          end
       else
          println("ReSHOP: unhandle case: solver stat $(m.solve_result); model stat $(m.model_result)")
          m.status = :Error
       end
    elseif solver_code == :Iteration || solver_code == :Resource
       m.status == :UserLimit
    elseif solver_code == :License || model_code == :LicenseError
       println("ReSHOP: License error. Check that you have a valid license")
       m.status == :Error
    elseif solver_code == :Capability
       if m.solver_name == ""
          sname = default
       else
          sname = m.solver_name
       end

       println("ReSHOP: solver $(sname) cannot solve the specific problem")
       m.status = :Error

    else
       println("ReSHOP: solver stat is $(m.solve_result) and model stat is $(m.model_result)")
       m.status = :Error
    end

    CONFIG[:debug] && println("status is $(m.status)")
 end

getsolvername(s::ReSHOPSolver) = basename(s.solver_name)

MPB.getreducedcosts(nlpm::ReSHOPNonlinearModel) = getreducedcosts(nlpm.inner)
MPB.getreducedcosts(quadm::ReSHOPLinearQuadraticModel) = getreducedcosts(quadm.inner)
MPB.getconstrduals(nlpm::ReSHOPNonlinearModel) = getconstrduals(nlpm.inner)
MPB.getconstrduals(quadm::ReSHOPLinearQuadraticModel) = getconstrduals(quadm.inner)

# Wrapper functions
status(m::ReSHOPMathProgBaseModel) = m.status
getsolution(m::ReSHOPMathProgBaseModel) = copy(m.solution)
getobjval(m::ReSHOPMathProgBaseModel) = m.objval
numvar(m::ReSHOPMathProgBaseModel) = m.nvar
numconstr(m::ReSHOPMathProgBaseModel) = m.ncon + length(m.quad_equs)
getsolvetime(m::ReSHOPMathProgBaseModel) = m.solve_time

# Access to solve results
get_solve_result(m::ReSHOPMathProgBaseModel) = m.solve_result
get_solve_result_num(m::ReSHOPMathProgBaseModel) = m.solve_result_num
get_model_result(m::ReSHOPMathProgBaseModel) = m.model_result
get_model_result_num(m::ReSHOPMathProgBaseModel) = m.model_result_num
get_solve_message(m::ReSHOPMathProgBaseModel) = m.solve_message
get_solve_exitcode(m::ReSHOPMathProgBaseModel) = m.solve_exitcode

for f in [:getvartype,:getsense,:optimize!,:status,:getsolution,:getobjval,:numvar,:numconstr,:getsolvetime]
    @eval MPB.$f(m::ReSHOPNonlinearModel) = $f(m.inner)
    @eval MPB.$f(m::ReSHOPLinearQuadraticModel) = $f(m.inner)
end

for f in [:get_solve_result,:get_solve_result_num,:get_solve_message,:get_solve_exitcode]
    @eval $f(m::ReSHOPNonlinearModel) = $f(m.inner)
    @eval $f(m::ReSHOPLinearQuadraticModel) = $f(m.inner)
end
for f in [:setvartype!,:setsense!,:setwarmstart!]
    @eval MPB.$f(m::ReSHOPNonlinearModel, x) = $f(m.inner, x)
    @eval MPB.$f(m::ReSHOPLinearQuadraticModel, x) = $f(m.inner, x)
end

# Deallocate the data
function reshop_cleanup(o::ReSHOPMathProgBaseModel)
    ctx_dealloc(o.reshop_ctx)
    ctx_dealloc(o.reshop_ctx_dest)
    reshop_options_dealloc(o.reshop_options)
    reshop_free(o.reshop_mdl)
    reshop_free(o.reshop_mdl_solver)
    if (!isempty(o.gams_dir))
        try
            rm(o.gams_dir, recursive=true, force=true)
        catch
            iswin && run(`cmd /C RMDIR /s /q $(o.gams_dir)`)
        end
    end
end


