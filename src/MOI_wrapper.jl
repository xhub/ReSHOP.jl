################################################################################
## MOI wrapper
#################################################################################
## This file is largely inspired from MOI wrappers existing in other
## JuliaOpt packages:
## https://www.juliaopt.org/
## The authors are indebted to the developers of JuliaOpt for
## the current MOI wrapper.
##
## In particular, KNITRO.jl and IpOPT.jl where used as an inspiration

RESHOP_OPTIONS = ["solver"]

# TODO(Xhub) is this useful?
const MOIU = MathOptInterface.Utilities

# Scalar Functions
const SF = Union{MOI.ScalarAffineFunction{Float64},
				 MOI.ScalarQuadraticFunction{Float64}}

# Vector Types
const VAF = MOI.VectorAffineFunction{Float64}
const VOV = MOI.VectorOfVariables

const ALLV = Union{MOI.SingleVariable,MOI.VectorOfVariables}

# Scalar Sets
const SS = Union{MOI.EqualTo{Float64},
					  MOI.GreaterThan{Float64},
					  MOI.LessThan{Float64},
					  MOI.Interval{Float64}}
# LinSets
const LS = Union{MOI.EqualTo{Float64},
					  MOI.GreaterThan{Float64},
					  MOI.LessThan{Float64}}
# VecLinSets
const VLS = Union{MOI.Nonnegatives,
				  MOI.Nonpositives,
				  MOI.Zeros}

# NonConvex Sets
const NCS = Union{MOI.ZeroOne,
									MOI.Integer,
									MOI.Semicontinuous,
									MOI.Semiinteger,
									MOI.SOS1,
									MOI.SOS2}

# NonPolyhedral Cones
const NPC = Union{MOI.SecondOrderCone,
MOI.RotatedSecondOrderCone,
MOI.ExponentialCone,
MOI.DualExponentialCone,
MOI.PowerCone,
MOI.DualPowerCone}

mutable struct Optimizer <: MOI.AbstractOptimizer
    # Get number of solve for restart.
    number_solved::Int
    sense::MOI.OptimizationSense
    status::Cint

    vov_mapping::Dict{MOI.ConstraintIndex{VOV, <: Union{VLS, MOI.SOS1{Float64}, MOI.SOS2{Float64}, MOI.Complements}}, Vector{Cuint}}
    vaf_mapping::Dict{MOI.ConstraintIndex{VAF, <: Union{VLS, MOI.Complements}}, Vector{Cuint}}
    quadfn_mapping::Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, <: SS}, Cuint}
    sos_sets::Dict{MOI.ConstraintIndex{VOV, <: Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}}, Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}}
    fake_cons_name::Dict{MOI.ConstraintIndex, String}
    delvars::Set{MOI.VariableIndex}

    ctx::Union{Ptr{context}, Nothing}
    ctx_dest::Ptr{context}
    mdl::Ptr{reshop_model}
    mdl_solver::Ptr{reshop_model}
    options::Ptr{reshop_options}
    rhp_options::Dict{String, Union{Cdouble, Cint, Bool, Cstring, String}}
    gams_dir::String
    avar_cache::Union{Ptr{abstract_var}, Nothing}
    solver_name::String
    solver_stack::String
    start_nl_cons::Int
    len_nl_cons::Cuint
end

function helper_options(ctx, options, reshop_opts::Ptr{reshop_options})
    solver_name = ""
    solver_stack = ""
    rhp_options = Dict{String, Union{Cdouble, Cint, Bool, Cstring, String}}()
    for (name, value) in options
        sname = string(name)
        if sname == "solver"
            solver_name = value
        elseif sname == "solver_stack"
            solver_stack = value
        else
            res = reshop_option_set(reshop_opts, sname, value)
            if (res != 0)
                 rhp_set_option(ctx, sname, value)
            end
            rhp_options[sname] = value
        end
    end
    return (solver_name, solver_stack, rhp_options)
end

function Optimizer(;options...)
    reshop_set_printops(stdout)
    # Create ReSHOP context.
    ctx = ctx_alloc()

    # TODO this is quite a hack just for the "output" option.
    # Refactoring option in ReSHOP will enable us to move on
    reshop_opts = reshop_options_alloc()
    solver_name, solver_stack, rhp_options = helper_options(ctx, options, reshop_opts)

    model = Optimizer(0, MOI.FEASIBILITY_SENSE, 0,
                      Dict{MOI.ConstraintIndex{VOV, <: VLS}, Cuint}(), Dict{MOI.ConstraintIndex{VAF, <: VLS}, Cuint}(),
                      Dict{MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, <: SS}, Cuint}(),
                      Dict{MOI.ConstraintIndex{VOV, <: Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}}, Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}}(),
                      Dict{MOI.ConstraintIndex,String}(), Set{MOI.VariableIndex}(),
                      ctx, Ptr{context}(C_NULL), Ptr{reshop_model}(C_NULL), Ptr{reshop_model}(C_NULL),
                      reshop_opts, rhp_options, "", nothing, solver_name, solver_stack, -1, 0)

    finalizer(MOI.empty!, model)
    return model
end

# Print Optimizer.
function Base.show(io::IO, model::Optimizer)
    println(io, "A MathOptInterface model with backend:")
    println(io, model.ctx)
    return
end

# copy
MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true
function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(model, src; kws...)
end

# Some utilities.
number_variables(model::Optimizer) = ctx_numvar(model.ctx) - length(model.delvars)
number_constraints(model::Optimizer) = ctx_m(model.ctx)

# Getter for solver's name.
function MOI.get(model::Optimizer, ::MOI.SolverName)
    if model.solver_name == ""
        return "ReSHOP"
    else
        return "ReSHOP ($(model.solver_name))"
    end
end

# MOI.Silent.
MOI.supports(model::Optimizer, ::MOI.Silent) = true
function MOI.get(model::Optimizer, ::MOI.Silent)
    return reshop_option_get(model.options, "output") == 0xF
end

function MOI.set(model::Optimizer, ::MOI.Silent, value)
    val = value ? 0xF : 0x8
    reshop_option_set(model.options, "output", val)
    return
end

# MOI.TimeLimitSec.
MOI.supports(model::Optimizer, ::MOI.TimeLimitSec) = false

# MOI.RawParameters
function MOI.supports(model::Optimizer, param::MOI.RawParameter)
    name = param.name
    if name in RESHOP_OPTIONS
        return true
    end
    return false
end

function MOI.set(model::Optimizer, p::MOI.RawParameter, value)
    if !MOI.supports(model, p)
        throw(MOI.UnsupportedAttribute)
    end
    if p.name == "solver"
         model.solver_name = value
    end
    #reshop_option_set(model.options, p.name, value)
    return
end

function MOI.get(model::Optimizer, p::MOI.RawParameter)
    if haskey(model.options, p.name)
        return model.options[p.name]
    end
    error("RawParameter with name $(p.name) is not set.")
end


##################################################
# Optimize
##################################################
function MOI.optimize!(model::Optimizer)
    model.number_solved += 1
    # If gams_dir already existed (from a past solve), get rid of it
    # TODO(xhub) model.gams_dir could be a list of directory to remove
    if (!isempty(model.gams_dir))
        try
            rm(model.gams_dir, recursive=true, force=true)
        catch
            iswin && run(`cmd /C RMDIR /s /q $(model.gams_dir)`)
        end
    end
    # TODO check if gams_dir and ctx_dest already exists, do not reallocate then.

    solver_stack = get_solverstack(model)

    if solver_stack == "GAMS"
      model.ctx_dest, model.gams_dir = reshop_setup_gams()
    elseif solver_stack == "RESHOP"
      model.ctx_dest = reshop_setup_ownsolver()
    else
      error("Unsupported solver stack $solver_stack")
    end

    # Calling from emp, we already have a mdl object
    if model.mdl == C_NULL
        model.mdl = reshop_alloc(model.ctx)
    end

    model.mdl_solver = reshop_alloc(model.ctx_dest)
    model.status = reshop_solve(model.mdl, model.mdl_solver, model.ctx_dest, model.solver_name)
    reshop_postprocess(model.mdl_solver)

    return
end

function MOI.empty!(model::Optimizer)
    ctx_dealloc(model.ctx_dest)
    model.ctx_dest = C_NULL
    reshop_free(model.mdl)
    model.mdl = C_NULL
    reshop_free(model.mdl_solver)
    model.mdl_solver = C_NULL

    if model.ctx != nothing
        ctx_dealloc(model.ctx)
        model.ctx = ctx_alloc()
        helper_options(model.ctx, model.rhp_options, model.options)
    end

    model.number_solved = 0
    model.vov_mapping = Dict()
    model.vaf_mapping = Dict()
    model.quadfn_mapping = Dict()
    model.sos_sets       = Dict()
    model.fake_cons_name = Dict()
    model.delvars    = Set()
    # TODO(xhub) cleanup avar?
    # TODO(xhub) check that it is not necessary to do anything here
    #set_options(model, model.options)
    if (!isempty(model.gams_dir))
        try
            rm(model.gams_dir, recursive=true, force=true)
        catch
            iswin && run(`cmd /C RMDIR /s /q $(model.gams_dir)`)
        end
    end
    model.start_nl_cons = -1
    model.len_nl_cons = 0
    return
end

function MOI.is_empty(model::Optimizer)
    return model.number_solved == 0 &&
           length(model.vov_mapping) == 0 &&
           length(model.vaf_mapping) == 0 &&
           length(model.quadfn_mapping) == 0 &&
           length(model.sos_sets) == 0 &&
           ctx_m(model.ctx) == 0 &&
           ctx_numvar(model.ctx) == 0
end

include(joinpath("MOI_wrapper", "getters.jl"))
include(joinpath("MOI_wrapper", "supports.jl"))
include(joinpath("MOI_wrapper", "variables.jl"))
include(joinpath("MOI_wrapper", "constraints.jl"))
include(joinpath("MOI_wrapper", "objective.jl"))
include(joinpath("MOI_wrapper", "results.jl"))
include(joinpath("MOI_wrapper", "utils.jl"))
include(joinpath("MOI_wrapper", "changes.jl"))
include(joinpath("MOI_wrapper", "nlp.jl"))
