mutable struct abstract_var
end

mutable struct context
end

mutable struct equtree
end

mutable struct equnode
end

mutable struct reshop_model
end

mutable struct reshop_var
end

mutable struct reshop_equ
end

mutable struct reshop_sp_matrix
end

mutable struct mathprgm
end

mutable struct ovf_definition
end

mutable struct reshop_options
end

mutable struct equil
end

mutable struct filter_ops
end


CONE_NONE = Cint(0)          # Unset/non-existent */
CONE_R_PLUS = Cuint(1)       # Non-negative real \f$\mathbb{R}_+\f$ */
CONE_R_MINUS = Cuint(2)      # Non-positive real \f$\mathbb{R}_-\f$  */
CONE_R = Cuint(3)            # Real \f$\mathbb{R}\f$ */
CONE_0 = Cuint(4)            # Zero singleton */
CONE_POLYHEDRAL = Cuint(5)   # Polyhedral cone */
CONE_SOC = Cuint(6)          # Second Order cone */
CONE_RSOC = Cuint(7)         # Rotated Second Order cone */
CONE_EXP = Cuint(8)          # Exponential cone */
CONE_DEXP = Cuint(9)         # Dual Exponential cone */
CONE_POWER = Cuint(10)       # Power cone */
CONE_DPOWER = Cuint(11)      # Dual Power cone */
__CONES_LEN = Cuint(12)




mutable struct option
	name::Vector{Cchar}
	typ::Cint
	value::Cdouble
end

# TODO with latest Julia
if VERSION < v"0.7"
	iswin = is_windows()
else
	iswin = Sys.iswindows()
end

#const libreshop = iswin ? "reshop" : "libreshop"
const term_str = iswin ? "\r\n" : "\n"

macro chk_index_or_size(idx)
quote
    if $(esc(idx)) >= reshop_valid_index_max
        error("the index/size returned by ReSHOP is not valid: upper limit is $(reshop_valid_index_max), value is ", $(esc(idx)))
    end
end
end

function ctx_add_lin_var(ctx::Ptr{context}, eidx, vidx, coeff::Cdouble)
	res = ccall((:rhp_equ_addvar, libreshop), Cint, (Ptr{context}, Cint, Cint, Cdouble), ctx, eidx, vidx, coeff)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_equ_add_linear(ctx::Ptr{context}, eidx, avar::Ptr{abstract_var}, coeffs::Vector{Cdouble})
	res = ccall((:rhp_equ_addlin, libreshop), Cint,
		    (Ptr{context}, Cint, Ptr{abstract_var}, Ref{Cdouble}),
			ctx, eidx, avar, coeffs)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_equ_add_linear_chk(ctx::Ptr{context}, eidx, avar::Ptr{abstract_var}, coeffs::Vector{Cdouble})
	res = ccall((:rhp_equ_addlinchk, libreshop), Cint,
		    (Ptr{context}, Cint, Ptr{abstract_var}, Ref{Cdouble}),
			ctx, eidx, avar, coeffs)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_alloc()
	ctx = ccall((:ctx_alloc, libreshop), Ptr{context}, (Cuint,), 2)
	if ctx == Ptr{context}(C_NULL)
		error("Could not allocate a context")
	end
	return ctx
end
function ctx_create(n, m)
  ctx = ctx_alloc()
	res = ccall((:model_reserve_eqns, libreshop), Cint, (Ptr{context}, Ptr{Cvoid}, Cuint), ctx, C_NULL, m);
	res != 0 && error("return code $res from ReSHOP")
	res = ccall((:model_reserve_vars, libreshop), Cint, (Ptr{context}, Ptr{Cvoid}, Cuint), ctx, C_NULL, n);
	res != 0 && error("return code $res from ReSHOP")
	res = ccall((:ctx_resize, libreshop), Cint, (Ptr{context}, Cuint, Cuint), ctx, n, m)
	res != 0 && error("return code $res from ReSHOP")
	return ctx
end

function ctx_getvar(ctx, idx)
	return ccall((:ctx_getvar, libreshop), Ptr{reshop_var}, (Ptr{context}, Cint), ctx, idx)
end

function ctx_m(ctx)
	return ccall((:ctx_m, libreshop), Cuint, (Ptr{context},), ctx)
end

function ctx_n(ctx)
	return ccall((:ctx_n, libreshop), Cuint, (Ptr{context},), ctx)
end

function ctx_numvar(ctx::Ptr{context})
	if ctx == Ptr{context}(C_NULL)
		return 0
	end
	res = ccall((:model_total_n, libreshop), Csize_t, (Ptr{context},), ctx)
	return res
end

function ctx_numequ(ctx::Ptr{context})
	if ctx == Ptr{context}(C_NULL)
		return 0
	end
	res = ccall((:model_total_m, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function ctx_setequname(ctx, eidx, name::String)
	res = ccall((:myo_set_equname, libreshop), Cint, (Ptr{context}, Cint, Cstring), ctx, eidx, name)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_setvarnames(ctx, names::Vector{String})
	res = ccall((:myo_set_varnames, libreshop), Cint, (Ptr{context}, Ptr{Ptr{Cchar}}, Cuint), ctx, names, length(names))
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_setvarname(ctx, vidx, name::String)
	res = ccall((:myo_set_varname, libreshop), Cint, (Ptr{context}, Cint, Cstring), ctx, vidx, name)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_getvarname(ctx, vidx)
	res = ccall((:myo_get_varname, libreshop), Cstring, (Ptr{context}, Cint), ctx, vidx)
	res != C_NULL || return ""
	return unsafe_string(res)
end

function hack_last_vidx(ctx::Ptr{context})
	res = ccall((:model_total_n, libreshop), Csize_t, (Ptr{context},), ctx) - 1
	@chk_index_or_size(res)
	return res
end

function hack_exportempinfo(ctx_gms::Ptr{context}, mdl::Ptr{reshop_model})
	res = ccall((:hack_exportempinfo, libreshop), Cint, (Ptr{context}, Ptr{reshop_model}), ctx_gms, mdl)
	res != 0 && error("return code $res from ReSHOP")
end

function hack_solver_log()
	ccall((:hack_solver_log, libreshop), Cvoid, ())
end

function ctx_setvarval(ctx::Ptr{context}, idx, val::Cdouble)
	res = ccall((:ctx_setvarlone, libreshop), Cint, (Ptr{context}, Cint, Cdouble), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_setvarlb(ctx::Ptr{context}, idx, val::Cdouble)
	res = ccall((:ctx_setvarlb, libreshop), Cint, (Ptr{context}, Cint, Cdouble), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_setvarub(ctx::Ptr{context}, idx, val::Cdouble)
	res = ccall((:ctx_setvarub, libreshop), Cint, (Ptr{context}, Cint, Cdouble), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_setvarfx(ctx::Ptr{context}, idx, val::Cdouble)
	ctx_setvarub(ctx, idx, val)
	ctx_setvarlb(ctx, idx, val)
end

function ctx_getvarbounds(ctx::Ptr{context}, idx)
	lb = Ref{Cdouble}(NaN)
	ub = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getvarbounds, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}, Ref{Cdouble}), ctx, idx, lb, ub)
	res != 0 && error("return code $res from ReSHOP")
	return (lb.x, ub.x)
end

function ctx_getvarbstat(ctx::Ptr{context}, idx)
	val = Ref{Cint}(-1)
	res = ccall((:ctx_getvarstatone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cint}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function ctx_getvarmult(ctx::Ptr{context}, idx)
	val = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getvarmone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

#function ctx_getvarname(ctx::Ptr{context}, idx)
#	len = 1024
#	buf = Vector{UInt8}(undef, len)
#	res = ccall((:ctx_getvarname, libreshop), Cint, (Ptr{context}, Cint, Ptr{UInt8}, Cuint), ctx, idx, buf, len)
#	res != 0 && error("return code $res from ReSHOP")
#	e = findfirst(isequal(0), buf)
#	return String(buf[1:e])
#end

function ctx_getvarbyname(ctx::Ptr{context}, name::String)
	val = Ref{Cint}(-1)
	res = ccall((:ctx_getvarbyname, libreshop), Cint, (Ptr{context}, Cstring, Ptr{Cint}), ctx, name, val)
	if res == 5
	  return -2
	end
	return val.x
end

function ctx_getequbyname(ctx::Ptr{context}, name::String)
	val = Ref{Cint}(-1)
	res = ccall((:ctx_getequbyname, libreshop), Cint, (Ptr{context}, Cstring, Ptr{Cint}), ctx, name, val)
	if res == 5
	  return -2
	end
	return val.x
end

function reshop_get_equtype(ctx::Ptr{context}, idx)
        typ = Ref{Cuint}(1000)
	cone = Ref{Cuint}(1000)
	res = ccall((:ctx_getequtype, libreshop), Cint, (Ptr{context}, Cint, Ref{Cuint}, Ref{Cuint}), ctx, idx, typ, cone)
	res != 0 && error("return code $res from ReSHOP")
	return (typ.x, cone.x)
end

function reshop_getvartype(ctx::Ptr{context}, idx)
	typ = Ref{Cuint}(1000)
	res = ccall((:ctx_getvartype, libreshop), Cint, (Ptr{context}, Cint, Ptr{Cuint}), ctx, idx, typ)
	res != 0 && error("return code $res from ReSHOP")
	return typ.x
end

function ctx_getvarval(ctx::Ptr{context}, idx)
	val = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getvarlone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function ctx_getequ(ctx, idx)
	equ = ccall((:ctx_getequ, libreshop), Ptr{reshop_equ}, (Ptr{context}, Cint), ctx, idx)
end

function ctx_getequval(ctx, idx)
	val = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getequlone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function ctx_getrhs(ctx::Ptr{context}, idx)
	val = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getrhs, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end


function ctx_getequmult(ctx::Ptr{context}, idx)
	val = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getequmone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function ctx_getequbstat(ctx::Ptr{context}, idx)
	val = Ref{Cint}(-1)
	res = ccall((:ctx_getequstatone, libreshop), Cint, (Ptr{context}, Cint, Ref{Cint}), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function ctx_getequname(ctx::Ptr{context}, eidx)
	res = ccall((:myo_get_equname, libreshop), Cstring, (Ptr{context}, Cint), ctx, eidx)
	res != C_NULL || return ""
	return unsafe_string(res)
end

function ctx_getequbounds(ctx::Ptr{context}, idx)
	lb = Ref{Cdouble}(NaN)
	ub = Ref{Cdouble}(NaN)
	res = ccall((:ctx_getequbounds, libreshop), Cint, (Ptr{context}, Cint, Ref{Cdouble}, Ref{Cdouble}), ctx, idx, lb, ub)
	res != 0 && error("return code $res from ReSHOP")
	return (lb.x, ub.x)
end

function emp_init(mdl::Ptr{reshop_model})
	res = ccall((:reshop_alloc_emp, libreshop), Cint, (Ptr{reshop_model},), mdl)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_create_equil(max_mp)
	return ccall((:mp_equil_alloc, libreshop), Ptr{equil}, (Cuint,), max_mp)
end

function emp_add_mp_mp(mp_parent, mp)
	res = ccall((:mathprgm_addmp, libreshop), Cint, (Ptr{mathprgm}, Ptr{mathprgm}), mp_parent, mp)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_add_mp_equil(mp_parent, mpe)
	res = ccall((:mathprgm_addequil, libreshop), Cint, (Ptr{mathprgm}, Ptr{equil}), mp_parent, mpe)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_add_equil(mdl::Ptr{reshop_model})
	mpe = Ref{Ptr{equil}}(C_NULL)
	res = ccall((:reshop_add_equil, libreshop), Cint, (Ptr{reshop_model}, Ref{Ptr{equil}}), mdl, mpe)
	res != 0 && error("return code $res from ReSHOP")

	return mpe.x
end

function emp_equil_add(mpe, mp)
	res = ccall((:mp_equil_add, libreshop), Cint, (Ptr{equil}, Ptr{mathprgm}), mpe, mp)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_ensure(mdl::Ptr{reshop_model}, nb)
	res = ccall((:reshop_ensure_mp, libreshop), Cint, (Ptr{reshop_model}, Cuint), mdl, nb)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_alloc(mdl::Ptr{reshop_model})
	mp = ccall((:mathprgm_alloc, libreshop), Ptr{mathprgm}, (Ptr{reshop_model},), mdl)
	return mp
end

function emp_mp_start(mp, typ)
	res = ccall((:mathprgm_addstart, libreshop), Cint, (Ptr{mathprgm}, Cuint), mp, typ)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_objdir(mp, sense)
	res = ccall((:mathprgm_addobjdir, libreshop), Cint, (Ptr{mathprgm}, Cint), mp, sense_to_reshop[sense])
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_objequ(mp, idx)
	res = ccall((:mathprgm_addobjequ, libreshop), Cint, (Ptr{mathprgm}, Cint), mp, idx)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_objvar(mp, idx)
	res = ccall((:mathprgm_addobjvar, libreshop), Cint, (Ptr{mathprgm}, Cint), mp, idx)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_var(mp, idx)
	res = ccall((:mathprgm_addvar, libreshop), Cint, (Ptr{mathprgm}, Cint), mp, idx)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_constraint(mp, idx)
	res = ccall((:mathprgm_addconstraint, libreshop), Cint, (Ptr{mathprgm}, Cint), mp, idx)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_vipair(mp, eidx, vidx)
	res = ccall((:mathprgm_addvipair, libreshop), Cint, (Ptr{mathprgm}, Cint, Cint), mp, eidx, vidx)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_mp_getobjvar(mp)
	return ccall((:mathprgm_getobjvar, libreshop), Cint, (Ptr{mathprgm},), mp)
end

function emp_mp_getobjdir(mp)
	return ccall((:mathprgm_getobjvar, libreshop), Cint, (Ptr{mathprgm},), mp)
end

function emp_mp_print(mp, ctx)
	ccall((:mathprgm_print, libreshop), Cint, (Ptr{mathprgm}, Ptr{context}), mp, ctx)
end

function reshop_postprocess(mdl_solver::Ptr{reshop_model})
	res = ccall((:rhp_postprocess, libreshop), Cint, (Ptr{reshop_model},), mdl_solver)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_set_root(mdl::Ptr{reshop_model}, mpe::Ptr{equil})
	res = ccall((:rhp_emproot_setmpe, libreshop), Cint, (Ptr{reshop_model}, Ptr{equil}), mdl, mpe)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_set_root(mdl::Ptr{reshop_model}, mp::Ptr{mathprgm})
	res = ccall((:rhp_emproot_setmp, libreshop), Cint, (Ptr{reshop_model}, Ptr{mathprgm}), mdl, mp)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_ovf(mdl, name, ovf_vidx, args)
	ovf_def = Ref{Ptr{ovf_definition}}(C_NULL)
	args_var = reshop_avar(length(args), args)
	res = ccall((:rhp_ovf_add, libreshop), Cint, (Ptr{reshop_model}, Cstring, Cint, Ptr{abstract_var}, Ref{Ptr{ovf_definition}}),
							                              mdl, name, ovf_vidx, args_var, ovf_def)
	res != 0 && error("return code $res from ReSHOP")
	reshop_avar_free(args_var)
	return ovf_def.x
end

function emp_ovf_param(ovf_def, param_name, scalar::Number)
	res = ccall((:rhp_ovf_param_add_scalar, libreshop), Cint, (Ptr{ovf_definition}, Cstring, Cdouble), ovf_def, param_name, scalar)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_ovf_param(ovf_def, param_name, arr::Vector)
  arrC = Vector{Cdouble}(arr)
	res = ccall((:rhp_ovf_param_add_vector, libreshop), Cint, (Ptr{ovf_definition}, Cstring, Cuint, Ptr{Cdouble}), ovf_def, param_name, length(arr), arrC)
	res != 0 && error("return code $res from ReSHOP")
end

function emp_ovf_check(ovf_def)
	res = ccall((:rhp_ovf_check, libreshop), Cint, (Ptr{ovf_definition},), ovf_def)
end

function equtree_var(ctx, tree, node, idx, coeff)
	res = ccall((:equtree_var, libreshop), Cint, (Ptr{context}, Ptr{equtree}, Ref{Ref{Ptr{equnode}}}, Cint, Cdouble),
		ctx,
		tree,
		node,
		idx,
		coeff)
	res != 0 && error("return code $res from ReSHOP")
end

function equtree_cst(ctx, tree, node, value)
	res = ccall((:equtree_cst, libreshop), Cint, (Ptr{context}, Ptr{equtree}, Ref{Ref{Ptr{equnode}}}, Cdouble),
		ctx,
		tree,
		node,
		value)
	res != 0 && error("return code $res from ReSHOP")
end

function equtree_arithm(tree, node, opcode, nb)
	res = ccall((:equtree_arithm, libreshop), Cint, (Ptr{equtree}, Ref{Ref{Ptr{equnode}}}, Cuint, Cuint),
		tree,
		node,
		opcode,
		nb)
	res != 0 && error("return code $res from ReSHOP")
end

function equtree_call(ctx, tree, node, fndata)
	res = ccall((:equtree_call, libreshop), Cint, (Ptr{context}, Ptr{equtree}, Ref{Ref{Ptr{equnode}}}, Cuint, Cuint),
		ctx,
		tree,
		node,
		fndata[1],
		fndata[2]
		)
	res != 0 && error("return code $res from ReSHOP")
end

function equtree_get_root_addr(tree::Ptr{equtree}, node::Ref{Ref{Ptr{equnode}}})
	res = ccall((:equtree_get_root_addr, libreshop), Cint, (Ptr{equtree}, Ref{Ref{Ptr{equnode}}}),
		tree,
		node)
	res != 0 && error("return code $res from ReSHOP")
end

function equtree_umin(ctx, tree, node)
	res = ccall((:equtree_umin, libreshop), Cint, (Ptr{equtree}, Ref{Ref{Ptr{equnode}}}),
		tree,
		node)
	res != 0 && error("return code $res from ReSHOP")
end

function equnode_get_child_addr(node::Ptr{equnode}, i::Int)
	child = Ref{Ref{Ptr{equnode}}}(C_NULL)
	res = ccall((:equnode_get_child_addr, libreshop), Cint, (Ptr{equnode}, Ref{Ref{Ptr{equnode}}}, Cuint),
		node,
		child,
		i)
	res != 0 && error("return code $res from ReSHOP")
	return child
end

function equnode_deref(node::Ref{Ref{Ptr{equnode}}})
	return ccall((:p2deref, libreshop), Ptr{equnode}, (Ref{Ref{Ptr{equnode}}},), node)
end

function ctx_get_solvername(ctx::Ptr{context})
	str = "?"^256
	res = ccall((:gams_getsolverstr, libreshop), Cint, (Ptr{context}, Cstring, Cuint), ctx, str, length(str))
	res != 0 && error("return code $res from ReSHOP")
	return str
end

function rhp_add_free_var(ctx::Ptr{context}, nb, avar::Ptr{abstract_var}=Ptr{abstract_var}(C_NULL))
	res = ccall((:model_add_free_vars, libreshop), Cint, (Ptr{context}, Cuint, Ptr{abstract_var}), ctx, nb, avar)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_add_neg_var(ctx::Ptr{context}, nb, avar::Ptr{abstract_var}=Ptr{abstract_var}(C_NULL))
	res = ccall((:model_add_neg_vars, libreshop), Cint, (Ptr{context}, Cuint, Ptr{abstract_var}), ctx, nb, avar)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_add_pos_var(ctx::Ptr{context}, nb, avar::Ptr{abstract_var}=Ptr{abstract_var}(C_NULL))
	res = ccall((:model_add_pos_vars, libreshop), Cint, (Ptr{context}, Cuint, Ptr{abstract_var}), ctx, nb, avar)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_alloc(ctx::Ptr{context})
	return ccall((:reshop_alloc, libreshop), Ptr{reshop_model}, (Ptr{context},), ctx)
end

function reshop_avar()
  return ccall((:avar_alloc, libreshop), Ptr{abstract_var}, ())
end

function reshop_avar(size, indices)
	indicesC = Vector{Cint}(indices)
	return ccall((:avar_alloc_list, libreshop), Ptr{abstract_var}, (Cuint, Ptr{Cint}), size, indicesC)
end

function rhp_avar_set(avar::Ptr{abstract_var}, indices::Integer)
	indicesC = Vector{Cint}([indices])
	res = ccall((:avar_set_list, libreshop), Cint, (Ptr{abstract_var}, Cuint, Ptr{Cint}), avar, length(indicesC), indicesC)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_avar_set(avar::Ptr{abstract_var}, indices::AbstractVector{Int32})
	indicesC = Vector{Cint}(indices)
	res = ccall((:avar_set_list, libreshop), Cint, (Ptr{abstract_var}, Cuint, Ptr{Cint}), avar, length(indices), indicesC)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_avar_free(avar::Ptr{abstract_var})
	return ccall((:avar_free, libreshop), Cvoid, (Ptr{abstract_var},), avar)
end

function reshop_avar_get(avar::Ptr{abstract_var}, i::UInt32)
	vidx = Ref{Cint}(-1)
	res = ccall((:avar_get, libreshop), Cint, (Ptr{abstract_var}, Cuint, Ref{Cint}), avar, i, vidx)
	res != 0 && error("return code $res from ReSHOP")
	return vidx.x
end

function reshop_avar_size(avar::Ptr{abstract_var})
	res = ccall((:avar_size, libreshop), Csize_t, (Ptr{abstract_var},), avar)
	@chk_index_or_size(res)
	return res
end

function reshop_decl_eqn(ctx::Ptr{context})
	minn = Ref{Cint}(-1)
	res = ccall((:rhp_add_equation, libreshop), Cint, (Ptr{context}, Ref{Cint}), ctx, minn)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_equ_add_lin_tree(ctx::Ptr{context}, eidx, qvals, avar::Ptr{abstract_var}, coeff)
	equ = ctx_getequ(ctx, eidx)
	res = ccall((:equ_add_lin_tree, libreshop), Cint, (Ptr{context}, Ptr{reshop_equ}, Ref{Cdouble}, Ptr{abstract_var}, Cdouble),
							ctx, equ, qvals, avar, coeff)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_equ_add_quadratic(ctx::Ptr{context}, eidx, mat::Ptr{reshop_sp_matrix}, avar::Ptr{abstract_var}, coeff)
	equ = ctx_getequ(ctx, eidx)
	res = ccall((:equ_add_quadratic, libreshop), Cint, (Ptr{context}, Ptr{reshop_equ}, Ptr{reshop_sp_matrix}, Ptr{abstract_var}, Cdouble),
							ctx, equ, mat, avar, coeff)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_equ_add_quadratic(ctx::Ptr{context}, eidx, vidx1, vidx2, coeffs)
	vidx1C = Vector{Cint}(vidx1)
	vidx2C = Vector{Cint}(vidx2)
	@assert length(vidx1) == length(vidx2) == length(coeffs)
	res = ccall((:rhp_equ_addquadabsolute, libreshop), Cint,
		    (Ptr{context}, Cint, Csize_t, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cdouble),
		    ctx, eidx, length(vidx1C), vidx1C, vidx2C, coeffs, 1.)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_free(mdl::Ptr{reshop_model})
	return ccall((:reshop_free, libreshop), Cvoid, (Ref{Ptr{reshop_model}},), mdl)
end

function reshop_mat_coo(ridx, cidx, vals)
	# beaware of dragons! --xhub
	ridxC = Vector{Cint}(ridx)
	cidxC = Vector{Cint}(cidx)
	mat = ccall((:rhp_mat_triplet, libreshop), Ptr{reshop_sp_matrix}, (Cuint, Cuint, Csize_t, Ptr{Cint}, Ptr{Cint}, Ref{Cdouble}),
							0, 0, length(vals), ridxC, cidxC, vals)
	return mat
end

function reshop_mat_free(mat)
	ccall((:rhp_mat_free, libreshop), Cvoid, (Ptr{reshop_sp_matrix},), mat)
end

function rhp_set_objeqn(ctx::Ptr{context}, eidx)
	res = ccall((:model_setobjequ, libreshop), Cint, (Ptr{context}, Cint), ctx, eidx)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_objvar(ctx::Ptr{context}, vidx)
	res = ccall((:ctx_setobjvar, libreshop), Cint, (Ptr{context}, Cint), ctx, vidx)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_getobjequ(ctx::Ptr{context})
	tmpCint = Ref{Cint}(0)
	res = ccall((:ctx_getobjequ, libreshop), Cint, (Ptr{context}, Ptr{Cint}), ctx, tmpCint)
	res != 0 && error("return code $res from ReSHOP")
	return tmpCint.x
end

function ctx_getobjvar(ctx::Ptr{context})
	tmpCint = Ref{Cint}(0)
	res = ccall((:ctx_getobjvar, libreshop), Cint, (Ptr{context}, Ptr{Cint}), ctx, tmpCint)
	res != 0 && error("return code $res from ReSHOP")
	return tmpCint.x
end

function reshop_set_rhs(ctx::Ptr{context}, idx, val)
	res = ccall((:ctx_setrhs, libreshop), Cint, (Ptr{context}, Cint, Cdouble), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_cst(ctx::Ptr{context}, idx, val)
	res = ccall((:rhp_equ_setcst, libreshop), Cint, (Ptr{context}, Cint, Cdouble), ctx, idx, val)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_perp(ctx::Ptr{context}, ei, vi)
	res = ccall((:ctx_setperp, libreshop), Cint, (Ptr{context}, RHP_IDXT, RHP_IDXT), ctx, ei, vi)
	res != 0 && error("return code $res from ReSHOP")
end


function reshop_set_equtype(ctx::Ptr{context}, idx, cone)
	res = ccall((:ctx_setequtype, libreshop), Cint, (Ptr{context}, Cint, Cuint, Cuint), ctx, idx, 2, cone)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_vartype(ctx::Ptr{context}, idx, typ)
	res = ccall((:ctx_setvartype, libreshop), Cint, (Ptr{context}, Cint, Cuint), ctx, idx, typ)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_vartype(ctx::Ptr{context}, avar::Ptr{abstract_var}, ::MOI.SOS1{Float64}, weights::AbstractVector{Cdouble})
  res = ccall((:rhp_set_var_sos1, libreshop), Cint, (Ptr{context}, Ptr{abstract_var}, Ptr{Cdouble}), ctx, avar, weights)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_vartype(ctx::Ptr{context}, avar::Ptr{abstract_var}, ::MOI.SOS2{Float64}, weights::AbstractVector{Cdouble})
  res = ccall((:rhp_set_var_sos2, libreshop), Cint, (Ptr{context}, Ptr{abstract_var}, Ptr{Cdouble}), ctx, avar, weights)
	res != 0 && error("return code $res from ReSHOP")
end

function reshop_get_treedata(ctx, eidx::Int)
	tree = ccall((:myo_getequtree, libreshop), Ptr{equtree}, (Ptr{context}, Cint), ctx, eidx)
	node = Ref{Ref{Ptr{equnode}}}(C_NULL)
	res = equtree_get_root_addr(tree, node)
	res != 0 && error("return code $res from ReSHOP")
	return (tree, node)
end

function rhp_add_equ(ctx::Ptr{context})
	return rhp_add_equs(ctx, 1)
end

function rhp_add_equs(ctx::Ptr{context}, nb)
	eidx_start = Int(ctx_numequ(ctx))
	res = ccall((:myo_addinit_equs, libreshop), Cint, (Ptr{context}, Cuint), ctx, nb)
	res != 0 && error("return code $res from ReSHOP")
	return eidx_start
end

function rhp_add_var(ctx::Ptr{context}, avar::Ptr{abstract_var})
	return rhp_add_var(ctx, 1, avar)
end

function rhp_add_var(ctx::Ptr{context}, nb, avar::Ptr{abstract_var})
	res = ccall((:rhp_add_vars, libreshop), Cint, (Ptr{context}, Cuint, Ptr{abstract_var}), ctx, nb, avar)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_get_solvestat(ctx::Ptr{context})
	tmpCint = Ref{Cint}(0)
	res = ccall((:ctx_getsolvestat, libreshop), Cint, (Ptr{context}, Ref{Cint}), ctx, tmpCint)
	res != 0 && error("return code $res from ReSHOP")
	return tmpCint.x
end

function rhp_get_modelstat(ctx::Ptr{context})
	tmpCint = Ref{Cint}(0)
	res = ccall((:ctx_getmodelstat, libreshop), Cint, (Ptr{context}, Ref{Cint}), ctx, tmpCint)
	res != 0 && error("return code $res from ReSHOP")
	return tmpCint.x
end

function reshop_options_alloc()
	return ccall((:hack_options_alloc, libreshop), Ptr{reshop_options}, ())
end

function reshop_options_dealloc(o::Ptr{reshop_options})
	ccall((:hack_options_dealloc, libreshop), Cvoid, (Ptr{reshop_options},), o)
end

function reshop_option_get(opt::Ptr{reshop_options}, k::String, val::Ref{Cdouble})
	res = ccall((:option_get_d, libreshop), Cint, (Ptr{reshop_options}, Cstring, Ptr{Cdouble}), opt, k, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function reshop_option_get(opt::Ptr{reshop_options}, k::String, val::Ref{Cstring})
	res = ccall((:option_get_i, libreshop), Cint, (Ptr{reshop_options}, Cstring, Ptr{Cstring}), opt, k, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function reshop_option_get(opt::Ptr{reshop_options}, k::String, val::Ref{Cint})
	res = ccall((:option_get_i, libreshop), Cint, (Ptr{reshop_options}, Cstring, Ptr{Cint}), opt, k, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function reshop_option_get(opt::Ptr{reshop_options}, k::String, val::Ref{Cuint})
	res = ccall((:option_get_i, libreshop), Cint, (Ptr{reshop_options}, Cstring, Ptr{Cuint}), opt, k, val)
	res != 0 && error("return code $res from ReSHOP")
	return val.x
end

function reshop_option_get(opt::Ptr{reshop_options}, k::String)
	res = ccall((:option_get_type, libreshop), Cint, (Ptr{reshop_options}, Cstring), opt, k)
	if haskey(reshop_option_type, res)
		fn, type = reshop_option_type[res]
		val = Ref{type}()
		reshop_option_get(opt, k, val)
		return val.x
	else
		error("unknown option named $k, ReSHOP code is $res")
	end
end

function reshop_option_set(opt::Ptr{reshop_options}, k::String, v::Bool)
	return ccall((:option_set_b, libreshop), Cint, (Ptr{reshop_options}, Cstring, Cint), opt, k, v)
end

function reshop_option_set(opt::Ptr{reshop_options}, k::String, v::Integer)
	return ccall((:option_set_i, libreshop), Cint, (Ptr{reshop_options}, Cstring, Cint), opt, k, v)
end

function reshop_option_set(opt::Ptr{reshop_options}, k::String, v::Cdouble)
	return ccall((:option_set_d, libreshop), Cint, (Ptr{reshop_options}, Cstring, Cdouble), opt, k, v)
end

function reshop_option_set(opt::Ptr{reshop_options}, k::String, v::String)
	return ccall((:option_set_str, libreshop), Cint, (Ptr{reshop_options}, Cstring, Cstring), opt, k, v)
end

function reshop_set_modeltype(ctx::Ptr{context}, model_type)
    res = ccall((:ctx_setmodeltype, libreshop), Cint, (Ptr{context}, Cint), ctx, model_type)
    res != 0 && error("return code $res from ReSHOP")
end

function reshop_set_objsense(ctx::Ptr{context}, sense)
    res = ccall((:ctx_setobjsense, libreshop), Cint, (Ptr{context}, Cint), ctx, sense_to_reshop[sense])
    res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_option(ctx::Ptr{context}, k::String, v::Bool)
	res = ccall((:rhp_set_option_b, libreshop), Cint, (Ptr{context}, Cstring, Cuchar), ctx, k, v)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_option(ctx::Ptr{context}, k::String, v::Cint)
	res = ccall((:rhp_set_option_i, libreshop), Cint, (Ptr{context}, Cstring, Cint), ctx, k, v)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_option(ctx::Ptr{context}, k::String, v::Cdouble)
	res = ccall((:rhp_set_option_d, libreshop), Cint, (Ptr{context}, Cstring, Cdouble), ctx, k, v)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_set_option(ctx::Ptr{context}, k::String, v::Union{Cstring, String})
	res = ccall((:rhp_set_option_s, libreshop), Cint, (Ptr{context}, Cstring, Cstring), ctx, k, v)
	res != 0 && error("return code $res from ReSHOP")
end

function ctx_dealloc(o::Ptr{context})
	ccall((:ctx_dealloc,  libreshop), Cvoid, (Ptr{context},), o)
end

function rhp_delete_var(ctx::Ptr{context}, vidx)
	res = ccall((:rhp_delete_var,  libreshop), Cint, (Ptr{context}, Cint), ctx, vidx)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_delete_eqn(ctx::Ptr{context}, eidx)
	res = ccall((:rhp_delete_eqn,  libreshop), Cint, (Ptr{context}, Cint), ctx, eidx)
	res != 0 && error("return code $res from ReSHOP")
end

function rhp_get_nb_bin_var(ctx::Ptr{context})
	res = ccall((:rhp_get_nb_var_bin, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_int_var(ctx::Ptr{context})
	res = ccall((:rhp_get_nb_var_int, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.GreaterThan{Float64}})
	res = ccall((:rhp_get_nb_var_lb, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.LessThan{Float64}})
	res = ccall((:rhp_get_nb_var_ub, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.EqualTo{Float64}})
	res = ccall((:rhp_get_nb_var_fx, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.Interval{Float64}})
	res = ccall((:rhp_get_nb_var_interval, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.SOS1{Float64}})
  res = ccall((:rhp_get_nb_var_sos1, libreshop), Csize_t, (Ptr{context},), ctx)
  @chk_index_or_size(res)
  return res
end

function rhp_get_nb_var(ctx::Ptr{context}, ::Type{MOI.SOS2{Float64}})
  res = ccall((:rhp_get_nb_var_sos2, libreshop), Csize_t, (Ptr{context},), ctx)
  @chk_index_or_size(res)
  return res
end

function rhp_get_nb_lequ(ctx::Ptr{context}, ::Type{MOI.GreaterThan{Float64}})
	res = ccall((:rhp_get_nb_lequ_ge, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_lequ(ctx::Ptr{context}, ::Type{MOI.LessThan{Float64}})
	res = ccall((:rhp_get_nb_lequ_le, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_nb_lequ(ctx::Ptr{context}, ::Type{MOI.EqualTo{Float64}})
	res = ccall((:rhp_get_nb_lequ_eq, libreshop), Csize_t, (Ptr{context},), ctx)
	@chk_index_or_size(res)
	return res
end

function rhp_get_sos_group(ctx::Ptr{context}, vidx, ::Type{MOI.SOS1{Float64}})
	grps = Ref{Ref{Cuint}}(C_NULL)
	res = ccall((:rhp_get_var_sos1, libreshop), Cint, (Ptr{context}, Cint, Ref{Ptr{Cuint}}), ctx, vidx, grps)
	res != 0 && error("return code $res from ReSHOP")
	return grps.x
end

function rhp_get_sos_group(ctx::Ptr{context}, vidx, ::Type{MOI.SOS2{Float64}})
	grps = Ref{Ptr{Cuint}}(C_NULL)
	res = ccall((:rhp_get_var_sos2, libreshop), Cint, (Ptr{context}, Cint, Ref{Ptr{Cuint}}), ctx, vidx, grps)
	res != 0 && error("return code $res from ReSHOP")
	return grps.x
end

function rhp_is_var_valid(ctx::Ptr{context}, vidx)
	return ccall((:rhp_is_var_valid, libreshop), Cint, (Ptr{context}, Cint,), ctx, vidx) == 1 ? true : false
end

function rhp_is_equ_valid(ctx::Ptr{context}, eidx)
	return ccall((:rhp_is_equ_valid, libreshop), Cint, (Ptr{context}, Cint,), ctx, eidx)  == 1 ? true : false
end
