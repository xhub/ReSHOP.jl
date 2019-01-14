function reshop_ovf(emp, ovf)
    # 0-based vs 1-based
    ovf_vidx = ovf.vidx-1
    argsC = Vector{Cint}(ovf.args .- 1)

    ovf_def = emp_ovf(emp, ovf.name, ovf_vidx, argsC)

    for (k, v) in ovf.params
        emp_ovf_param(ovf_def, k, v)
    end

    emp_ovf_check(ovf_def)
end
