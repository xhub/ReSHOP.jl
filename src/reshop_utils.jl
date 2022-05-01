function is_valid_idx(idx)
    if idx == RHP_INDEX_NA || idx == RHP_INDEX_INVALID
        return false
    end

    return true
end

function _c_print(data::Ptr{Cvoid}, mode::Cuint, msg::Ptr{Cchar})
    io = unsafe_pointer_to_objref(data)::IO
   print(io, unsafe_string(msg))
    return
end

function _c_print_stdout(data::Ptr{Cvoid}, mode::Cuint, msg::Ptr{Cchar})
  print(stdout, unsafe_string(msg))
  return
end

function reshop_set_printops(io)
  _C_PRINT = @cfunction(_c_print, Cvoid, (Ptr{Cvoid}, Cuint, Ptr{Cchar}))
  ccall((:reshop_set_printops, libreshop), Cvoid, (Ref{Cvoid}, Ref{Cvoid}), pointer_from_objref(io), _C_PRINT)
  return
end

function reshop_set_printops()
  _C_PRINT = @cfunction(_c_print_stdout, Cvoid, (Ptr{Cvoid}, Cuint, Ptr{Cchar}))
  ccall((:reshop_set_printops, libreshop), Cvoid, (Ref{Cvoid}, Ref{Cvoid}), C_NULL, _C_PRINT)
  return
end

