module Libcint

using Libdl
using libcint_jll

const lib = Ref{Ptr{Cvoid}}(0)

function __init__()
    lib[] = Libdl.dlopen(libcint)
    return nothing
end

function func_by_name(name::String)
    func = Libdl.dlsym(lib.x, name)
end

function eval(
    func::Ptr{Nothing},
    buf::AbstractArray{Cdouble},
    shls::AbstractVector{Cint},
    atm::AbstractVector{Cint},
    natm::Cint,
    bas::AbstractVector{Cint},
    nbas::Cint,
    env::AbstractVector{Cdouble},
    )
    not0 = ccall(
        func,
        Cint,
        (Ref{Cdouble}, Ref{Cint}, Ref{Cint}, Cint, Ref{Cint}, Cint, Ref{Cdouble}, Ptr{Cvoid}),
        buf, shls, atm, natm, bas, nbas, env, C_NULL,
        )
end

function eval(
    func::Ptr{Nothing},
    i::Int,
    bas::AbstractVector{Cint}
    )
    n = ccall(
        func,
        Cint,
        (Cint, Ref{Cint}),
        i, bas,
        )
end

end # module