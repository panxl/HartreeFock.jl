# for _atm, _bas, _env
const CHARGE_OF  = 1
const PTR_COORD  = 2
const NUC_MOD_OF = 3
const PTR_ZETA   = PTR_FRAC_CHARGE = 4
const ATM_SLOTS  = 6
const ATOM_OF    = 1
const ANG_OF     = 2
const NPRIM_OF   = 3
const NCTR_OF    = 4
const RADI_POWER = 4 # for ECP
const KAPPA_OF   = 5
const SO_TYPE_OF = 5 # for ECP
const PTR_EXP    = 6
const PTR_COEFF  = 7
const BAS_SLOTS  = 8
# pointer to env
const PTR_LIGHT_SPEED = 1
const PTR_COMMON_ORIG = 2
const PTR_RINV_ORIG   = 5
const PTR_RINV_ZETA   = 8
const PTR_RANGE_OMEGA = 9
const PTR_F12_ZETA    = 10
const PTR_GTG_ZETA    = 11
const AS_RINV_ORIG_ATOM = 18
const AS_ECPBAS_OFFSET = 19
const AS_NECPBAS      = 20
const PTR_ENV_START   = 21
# parameters from libcint
const NUC_POINT = 1
const NUC_GAUSS = 2
# nucleus with fractional charges. It can be used to mimic MM particles
const NUC_FRAC_CHARGE = 3
const NUC_ECP = 4  # atoms with pseudo potential

abstract type AbstractIntor end

struct CIntor{N, M} <: AbstractIntor
    atm::SVector{N, SVector{ATM_SLOTS, Cint}}
    bas::SVector{M, SVector{BAS_SLOTS, Cint}}
    env::Vector{Cdouble}
    shl::SVector{M, Int}
end

function CIntor(nuclei::Nuclei, basis::Basis)
    env = zeros(Cdouble, PTR_ENV_START-1)
    atm = make_atm!(nuclei, env)
    bas = make_bas!(nuclei, basis, env)
    shl = shell_size(nuclei, basis)
    N = length(atm)
    M = length(bas)
    return CIntor{N, M}(atm, bas, env, shl)
end

CIntor(mole::Mole) = CIntor(mole.nuclei, mole.basis)

function make_atm!(nuclei::Nuclei, env::Vector{Float64})
    N = length(nuclei)
    zeta = 0.0
    atm = Vector{SVector{ATM_SLOTS, Cint}}(undef, N)
    for i = 1:N
        nuc_charge = nuclei.numbers[i]
        ptr = length(env)
        _atm = zeros(Cint, ATM_SLOTS)
        _atm[CHARGE_OF] = nuc_charge
        _atm[PTR_COORD] = ptr
        _atm[NUC_MOD_OF] = NUC_POINT
        _atm[PTR_ZETA] = ptr + 3
        atm[i] = _atm
        _env = [nuclei.positions[i]; zeta]
        append!(env, _env)
    end
    return atm
end

function make_bas!(nuclei::Nuclei, basis::Basis, env::Vector{Float64})
    ptr_exp = []
    ptr_coeff = []
    for i = 1:length(basis.shells)
        _ptr_exp = []
        _ptr_coeff = []
        for shell in basis.shells[i]
            ptr = length(env)
            append!(env, shell.exponents)
            push!(_ptr_exp, ptr)
            for coeffs in shell.coefficients
                coeffs = normalizaton_sph(shell.angular_momentum, shell.exponents, coeffs)
                append!(env, coeffs)
            end
            push!(_ptr_coeff, ptr + length(shell.exponents))
        end
        push!(ptr_exp, _ptr_exp)
        push!(ptr_coeff, _ptr_coeff)
    end

    bas = SVector{BAS_SLOTS, Cint}[]
    for i = 1:length(nuclei)
        number = nuclei.numbers[i]
        shells_idx = findfirst(==(number), basis.numbers)
        shells = basis.shells[shells_idx]
        for (j, shell) in enumerate(shells)
            _bas = zeros(Cint, BAS_SLOTS)
            _bas[ATOM_OF] = i - 1
            _bas[ANG_OF] = shell.angular_momentum
            _bas[NPRIM_OF] = length(shell.exponents)
            _bas[NCTR_OF] = length(shell.coefficients)
            _bas[KAPPA_OF] = 0
            _bas[PTR_EXP] = ptr_exp[shells_idx][j]
            _bas[PTR_COEFF] = ptr_coeff[shells_idx][j]
            push!(bas, _bas)
        end
    end
    return bas
end

function shell_size(nuclei::Nuclei, basis::Basis)
    size = Int[]
    for i = 1:length(nuclei)
        number = nuclei.numbers[i]
        shells_idx = findfirst(==(number), basis.numbers)
        shells = basis.shells[shells_idx]
        for (j, shell) in enumerate(shells)
            push!(size, 2*shell.angular_momentum + 1)
        end
    end
    return size
end

function (intor::CIntor)(intor_name::String)
    cint_name = 'c'*intor_name*"_sph"
    func = Libcint.func_by_name(cint_name)
    atm = reinterpret(Cint, Array(intor.atm))
    natm = Cint(length(intor.atm))
    bas = reinterpret(Cint, Array(intor.bas))
    nbas = Cint(length(intor.bas))
    if startswith(cint_name, "cint1e")
        return getints2c(func, atm, natm, bas, nbas, intor.env, intor.shl)
    elseif startswith(cint_name, "cint2e")
        return getints4c(func, atm, natm, bas, nbas, intor.env, intor.shl)
    end
end

function getints2c(
    func::Ptr{Nothing},
    atm::AbstractVector{Cint},
    natm::Cint,
    bas::AbstractVector{Cint},
    nbas::Cint,
    env::AbstractVector{Cdouble},
    shl::AbstractVector{Int},
    )
    M = BlockArray{Cdouble}(undef, shl, shl)
    n = length(shl)
    shls = Vector{Cint}(undef, 2)
    for i = 1:n, j = i:n
        shls .= (i-1,j-1)
        Libcint.eval(func, M[Block(i,j)], shls, atm, natm, bas, nbas, env)
    end
    symmetrize!(M)
    return Array(M)
end

function getints4c(
    func::Ptr{Nothing},
    atm::AbstractVector{Cint},
    natm::Cint,
    bas::AbstractVector{Cint},
    nbas::Cint,
    env::AbstractVector{Cdouble},
    shl::AbstractVector{Int},
    )
    T = BlockArray{Cdouble}(undef, shl, shl, shl, shl)
    n = length(shl)
    shls = Vector{Cint}(undef, 4)
    for i = 1:n, j = i:n
        for k = 1:n, l = k:n
            shls .= (i-1,j-1,k-1,l-1)
            Libcint.eval(func, T[Block(i,j,k,l)], shls, atm, natm, bas, nbas, env)
        end
    end
    symmetrize!(T)
    return T
end

function symmetrize!(M::AbstractMatrix{Float64})
    n = LinearAlgebra.checksquare(M)
    for i = 1:n, j = i+1:n
        M[j, i] = M[i, j]
    end
end

function symmetrize!(T::AbstractArray{Float64,4})
    n = LinearAlgebra.checksquare(T)
    for i = 1:n, j = i:n, k = 1:n, l = k:n
        T[i, j, l, k] = T[i, j, k, l]
        T[j, i, k, l] = T[i, j, k, l]
        T[j, i, l, k] = T[i, j, k, l]
    end
end

function normalizaton_sph(n::Int, es::AbstractVector{Float64}, d::AbstractVector{Float64})
    N = length(d)
    cs = Vector{Float64}(undef, N)
    for i=1:N
        cs[i] = d[i] * normalizaton_sph(n, es[i])
    end

    val = 0.0
    for i = 1:N
        val += (cs[i] / normalizaton_sph(n, es[i]))^2
        for j = i+1:N
            val += 2 * cs[i] * cs[j] / normalizaton_sph(n, (es[i] + es[j])/2)^2
        end
    end
    cs ./= val
    return cs
end

"""Normalized factor for GTO radial part"""
function normalizaton_sph(n::Int, a::Float64)
    return sqrt(2^(2*n+3) * factorial(n+1) * (2*a)^(n+1.5) /
                (factorial(2*n+2) * sqrt(π)))
end