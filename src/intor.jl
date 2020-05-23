struct GTO
    ζ::Float64
    L::Vec3{Int}
    norm::Float64
end

function GTO(
    ζ::Float64,
    L::AbstractVector{Int},
)
    norm = normalizaton(ζ, Vec3(L))
    return GTO(ζ, L, norm)
end

struct CGTO{N}
    d::SVector{N, Float64}
    gtos::SVector{N, GTO}
    norm::Float64
end

function CGTO(
    d::AbstractVector{Float64},
    ζ::AbstractVector{Float64},
    L::AbstractVector{Int},
)
    N = length(d)
    gtos = Vector{GTO}(undef, N)
    for i = 1:N
        gtos[i] = GTO(ζ[i], L)
    end
    norm = normalizaton(d, gtos)
    return CGTO{N}(d, gtos, norm)
end

Base.iterate(cgto::CGTO, state=1) = state > length(cgto) ? nothing : (cgto[state], state+1)
Base.eltype(::CGTO) = GTO
Base.length(cgto::CGTO) = length(cgto.d)
Base.size(cgto::CGTO) = length(cgto.d)

Base.getindex(cgto::CGTO, i::Int) = cgto.gtos[i]
Base.firstindex(cgto::CGTO) = 1
Base.lastindex(cgto::CGTO) = length(cgto)

struct Intor <: AbstractIntor
    nuclei::Nuclei
    cgtos::Vector{CGTO}
    idx::Vector{Int}
end

function Intor(nuclei::Nuclei, basis::Basis)
    cgtos = CGTO[]
    idx = Int[]
    for i = 1:length(nuclei)
        number = nuclei.numbers[i]
        shells_idx = findfirst(==(number), basis.numbers)
        shells = basis.shells[shells_idx]
        for (j, shell) in enumerate(shells)
            exponents = shell.exponents
            angular_momentum = shell.angular_momentum
            for coefficients in shell.coefficients
                for L in ANGULAR_MOMENTUM[angular_momentum + 1]
                    push!(cgtos, CGTO(coefficients, exponents, collect(L)))
                    push!(idx, i)
                end
            end
        end
    end
    return Intor(nuclei, cgtos, idx)
end

Intor(mole::Mole) = Intor(mole.nuclei, mole.basis)

function (intor::Intor)(intor_name::String)
    func = getfield(@__MODULE__, Symbol(intor_name))
    if startswith(intor_name, "int1e")
        return getints2c(func, intor.nuclei, intor.cgtos, intor.idx)
    elseif startswith(intor_name, "int2e")
        return getints4c(func, intor.nuclei, intor.cgtos, intor.idx)
    end
end

function (intor::Intor)(intor_name::String, args...)
    func = getfield(@__MODULE__, Symbol(intor_name))
    return getints2c(func, intor.nuclei, intor.cgtos, intor.idx, args...)
end

function normalizaton(
    a::Float64,
    L::Int,
    )
    return sqrt(sqrt(2*a/π) * (8*a)^L * factorial(L) / factorial(2*L))
end

function normalizaton(
    a::Float64,
    L::Vec3{Int},
    )
    return normalizaton(a, L[1]) * normalizaton(a, L[2]) * normalizaton(a, L[3])
end

normalizaton(gto::GTO) = normalizaton(gto.ζ, gto.L)

function normalizaton(
    d::AbstractVector{Float64},
    gtos::AbstractVector{GTO},
    )
    N = length(d)
    ret = 0.0
    for i = 1:N, j = 1:N
        ret += d[i] * d[j] * int1e_ovlp(gtos[i], gtos[j], Vec3(zeros(3)))
    end
    return 1 / sqrt(ret)
end

function int1e_ovlp(
    A::GTO,
    B::GTO,
    RAB::Vec3{Float64},
    )
    return A.norm * B.norm * int1e_ovlp(A.ζ, A.L, B.ζ, B.L, RAB)
end

int1e_ovlp(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}) = contract(int1e_ovlp, A, RA, B, RB)

function int1e_kin(
    A::GTO,
    B::GTO,
    RAB::Vec3{Float64},
    )
    return A.norm * B.norm * int1e_kin(A.ζ, A.L, B.ζ, B.L, RAB)
end

int1e_kin(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}) = contract(int1e_kin, A, RA, B, RB)

function int1e_rinv(
    A::GTO,
    RA::Vec3{Float64},
    B::GTO,
    RB::Vec3{Float64},
    RC::Vec3{Float64},
    )
    RAB = RA - RB
    RP = gaussian_product_center(A.ζ, RA, B.ζ, RB)
    RPC = RP - RC
    return A.norm * B.norm * int1e_rinv(A.ζ, A.L, B.ζ, B.L, RAB, RPC)
end

int1e_rinv(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, RC::Vec3{Float64}) = contract(int1e_rinv, A, RA, B, RB, RC)

int1e_nuc(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, Q::Float64, RC::Vec3{Float64}) = -Q * contract(int1e_rinv, A, RA, B, RB, RC)

function int2e(
    A::GTO,
    RA::Vec3{Float64},
    B::GTO,
    RB::Vec3{Float64},
    C::GTO,
    RC::Vec3{Float64},
    D::GTO,
    RD::Vec3{Float64},
    )
    RAB = RA - RB
    RCD = RC - RD
    RP = gaussian_product_center(A.ζ, RA, B.ζ, RB)
    RQ = gaussian_product_center(C.ζ, RC, D.ζ, RD)
    RPQ = RP - RQ
    return A.norm * B.norm * C.norm * D.norm * int2e(A.ζ, A.L, B.ζ, B.L, C.ζ, C.L, D.ζ, D.L, RAB, RCD, RPQ)
end

function int2e(
    A::CGTO,
    RA::Vec3{Float64},
    B::CGTO,
    RB::Vec3{Float64},
    C::CGTO,
    RC::Vec3{Float64},
    D::CGTO,
    RD::Vec3{Float64},
    )
    contract(int2e, A, RA, B, RB, C, RC, D, RD)
end

function contract(f::Function, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64})
    RAB = RA - RB
    ret = 0.0
    for i = 1:length(A), j = 1:length(B)
        ret += A.d[i] * B.d[j] * f(A[i], B[j], RAB)
    end
    return A.norm * B.norm * ret
end

function contract(f::Function, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, RC::Vec3{Float64})
    ret = 0.0
    for i = 1:length(A), j = 1:length(B)
        ret += A.d[i] * B.d[j] * f(A[i], RA, B[j], RB, RC)
    end
    return A.norm * B.norm * ret
end

function contract(f::Function, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, C::CGTO, RC::Vec3{Float64}, D::CGTO, RD::Vec3{Float64})
    ret = 0.0
    for i = 1:length(A), j = 1:length(B), k = 1:length(C), l = 1:length(D)
        ret += A.d[i] * B.d[j] * C.d[k] * D.d[l] * 
            f(A[i], RA, B[j], RB, C[k], RC, D[l], RD)
    end
    return A.norm * B.norm * C.norm * D.norm * ret
end

function getints2c(
    func::Function,
    nuclei::Nuclei,
    cgtos::Vector{CGTO},
    idx::Vector{Int},
    )
    n = length(cgtos)
    M = zeros(n, n)
    for i=1:n, j=i:n
        A = cgtos[i]
        RA = nuclei.positions[idx[i]]
        B = cgtos[j]
        RB = nuclei.positions[idx[j]]
        M[i, j] = func(A, RA, B, RB)
    end
    symmetrize!(M)
    return M
end

function getints2c(
    func::typeof(int1e_nuc),
    nuclei::Nuclei,
    cgtos::Vector{CGTO},
    idx::Vector{Int},
    )
    n = length(cgtos)
    M = zeros(n, n)
    for i=1:n, j=i:n
        A = cgtos[i]
        RA = nuclei.positions[idx[i]]
        B = cgtos[j]
        RB = nuclei.positions[idx[j]]
        for (Q, RC) in zip(nuclei.charges, nuclei.positions)
            M[i,j] += func(A, RA, B, RB, Q, RC)
        end
    end
    symmetrize!(M)
    return M
end

function getints2c(
    func::typeof(int1e_nuc),
    nuclei::Nuclei,
    cgtos::Vector{CGTO},
    idx::Vector{Int},
    pg::AbstractParticleGroup,
    )
    n = length(cgtos)
    M = zeros(n, n)
    for i=1:n, j=i:n
        A = cgtos[i]
        RA = nuclei.positions[idx[i]]
        B = cgtos[j]
        RB = nuclei.positions[idx[j]]
        for (Q, RC) in zip(pg.charges, pg.positions)
            M[i,j] += func(A, RA, B, RB, Q, RC)
        end
    end
    symmetrize!(M)
    return M
end

function getints2c(
    func::typeof(int1e_rinv),
    nuclei::Nuclei,
    cgtos::Vector{CGTO},
    idx::Vector{Int},
    point::Vec3{Float64},
    )
    n = length(cgtos)
    M = zeros(n, n)
    for i=1:n, j=i:n
        A = cgtos[i]
        RA = nuclei.positions[idx[i]]
        B = cgtos[j]
        RB = nuclei.positions[idx[j]]
        M[i,j] += func(A, RA, B, RB, point)
    end
    symmetrize!(M)
    return M
end

function getints4c(
    func::Function,
    nuclei::Nuclei,
    cgtos::Vector{CGTO},
    idx::Vector{Int},
    )
    n = length(cgtos)
    T = zeros(n, n, n, n)
    for i = 1:n, j = i:n
        for k = 1:n, l = k:n
            A = cgtos[i]
            RA = nuclei.positions[idx[i]]
            B = cgtos[j]
            RB = nuclei.positions[idx[j]]
            C = cgtos[k]
            RC = nuclei.positions[idx[k]]
            D = cgtos[l]
            RD = nuclei.positions[idx[l]]
            T[i,j,k,l] = func(A, RA, B, RB, C, RC, D, RD)
        end
    end
    symmetrize!(T)
    return T
end