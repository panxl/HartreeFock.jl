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
    R::Vec3{Float64}
    norm::Float64
end

function CGTO(
    d::AbstractVector{Float64},
    ζ::AbstractVector{Float64},
    L::AbstractVector{Int},
    R::Vec3{Float64},
)
    N = length(d)
    gtos = Vector{GTO}(undef, N)
    for i = 1:N
        gtos[i] = GTO(ζ[i], L)
    end
    norm = normalizaton(d, gtos)
    return CGTO{N}(d, gtos, R, norm)
end

Base.iterate(cgto::CGTO, state=1) = state > length(cgto) ? nothing : (cgto[state], state+1)
Base.eltype(cgto::CGTO) = GTO
Base.length(cgto::CGTO) = length(cgto.d)
Base.size(cgto::CGTO) = length(cgto.d)

Base.getindex(cgto::CGTO, i::Int) = cgto.gtos[i]
Base.firstindex(cgto::CGTO) = 1
Base.lastindex(cgto::CGTO) = length(cgto)

struct Basis
    cgtos::Vector{CGTO}
end

Basis() = Basis(Vector{CGTO}[])

function Base.push!(basis::Basis, cgto::CGTO)
    push!(basis.cgtos, cgto)
end

Base.iterate(basis::Basis, state=1) = state > length(basis) ? nothing : (basis[state], state+1)
Base.eltype(basis::Basis) = CGTO
Base.length(basis::Basis) = length(basis.cgtos)
Base.size(basis::Basis) = length(basis.cgtos)

Base.getindex(basis::Basis, i::Int) = basis.cgtos[i]
Base.firstindex(basis::Basis) = 1
Base.lastindex(basis::Basis) = length(basis)

function Basis(nuclei::Nuclei, basis_dict::Dict)
    basis = Basis()
    for nucleus in nuclei
        for shell in basis_dict[string(nucleus.number)]["electron_shells"]
            exponents = parse.(Float64, shell["exponents"])
            for angular in shell["angular_momentum"]
                coefficients = parse.(Float64, shell["coefficients"][angular + 1])
                for L in ANGULAR_MOMENTUM[angular + 1]
                    push!(basis, CGTO(coefficients, exponents, collect(L), nucleus.position))
                end
            end
        end
    end
    return basis
end

function Basis(nuclei::Nuclei, basis_str::String)
    file = joinpath(@__DIR__, "data/basis/$(basis_str).json")
    basis_dict = JSON.parsefile(file)["elements"]
    return Basis(nuclei, basis_dict)
end