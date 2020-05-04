struct Basis
    cgto::Vector{CGTO}
end

Basis() = Basis(Vector{CGTO}[])

function Basis(atoms::Atoms, basis_dict::Dict)
    basis = Basis()
    for atom in atoms
        for shell in basis_dict[string(atom.number)]["electron_shells"]
            exponents = parse.(Float64, shell["exponents"])
            for angular in shell["angular_momentum"]
                coefficients = parse.(Float64, shell["coefficients"][angular + 1])
                for L in ANGULAR[angular + 1]
                    push!(basis, CGTO(coefficients, exponents, atom.position, collect(L)))
                end
            end
        end
    end
    return basis
end

function Basis(atoms::Atoms, basis_str::String)
    file = joinpath(@__DIR__, "data/basis/$(basis_str).json")
    basis_dict = JSON.parsefile(file)["elements"]
    return Basis(atoms, basis_dict)
end

Base.iterate(basis::Basis, state=1) = state > length(basis) ? nothing : (basis[state], state+1)
Base.eltype(basis::Basis) = CGTO
Base.length(basis::Basis) = length(basis.cgto)
Base.size(basis::Basis) = length(basis)

Base.getindex(basis::Basis, i::Int) = basis.cgto[i]
Base.firstindex(basis::Basis) = 1
Base.lastindex(basis::Basis) = length(basis)

function Base.push!(basis::Basis, cgto::CGTO)
    push!(basis.cgto, cgto)
end