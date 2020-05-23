abstract type AbstractMole end

struct Mole <: AbstractMole
    nuclei::Nuclei
    basis::Basis
    net_charge::Int
    multiplicity::Int
end

function Mole(
    nuclei::Nuclei,
    basis::Union{Dict,String},
    net_charge::Int,
    multiplicity::Int,
    )
    basis = Basis(nuclei, basis)
    return Mole(nuclei, basis, net_charge, multiplicity)
end

function Mole(
    nuclei::String,
    basis::Union{Basis,Dict,String},
    net_charge::Int,
    multiplicity::Int,
    )
    nuclei = Nuclei(nuclei)
    return Mole(nuclei, basis, net_charge, multiplicity)
end

function Mole(atoms, basis)
    return Mole(atoms, basis, 0, 1)
end

abstract type AbstractEnv end

struct Env <: AbstractEnv
    pointcharges::PointCharges
end

function get_basis_idx(nuclei::Nuclei, basis::Basis)
    idx = Int[]
    for i = 1:length(nuclei)
        number = nuclei.numbers[i]
        shells_idx = findfirst(==(number), basis.numbers)
        shells = basis.shells[shells_idx]
        for (j, shell) in enumerate(shells)
            angular_momentum = shell.angular_momentum
            for coefficients in shell.coefficients
                for _ = 1:2*angular_momentum+1
                    push!(idx, i)
                end
            end
        end
    end
    return idx
end

get_basis_idx(mole::Mole) = get_basis_idx(mole.nuclei, mole.basis)