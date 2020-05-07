abstract type AbstractMole end

mutable struct Mole <: AbstractMole
    nuclei::Nuclei
    basis::Basis
    net_charge::Int
    multiplicity::Int
    density_matrix::Matrix{Float64}
    function Mole(
        nuclei::Nuclei,
        basis::Basis,
        net_charge::Int,
        multiplicity::Int,
        )
        n = length(basis)
        density_matrix = zeros(n, n)
        new(nuclei, basis, net_charge, multiplicity, density_matrix)
    end
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