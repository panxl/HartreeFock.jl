abstract type AbstractMole end

struct Mole <: AbstractMole
    atoms::Atoms
    basis::Basis
    charge::Int
    multiplicity::Int
end

function Mole(atoms::Atoms, basis_dict::Dict, charge::Int, multiplicity::Int)
    basis = Basis(atoms, basis_dict)
    return Mole(atoms, basis, charge, multiplicity)
end

function Mole(atoms::Atoms, basis_str::String, charge::Int, multiplicity::Int)
    basis = Basis(atoms, basis_str)
    return Mole(atoms, basis, charge, multiplicity)
end

function Mole(atoms, basis)
    return Mole(atoms, basis, 0, 1)
end

abstract type AbstractEnv end

struct Env <: AbstractEnv
    pcharges::PCharges
end