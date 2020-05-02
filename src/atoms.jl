abstract type AbstractParticle end
abstract type AbstractParticles end

Base.iterate(particles::AbstractParticles, state=1) = state > length(particles) ? nothing : (particles[state], state+1)
Base.length(particles::AbstractParticles) = length(particles.positions)
Base.size(particles::AbstractParticles) = length(particles)

Base.firstindex(particles::AbstractParticles) = 1
Base.lastindex(particles::AbstractParticles) = length(particles)

struct Atom <: AbstractParticle
    position::Vec3{Float64}
    number::Int
    charge::Float64
end

Atom(position::AbstractVector{Float64}, number::Int) = Atom(position, number, Float64(number))

struct Atoms <: AbstractParticles
    positions::Vector{Vec3{Float64}}
    numbers::Vector{Int}
    charges::Vector{Float64}
end

Atoms() = Atoms(Vector{Vec3{Float64}}[], Vector{Int}[], Vector{Float64}[])
Atoms(file::String) = Atoms(parser(file)...)

function Base.push!(atoms::Atoms, atom::Atom)
    push!(atoms.positions, atom.position)
    push!(atoms.numbers, atom.number)
    push!(atoms.charges, atom.charge)
end

Base.getindex(atoms::Atoms, i::Int) = Atom(atoms.positions[i], atoms.numbers[i], atoms.charges[i])

struct PCharge <: AbstractParticle
    position::Vec3{Float64}
    charge::Float64
end

struct PCharges <: AbstractParticles
    positions::Vector{Vec3{Float64}}
    charges::Vector{Float64}
end

PCharges() = PCharges(Vector{Vec3{Float64}}[], Vector{Float64}[])

function Base.push!(pcharges::PCharges, pcharge::PCharge)
    push!(pcharges.positions, pcharge.position)
    push!(pcharges.charges, pcharge.charge)
end

Base.getindex(pcharges::PCharges, i::Int) = PCharge(pcharges.positions[i], pcharges.charges[i])

struct Point <: AbstractParticle
    position::Vec3{Float64}
end

struct Points <: AbstractParticles
    positions::Vector{Vec3{Float64}}
end

Points() = Points(Vector{Vec3{Float64}}[])

function Base.push!(points::Points, point::Point)
    push!(points.positions, point.position)
end

Base.getindex(points::Points, i::Int) = PCharge(points.positions[i])