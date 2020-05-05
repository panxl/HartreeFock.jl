#####
##### Abstract Particle and Particle Group
#####

abstract type AbstractParticle end

position(p::AbstractParticle) = p.position

abstract type AbstractParticleGroup end

position(g::AbstractParticleGroup) = g.positions

Base.iterate(g::AbstractParticleGroup, state=1) = state > length(g) ? nothing : (g[state], state+1)
Base.length(g::AbstractParticleGroup) = length(g.positions)
Base.size(g::AbstractParticleGroup) = length(g)

Base.firstindex(g::AbstractParticleGroup) = 1
Base.lastindex(g::AbstractParticleGroup) = length(g)

#####
##### Nucleus and Nuclei
#####

struct Nucleus <: AbstractParticle
    position::Vec3{Float64}
    charge::Float64
    number::Int
end

function Nucleus(
    position::AbstractVector{Float64},
    number::Int,
    )
    Nucleus(position, Float64(number), number)
end

function Nucleus(
    p::AbstractParticle,
    )
    Nucleus(p.position, p.number)
end

struct Nuclei <: AbstractParticleGroup
    positions::Vector{Vec3{Float64}}
    charges::Vector{Float64}
    numbers::Vector{Int}
end

function Nuclei(
    positions::Vector{Vec3{Float64}},
    numbers::Vector{Int}
    )
    return Nuclei(positions, Float64.(numbers), numbers)
end

Nuclei() = Nuclei(Vector{Vec3{Float64}}[], Vector{Float64}[], Vector{Int}[])
Nuclei(file::String) = Nuclei(parse_geom(file)...)

function Base.push!(g::Nuclei, p::Nucleus)
    push!(g.positions, p.position)
    push!(g.charges, p.charge)
    push!(g.numbers, p.number)
end

Base.getindex(g::Nuclei, i::Int) = Nucleus(g.positions[i], g.charges[i], g.numbers[i])

#####
##### Point Charge and Point Charges
#####

struct PointCharge <: AbstractParticle
    position::Vec3{Float64}
    charge::Float64
end

function PointCharge(
    p::AbstractParticle,
    )
    return PointCharge(p.position, p.charge)
end

function PointCharge(
    p::AbstractParticle,
    charge::Float64
    )
    return PointCharge(p.position, charge)
end

struct PointCharges <: AbstractParticleGroup
    positions::Vector{Vec3{Float64}}
    charges::Vector{Float64}
end

PointCharges(g::AbstractParticleGroup) = PointCharges(g.positions, g.charges)
PointCharges() = PointCharges(Vector{Vec3{Float64}}[], Vector{Float64}[])
PointCharges(file::String) = PointCharges(parse_pc(file)...)

function Base.push!(g::PointCharges, p::PointCharge)
    push!(g.positions, p.position)
    push!(g.charges, p.charge)
end

Base.getindex(g::PointCharges, i::Int) = PointCharge(g.positions[i], g.charges[i])

#####
##### Gaussian Charge and Gaussian Charges
#####

struct GaussianCharge <: AbstractParticle
    position::Vec3{Float64}
    charge::Float64
    sigma::Float64
end

function GaussianCharge(
    p::AbstractParticle,
    sigma::Float64,
    )
    return GaussianCharge(p.position, p.charge, sigma)
end

struct GaussianCharges <: AbstractParticleGroup
    positions::Vector{Vec3{Float64}}
    charges::Vector{Float64}
    sigmas::Vector{Float64}
end

GaussianCharges(g::AbstractParticleGroup, sigmas::Vector{Float64}) = GaussianCharges(g.positions, g.charges, sigmas)
GaussianCharges() = GaussianCharges(Vector{Vec3{Float64}}[], Vector{Float64}[], Vector{Float64}[])

function Base.push!(g::GaussianCharges, p::GaussianCharge)
    push!(g.positions, p.position)
    push!(g.charges, p.charge)
    push!(g.sigmas, p.sigma)
end

Base.getindex(g::GaussianCharges, i::Int) = GaussianCharge(g.positions[i], g.charges[i], g.sigmas[i])

#####
##### Point and Points
#####

struct Point <: AbstractParticle
    position::Vec3{Float64}
end

struct Points <: AbstractParticleGroup
    positions::Vector{Vec3{Float64}}
end

Points(g::AbstractParticleGroup) = Points(g.positions)
Points() = Points(Vector{Vec3{Float64}}[])

function Base.push!(g::Points, p::Point)
    push!(g.positions, p.position)
end

Base.getindex(g::Points, i::Int) = Point(g.positions[i])

#####
##### Typed Atom and Typed Atoms
#####

struct TypedAtom <: AbstractParticle
    position::Vec3{Float64}
    charge::Float64
    type::Int
end

struct TypedAtoms <: AbstractParticleGroup
    positions::Vector{Vec3{Float64}}
    charges::Vector{Float64}
    types::Vector{Int}
end

TypedAtoms() = TypedAtom(Vector{Vec3{Float64}}[], Vector{Float64}[], Vector{Int}[])

function Base.push!(g::TypedAtoms, p::TypedAtom)
    push!(g.positions, p.position)
    push!(g.charges, p.charge)
    push!(g.types, p.type)
end

Base.getindex(g::TypedAtoms, i::Int) = TypedAtom(g.positions[i], g.charges[i], g.types[i])