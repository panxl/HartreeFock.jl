#####
##### Electrostatic Energy
#####

function electrostatic_energy(
    qa::Float64,
    qb::Float64,
    d::Float64,
    )
    return qa * qb / d
end

function electrostatic_energy(
    qa::Float64,
    ra::Vec3{Float64},
    qb::Float64,
    rb::Vec3{Float64},
    )
    d = norm(ra - rb)
    return electrostatic_energy(qa, qb, d)
end

function electrostatic_energy(
    a::AbstractParticle,
    b::AbstractParticle,
    )
    return electrostatic_energy(a.charge, a.position, b.charge, b.position)
end

function electrostatic_energy(
    g::AbstractParticleGroup,
    )
    ret = 0.0
    for i=1:length(g)-1, j=i+1:length(g)
        @inbounds ret += electrostatic_energy(g[i], g[j])
    end
    return ret
end

function electrostatic_energy(
    ga::AbstractParticleGroup,
    gb::AbstractParticleGroup,    
    )
    ret = 0.0
    for i=1:length(ga), j=1:length(gb)
        @inbounds ret += electrostatic_energy(ga[i], gb[j])
    end
    return ret
end

#####
##### Electrostatic Potential
#####

function electrostatic_potential(
    q::Float64,
    d::Float64,
    )
    return q / d
end

function electrostatic_potential(
    q::Float64,
    ra::Vec3{Float64},
    rb::Vec3{Float64},
    )
    d = norm(ra - rb)
    return electrostatic_potential(q, d)
end

function electrostatic_potential(
    q::Float64,
    σ::Float64,
    d::Float64,
    )
    return q / d * erf(d / (sqrt(2) * σ))
end

function electrostatic_potential(
    q::Float64,
    σ::Float64,
    ra::Vec3{Float64},
    rb::Vec3{Float64},
    )
    d = norm(ra - rb)
    return electrostatic_potential(q, σ, d)
end

function electrostatic_potential(
    particle::AbstractParticle,
    point::Vec3{Float64},
    )
    return electrostatic_potential(particle.charge, particle.position, point)
end

function electrostatic_potential(
    particle::GaussianCharge,
    point::Vec3{Float64},
    )
    return electrostatic_potential(particle.charge, particle.sigma, particle.position, point)
end

function electrostatic_potential(
    particles::AbstractParticleGroup,
    point::Vec3{Float64},
    )
    ret = 0.0
    for i=1:length(particles)
        @inbounds ret += electrostatic_potential(particles[i], point)
    end
    return ret
end

function electrostatic_potential(
    scf::SCF,
    point::Vec3{Float64},
    )
    Vnuc = electrostatic_potential(scf.mole.nuclei, point)
    intor = Intor(scf.mole)
    Vel = intor("int1e_rinv", point) ⋅ scf.P * 2
    return Vnuc - Vel
end

function electrostatic_potential(
    obj::Union{AbstractParticle,AbstractParticleGroup,SCF},
    points::AbstractVector{Vec3{Float64}},
    )
    ret = zeros(length(points))
    for (i, point) in enumerate(points)
        ret[i] = electrostatic_potential(obj, point)
    end
    return ret
end