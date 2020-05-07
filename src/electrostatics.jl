
#####
##### Electrostatics
#####

function electrostatic_potential(
    Q::Float64,
    d::Float64,
    )
    return Q / d
end

function electrostatic_potential(
    Q::Float64,
    ra::Vec3{Float64},
    rb::Vec3{Float64},
    )
    d = norm(ra - rb)
    return electrostatic_potential(Q, d)
end

function electrostatic_potential(
    Q::Float64,
    σ::Float64,
    d::Float64,
    )
    return Q  / d * erf(d / (sqrt(2) * σ))
end

function electrostatic_potential(
    Q::Float64,
    σ::Float64,
    ra::Vec3{Float64},
    rb::Vec3{Float64},
    )
    d = norm(ra - rb)
    return electrostatic_potential(Q, σ, d)
end

function electrostatic_potential(
    pa::AbstractParticle,
    pb::AbstractParticle,
    )
    return electrostatic_potential(pa.charge, pa.position, pb.position)
end

function electrostatic_potential(
    pa::GaussianCharge,
    pb::AbstractParticle,
    )
    return electrostatic_potential(pa.charge, pa.sigma, pa.position, pb.position)
end

function electrostatic_potential(
    ga::AbstractParticleGroup,
    pb::AbstractParticle,
    )
    ret = 0.0
    for pa in ga
        ret += electrostatic_potential(pa, pb)
    end
    return ret
end

function electrostatic_potential(
    mole::Mole,
    pb::AbstractParticle,
    )
    Vnuc = electrostatic_potential(mole.nuclei, pb)
    Vel = sum(rinv_matrix(mole.basis, pb) .* mole.density_matrix)
    return Vnuc - Vel
end

function electrostatic_potential(
    mole::Union{AbstractParticle,AbstractParticleGroup,Mole},
    gb::AbstractParticleGroup,
    )
    ret = zeros(length(gb))
    for (i, pb) in enumerate(gb)
        ret[i] = electrostatic_potential(mole, pb)
    end
    return ret
end