nuclear_repulsion(ZA::Float64, RA::Vec3{Float64}, ZB::Float64, RB::Vec3{Float64}) = ZA * ZB / norm(RA - RB)
nuclear_repulsion(A::AbstractParticle, B::AbstractParticle) = nuclear_repulsion(A.charge, A.position, B.charge, B.position)

function nuclear_repulsion(particles::AbstractParticleGroup)
    ret = 0.0
    for i = 1:length(particles)-1, j = i+1:length(particles)
        ret += nuclear_repulsion(particles[i], particles[j])
    end
    return ret
end

function nuclear_repulsion(particle_a::AbstractParticleGroup, particle_b::AbstractParticleGroup)
    ret = 0.0
    for i = 1:length(particle_a), j = 1:length(particle_b)
        ret += nuclear_repulsion(particle_a[i], particle_b[j])
    end
    return ret
end

nuclear_repulsion(mole::Mole) = nuclear_repulsion(mole.nuclei)
nuclear_repulsion(mole::Mole, env::Env) = nuclear_repulsion(mole.nuclei, env.pointcharges)

function overlap_matrix(basis::Basis)
    n = length(basis)
    M = zeros(n, n)
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = overlap(basis[i], basis[j])
    end
    return M
end

overlap_matrix(mole::Mole) = overlap_matrix(mole.basis)

function kinetic_matrix(basis::Basis)
    n = length(basis)
    M = zeros(n, n)
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = kinetic(basis[i], basis[j])
    end
    return M
end

kinetic_matrix(mole::Mole) = kinetic_matrix(mole.basis)

function nuclear_attraction_matrix(
    basis::Basis,
    RC::Vec3{Float64},
    Z::Float64
    )
    n = length(basis)
    M = zeros(n, n)
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = -Z * rinv(basis[i], basis[j], RC)
    end
    return M
end

function nuclear_attraction_matrix(
    basis::Basis,
    p::AbstractParticle,
    )
    return nuclear_attraction_matrix(basis, p.position, p.charge)
end

function nuclear_attraction_matrix(
    basis::Basis,
    RC::Vector{Vec3{Float64}},
    Z::Vector{Float64}
    )
    n = length(basis)
    M = zeros(n, n)
    for (i,j) in pairs(n)
        for (rc, z) in zip(RC, Z)
            M[i, j] += -z * rinv(basis[i], basis[j], rc)
        end
        M[j, i] = M[i, j]
    end
    return M
end

function nuclear_attraction_matrix(
    basis::Basis,
    g::AbstractParticleGroup,
    )
    return nuclear_attraction_matrix(basis, g.positions, g.charges)
end

nuclear_attraction_matrix(mole::Mole) = nuclear_attraction_matrix(mole.basis, mole.nuclei)

nuclear_attraction_matrix(mole::Mole, env::Env) = nuclear_attraction_matrix(mole.basis, env.pointcharges)

function electron_repulsion_tensor(basis::Basis)
    n = length(basis)
    N = binomial(binomial(n + 1, 2) + 1, 2)
    tensor = zeros(N)
    for (i,j,k,l) in basis_iterator(n)
        tensor[basis_index(i,j,k,l)] = electron_repulsion(basis[i], basis[j], basis[k], basis[l])
    end
    return tensor
end

electron_repulsion_tensor(mole::Mole) = electron_repulsion_tensor(mole.basis)

function twoe_fock_matrix(P::Matrix{Float64}, T::Vector{Float64})
    n = size(P, 1)
    G = zeros(n, n)
    for (i,j) in pairs(n)
        for k=1:n, l=1:n
            G[i,j] = G[j,i] += (2 * T[basis_index(i,j,k,l)] - T[basis_index(i,k,j,l)]) * P[k, l]
        end
    end
    return G
end

twoe_fock_matrix(P::Matrix{Float64}, basis::Basis) = twoe_fock_matrix(P, electron_repulsion_tensor(basis))
twoe_fock_matrix(P::Matrix{Float64}, mole::Mole) = twoe_fock_matrix(P, mole.basis)

function density_matrix(C::Matrix{Float64}, N::Int)
    n = size(C, 1)
    P = zeros(n, n)
    for i = 1:n, j = 1:n
        for a = 1:N
            P[i,j] += C[i,a] * C[j,a]
        end
    end
    return P
end

function scf(
        S::Matrix{Float64},
        T::Matrix{Float64},
        V::Matrix{Float64},
        int2e::Vector{Float64},
        P0::Matrix{Float64},
        N::Int,
        conv_tol::Float64=1e-9,
        max_cycle::Int=64,
    )
    n = size(S, 1)
    Hcore = T + V
    P = copy(P0)
    Eel = 0.0

    for cycle = 1:max_cycle
        G = twoe_fock_matrix(P, int2e)
        F = Hcore + G
        E, C = eigen(F, S)
        P = density_matrix(C, N)

        Eold = Eel
        Eel = 0.0
        for i=1:n, j=1:n
            Eel += P[i,j] * (Hcore[i,j] + F[i,j])
        end

        println("Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel)")

        if abs(Eel - Eold) < conv_tol
            println("SCF converged: Eel = $(Eel)")
            return Eel, 2 .* P
        end
    end
    println("SCF failed after $(max_cycle) steps")
end

function scf(mole::Mole)
    S = overlap_matrix(mole)
    T = kinetic_matrix(mole)
    V = nuclear_attraction_matrix(mole)
    int2e = electron_repulsion_tensor(mole)
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P = scf(S, T, V, int2e, P0, N)
    Enuc = nuclear_repulsion(mole)
    return Eel + Enuc, P
end

function scf(mole::Mole, env::Env)
    S = overlap_matrix(mole)
    T = kinetic_matrix(mole)
    V = nuclear_attraction_matrix(mole) + nuclear_attraction_matrix(mole, env)
    int2e = electron_repulsion_tensor(mole)
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P = scf(S, T, V, int2e, P0, N)
    Enuc = nuclear_repulsion(mole) + nuclear_repulsion(mole, env)
    return Eel + Enuc, P
end