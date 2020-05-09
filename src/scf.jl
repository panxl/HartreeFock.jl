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
    M = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = overlap(basis[i], basis[j])
    end
    return copy(M)
end

overlap_matrix(mole::Mole) = overlap_matrix(mole.basis)

function kinetic_matrix(basis::Basis)
    n = length(basis)
    M = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = kinetic(basis[i], basis[j])
    end
    return copy(M)
end

kinetic_matrix(mole::Mole) = kinetic_matrix(mole.basis)

function rinv_matrix(
    basis::Basis,
    RC::Vec3{Float64},
    )
    n = length(basis)
    M = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = rinv(basis[i], basis[j], RC)
    end
    return copy(M)
end

function rinv_matrix(
    basis::Basis,
    p::AbstractParticle,
    )
    return rinv_matrix(basis, p.position)
end

function nuclear_attraction_matrix(
    basis::Basis,
    RC::Vec3{Float64},
    Z::Float64
    )
    n = length(basis)
    M = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for (i,j) in pairs(n)
        M[i, j] = M[j, i] = -Z * rinv(basis[i], basis[j], RC)
    end
    return copy(M)
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
    M = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for i in eachindex(M)
        M[i] = 0.0
    end
    for (i,j) in pairs(n)
        for (rc, z) in zip(RC, Z)
            M[i, j] += -z * rinv(basis[i], basis[j], rc)
        end
        M[j, i] = M[i, j]
    end
    return copy(M)
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
    tensor = Zygote.Buffer(Vector{Float64}(undef, N))
    for (i,j,k,l) in basis_iterator(n)
        tensor[basis_index(i,j,k,l)] = electron_repulsion(basis[i], basis[j], basis[k], basis[l])
    end
    return copy(tensor)
end

electron_repulsion_tensor(mole::Mole) = electron_repulsion_tensor(mole.basis)

function twoe_fock_matrix(P::Matrix{Float64}, T::Vector{Float64})
    n = size(P, 1)
    G = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for i in eachindex(G)
        G[i] = 0.0
    end
    for (i,j) in pairs(n)
        for k=1:n, l=1:n
            G[i,j] = G[j,i] += (2 * T[basis_index(i,j,k,l)] - T[basis_index(i,k,j,l)]) * P[k, l]
        end
    end
    return copy(G)
end

twoe_fock_matrix(P::Matrix{Float64}, basis::Basis) = twoe_fock_matrix(P, electron_repulsion_tensor(basis))
twoe_fock_matrix(P::Matrix{Float64}, mole::Mole) = twoe_fock_matrix(P, mole.basis)

function density_matrix(C::Matrix{Float64}, N::Int)
    n = size(C, 1)
    P = Zygote.Buffer(Matrix{Float64}(undef, n, n))
    for i in eachindex(P)
        P[i] = 0.0
    end
    for i = 1:n, j = 1:n
        for a = 1:N
            P[i,j] += C[i,a] * C[j,a]
        end
    end
    return copy(P)
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
    e = zeros(n)
    C = zeros(n, n)

    λ, U = eigen(Symmetric(S))

    v = Zygote.Buffer(λ)
    for i in eachindex(v)
        v[i] = sqrt(1 / λ[i])
    end
    v = copy(v)

    Sp = U * Diagonal(v) * U

    for cycle = 1:max_cycle
        G = twoe_fock_matrix(P, int2e)
        F = Hcore + G
        # e, C = eigen(F, S)

        Fp = Symmetric(Sp * F * Sp)
        e, Cp = eigen(Fp)
        C = Sp * Cp

        P = density_matrix(C, N)

        Eold = Eel
        Eel = 0.0
        for i=1:n, j=1:n
            Eel += P[i,j] * (Hcore[i,j] + F[i,j])
        end

        println("Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel)")

        if abs(Eel - Eold) < dropgrad(conv_tol)
            println("SCF converged: Eel = $(Eel)")
            return Eel, 2 .* copy(P), copy(e), copy(C)
        end
    end
    println("SCF failed after $(max_cycle) steps")
end

struct SCF
    mole::Mole
    S::Matrix{Float64}
    T::Matrix{Float64}
    V::Matrix{Float64}
    int2e::Vector{Float64}
    P::Matrix{Float64}
    C::Matrix{Float64}
    e::Vector{Float64}
    Eel::Float64
    Enuc::Float64
end

function SCF(mole::Mole)
    S = overlap_matrix(mole)
    T = kinetic_matrix(mole)
    V = nuclear_attraction_matrix(mole)
    int2e = electron_repulsion_tensor(mole)
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P, e, C = scf(S, T, V, int2e, P0, N)
    Enuc = nuclear_repulsion(mole)
    return SCF(mole, S, T, V, int2e, P, C, e, Eel, Enuc)
end

function SCF(mole::Mole, env::Env)
    S = overlap_matrix(mole)
    T = kinetic_matrix(mole)
    V = nuclear_attraction_matrix(mole) + nuclear_attraction_matrix(mole, env)
    int2e = electron_repulsion_tensor(mole)
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P, e, C = scf(S, T, V, int2e, P0, N)
    Enuc = nuclear_repulsion(mole) + nuclear_repulsion(mole, env)
    return SCF(mole, S, T, V, int2e, P, C, e, Eel, Enuc)
end

function total_energy(scf::SCF)
    return scf.Eel + scf.Enuc
end

function total_energy(mole::Mole)
    scf = SCF(mole)
    return scf.Eel + scf.Enuc
end

function total_energy(mole::Mole, env::Env)
    scf = SCF(mole, env)
    return scf.Eel + scf.Enuc
end

function electron_energy(mole::Mole)
    scf = SCF(mole)
    return scf.Eel
end

function nuclear_repulsion_energy(mole::Mole)
    scf = SCF(mole)
    return scf.Enuc
end