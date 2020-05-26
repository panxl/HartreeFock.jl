function two_electron_fock_matrix!(G::Matrix{Float64}, P::Matrix{Float64}, T::AbstractArray{Float64,4})
    n = LinearAlgebra.checksquare(G)
    for i=1:n, j=i:n
        G[i,j] = 0.0
        for k=1:n, l=1:n
            G[i,j] += (2 * T[i,j,k,l] - T[i,k,j,l]) * P[k,l]
        end
        G[j,i] = G[i,j]
    end
    return nothing
end

function two_electron_fock_matrix(P::Matrix{Float64}, T::AbstractArray{Float64,4})
    G = similar(P)
    two_electron_fock_matrix!(G, P, T)
    return G
end

function density_matrix!(P::Matrix{Float64}, C::Matrix{Float64}, N::Int)
    n = LinearAlgebra.checksquare(P)
    for i=1:n, j=1:n
        P[i,j] = 0.0
        for a=1:N
            P[i,j] += C[i,a] * C[j,a]
        end
    end
    return nothing
end

function density_matrix(C::Matrix{Float64}, N::Int)
    P = similar(C)
    density_matrix!(P, C, N)
    return P
end

function scf(
        S::Matrix{Float64},
        T::Matrix{Float64},
        V::Matrix{Float64},
        int2e::AbstractArray{Float64,4},
        P0::Matrix{Float64},
        N::Int,
        e_tol::Float64=1e-7,
        d_tol::Float64=1e-7,
        max_cycle::Int=64,
    )
    n = size(S, 1)
    Hcore = T + V

    # Guess Fock matrix
    P = copy(P0)
    G = two_electron_fock_matrix(P, int2e)
    F = Hcore + G

    # Guess electron energy
    Eel = P ⋅ (Hcore + F)
    @info "Cycle 0: Eel = $(Eel)"

    diis = DIIS(n)

    for cycle = 1:max_cycle
        # Update density matrix
        e, C = eigen(F, S)
        density_matrix!(P, C, N)

        # Update Fock matrix
        two_electron_fock_matrix!(G, P, int2e)
        F = Hcore + G

        # Calculate electron energy
        Eold = Eel
        Eel = P ⋅ (Hcore + F)

        # Calculate DIIS error
        SPF = S * P * F
        residual = SPF' - SPF
        diis_err = norm(residual) / sqrt(length(residual))

        # Test convergence
        if abs(Eel - Eold) < e_tol && diis_err < d_tol
            @info "Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err) Converged!"
            return Eel, P, e, C
        else
            @info "Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err)"
        end

        # Build DIIS Fock matrix
        diis.add_trial(F)
        diis.add_residual(residual)
        F = diis.get_F()
    end
    @error "SCF failed after $(max_cycle) cycles"
end

struct SCF
    mole::Mole
    S::Matrix{Float64}
    T::Matrix{Float64}
    V::Matrix{Float64}
    int2e::Array{Float64,4}
    P::Matrix{Float64}
    C::Matrix{Float64}
    e::Vector{Float64}
    Eel::Float64
    Enuc::Float64
end

function SCF(mole::Mole)
    intor = CIntor(mole)
    return SCF(mole, intor)
end

function SCF(mole::Mole, intor::AbstractIntor)
    S = intor("int1e_ovlp")
    T = intor("int1e_kin")
    V = intor("int1e_nuc")
    int2e = intor("int2e")
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P, e, C = scf(S, T, V, int2e, P0, N)
    Enuc = electrostatic_energy(mole.nuclei)
    return SCF(mole, S, T, V, int2e, P, C, e, Eel, Enuc)
end

function SCF(mole::Mole, env::Env)
    intor = Intor(mole)
    S = intor("int1e_ovlp")
    T = intor("int1e_kin")
    V = intor("int1e_nuc") + intor("int1e_nuc", env.pointcharges)
    int2e = intor("int2e")
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    Eel, P, e, C = scf(S, T, V, int2e, P0, N)
    Enuc = electrostatic_energy(mole.nuclei) + electrostatic_energy(mole.nuclei, env.pointcharges)
    return SCF(mole, S, T, V, int2e, P, C, e, Eel, Enuc)
end

function total_energy(scf::SCF)
    return scf.Eel + scf.Enuc
end

function total_energy(mole::Mole)
    Eel = electron_energy(mole)
    Enuc = nuclear_repulsion_energy(mole)
    return Eel + Enuc
end

function total_energy(mole::Mole, env::Env)
    scf = SCF(mole, env)
    return scf.Eel + scf.Enuc
end

function electron_energy(mole::Mole)
    intor = CIntor(mole)
    S = intor("int1e_ovlp")
    T = intor("int1e_kin")
    V = intor("int1e_nuc")
    int2e = intor("int2e")
    P0 = zeros(size(S))
    N = div(sum(mole.nuclei.numbers) + mole.net_charge, 2)
    return electron_energy(S, T, V, int2e, P0, N)
end

function electron_energy(
    S::Matrix{Float64},
    T::Matrix{Float64},
    V::Matrix{Float64},
    int2e::AbstractArray{Float64,4},
    P0::Matrix{Float64},
    N::Int,
    )
    Eel, P, e, C = scf(S, T, V, int2e, P0, N)
    return Eel
end

function nuclear_repulsion_energy(mole::Mole)
    return electrostatic_energy(mole.nuclei)
end