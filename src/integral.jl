function E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    if i < 0 || j < 0 || t < 0 || t > (i + j)
        # out of bounds for i, j, or t
        return 0.0
    end

    p = a + b

    if i == j == t == 0
        return exp(-a*b/p*Rx*Rx)
    end

    if j == 0
        return (0.5/p)*E(a,b,Rx,i-1,j,t-1) -
               (b/p*Rx)*E(a,b,Rx,i-1,j,t)  +
               (t+1)*E(a,b,Rx,i-1,j,t+1)
    else
        return (0.5/p)*E(a,b,Rx,i,j-1,t-1) +
               (a/p*Rx)*E(a,b,Rx,i,j-1,t)  +
               (t+1)*E(a,b,Rx,i,j-1,t+1)
    end
end

@adjoint function E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    return E(a, b, Rx, i, j, t), function (Δ)
        x̄ = a*E(a,b,Rx,i+1,j)     - 
            0.5*i*E(a,b,Rx,i-1,j) -
            b*E(a,b,Rx,i,j+1)     +
            0.5*j*E(a,b,Rx,i,j-1)
        return nothing, nothing, x̄*Δ, nothing, nothing, nothing
    end
end

# ∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int) = gradient(E, a, b, Rx, i, j)[3]

function ∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int)
    return a*E(a,b,Rx,i+1,j)     - 
           0.5*i*E(a,b,Rx,i-1,j) -
           b*E(a,b,Rx,i,j+1)     +
           0.5*j*E(a,b,Rx,i,j-1)
end

# ∇∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int) = gradient(∇E, a, b, Rx, i, j)[3]

function ∇∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int)
    return -i*j*E(a,b,Rx,i-1,j-1,0)   +
            2*i*b*E(a,b,Rx,i-1,j+1,0) +
            2*a*j*E(a,b,Rx,i+1,j-1,0) -
            4*a*b*E(a,b,Rx,i+1,j+1,0)
end

∇∇∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int) = gradient(∇∇E, a, b, Rx, i, j)[3]

function normalizaton(
    a::Float64,
    L::Vec3{Int},
    )
    return 1 / sqrt(overlap(a, L, a, L, Vec3(zeros(3))))
end

normalizaton(gto::GTO) = normalizaton(gto.ζ, gto.L)

function normalizaton(
    d::Vector{Float64},
    gtos::Vector{GTO},
    )
    N = length(d)
    ret = 0.0
    for i = 1:N, j = 1:N
        ret += d[i] * d[j] * overlap(gtos[i], gtos[j], Vec3(zeros(3)))
    end
    return ret
end

normalizaton(cgto::CGTO) = normalizaton(cgto.d, cgto.gtos)

function overlap(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB[1], LA[1], LB[1])
    Sy = E(a, b, RAB[2], LA[2], LB[2])
    Sz = E(a, b, RAB[3], LA[3], LB[3])
    return (PI / p)^(1.5) * Sx * Sy * Sz
end

@adjoint function overlap(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB[1], LA[1], LB[1])
    Sy = E(a, b, RAB[2], LA[2], LB[2])
    Sz = E(a, b, RAB[3], LA[3], LB[3])
    return (PI / p)^(1.5) * Sx * Sy * Sz, function (Δ)
        ∇Sx = ∇E(a, b, RAB[1], LA[1], LB[1])
        ∇Sy = ∇E(a, b, RAB[2], LA[2], LB[2])
        ∇Sz = ∇E(a, b, RAB[3], LA[3], LB[3])
        gradient = zeros(3)
        gradient[1] = (PI / p)^(1.5) * ∇Sx * Sy * Sz
        gradient[2] = (PI / p)^(1.5) * Sx * ∇Sy * Sz
        gradient[3] = (PI / p)^(1.5) * Sx * Sy * ∇Sz
        return nothing, nothing, nothing, nothing, gradient*Δ
    end
end

function overlap(
    A::GTO,
    B::GTO,
    RAB::Vec3{Float64},
    )
    return A.norm * B.norm * overlap(A.ζ, A.L, B.ζ, B.L, RAB)
end

overlap(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}) = contract(overlap, A, RA, B, RB)

function kinetic(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB[1], LA[1], LB[1])
    Sy = E(a, b, RAB[2], LA[2], LB[2])
    Sz = E(a, b, RAB[3], LA[3], LB[3])
    ∇∇Sx = ∇∇E(a, b, RAB[1], LA[1], LB[1])
    ∇∇Sy = ∇∇E(a, b, RAB[2], LA[2], LB[2])
    ∇∇Sz = ∇∇E(a, b, RAB[3], LA[3], LB[3])
    return -0.5 * (PI / p)^(1.5) * (∇∇Sx * Sy * Sz + Sx * ∇∇Sy * Sz + Sx * Sy * ∇∇Sz)
end

@adjoint function kinetic(
    a::Float64,
    LA::Vec3{Int},
    b::Float64,
    LB::Vec3{Int},
    RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB[1], LA[1], LB[1])
    Sy = E(a, b, RAB[2], LA[2], LB[2])
    Sz = E(a, b, RAB[3], LA[3], LB[3])
    ∇∇Sx = ∇∇E(a, b, RAB[1], LA[1], LB[1])
    ∇∇Sy = ∇∇E(a, b, RAB[2], LA[2], LB[2])
    ∇∇Sz = ∇∇E(a, b, RAB[3], LA[3], LB[3])
    return -0.5 * (PI / p)^(1.5) * (∇∇Sx * Sy * Sz + Sx * ∇∇Sy * Sz + Sx * Sy * ∇∇Sz), function (Δ)
        ∇Sx = ∇E(a, b, RAB[1], LA[1], LB[1])
        ∇Sy = ∇E(a, b, RAB[2], LA[2], LB[2])
        ∇Sz = ∇E(a, b, RAB[3], LA[3], LB[3])
        ∇∇∇Sx = ∇∇∇E(a, b, RAB[1], LA[1], LB[1])
        ∇∇∇Sy = ∇∇∇E(a, b, RAB[2], LA[2], LB[2])
        ∇∇∇Sz = ∇∇∇E(a, b, RAB[3], LA[3], LB[3])
        gradient = zeros(3)
        gradient[1] = -0.5 * (PI / p)^(1.5) * (∇∇∇Sx * Sy * Sz + ∇Sx * ∇∇Sy * Sz + ∇Sx * Sy * ∇∇Sz)
        gradient[2] = -0.5 * (PI / p)^(1.5) * (∇∇Sx * ∇Sy * Sz + Sx * ∇∇∇Sy * Sz + Sx * ∇Sy * ∇∇Sz)
        gradient[3] = -0.5 * (PI / p)^(1.5) * (∇∇Sx * Sy * ∇Sz + Sx * ∇∇Sy * ∇Sz + Sx * Sy * ∇∇∇Sz)
        return nothing, nothing, nothing, nothing, gradient*Δ
    end
end

function kinetic(
    A::GTO,
    B::GTO,
    RAB::Vec3{Float64},
    )
    return A.norm * B.norm * kinetic(A.ζ, A.L, B.ζ, B.L, RAB)
end

kinetic(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}) = contract(kinetic, A, RA, B, RB)

function R(t::Int, u::Int, v::Int, n::Int, p::Float64, RPC::Vec3{Float64})
    if t == u == v == 0
        return (-2 * p)^n * boys(n, p*norm2(RPC))
    elseif u == v == 0
        if t > 1
            return (t-1) * R(t-2, 0, 0, n+1, p, RPC) +
                   RPC[1] * R(t-1, 0, 0, n+1, p, RPC)
        else
            return RPC[1] * R(0, 0, 0, n+1, p, RPC)
        end
    elseif v == 0
        if u > 1
            return (u-1) * R(t, u-2, 0, n+1, p, RPC) +
                   RPC[2] * R(t, u-1, 0, n+1, p, RPC)
        else
            return RPC[2] * R(t, 0, 0, n+1, p ,RPC)
        end
    else
        if v > 1
            return (v-1) * R(t, u, v-2, n+1, p, RPC) +
                   RPC[3] * R(t, u, v-1, n+1, p, RPC)
        else
            return RPC[3] * R(t, u, 0, n+1, p, RPC)
        end
    end
end

function rinv(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
        RPC::Vec3{Float64},
    )
    p = a + b
    val = 0.0
    for t = 0:LA[1] + LB[1], u = 0:LA[2] + LB[2], v = 0:LA[3] + LB[3]
        val += E(a, b, RAB[1], LA[1], LB[1], t) *
               E(a, b, RAB[2], LA[2], LB[2], u) *
               E(a, b, RAB[3], LA[3], LB[3], v) *
               R(t, u, v, 0, p, RPC)
    end
    return 2 * PI / p * val
end

function rinv(
    A::GTO,
    RA::Vec3{Float64},
    B::GTO,
    RB::Vec3{Float64},
    RC::Vec3{Float64},
    )
    RAB = RA - RB
    RP = gaussian_product_center(A.ζ, RA, B.ζ, RB)
    RPC = RP - RC
    return A.norm * B.norm * rinv(A.ζ, A.L, B.ζ, B.L, RAB, RPC)
end

rinv(A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, RC::Vec3{Float64}) = contract(rinv, A, RA, B, RB, RC)

function electron_repulsion(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        c::Float64,
        LC::Vec3{Int},
        d::Float64,
        LD::Vec3{Int},
        RAB::Vec3{Float64},
        RCD::Vec3{Float64},
        RPQ::Vec3{Float64},
    )
    p = a + b
    q = c + d
    ω = p * q / (p + q)
    val = 0.0
    for t = 0:LA[1] + LB[1], u = 0:LA[2] + LB[2], v = 0:LA[3] + LB[3]
        for τ = 0:LC[1] + LD[1], μ = 0:LC[2] + LD[2], υ = 0:LC[3] + LD[3]
            val += E(a, b, RAB[1], LA[1], LB[1], t) *
                   E(a, b, RAB[2], LA[2], LB[2], u) *
                   E(a, b, RAB[3], LA[3], LB[3], v) *
                   E(c, d, RCD[1], LC[1], LD[1], τ) *
                   E(c, d, RCD[2], LC[2], LD[2], μ) *
                   E(c, d, RCD[3], LC[3], LD[3], υ) *
                   R(t+τ, u+μ, v+υ, 0, ω, RPQ) *
                   (-1)^(τ+μ+υ)
        end
    end
    val *= 2 * PI^2.5 / (p * q * sqrt(p + q))
    return val
end

function electron_repulsion(
    A::GTO,
    RA::Vec3{Float64},
    B::GTO,
    RB::Vec3{Float64},
    C::GTO,
    RC::Vec3{Float64},
    D::GTO,
    RD::Vec3{Float64},
    )
    RAB = RA - RB
    RCD = RC - RD
    RP = gaussian_product_center(A.ζ, RA, B.ζ, RB)
    RQ = gaussian_product_center(C.ζ, RC, D.ζ, RD)
    RPQ = RP - RQ
    return A.norm * B.norm * C.norm * D.norm *
        electron_repulsion(A.ζ, A.L, B.ζ, B.L, C.ζ, C.L, D.ζ, D.L, RAB, RCD, RPQ)
end

function electron_repulsion(
    A::CGTO,
    RA::Vec3{Float64},
    B::CGTO,
    RB::Vec3{Float64},
    C::CGTO,
    RC::Vec3{Float64},
    D::CGTO,
    RD::Vec3{Float64},
    )
    contract(electron_repulsion, A, RA, B, RB, C, RC, D, RD)
end

function contract(f, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64})
    RAB = RA - RB
    ret = 0.0
    for i = 1:length(A), j = 1:length(B)
        ret += A.d[i] * B.d[j] * f(A[i], B[j], RAB)
    end
    return A.norm * B.norm * ret
end

function contract(f, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, RC::Vec3{Float64})
    ret = 0.0
    for i = 1:length(A), j = 1:length(B)
        ret += A.d[i] * B.d[j] * f(A[i], RA, B[j], RB, RC)
    end
    return A.norm * B.norm * ret
end

function contract(f, A::CGTO, RA::Vec3{Float64}, B::CGTO, RB::Vec3{Float64}, C::CGTO, RC::Vec3{Float64}, D::CGTO, RD::Vec3{Float64})
    ret = 0.0
    for i = 1:length(A), j = 1:length(B), k = 1:length(C), l = 1:length(D)
        ret += A.d[i] * B.d[j] * C.d[k] * D.d[l] * 
            f(A[i], RA, B[j], RB, C[k], RC, D[l], RD)
    end
    return A.norm * B.norm * C.norm * D.norm * ret
end