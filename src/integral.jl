function E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    if t < 0 || t > (i + j)
        # out of bounds for t
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

function overlap(
        a::Float64,
        RA::Vec3{Float64},
        LA::Vec3{Int},
        b::Float64,
        RB::Vec3{Float64},
        LB::Vec3{Int},
    )
    p = a + b
    RAB = RA - RB
    Sx = E(a, b, RAB.x, LA.x, LB.x, 0)
    Sy = E(a, b, RAB.y, LA.y, LB.y, 0)
    Sz = E(a, b, RAB.z, LA.z, LB.z, 0)
    return (π / p)^(1.5) * Sx * Sy * Sz
end

overlap(A::GTO, B::GTO) = A.norm * B.norm * overlap(A.ζ, A.R, A.L, B.ζ, B.R, B.L)
overlap(A::CGTO, B::CGTO) = contract(overlap, A, B)

function K(a::Float64, b::Float64, Rx::Float64, i::Int64, j::Int64)
    if i == j == 0
        return 2*a*b*E(a,b,Rx,1,1,0)
    end

    if j == 0
        return -i*b*E(a,b,Rx,i-1,1,0) + 2*a*b*E(a,b,Rx,i+1,1,0)
    elseif i == 0
        return -a*j*E(a,b,Rx,1,j-1,0) + 2*a*b*E(a,b,Rx,1,j+1,0)
    else
        return 0.5*i*j*E(a,b,Rx,i-1,j-1,0) -
               i*b*E(a,b,Rx,i-1,j+1,0)     -
               a*j*E(a,b,Rx,i+1,j-1,0)     +
               2*a*b*E(a,b,Rx,i+1,j+1,0)
    end
end

function kinetic(
        a::Float64,
        RA::Vec3{Float64},
        LA::Vec3{Int},
        b::Float64,
        RB::Vec3{Float64},
        LB::Vec3{Int},
    )
    p = a + b
    RAB = RA - RB
    Tx = K(a, b, RAB.x, LA.x, LB.x) * E(a, b, RAB.y, LA.y, LB.y) * E(a, b, RAB.z, LA.z, LB.z)
    Ty = E(a, b, RAB.x, LA.x, LB.x) * K(a, b, RAB.y, LA.y, LB.y) * E(a, b, RAB.z, LA.z, LB.z)
    Tz = E(a, b, RAB.x, LA.x, LB.x) * E(a, b, RAB.y, LA.y, LB.y) * K(a, b, RAB.z, LA.z, LB.z)
    return (π / p)^(1.5) * (Tx + Ty + Tz)
end

kinetic(A::GTO, B::GTO) = A.norm * B.norm * kinetic(A.ζ, A.R, A.L, B.ζ, B.R, B.L)
kinetic(A::CGTO, B::CGTO) = contract(kinetic, A, B)

function R(t::Int, u::Int, v::Int, n::Int, p::Float64, RPC::Vec3{Float64})
    if t == u == v == 0
        return (-2 * p)^n * boys(n, p*norm2(RPC))
    elseif u == v == 0
        if t > 1
            return (t-1) * R(t-2, u, v, n+1, p, RPC) +
                   RPC[1] * R(t-1, u, v, n+1, p, RPC)
        else
            return RPC[1] * R(t-1, u, v, n+1, p, RPC)
        end
    elseif v == 0
        if u > 1
            return (u-1) * R(t, u-2, v, n+1, p, RPC) +
                   RPC[2] * R(t, u-1, v, n+1, p, RPC)
        else
            return RPC[2] * R(t, u-1, v, n+1, p ,RPC)
        end
    else
        if v > 1
            return (v-1) * R(t, u, v-2, n+1, p, RPC) +
                   RPC[3] * R(t, u, v-1, n+1, p, RPC)
        else
            return RPC[3] * R(t, u, v-1, n+1, p, RPC)
        end
    end
end

function nuclear_attraction(
        a::Float64,
        RA::Vec3{Float64},
        LA::Vec3{Int},
        b::Float64,
        RB::Vec3{Float64},
        LB::Vec3{Int},
        RC::Vec3{Float64},
    )
    p = a + b
    RP = gaussian_product_center(a, RA, b, RB)
    RAB = RA - RB
    RPC = RP - RC

    val = 0.0
    for t = 0:LA.x + LB.x, u = 0:LA.y + LB.y, v = 0:LA.z + LB.z
        val += E(a, b, RAB.x, LA.x, LB.x, t) *
               E(a, b, RAB.y, LA.y, LB.y, u) *
               E(a, b, RAB.z, LA.z, LB.z, v) *
               R(t, u, v, 0, p, RPC)
    end

    return -2 * π / p * val
end

nuclear_attraction(A::GTO, B::GTO, RC::Vec3{Float64}) = A.norm * B.norm * nuclear_attraction(A.ζ, A.R, A.L, B.ζ, B.R, B.L, RC)
nuclear_attraction(A::GTO, B::GTO, RC::Vec3{Float64}, Z::Float64) = Z * A.norm * B.norm * nuclear_attraction(A.ζ, A.R, A.L, B.ζ, B.R, B.L, RC)
nuclear_attraction(A::CGTO, B::CGTO, RC::Vec3{Float64}) = contract(nuclear_attraction, A, B, RC)
nuclear_attraction(A::CGTO, B::CGTO, RC::Vec3{Float64}, Z::Float64) = Z * contract(nuclear_attraction, A, B, RC)

function electron_repulsion(
        a::Float64,
        RA::Vec3{Float64},
        LA::Vec3{Int},
        b::Float64,
        RB::Vec3{Float64},
        LB::Vec3{Int},
        c::Float64,
        RC::Vec3{Float64},
        LC::Vec3{Int},
        d::Float64,
        RD::Vec3{Float64},
        LD::Vec3{Int},
    )
    RAB = RA - RB
    RCD = RC - RD
    p = a + b
    q = c + d
    ω = p * q / (p + q)
    RP = gaussian_product_center(a, RA, b, RB)
    RQ = gaussian_product_center(c, RC, d, RD)
    RPQ = RP - RQ

    val = 0.0
    for t = 0:LA.x + LB.x, u = 0:LA.y + LB.y, v = 0:LA.z + LB.z
        for τ = 0:LC.x + LD.x, μ = 0:LC.y + LD.y, υ = 0:LC.z + LD.z
            val += E(a, b, RAB.x, LA.x, LB.x, t) *
                   E(a, b, RAB.y, LA.y, LB.y, u) *
                   E(a, b, RAB.z, LA.z, LB.z, v) *
                   E(c, d, RCD.x, LC.x, LD.x, τ) *
                   E(c, d, RCD.y, LC.y, LD.y, μ) *
                   E(c, d, RCD.z, LC.z, LD.z, υ) *
                   R(t+τ, u+μ, v+υ, 0, ω, RPQ) *
                   (-1)^(τ+μ+υ)
        end
    end
    val *= 2 * π^2.5 / (p * q * sqrt(p + q))
    return val
end

electron_repulsion(A::GTO, B::GTO, C::GTO, D::GTO) = A.norm * B.norm * C.norm * D.norm *
    electron_repulsion(A.ζ, A.R, A.L, B.ζ, B.R, B.L, C.ζ, C.R, C.L, D.ζ, D.R, D.L)
electron_repulsion(A::CGTO, B::CGTO, C::CGTO, D::CGTO) = contract(electron_repulsion, A, B, C, D)

function contract(f, A::CGTO, B::CGTO)
    s = 0.0
    for i = 1:length(A), j = 1:length(B)
        s += A.d[i] * B.d[j] * f(A[i], B[j])
    end
    return A.norm * B.norm * s
end

function contract(f, A::CGTO, B::CGTO, RC::Vec3{Float64})
    s = 0.0
    for i = 1:length(A), j = 1:length(B)
        s += A.d[i] * B.d[j] * f(A[i], B[j], RC)
    end
    return A.norm * B.norm * s
end

function contract(f, A::CGTO, B::CGTO, C::CGTO, D::CGTO)
    s = 0.0
    for i = 1:length(A), j = 1:length(B), k = 1:length(C), l = 1:length(D)
        s += A.d[i] * B.d[j] * C.d[k] * D.d[l] * f(A[i], B[j], C[k], D[l])
    end
    return A.norm * B.norm * C.norm * D.norm * s
end