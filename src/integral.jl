"Hermite expansion coefficients for Gaussian product"
function E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    if i < 0 || j < 0 || t < 0 || t > (i + j)
        # out of bounds for i, j, or t
        return 0.0
    end

    p = a + b

    if i == j == t == 0 # base case
        return exp(-a*b/p*Rx*Rx)
    end

    if j == 0 # decrement index i
        return (0.5/p)*E(a,b,Rx,i-1,0,t-1) -
               (b/p*Rx)*E(a,b,Rx,i-1,0,t)  +
               (t+1)*E(a,b,Rx,i-1,0,t+1)
    else # decrement index j
        return (0.5/p)*E(a,b,Rx,i,j-1,t-1) +
               (a/p*Rx)*E(a,b,Rx,i,j-1,t)  +
               (t+1)*E(a,b,Rx,i,j-1,t+1)
    end
end

function ∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    if Rx == 0.0
        return 0.0
    end

    return 2*a*E(a,b,Rx,i+1,j,t) - i*E(a,b,Rx,i-1,j,t)
end

function ∇∇E(a::Float64, b::Float64, Rx::Float64, i::Int, j::Int, t::Int=0)
    return -i*j*E(a,b,Rx,i-1,j-1,t)   +
            2*i*b*E(a,b,Rx,i-1,j+1,t) +
            2*a*j*E(a,b,Rx,i+1,j-1,t) -
            4*a*b*E(a,b,Rx,i+1,j+1,t)
end

function int1e_ovlp(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB.x, LA.x, LB.x, 0)
    Sy = E(a, b, RAB.y, LA.y, LB.y, 0)
    Sz = E(a, b, RAB.z, LA.z, LB.z, 0)
    return (π / p)^(1.5) * Sx * Sy * Sz
end

function int1e_kin(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
    )
    p = a + b
    Sx = E(a, b, RAB.x, LA.x, LB.x)
    Sy = E(a, b, RAB.y, LA.y, LB.y)
    Sz = E(a, b, RAB.z, LA.z, LB.z)
    ∇∇Sx = ∇∇E(a, b, RAB.x, LA.x, LB.x)
    ∇∇Sy = ∇∇E(a, b, RAB.y, LA.y, LB.y)
    ∇∇Sz = ∇∇E(a, b, RAB.z, LA.z, LB.z)
    return -0.5 * (π / p)^(1.5) *  (∇∇Sx * Sy * Sz + Sx * ∇∇Sy * Sz + Sx * Sy * ∇∇Sz)
end

"Auxiliary Hermite integrals"
function R(t::Int, u::Int, v::Int, n::Int, p::Float64, RPC::Vec3{Float64})
    if t == u == v == 0 # base case
        return (-2 * p)^n * boys(n, p * (RPC[1]*RPC[1] + RPC[2]*RPC[2] + RPC[3]*RPC[3]))
    elseif u == v == 0 # decrement index t
        if t == 1
            return RPC[1] * R(0, 0, 0, n+1, p, RPC)
        end

        return (t-1) * R(t-2, 0, 0, n+1, p, RPC) +
               RPC[1] * R(t-1, 0, 0, n+1, p, RPC)
    elseif v == 0 # decrement index u
        if u == 1
            return RPC[2] * R(t, 0, 0, n+1, p ,RPC)
        end

        return (u-1) * R(t, u-2, 0, n+1, p, RPC) +
               RPC[2] * R(t, u-1, 0, n+1, p, RPC)
    else # decrement index v
        if v == 1
            return RPC[3] * R(t, u, 0, n+1, p, RPC)
        end

        return (v-1) * R(t, u, v-2, n+1, p, RPC) +
               RPC[3] * R(t, u, v-1, n+1, p, RPC)
    end
end

function int1e_rinv(
        a::Float64,
        LA::Vec3{Int},
        b::Float64,
        LB::Vec3{Int},
        RAB::Vec3{Float64},
        RPC::Vec3{Float64},
    )
    p = a + b
    ret = 0.0
    for t = 0:LA.x + LB.x, u = 0:LA.y + LB.y, v = 0:LA.z + LB.z
        ret += E(a, b, RAB.x, LA.x, LB.x, t) *
               E(a, b, RAB.y, LA.y, LB.y, u) *
               E(a, b, RAB.z, LA.z, LB.z, v) *
               R(t, u, v, 0, p, RPC)
    end
    return 2 * π / p * ret
end

function int2e(
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
    ret = 0.0
    for t = 0:LA.x + LB.x, u = 0:LA.y + LB.y, v = 0:LA.z + LB.z
        for τ = 0:LC.x + LD.x, μ = 0:LC.y + LD.y, υ = 0:LC.z + LD.z
            ret += E(a, b, RAB.x, LA.x, LB.x, t) *
                   E(a, b, RAB.y, LA.y, LB.y, u) *
                   E(a, b, RAB.z, LA.z, LB.z, v) *
                   E(c, d, RCD.x, LC.x, LD.x, τ) *
                   E(c, d, RCD.y, LC.y, LD.y, μ) *
                   E(c, d, RCD.z, LC.z, LD.z, υ) *
                   R(t+τ, u+μ, v+υ, 0, ω, RPQ) *
                   (-1)^(τ+μ+υ)
        end
    end
    ret *= 2 * π^2.5 / (p * q * sqrt(p + q))
    return ret
end