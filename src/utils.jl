struct Vec3{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

StaticArrays.similar_type(::Type{<:Vec3}, ::Type{T}, ::Size{(3,)}) where {T} = Vec3{T}

function gaussian_product_center(
        a::Float64,
        RA::Vec3{Float64},
        b::Float64,
        RB::Vec3{Float64},
    )
    p = a + b
    return (a * RA + b * RB) / p
end

function norm2(R::Vec3{Float64})
    s = 0.0
    for r in R
        s += abs2(r)
    end
    return s
end

gamma_inc_p(a, x) = gamma_inc(a, x, 0)[1]

function boys(n::Int, x::Float64)
    if x < 1e-7
        return 1 / (2 * n + 1) - x / (2 * n + 3) + x^2 / (2 * n + 5) / 2
    else
        return 1 / (2 * x^(n + 0.5)) * gamma(n + 0.5) * gamma_inc_p(n + 0.5, x)
    end
end