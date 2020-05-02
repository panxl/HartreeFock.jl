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

function boys(n::Int, x::Float64)
    if x < 1e-7
        return 1 / (2 * n + 1) - x / (2 * n + 3)
    else
        return 1 / (2 * x^(n + 0.5)) * gamma(n + 0.5) * gamma_inc(n + 0.5, x, 0)[1]
    end
end

triangle(i::Int) = div(i*(i+1), 2)
triangle(i::Int, j::Int) = i<j ? triangle(j-1)+i : triangle(i-1)+j
pairs(n::Int) = ((i,j) for i = 1:n for j = i:n)

basis_iterator(n::Int) = ((i,j,k,l) for (i,j) in pairs(n) for (k,l) in pairs(n) if triangle(i,j) <= triangle(k,l))
basis_index(i::Int, j::Int, k::Int, l::Int) = triangle(triangle(i,j), triangle(k,l))