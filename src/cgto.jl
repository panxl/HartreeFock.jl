struct GTO
    ζ::Float64
    R::Vec3{Float64}
    L::Vec3{Int}
    norm::Float64
end

function GTO(
    ζ::Float64,
    R::AbstractVector{Float64},
    L::AbstractVector{Int},
)
    gto = GTO(ζ, R, L, 1.0)
    norm = 1 / sqrt(overlap(gto, gto))
    return GTO(ζ, R, L, norm)
end

gaussian_product_center(A::GTO, B::GTO) = gaussian_product_center(A.ζ, A.R, B.ζ, B.R)

struct CGTO{N}
    d::SVector{N, Float64}
    ζ::SVector{N, Float64}
    R::Vec3{Float64}
    L::Vec3{Int}
    pnorm::SVector{N, Float64}
    norm::Float64
end

function CGTO(
    d::AbstractVector{Float64},
    ζ::AbstractVector{Float64},
    R::AbstractVector{Float64},
    L::AbstractVector{Int},
    norm::Float64,
)
    pnorm = similar(d)
    for i = 1:length(d)
        gto = GTO(ζ[i], R, L)
        pnorm[i] = gto.norm
    end
    CGTO{length(d)}(d, ζ, R, L, pnorm, norm)
end

function CGTO(
    d::AbstractVector{Float64},
    ζ::AbstractVector{Float64},
    R::AbstractVector{Float64},
    L::AbstractVector{Int},
)
    cgto = CGTO(d, ζ, R, L, 1.0)
    norm = 1 / sqrt(overlap(cgto, cgto))
    CGTO{length(d)}(d, ζ, R, L, cgto.pnorm, norm)
end

Base.iterate(cgto::CGTO, state=1) = state > length(cgto) ? nothing : (cgto[state], state+1)
Base.eltype(cgto::CGTO) = GTO
Base.length(cgto::CGTO) = length(cgto.d)
Base.size(cgto::CGTO) = length(cgto)

Base.getindex(cgto::CGTO, i::Int) = GTO(cgto.ζ[i], cgto.R, cgto.L, cgto.pnorm[i])
Base.firstindex(cgto::CGTO) = 1
Base.lastindex(cgto::CGTO) = length(cgto)