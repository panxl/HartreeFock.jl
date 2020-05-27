struct DIIS
    trial_vector::Array{Float64,3}
    residual_vector::Array{Float64,3}
    L::Matrix{Float64}
    depth::Int
    cycle::MVector{1,Int}
    function DIIS(n::Integer, depth::Integer=10)
        trial_vector = zeros(n, n, depth)
        residual_vector =  zeros(n, n, depth)
        L = zeros(depth+1, depth+1)
        L[1, 2:end] .= 1
        L[2:end, 1] .= 1
        new(trial_vector, residual_vector, L, depth, [0])
    end
end

function update!(
    diis::DIIS,
    trial::AbstractMatrix{Float64},
    residual::AbstractMatrix{Float64},
    )
    ptr = mod(diis.cycle[], diis.depth) + 1
    diis.trial_vector[:, :, ptr] = trial
    diis.residual_vector[:, :, ptr] = residual

    n = diis.cycle[] < diis.depth ? ptr : diis.depth
    for i = 1:n
        diis.L[i+1, ptr+1] = diis.residual_vector[:, :, i] ⋅ residual
        diis.L[ptr+1, i+1] = diis.L[i+1, ptr+1]
    end

    coeff = inv(diis.L[1:n+1, 1:n+1])[2:n+1, 1]

    m = LinearAlgebra.checksquare(trial)
    for i = 1:m, j = i:m
        trial[i, j] = 0.0
        for k = 1:n
            trial[i, j] += coeff[k] * diis.trial_vector[i, j, k]
        end
        trial[j, i] = trial[i, j]
    end
    diis.cycle[] += 1
end