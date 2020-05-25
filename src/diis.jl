function DIIS(n::Int, depth::Int=10)
    trial_vector = SMatrix{n, n, Float64}[]
    residual_vector = SMatrix{n, n, Float64}[]

    L = zeros(depth+1, depth+1)
    L[1, :] .= 1
    L[:, 1] .= 1
    L[1, 1] = 0

    function add_trial(trial::Matrix{Float64})
        push!(trial_vector, trial)
        if length(trial_vector) > depth
            popfirst!(trial_vector)
        end
    end

    function add_residual(residual::Matrix{Float64})
        push!(residual_vector, residual)
        if length(residual_vector) > depth
            popfirst!(residual_vector)
            B = @view get_L()[2:end, 2:end]
            B = ShiftedArrays.circshift(B, (-1, -1))
        else
            B = @view get_L()[2:end, 2:end]   
        end
        for i = 1:length(residual_vector)
            B[i, end] = residual_vector[i] ⋅ residual
        end
    end

    function get_L()
        return @view L[1:length(residual_vector)+1, 1:length(residual_vector)+1]
    end

    function get_coeff()
        return inv(Symmetric(get_L()))[1, 2:end]
    end

    function get_F()
        coeff = get_coeff()
        return Array(sum(coeff .* trial_vector))
    end

    #Declare public:
    ()->(trial_vector, residual_vector, add_trial, add_residual, get_F)
end