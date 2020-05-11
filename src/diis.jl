using ShiftedArrays

function DIIS(n::Int, depth::Int=7)
    trial_vector = SMatrix{n, n, Float64}[]
    residual_vector = SMatrix{n, n, Float64}[]

    L = zeros(depth+1, depth+1)
    L[begin, :] .= -1
    L[:, begin] .= -1
    L[begin, begin] = 0

    function add_trial(trial)
        push!(trial_vector, trial)
        if length(trial_vector) > depth
            popfirst!(trial_vector)
        end
    end

    function add_residual(residual)
        push!(residual_vector, residual)
        if length(residual_vector) > depth
            popfirst!(residual_vector)
            B = @view get_L()[2:end, 2:end]
            B = ShiftedArrays.circshift(B, (-1, -1))
        else
            B = @view get_L()[2:end, 2:end]   
        end
        for i = 1:length(residual_vector)
            B[i, end] = sum(residual_vector[i] .* residual)
        end
    end

    function get_L()
        return @view L[1:length(residual_vector)+1, 1:length(residual_vector)+1]
    end

    function get_coeff()
        coeff = -inv(Symmetric(get_L()))[begin, begin+1:end]
    end

    #Declare public:
    ()->(trial_vector, residual_vector, add_trial, add_residual, get_coeff)
end