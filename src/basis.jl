struct Shell{M,N}
    exponents::SVector{N, Float64}
    angular_momentum::Int
    coefficients::SVector{M, SVector{N, Float64}}
end

struct Basis
    shells::Vector{Vector{Shell}}
    numbers::Vector{Int}
end

function Basis(numbers::AbstractVector{Int}, basis_dict::Dict)
    numbers = sort(unique(numbers))
    basis = Basis(Vector{Shell}[], Int[])
    for number in numbers
        shells = Shell[]
        for shell in basis_dict[string(number)]["electron_shells"]
            exponents = parse.(Float64, shell["exponents"])
            angular_momentum = shell["angular_momentum"]
            N = length(exponents)
            if length(angular_momentum) == length(shell["coefficients"])
                for i = 1:length(angular_momentum)
                    coefficients = [parse.(Float64, shell["coefficients"][i])]
                    push!(shells, Shell{1, N}(exponents, angular_momentum[i], coefficients))
                end
            elseif length(angular_momentum) == 1
                M = length(shell["coefficients"])
                coefficients = []
                for i = 1:M
                    push!(coefficients, parse.(Float64, shell["coefficients"][i]))
                end
                push!(shells, Shell{M, N}(exponents, angular_momentum[1], coefficients))
            end
        end
        push!(basis.shells, shells)
        push!(basis.numbers, number)
    end
    return basis
end

function Basis(numbers::AbstractVector{Int}, basis_str::String)
    file = joinpath(@__DIR__, "data/basis/$(basis_str).json")
    basis_dict = JSON.parsefile(file)["elements"]
    return Basis(numbers, basis_dict)
end

function Basis(nuclei::Nuclei, basis::Union{Dict,String})
    return Basis(nuclei.numbers, basis)
end