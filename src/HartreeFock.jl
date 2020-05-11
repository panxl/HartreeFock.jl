module HartreeFock

using JSON
using LinearAlgebra
using SpecialFunctions
using StaticArrays
using Zygote
using Zygote: @adjoint, dropgrad

include("constants.jl")
include("parser.jl")
include("utils.jl")
include("particles.jl")
include("basis.jl")
include("mole.jl")
include("integral.jl")
include("diis.jl")
include("scf.jl")
include("electrostatics.jl")

end # module