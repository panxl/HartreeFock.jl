module HartreeFock

using LinearAlgebra: norm, eigen
using SpecialFunctions: gamma, gamma_inc
using StaticArrays
using JSON

include("parser.jl")
include("utils.jl")
include("atoms.jl")
include("cgto.jl")
include("basis.jl")
include("mole.jl")
include("integral.jl")
include("scf.jl")

end # module