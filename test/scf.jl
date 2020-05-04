using HartreeFock: Mole, Env, scf, PointCharges
using Test

H2 = Mole("
H   0.000000  0.000000  0.000000
H   0.740848  0.000000  0.000000
", "STO-3G")
E, P = scf(H2)
@test E ≈ -1.116714325

H2O = Mole(joinpath(@__DIR__, "data/water.xyz"), "STO-3G")
E, P = scf(H2O)
@test E ≈ -74.962927947

pointcharges = PointCharges("
-0.834    0.00000  0.00000 -5.00000
 0.417    0.00000  0.00000 -4.04280
 0.417   -0.92663  0.00000 -5.23999
")
env = Env(pointcharges)
E, P = scf(H2O, env)
@test E ≈ -74.964115862