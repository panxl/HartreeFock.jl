using HartreeFock: Mole, Env, scf!, PointCharges, electrostatic_potential
using Test

H2 = Mole("
H   0.000000  0.000000  0.000000
H   0.740848  0.000000  0.000000
", "STO-3G")
E = scf!(H2)
@test E ≈ -1.116714325 atol=1e-7

H2O = Mole(joinpath(@__DIR__, "data/water.xyz"), "STO-3G")
E = scf!(H2O)
@test E ≈ -74.962927947 atol=1e-7

pointcharges = PointCharges("
-0.834    0.00000  0.00000 -5.00000
 0.417    0.00000  0.00000 -4.04280
 0.417   -0.92663  0.00000 -5.23999
")
env = Env(pointcharges)
E = scf!(H2O, env)
@test E ≈ -74.964115862

V = electrostatic_potential(H2O, pointcharges)
@test all(isapprox.(V, [-0.73578734E-02, -0.11128084E-01, -0.64441894E-02], atol=1e-7))