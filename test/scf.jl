using HartreeFock: Atom, Atoms, PCharge, PCharges, Mole, Env, scf

atoms = Atoms()
push!(atoms, Atom([0.0,0.0,0.0], 1))
push!(atoms, Atom([1.4,0.0,0.0], 1))
H2 = Mole(atoms, "STO-3G")
E, P = scf(H2)
@test E ≈ -1.116714325

H2O = Mole(Atoms("data/water.xyz"), "STO-3G")
E, P = scf(H2O)
@test E ≈ -74.962927947

CODATA08_BOHR_TO_A = 0.5291772085936
pcharges = PCharges()
push!(pcharges, PCharge([ 0.00000, 0.00000, -5.00000] / CODATA08_BOHR_TO_A, -0.834))
push!(pcharges, PCharge([ 0.00000, 0.00000, -4.04280] / CODATA08_BOHR_TO_A,  0.417))
push!(pcharges, PCharge([-0.92663, 0.00000, -5.23999] / CODATA08_BOHR_TO_A,  0.417))
env = Env(pcharges)
E, P = scf(H2O, env)
@test E ≈ -74.964115862
