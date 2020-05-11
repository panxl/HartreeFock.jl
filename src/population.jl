function mulliken_charges(scf::SCF)
    mulliken_charges = copy(scf.mole.nuclei.charges)
    orb_charges = sum(scf.P .* scf.S, dims=2)
    for (i, q) in enumerate(orb_charges)
        mulliken_charges[scf.mole.basis.ids[i]] -= q
    end
    return mulliken_charges
end