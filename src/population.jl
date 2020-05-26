function mulliken_charges(scf::SCF)
    mulliken_charges = copy(scf.mole.nuclei.charges)
    orb_charges = sum(2 .* scf.P .* scf.S, dims=2)
    basis_idx = get_basis_idx(scf.mole)
    for (i, q) in enumerate(orb_charges)
        mulliken_charges[basis_idx[i]] -= q
    end
    return mulliken_charges
end