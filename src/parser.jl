function parser(file::String)

    lines = readlines(file)

    try
        return read_xyz(lines)
    catch e
    end

    try
        return JSON.parse(join(lines))
    catch e
        println("Unrecognized file type")
    end
end

const ELEMENTS = [
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi"
    ]

const CODATA08_BOHR_TO_A = 0.5291772085936

function read_xyz(lines::Vector{String})
    n = parse(Int, lines[1])
    positions = Vector{Vec3{Float64}}(undef, n)
    numbers = Vector{Int}(undef, n)
    charges = Vector{Float64}(undef, n)
    for i = 1:n
        s, x, y, z = split(lines[i + 2])
        positions[i] = parse.(Float64, [x, y, z]) / CODATA08_BOHR_TO_A
        numbers[i] = findfirst(==(titlecase(s)), ELEMENTS)
        charges[i] = float(numbers[i])
    end
    return positions, numbers, charges
end

read_xyz(file::String) = read_xyz(readlines(file))
