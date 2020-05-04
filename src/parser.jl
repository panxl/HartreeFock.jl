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