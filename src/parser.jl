function parse_geom(string::AbstractString)
    if isfile(string)
        return read_xyz_file(string)
    else
        return read_xyz_string(string)
    end
end

function read_xyz_line(line::AbstractString)
    s, x, y, z = split(line)[1:4]
    position = parse.(Float64, [x, y, z]) / CODATA08_BOHR_TO_A
    number = findfirst(==(titlecase(s)), ELEMENTS)
    return position, number
end

function read_xyz_file(file::AbstractString)
    io = open(file, "r")
    n = parse(Int, readline(io))
    readline(io)
    positions = Vector{Vec3{Float64}}(undef, n)
    numbers = Vector{Int}(undef, n)
    for i = 1:n
        line = readline(io)
        positions[i], numbers[i] = read_xyz_line(line)
    end
    close(io)
    return positions, numbers
end

function read_xyz_string(string::AbstractString)
    lines = split(strip(string), "\n")
    n = length(lines)
    positions = Vector{Vec3{Float64}}(undef, n)
    numbers = Vector{Int}(undef, n)
    for (i, line) in enumerate(lines)
        positions[i], numbers[i] = read_xyz_line(line)
    end
    return positions, numbers
end

function parse_pc(string::AbstractString)
    if isfile(string)
        return read_pc_file(string)
    else
        return read_pc_string(string)
    end
end

function read_pc_line(line::AbstractString)
    q, x, y, z = split(line)[1:4]
    position = parse.(Float64, [x, y, z]) / CODATA08_BOHR_TO_A
    charge = parse(Float64, q)
    return position, charge
end

function read_pc_file(file::AbstractString)
    io = open(file, "r")
    n = parse(Int, readline(io))
    positions = Vector{Vec3{Float64}}(undef, n)
    charges = Vector{Float64}(undef, n)
    for i = 1:n
        line = readline(io)
        positions[i], charges[i] = read_pc_line(line)
    end
    close(io)
    return positions, charges
end

function read_pc_string(string::AbstractString)
    lines = split(strip(string), "\n")
    n = length(lines)
    positions = Vector{Vec3{Float64}}(undef, n)
    charges = Vector{Float64}(undef, n)
    for (i, line) in enumerate(lines)
        positions[i], charges[i] = read_pc_line(line)
    end
    return positions, charges
end