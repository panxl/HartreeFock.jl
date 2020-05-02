using PyCall

function __init__()
    py"""
    import basis_set_exchange as bse

    def dump_json(basis, path=None):
        if path is None:
            path = "."
        bse.fileio.write_json_basis(f"{path}/{basis}.json", bse.get_basis(basis))
    """
end

dump_json(basis::String, path::String) = py"dump_json"(basis, path)
dump_json(basis::String) = py"dump_json"(basis)