from phonopy.structure.atoms import symbol_map, atom_data


def mass_from_symbol(symbol: str):
    index = symbol_map[symbol.title()]
    data = atom_data[index]
    return data[3]
