
amino_acids = ["ALA", "ARG", "ASP", "ASH", "ASN", "CYS", "CYM", "CYX",
               "GLN", "GLU", "GLH", "GLY", "HIS", "HIE", "HID", "HIP",
               "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER",
               "THR", "TRP", "TYR", "VAL"]

non_polar = ["ALA", "GLY", "ILE", "LEU", "PHE", "PRO", "VAL"]

sc_atoms = {"ARG": ["NE", "NH1", "NH2"],
            "ASP": ["OD1", "OD2"], "ASH": ["OD1", "OD2"],
            "ASN": ["OD1", "ND2"],
            "CYS": ["SG"], "CYM": ["SG"], "CYX": ["SG"],
            "GLN": ["OE1", "NE2"],
            "GLU": ["OE1", "OE2"], "GLH": ["OE1", "OE2"],
            "HID": ["ND1", "NE2"], "HIE": ["ND1", "NE2"], "HIP": ["ND1", "NE2"], "HIS": ["ND1", "NE2"],
            "LYS": ["NZ"], "LYN": ["NZ"],
            "MET": ["SD"],
            "SER": ["OG"],
            "THR": ["OG1"],
            "TRP": ["NE1"],
            "TYR": ["OH"]}

def atom_search(resn):
    atoms = ["O", "N"]

    if resn not in amino_acids:
        raise ValueError(f"'{resn}' is not a valid amino acid.")

    if resn in amino_acids:
        if resn in sc_atoms.keys():
            atoms += sc_atoms[resn]

        return atoms

def sc_atom_search(resn):
    if resn not in amino_acids:
        raise ValueError(f"'{resn}' is not a valid amino acid.")

    if resn in amino_acids:
        if resn in sc_atoms.keys():

            return sc_atoms[resn]

