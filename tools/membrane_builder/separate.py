
def SplitComponents(input_file, ligand_name, output_ligand='ligand.pdb'):

    file = open(input_file  )
    lines = file.readlines()

    # Standard amino acid residues
    # AMINO_ACIDS = {
    #     'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    #     'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    #     'SEC', 'PYL' , 'HID' # Selenocysteine and Pyrrolysine
    # }

    AMINO_ACIDS = [

        # =========================
        # Canonical amino acids
        # =========================
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",

        # =========================
        # N-terminal variants
        # =========================
        "NALA", "NARG", "NASN", "NASP", "NCYS",
        "NGLN", "NGLU", "NGLY", "NHIS", "NILE",
        "NLEU", "NLYS", "NMET", "NPHE", "NPRO",
        "NSER", "NTHR", "NTRP", "NTYR", "NVAL",

        # =========================
        # C-terminal variants
        # =========================
        "CALA", "CARG", "CASN", "CASP", "CCYS",
        "CGLN", "CGLU", "CGLY", "CHIS", "CILE",
        "CLEU", "CLYS", "CMET", "CPHE", "CPRO",
        "CSER", "CTHR", "CTRP", "CTYR", "CVAL",

        # =========================
        # Histidine protonation states
        # =========================
        "HID", "HIE", "HIP",

        # =========================
        # Acidic residue protonation
        # =========================
        "ASH", "GLH",

        # =========================
        # Cysteine variants
        # =========================
        "CYM",   # deprotonated
        "CYX",   # disulfide bonded

        # =========================
        # Lysine neutral variant
        # =========================
        "LYN",
    ]

    ligand_lines = []

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # Extract residue name from columns 18-20 (0-indexed: 17:20)
            resname = line[17:20].strip()
            if resname == ligand_name:
                ligand_lines.append(line)

    if len(ligand_lines) > 0:
        if ligand_lines[0] == "TER\n":
            ligand_lines.pop(0)

    # Always write output file
    with open(output_ligand, 'w') as l_file:
        if len(ligand_lines) > 0:
            l_file.writelines(ligand_lines)
            l_file.write("TER\n")
        else:
            # Write a comment indicating no ligand found
            l_file.write(f"REMARK  No ligand '{ligand_name}' found in input file\n")       

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Separate PDB file into protein, ligand, and membrane/solvent components.")
    parser.add_argument("--input-file", help="Input PDB file")
    parser.add_argument("--ligand-name", "-l", required=False, default="XXX",  help="Ligand residue name to separate")
    parser.add_argument("--output-ligand", default="ligand.pdb", help="Output PDB file for ligand")

    args = parser.parse_args()

    SplitComponents(args.input_file, args.ligand_name, args.output_ligand)