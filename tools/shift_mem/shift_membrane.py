#!/usr/bin/env python3
import argparse
import os
import re
from dataclasses import dataclass

import numpy as np

# --- Basic residue classification ---

AA_RESNAMES = {
    # Standard amino acids
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    # Common AMBER protonation/tautomer variants
    "HID","HIE","HIP","ASH","GLH","LYN","CYM","CYX"
}

WATER_RESNAMES = {"WAT", "TIP3", "HOH", "SOL"}

# --- Mass table (extend as needed) ---
ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
    "F": 18.998,
    "CL": 35.45,
    "NA": 22.990,
    "K": 39.0983,
    "BR": 79.904,
    "I": 126.90447,
    "MG": 24.305,
    "CA": 40.078,
    "ZN": 65.38,
}

@dataclass
class Coord:
    x: float
    y: float
    z: float

def parse_pdb_atom_line(line: str):
    """
    Returns (record, atom_name, resname, x, y, z) or None if not ATOM/HETATM.
    Uses fixed-column PDB parsing (more robust than split()).
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    atom_name = line[12:16].strip()
    resname = line[17:20].strip()
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    record = line[0:6].strip()
    return record, atom_name, resname, x, y, z

def element_from_atom_name(atom_name: str) -> str:
    """
    Extract element from an atom name like:
      'CA' (alpha carbon) -> 'C' (not calcium)
      'Cl-' -> 'CL'
      'K+'  -> 'K'
      '1H'  -> 'H'
      'OW'  -> 'O'
    Strategy:
      - keep letters only
      - if starts with CL/BR/NA/CA/MG/ZN etc, use 2-letter element
      - else use first letter
    """
    letters = re.sub(r"[^A-Za-z]", "", atom_name).upper()
    if not letters:
        return ""
    # Prefer known 2-letter elements if present
    if len(letters) >= 2 and letters[:2] in ATOMIC_MASS:
        return letters[:2]
    # Special case: "CA" in proteins is usually Carbon alpha, not Calcium
    if letters.startswith("CA") and atom_name.strip().upper() == "CA":
        return "C"
    return letters[0]

def atom_mass_from_name(atom_name: str) -> float:
    elem = element_from_atom_name(atom_name)
    if not elem or elem not in ATOMIC_MASS:
        raise KeyError(f"Unknown element parsed from atom name '{atom_name}' -> '{elem}'. "
                       f"Add it to ATOMIC_MASS if needed.")
    return ATOMIC_MASS[elem]

def protein_com_pdb(pdb_path: str) -> Coord:
    masses = []
    coords = []
    with open(pdb_path, "r") as f:
        for line in f:
            parsed = parse_pdb_atom_line(line)
            if not parsed:
                continue
            _, atom_name, resname, x, y, z = parsed
            if resname not in AA_RESNAMES:
                continue
            m = atom_mass_from_name(atom_name)
            masses.append(m)
            coords.append((x, y, z))

    if not masses:
        raise RuntimeError(f"No protein atoms found (resname in AA list) in: {pdb_path}")

    masses_np = np.array(masses, dtype=float)
    coords_np = np.array(coords, dtype=float)
    com = (coords_np * masses_np[:, None]).sum(axis=0) / masses_np.sum()
    return Coord(float(com[0]), float(com[1]), float(com[2]))

def water_box_size(pdb_path: str) -> Coord:
    pts = []
    with open(pdb_path, "r") as f:
        for line in f:
            parsed = parse_pdb_atom_line(line)
            if not parsed:
                continue
            _, _, resname, x, y, z = parsed
            if resname in WATER_RESNAMES:
                pts.append((x, y, z))
    if not pts:
        return Coord(0.0, 0.0, 0.0)
    pts = np.array(pts, dtype=float)
    return Coord(float(pts[:,0].max() - pts[:,0].min()),
                 float(pts[:,1].max() - pts[:,1].min()),
                 float(pts[:,2].max() - pts[:,2].min()))

def shift_line_xyz(line: str, dx: float, dy: float, dz: float) -> str:
    # Preserve original formatting as much as possible; only overwrite XYZ fields.
    x = float(line[30:38]) + dx
    y = float(line[38:46]) + dy
    z = float(line[46:54]) + dz
    return f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", required=True, help="Original input receptor PDB file")
    ap.add_argument("-m", required=True, help="Packmol-memgen generated protein+membrane PDB to shift")
    ap.add_argument("-o", required=True, help="Output name of shifted membrane/solvent/ions PDB")
    args = ap.parse_args()

    if not (os.path.isfile(args.i) and os.path.isfile(args.m)):
        raise FileNotFoundError("Input file not found.")

    prot_com = protein_com_pdb(args.i)
    mem_prot_com = protein_com_pdb(args.m)

    print(f"Original protein COM: {prot_com.x:.3f} {prot_com.y:.3f} {prot_com.z:.3f}\n")
    print(f"Protein COM in membrane system: {mem_prot_com.x:.3f} {mem_prot_com.y:.3f} {mem_prot_com.z:.3f}\n")

    dx = prot_com.x - mem_prot_com.x
    dy = prot_com.y - mem_prot_com.y
    dz = prot_com.z - mem_prot_com.z

    print(f"Writing {args.o} (non-protein only), shifted by {dx:.3f} {dy:.3f} {dz:.3f}\n")

    with open(args.o, "w") as out, open(args.m, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                # Write out membrane/water/ions/etc — everything that is NOT protein
                if resname not in AA_RESNAMES:
                    out.write(shift_line_xyz(line, dx, dy, dz))
            elif line.startswith(("TER", "END", "MODEL", "ENDMDL", "REMARK", "CRYST1")):
                out.write(line)

    box = water_box_size(args.m)
    if box.x > 0:
        print(f"Water box size (from water atoms): X,Y,Z = {box.x:.3f} {box.y:.3f} {box.z:.3f}\n")

if __name__ == "__main__":
    main()
