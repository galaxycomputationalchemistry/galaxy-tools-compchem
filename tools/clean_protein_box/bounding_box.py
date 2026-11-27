import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

TLEAP_TEMPLATE = """set default nocenter on
set default reorder_residues off

source leaprc.protein.ff14SB
source leaprc.gaff2
{ligand_block}
# Load protein system
mol = loadpdb {pdb}
{combine_block}
# Automatically create VDW bounding box + 10 A buffer
setBox {system_var} vdw

check {system_var}

saveamberparm {system_var} {parm7} {rst7}
savepdb {system_var} {leap_pdb}

quit
"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a solvated Amber system: run tleap, size a VDW "
        "bounding box, and box the restart file with ChBox."
    )
    parser.add_argument(
        "--pdb",
        default="2RH1.pdb",
        help="Input PDB file. May contain a co-crystallized ligand, other "
        "HETATM hetero groups, and/or non-protein residues written as "
        "plain ATOM records (e.g. lipids/waters/ions in a membrane-system "
        "PDB). All HETATM records are stripped, then any remaining ATOM "
        "residue that isn't a standard amino acid is stripped too, to "
        "build a protein-only PDB for tleap; pass --frcmod/--mol2 to "
        "merge a ligand back in from its own coordinates.",
    )
    parser.add_argument(
        "--clean-pdb",
        default="protein_clean.pdb",
        help="Where to write the protein-only PDB used to build the tleap system",
    )
    parser.add_argument(
        "--buffer", type=float, default=10.0, help="Buffer added on each side (A)"
    )

    parser.add_argument(
        "--manual-box-x", type=float, help="Manual box X dimension (A). If provided with Y and Z, skip auto-calculation."
    )
    parser.add_argument(
        "--manual-box-y", type=float, help="Manual box Y dimension (A). If provided with X and Z, skip auto-calculation."
    )
    parser.add_argument(
        "--manual-box-z", type=float, help="Manual box Z dimension (A). If provided with X and Y, skip auto-calculation."
    )
    parser.add_argument(
        "--manual-box-alpha", type=float, default=90.0, help="Manual box alpha angle (degrees). Default: 90.0"
    )
    parser.add_argument(
        "--manual-box-beta", type=float, default=90.0, help="Manual box beta angle (degrees). Default: 90.0"
    )
    parser.add_argument(
        "--manual-box-gamma", type=float, default=90.0, help="Manual box gamma angle (degrees). Default: 90.0"
    )

    parser.add_argument(
        "--tleap-file",
        help="Existing tleap input file to run. If omitted, one is generated "
        "from the built-in template.",
    )
    parser.add_argument(
        "--tleap-out",
        default="tleap.in",
        help="Where to write the generated tleap input file (ignored if "
        "--tleap-file is given)",
    )
    parser.add_argument(
        "--frcmod",
        help="Ligand frcmod file. If given together with --mol2, the "
        "generated tleap input loads the ligand; if omitted, the tleap "
        "input is built for the protein only.",
    )
    parser.add_argument(
        "--mol2",
        help="Ligand mol2 file. If given together with --frcmod, the "
        "generated tleap input loads the ligand; if omitted, the tleap "
        "input is built for the protein only.",
    )
    parser.add_argument(
        "--parm7", default="protein.parm7", help="Topology written by tleap"
    )
    parser.add_argument(
        "--rst7", default="protein.rst7", help="Restart file written by tleap"
    )
    parser.add_argument(
        "--leap-pdb",
        default="protein_leap.pdb",
        help="Combined PDB written by tleap",
    )
    parser.add_argument(
        "--box-rst7",
        default="protein_box.rst7",
        help="Restart file written by ChBox",
    )

    parser.add_argument("--tleap-exe", default="tleap")
    parser.add_argument("--chbox-exe", default="ChBox")

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would run without executing them",
    )

    return parser.parse_args()


def compute_bounding_box(pdb_file, buffer_angstrom):
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    coords = np.array(coords)
    min_xyz = coords.min(axis=0)
    max_xyz = coords.max(axis=0)
    size = max_xyz - min_xyz
    buffered_size = size + 2.0 * buffer_angstrom
    return size, buffered_size


# Standard amino acid residue names tleap/ff14SB recognize, including
# alternate protonation states (HID/HIE/HIP, ASH, GLH, LYN, CYX/CYM) and
# the common N-/C-terminal caps.
STANDARD_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH",
    "CYS", "CYM", "CYX",
    "GLN", "GLU", "GLH", "GLY",
    "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "LYN",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "ACE", "NME", "NHE", "NH2",
}


def write_clean_pdb(pdb_file, clean_pdb_out):
    """Two-step clean: (1) drop every HETATM record outright, then (2) of
    the remaining ATOM records, keep only those whose residue is a
    standard amino acid. Step 2 is what catches non-protein residues
    (lipids, water, ions, ...) that Amber-prepped/membrane-system PDBs
    (e.g. packmol-memgen/charmmlipid2amber output) commonly write as
    plain ATOM records rather than HETATM -- a HETATM-only filter would
    miss them. CONECT records are dropped too, since their atom-serial
    references may no longer line up once atoms are removed. TER records
    are kept only when they actually terminate a kept chain.
    """
    kept_lines = []
    dropped_resnames = set()
    dropped_count = 0
    last_kept_was_atom = False

    with open(pdb_file) as f:
        for line in f:
            record = line[:6]
            if record == "HETATM":
                continue
            if record == "ATOM  ":
                resname = line[17:20].strip().upper()
                if resname in STANDARD_PROTEIN_RESNAMES:
                    kept_lines.append(line)
                    last_kept_was_atom = True
                else:
                    dropped_resnames.add(resname)
                    dropped_count += 1
                continue
            if record == "CONECT":
                continue
            if record.startswith("TER"):
                if last_kept_was_atom:
                    kept_lines.append(line)
                    last_kept_was_atom = False
                continue
            kept_lines.append(line)

    Path(clean_pdb_out).write_text("".join(kept_lines))

    if dropped_count:
        print(
            f"Cleaned PDB: removed {dropped_count} non-protein ATOM "
            f"atom(s) from residue(s) {sorted(dropped_resnames)} (plus "
            "all HETATM records)"
        )
    else:
        print("Cleaned PDB: removed all HETATM records; no non-protein "
              "ATOM residues found")

    return dropped_count, dropped_resnames


def resolve_tleap_file(tleap_file, tleap_out, frcmod, mol2, pdb, parm7, rst7, leap_pdb):
    if tleap_file:
        path = Path(tleap_file).resolve()
        if not path.exists():
            raise FileNotFoundError(f"--tleap-file not found: {path}")
        return path

    if bool(frcmod) != bool(mol2):
        raise ValueError("--frcmod and --mol2 must be given together, or not at all")

    if frcmod and mol2:
        ligand_block = (
            f"\n# Ligand parameters (residue name comes from the mol2 file)"
            f"\nloadamberparams {frcmod}"
            f"\nligand = loadmol2 {mol2}\n"
        )
        combine_block = "\n# Merge the ligand (its own mol2 coordinates) into the protein\ncomplex = combine { mol ligand }\n"
        system_var = "complex"
    else:
        ligand_block = ""
        combine_block = ""
        system_var = "mol"

    path = Path(tleap_out).resolve()
    path.write_text(
        TLEAP_TEMPLATE.format(
            ligand_block=ligand_block,
            combine_block=combine_block,
            system_var=system_var,
            pdb=pdb,
            parm7=parm7,
            rst7=rst7,
            leap_pdb=leap_pdb,
        )
    )
    return path


def run_command(cmd, cwd, dry_run):
    print(f"$ {' '.join(str(c) for c in cmd)}")
    if dry_run:
        return
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except FileNotFoundError:
        sys.exit(f"Error: '{cmd[0]}' not found on PATH. Is AmberTools loaded?")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error: '{cmd[0]}' failed with exit code {e.returncode}")


def main():
    args = parse_args()
    work_dir = Path.cwd()
    pdb_abs = Path(args.pdb).resolve()
    frcmod_abs = Path(args.frcmod).resolve() if args.frcmod else None
    mol2_abs = Path(args.mol2).resolve() if args.mol2 else None

    # 1. Strip all HETATM records, then strip any remaining ATOM residue
    #    that isn't a standard amino acid (lipids/waters/ions written as
    #    ATOM), to get a protein-only PDB.
    clean_pdb_path = (work_dir / args.clean_pdb).resolve()
    write_clean_pdb(pdb_abs, clean_pdb_path)
    print(f"clean protein PDB: {clean_pdb_path}")

    # 2. Generate (or reuse) the tleap input, then run tleap. The ligand (if
    #    --frcmod/--mol2 are given) is merged in from its own mol2 coordinates.
    tleap_path = resolve_tleap_file(
        args.tleap_file,
        args.tleap_out,
        frcmod_abs,
        mol2_abs,
        clean_pdb_path,
        args.parm7,
        args.rst7,
        args.leap_pdb,
    )
    print(f"tleap input file:   {tleap_path}")
    run_command([args.tleap_exe, "-f", str(tleap_path)], work_dir, args.dry_run)

    # 3. Compute the VDW bounding box + buffer from the CLEANED protein PDB
    #    (after removing HETATM and non-protein residues), optionally with ligand.
    #    OR use manually provided box dimensions if given.
    box_alpha = 90.0
    box_beta = 90.0
    box_gamma = 90.0
    
    if args.manual_box_x and args.manual_box_y and args.manual_box_z:
        # Use manual box dimensions
        buffered_size = np.array([args.manual_box_x, args.manual_box_y, args.manual_box_z])
        box_alpha = args.manual_box_alpha
        box_beta = args.manual_box_beta
        box_gamma = args.manual_box_gamma
        print("Using manual box dimensions (from crystallography data):")
        print(f"Box size (A): X={buffered_size[0]:.3f}, Y={buffered_size[1]:.3f}, Z={buffered_size[2]:.3f}")
        print(f"Angles (degrees): alpha={box_alpha:.2f}, beta={box_beta:.2f}, gamma={box_gamma:.2f}")
    elif args.manual_box_x or args.manual_box_y or args.manual_box_z:
        # Partial manual dimensions provided
        sys.exit("Error: --manual-box-x, --manual-box-y, and --manual-box-z must all be provided together, or not at all")
    else:
        # Auto-calculate from protein
        size, buffered_size = compute_bounding_box(clean_pdb_path, args.buffer)
        print("Raw box size (A):        ", size)
        print(f"Box size with {args.buffer} A buffer:", buffered_size)
        print(f"Angles (degrees): alpha={box_alpha:.2f}, beta={box_beta:.2f}, gamma={box_gamma:.2f}")

    # 4. Box the tleap restart file with ChBox.
    rst7_path = work_dir / args.rst7
    box_rst7_path = work_dir / args.box_rst7
    run_command(
        [
            args.chbox_exe,
            "-c", str(rst7_path),
            "-o", str(box_rst7_path),
            "-X", str(buffered_size[0]),
            "-Y", str(buffered_size[1]),
            "-Z", str(buffered_size[2]),
            "-al", str(box_alpha),
            "-bt", str(box_beta),
            "-gm", str(box_gamma),
        ],
        work_dir,
        args.dry_run,
    )

    print(f"parm7:    {work_dir / args.parm7}")
    print(f"box rst7: {box_rst7_path}")


if __name__ == "__main__":
    main()
