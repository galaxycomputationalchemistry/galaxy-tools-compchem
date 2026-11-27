# #!/usr/bin/env python3
# """
# Configuration file generator for membrane_builder tool.
# Generates configuration files with system components (lipids, ligands, ions).
# """

# import argparse
# import sys
# import json
# from pathlib import Path
# from typing import Dict, List, Optional, Tuple

# # Get the directory of this script for finding JSON files
# SCRIPT_DIR = Path(__file__).parent
# LIPIDS_JSON_PATH = SCRIPT_DIR / "lipids.json"
# LIPID_FRAGMENTS_PATH = SCRIPT_DIR / "lipid_fragments.json"

# def load_json(file_path: Path) -> Dict:
#     """Load a JSON file and return its content as a dictionary."""
#     with open(file_path, 'r') as file:
#         return json.load(file)
    

# def generate_config(lipids: List[str], ligands: List[str], ions: List[str], output_path: Path) -> None:
#     frgs_dict = load_json(LIPID_FRAGMENTS_PATH)
#     lipid_frgs = list(set(i for l in lipids[0].split(',') for i in frgs_dict[l]))

#     lipid_line, ligand_line, ion_line = None, None, None

#     if len(lipid_frgs)  > 1:
#         lipid_line = "lipid:" + ",".join(lipid_frgs)
#     else:
#         lipid_line = "lipid:"+lipid_frgs[0]

#     if ligands:
#         ligand_line = "ligand:" + ligands

#     if ions:
#         if len(ions[0].split(',')) > 1:
#             ion_line = "ions:" + ",".join(ions[0].split(','))
#             print("Line 41")
#         else:
#             print("Line 43")
#             ion_line = "ions:" + ions[0]
#     print(ion_line)
        
#     with open(output_path, 'w') as outfile:
#         if lipid_line:
#             outfile.write(lipid_line + "\n")
#         if ligand_line:
#             outfile.write(ligand_line + "\n")     
#         if ion_line:
#             outfile.write(ion_line + "\n")  



# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Generate configuration file for membrane_builder.")
#     parser.add_argument('--lipids', nargs='+', required=True, help='List of lipids to include in the membrane.')
#     parser.add_argument('--ligands', type=str, help='List of ligands to include in the system.')
#     parser.add_argument('--ions', nargs='*', default=[], help='List of ions to include in the system.')
#     parser.add_argument('--output', type=Path, required=True, help='Output path for the configuration file.')

#     args = parser.parse_args()

#     generate_config(args.lipids, args.ligands, args.ions, args.output)

#!/usr/bin/env python3
"""
Configuration file generator for membrane_builder tool.
Generates configuration files with system components (lipids, ligands, ions).
"""

import argparse
import sys
import json
from pathlib import Path
from typing import Dict, List

# Get the directory of this script for finding JSON files
SCRIPT_DIR = Path(__file__).parent
LIPIDS_JSON_PATH = SCRIPT_DIR / "lipids.json"
LIPID_FRAGMENTS_PATH = SCRIPT_DIR / "lipid_fragments.json"

# 🔑 Ion normalization map → GROMACS index names
ION_MAP = {
    # Sodium
    "NA": "NA+",
    "NA+": "NA+",
    "SOD": "NA+",
    "SODIUM": "NA+",

    # Potassium
    "K": "K+",
    "K+": "K+",
    "POTASSIUM": "K+",

    # Calcium
    "CA": "CA2+",
    "CA2+": "CA2+",
    "CALCIUM": "CA2+",

    # Magnesium
    "MG": "MG2+",
    "MG2+": "MG2+",
    "MAGNESIUM": "MG2+",

    # Chloride
    "CL": "CL-",
    "CL-": "CL-",
    "CLA": "CL-",
    "CHLORIDE": "CL-"
}

def load_json(file_path: Path) -> Dict:
    """Load a JSON file and return its content as a dictionary."""
    with open(file_path, 'r') as file:
        return json.load(file)


def normalize_ion(ion: str) -> str:
    """
    Normalize ion name to GROMACS index group naming.
    Example: Na+, na, sodium → NA+
    """
    ion_clean = ion.strip().upper()
    return ION_MAP.get(ion_clean, ion_clean)  # fallback if unknown


def process_ions(ions: List[str]) -> List[str]:
    """
    Flatten, split, and normalize ion list.
    """
    if not ions:
        return []

    split_ions = []
    for item in ions:
        split_ions.extend(item.split(','))

    # Normalize and deduplicate
    normalized = [normalize_ion(i) for i in split_ions if i.strip()]
    return list(dict.fromkeys(normalized))  # preserve order, remove duplicates


def generate_config(lipids: List[str], ligands: str, ions: List[str], output_path: Path) -> None:
    frgs_dict = load_json(LIPID_FRAGMENTS_PATH)

    lipid_frgs = list(set(i for l in lipids[0].split(',') for i in frgs_dict[l]))

    # Lipid line
    if len(lipid_frgs) > 1:
        lipid_line = "lipid:" + ",".join(lipid_frgs)
    else:
        lipid_line = "lipid:" + lipid_frgs[0]

    # Ligand line
    ligand_line = None
    if ligands:
        ligand_line = "ligand:" + ligands

    # Ion line (🔥 FIXED HERE)
    ion_line = None
    normalized_ions = process_ions(ions)

    if normalized_ions:
        ion_line = "ions:" + ",".join(normalized_ions)

    # Debug print
    print("Final ion line:", ion_line)

    # Write output
    with open(output_path, 'w') as outfile:
        if lipid_line:
            outfile.write(lipid_line + "\n")
        if ligand_line:
            outfile.write(ligand_line + "\n")
        if ion_line:
            outfile.write(ion_line + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate configuration file for membrane_builder.")
    parser.add_argument('--lipids', nargs='+', required=True, help='List of lipids to include in the membrane.')
    parser.add_argument('--ligands', type=str, help='List of ligands to include in the system.')
    parser.add_argument('--ions', nargs='*', default=[], help='List of ions to include in the system.')
    parser.add_argument('--output', type=Path, required=True, help='Output path for the configuration file.')

    args = parser.parse_args()

    generate_config(args.lipids, args.ligands, args.ions, args.output)