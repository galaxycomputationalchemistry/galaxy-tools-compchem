#!/usr/bin/env python3

import subprocess
import sys
import argparse
import tempfile
import os
import re


def generate_command_for_index_file(args):
    """
    Generate GROMACS index file with custom groups based on configuration.
    
    Creates two main combined groups:
    1. Protein + all lipids + all ligands (e.g., Protein_CHL_PA_PC_CAU)
    2. Water + all ions (e.g., Water_NA+_CL-)
    
    Args:
        args: Argparse namespace with config, gro_file, and output attributes
        
    Returns:
        Tuple of (tc_grps, comm_grps) as space-separated strings for temperature and COM coupling
    """
    
    monovalent_ions = [
        "Na+", "NA+",
        "K+",  "K+",
        "Li+", "LI+",
        "Rb+", "RB+",
        "Cs+", "CS+",
        "Cl-", "CL-",
        "F-",  "F-",
        "Br-", "BR-",
        "I-",  "I-",
    ]

    tc_grps = []
    comm_grps = []

    index_dict = {}
    extracted_lines = []

    lipid_frgs_conf = []
    ligand_conf = []
    ions = []

    # Read configuration file for lipid, ligand, and ion definitions
    with open(args.config, 'r') as f:
        conf_lines = f.readlines()

    # Run gmx make_ndx to get the list of available groups
    try:
        with tempfile.NamedTemporaryFile(suffix='.ndx', delete=False) as tmp:
            tmp_ndx = tmp.name
        
        result = subprocess.run(
            f"echo 'q' | gmx make_ndx -f {args.gro_file} -o {tmp_ndx}",
            shell=True,
            capture_output=True,
            text=True
        )
        index_log_lines = result.stdout.splitlines() + result.stderr.splitlines()
        
        # Clean up temporary file
        if os.path.exists(tmp_ndx):
            os.remove(tmp_ndx)
    except subprocess.CalledProcessError as e:
        print(f"Error running initial gmx make_ndx: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse configuration file
    # Format expected:
    #   lipid:PA,PC,OL,CHL
    #   ligand:CAU
    for l in conf_lines:
        line = l.strip()
        # Skip empty lines and comments
        if not line or line.startswith(';') or line.startswith('#'):
            continue
        # Remove leading "-" if present (for backward compatibility)
        line = line.lstrip('- ')
        
        if ":" in line:
            key, value = line.split(":", 1)
            key = key.strip().lower()
            if key == "lipid":
                lipid_frgs_conf = [x.strip() for x in value.split(',') if x.strip()]
            elif key == "ligand":
                ligand_conf = [x.strip() for x in value.split(',') if x.strip()]

    # Parse index_log to build index_dict mapping names to index numbers
    for i in index_log_lines:
        parts = i.split()
        if len(parts) == 5 and parts[2] == ":" and parts[4] == 'atoms':
            extracted_lines.append(parts)
            index_dict[parts[1]] = parts[0]
            
            # Identify ions by matching against monovalent_ions list
            if parts[1].strip() in monovalent_ions:
                ions.append(parts[1].strip())

    # Check if we found any groups
    if not extracted_lines:
        print("Error: Could not parse any index groups from gmx make_ndx output", file=sys.stderr)
        sys.exit(1)

    # Get the highest index number to calculate new indices
    end_ndx = int(extracted_lines[-1][0])

    # Filter to only include groups that exist in the system
    lipid_frgs = [f for f in lipid_frgs_conf if f in index_dict]
    ligands = [l for l in ligand_conf if l in index_dict]
    ions = [ion for ion in ions if ion in index_dict]

    print(f"Available index groups: {list(index_dict.keys())}")
    if lipid_frgs:
        print(f"Found lipid groups: {lipid_frgs}")
    if ligands:
        print(f"Found ligand groups: {ligands}")
    if ions:
        print(f"Found ion groups: {ions}")

    # Build commands list
    command_list = []
    
    # =========================================================================
    # GROUP 1: Protein + Lipids + Ligands combined (e.g., Protein_CHL_PA_PC_CAU)
    # =========================================================================
    all_membrane_components = lipid_frgs + ligands
    
    if all_membrane_components:
        # Get indices for all components
        component_indices = [index_dict[c] for c in all_membrane_components]
        
        if len(all_membrane_components) == 1:
            # Single component: Protein | component
            command_list.append(f"echo '{index_dict['Protein']} | {component_indices[0]}'")
            end_ndx += 1
            group_name = f"Protein_{all_membrane_components[0]}"
        else:
            # Multiple components: first combine all components, then combine with Protein
            # Step 1: Combine all lipids and ligands
            command_list.append(f"echo '{' | '.join(component_indices)}'")
            end_ndx += 1
            combined_components_idx = end_ndx
            
            # Step 2: Combine Protein with combined components
            command_list.append(f"echo '{index_dict['Protein']} | {combined_components_idx}'")
            end_ndx += 1
            group_name = f"Protein_{'_'.join(all_membrane_components)}"
        
        tc_grps.append(group_name)
        comm_grps.append(group_name)
        print(f"Creating group: {group_name}")
    else:
        # No lipids or ligands, just use Protein
        tc_grps.append("Protein")
        comm_grps.append("Protein")

    # =========================================================================
    # GROUP 2: Water + Ions combined (e.g., Water_NA+_CL-)
    # =========================================================================
    if ions:
        ion_indices = [index_dict[ion] for ion in ions]
        
        if len(ions) == 1:
            # Single ion: Water | ion
            command_list.append(f"echo '{index_dict['Water']} | {ion_indices[0]}'")
            end_ndx += 1
            group_name = f"Water_{ions[0]}"
        else:
            # Multiple ions: first combine all ions, then combine with Water
            # Step 1: Combine all ions
            command_list.append(f"echo '{' | '.join(ion_indices)}'")
            end_ndx += 1
            combined_ions_idx = end_ndx
            
            # Step 2: Combine Water with combined ions
            command_list.append(f"echo '{index_dict['Water']} | {combined_ions_idx}'")
            end_ndx += 1
            group_name = f"Water_{'_'.join(ions)}"
        
        tc_grps.append(group_name)
        comm_grps.append(group_name)
        print(f"Creating group: {group_name}")
    else:
        # No ions, just use Water if it exists
        if "Water" in index_dict:
            tc_grps.append("Water")
            comm_grps.append("Water")

    # Build final gmx make_ndx command
    if command_list:
        cmd_string = f"({' ; '.join(command_list)} ; echo 'q') | gmx make_ndx -f {args.gro_file} -o {args.output}"
    else:
        # No custom groups to create, just generate the default index file
        cmd_string = f"echo 'q' | gmx make_ndx -f {args.gro_file} -o {args.output}"

    print("\n" + "=" * 80)
    print("Generated command:")
    print(cmd_string)
    print("=" * 80 + "\n")

    # Execute the command
    try:
        result = subprocess.run(cmd_string, shell=True, check=True, capture_output=True, text=True)
        print(f"\nIndex file successfully created: {args.output}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing gmx make_ndx: {e}", file=sys.stderr)
        sys.exit(1)

    # Now show the final index groups by reading the created index file
    print("\n" + "=" * 80)
    print("FINAL INDEX GROUPS IN GENERATED FILE:")
    print("=" * 80)
    
    try:
        # Run gmx make_ndx with the new index file to list all groups
        result = subprocess.run(
            f"echo 'q' | gmx make_ndx -f {args.gro_file} -n {args.output} -o {args.output}",
            shell=True,
            capture_output=True,
            text=True
        )
        
        # Parse and display the group list using regex for reliable parsing
        all_output = result.stdout + result.stderr
        # Pattern matches: "  0 System              : 43365 atoms" or "24 Protein_PA_PC_OL_CHL_CAU: 22160 atoms"
        pattern = r'^\s*(\d+)\s+(\S+)\s*:\s*(\d+)\s+atoms'
        for line in all_output.splitlines():
            match = re.match(pattern, line)
            if match:
                idx, name, atoms = match.groups()
                print(f"  {idx:>3} {name:<30} : {atoms:>6} atoms")
                
    except subprocess.CalledProcessError as e:
        print(f"Warning: Could not list final index groups: {e}")

    print("=" * 80)
    print(f"\nTemperature coupling groups (tc_grps): {' '.join(tc_grps)}")
    print(f"COM motion removal groups (comm_grps): {' '.join(comm_grps)}")

    return " ".join(tc_grps), " ".join(comm_grps)


def main(args):
    """Main entry point."""
    generate_command_for_index_file(args)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automate GROMACS index file creation.")
    parser.add_argument("-c", "--config", default="index_groups.txt", help="Configuration file (default: index_groups.txt)")
    parser.add_argument("-f", "--gro-file", default="system.gro", help="Input GRO file (default: system.gro)")
    parser.add_argument("-o", "--output", default="ndx.ndx", help="Output index file (default: ndx.ndx)")
    args = parser.parse_args()
    main(args)