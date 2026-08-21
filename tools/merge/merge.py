#!/usr/bin/env python3
"""
Merge PDB files by extracting matching residues and appending to System.pdb.
"""

import argparse
import sys


def parse_pdb_line(line):
    """Parse a PDB ATOM/HETATM line and extract key fields."""
    if len(line) < 54:
        return None
    
    try:
        record_name = line[0:6].strip()      # ATOM or HETATM
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21].strip()
        res_num = int(line[22:26].strip())
        
        return {
            'record': record_name,
            'atom_num': atom_num,
            'atom_name': atom_name,
            'res_name': res_name,
            'chain': chain,
            'res_num': res_num,
            'line': line
        }
    except:
        return None


def format_water_line(line):
    """Format a water residue line: HOH→WAT and add proper TIP3 element column."""
    if len(line) < 80:
        # Pad line to at least 80 characters
        line = line.ljust(80)
    
    # Replace HOH with WAT at columns 17-20
    line_list = list(line)
    if len(line) >= 20:
        line_list[17:20] = list("WAT")
    
    # Replace element column (66-80) with "      TIP3    " (proper spacing for water)
    if len(line_list) >= 76:
        element_str = "      TIP3    "
        line_list[66:80] = list(element_str.ljust(14))
    
    return ''.join(line_list).rstrip() + '\n'


def merge_pdb_files(water_file, output_file=None, 
                    residue_nums=None, residue_names=None, system_file=None):
    """
    Extract matching ATOM/HETATM lines from water file and append to System.pdb.
    Converts HOH to WAT and formats element column with TIP3.
    Adds TER record after each residue.
    
    Args:
        water_file: Path to water/ligand PDB file to extract from
        output_file: Path to save merged PDB file
        residue_nums: List of residue numbers to match
        residue_names: List of residue names to match
        system_file: Path to System.pdb file to append to
    """
    
    if residue_nums is None:
        residue_nums = []
    if residue_names is None:
        residue_names = []
    
    # Ensure residue numbers are integers
    residue_nums = [int(r) if isinstance(r, str) else r for r in residue_nums]
    
    print(f"Extracting from water file: {water_file}")
    print(f"Searching for residues: {residue_nums} with names: {residue_names}")
    
    # Read water file
    try:
        with open(water_file, 'r') as f:
            water_lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found: {water_file}")
        return False
    
    # Read system file if provided
    system_lines = []
    if system_file:
        try:
            with open(system_file, 'r') as f:
                system_lines = f.readlines()
            print(f"Reading system file: {system_file}")
        except FileNotFoundError:
            print(f"Error: File not found: {system_file}")
            return False
    
    # Extract matching lines from water file
    extracted_lines = []
    matched_count = 0
    previous_res_num = None
    
    for line in water_lines:
        # Only process ATOM and HETATM lines
        if line.startswith('ATOM') or line.startswith('HETATM'):
            parsed = parse_pdb_line(line)
            
            if parsed:
                # Check if this line matches our criteria
                res_num_match = parsed['res_num'] in residue_nums
                res_name_match = parsed['res_name'] in residue_names
                
                # If residue_nums or residue_names are empty, match all
                if (not residue_nums or res_num_match) and \
                   (not residue_names or res_name_match):
                    
                    # Add TER if residue number changed
                    if previous_res_num is not None and previous_res_num != parsed['res_num']:
                        extracted_lines.append("TER\n")
                    
                    # Format the line: HOH→WAT and TIP3 element column
                    formatted_line = format_water_line(line)
                    extracted_lines.append(formatted_line)
                    matched_count += 1
                    previous_res_num = parsed['res_num']
                    print(f"✓ Matched: Res {parsed['res_num']} ({parsed['res_name']}) - {parsed['atom_name']}")
    
    # Add final TER after extracted atoms
    if matched_count > 0:
        extracted_lines.append("TER\n")
    
    if matched_count == 0:
        print("No matching residues found!")
        return False
    
    # Merge with system file
    if system_lines:
        # Remove END line from system file if present
        output_lines = [line for line in system_lines if not line.startswith('END')]
        # Append extracted lines
        output_lines.extend(extracted_lines)
        # Add END at the end
        output_lines.append("END\n")
    else:
        # Just use extracted lines with END
        output_lines = extracted_lines + ["END\n"]
    
    # Write output file
    if output_file is None:
        output_file = 'merged.pdb'
    
    try:
        with open(output_file, 'w') as f:
            f.writelines(output_lines)
        print(f"\n✓ Extracted {matched_count} matching atoms")
        print(f"✓ Saved to: {output_file}")
        return True
    except Exception as e:
        print(f"Error writing output file: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Extract matching residues and append to System.pdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract residues 293 and 309 with name HOH
  python merge.py in.wat.pdb -n 293 309 -r HOH
  
  # Extract and append to System.pdb
  python merge.py in.wat.pdb -n 293 309 -r HOH -s System.pdb -o System_with_water.pdb
  
  # Extract and save to extracted.pdb
  python merge.py in.wat.pdb -n 293 309 -r HOH -o extracted.pdb
        """
    )
    
    parser.add_argument('water_file', help='Water/ligand PDB file to extract from')
    parser.add_argument('-n', '--residues', type=int, nargs='+', 
                        help='Residue numbers to match (e.g., 293 309)')
    parser.add_argument('-r', '--resnames', nargs='+', 
                        help='Residue names to match (e.g., HOH WAT)')
    parser.add_argument('-o', '--output', help='Output file (default: extracted.pdb)')
    parser.add_argument('-s', '--system', help='System.pdb file to append extracted residues to')
    
    args = parser.parse_args()
    
    success = merge_pdb_files(
        args.water_file,
        args.output,
        residue_nums=args.residues,
        residue_names=args.resnames,
        system_file=args.system
    )
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()




