#!/usr/bin/env python3
import sys

# Check if arguments are provided
if len(sys.argv) < 3:
    print("Usage: python fixed.py <input_file> <output_file>")
    print("Example: python fixed.py wats.pdb fixed_pdb.pdb")
    sys.exit(1)

# Get arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

file = open(input_file, "r")
out_file = open(output_file, "w")

lines = file.readlines()
file.close()

for line in lines:
    splitted = line.split('WAT')
    fixed_line = splitted[0]+splitted[1]+"WAT"+splitted[2]+"\n"
    out_file.write(fixed_line )

out_file.close()
print(f"✓ Processing complete. Output saved to: {output_file}")




