#!/usr/bin/env python3
import argparse


def __main__():
    parser = argparse.ArgumentParser(
        description='Adds New topologies to gromacs topology file')
    parser.add_argument(
                        '--top_file', default=None,
                        help="Topology file input")
    parser.add_argument(
                        '--mol_file', default=None,
                        help='molecule type and bonded parameters input')
    parser.add_argument(
                        '--atom_file', default=None,
                        help='atomtype and nonbonded parameters input')
    parser.add_argument(
                        '--out', default=None,
                        help='Path to output')
    args = parser.parse_args()
    with open(args.out, 'w') as fh_out:
        with open(args.top_file, 'r') as fh:
            # these two short loop takes care of
            # adding the atom types and molecule types.
            for line in fh:
                fh_out.write(line)
                if ';name   bond_type' in line:
                    for contents in open(args.atom_file):
                        fh_out.write(contents)
                    break
            for line in fh:
                if '[ system ]' in line:
                    fh_out.write("\n; Begin NewTopologyInfo\n")
                    for contents in open(args.mol_file):
                        fh_out.write(contents)
                    fh_out.write("; end NewTopologyInfo\n\n")
                    fh_out.write(line)
                    break
                fh_out.write(line)
            for line in fh:
                fh_out.write(line)


if __name__ == "__main__":
    __main__()
