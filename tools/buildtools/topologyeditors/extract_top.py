#!/usr/bin/env python3
import argparse


def __main__():
    parser = argparse.ArgumentParser(
        description='Adds New topologies to gromacs topology file')
    parser.add_argument(
                        '--top_file', default=None,
                        help="Topologies input")
    parser.add_argument(
                        '--out_bondparam', default=None,
                        help='moleculetype section')
    parser.add_argument(
                        '--out_nonbondparam', default=None,
                        help='atomtypes section')

    args = parser.parse_args()
    # extracts the atom types with nonbonded terms from
    # the new molecules and puts them in a new file
    inFile = open(args.top_file)
    outFile = open(args.out_nonbondparam, "w")
    buffer = []
    for line in inFile:
        if line.startswith(";name   bond_type"):
            buffer = ['']
        elif line.startswith("[ moleculetype ]"):
            outFile.write("".join(buffer))
            buffer = []
        elif buffer:
            buffer.append(line)
    inFile.close()
    outFile.close()

    # extracts the molecule types (rest of the force field parameters)
    # with bonded terms and puts them in a new file
    inFile = open(args.top_file)
    outFile = open(args.out_bondparam, "w")
    buffer = []
    for line in inFile:
        if line.startswith("[ moleculetype ]"):
            buffer = ["\n[ moleculetype ]\n"]
        elif line.startswith("[ system ]"):
            outFile.write("".join(buffer))
            buffer = []
        elif buffer:
            buffer.append(line)
    inFile.close()
    outFile.close()


if __name__ == "__main__":
    __main__()
