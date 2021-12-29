#!/usr/bin/env python3
import argparse
END_OF_MOLECULE = ('[ moleculetype ]', '[ system ]')


def __main__():
    parser = argparse.ArgumentParser(
        description='Add restriction to gromacs topology file')
    parser.add_argument(
                        '--top_file', default=None,
                        help="Topology file input")
    parser.add_argument(
                        '--res_file', default=None,
                        help='Restraint input')
    parser.add_argument(
                        '--molecule', default=None,
                        help='Target Molecule Name you restrained')
    parser.add_argument(
                        '--out', default=None,
                        help='Path to output')
    args = parser.parse_args()
    with open(args.out, 'w') as fh_out:
        with open(args.top_file, 'r') as fh:
            # for now, we will avoid using 'for line in fh:',
            # since we have multiple places where we might want
            # to read the next line
            while True:
                line = fh.readline()
                if not line:
                    # eof
                    break
                # always write out the line
                fh_out.write(line)
                # check if line matches molecule, then check if
                # molecule name matches args.molecule
                if line.strip().startswith('[ moleculetype ]'):
                    not_found_molecule = True
                    while not_found_molecule:
                        line = fh.readline()
                        if not line:
                            # eof
                            break
                        # always write this line
                        fh_out.write(line)
                        if not line.strip().startswith(';') or (line.strip() and not line.strip().startswith(';')):
                            # this line should be the name line,
                            fields = line.strip().split()
                            if fields[0] == args.molecule:
                                # found our molecule!
                                while True:
                                    line = fh.readline()
                                    if not line:
                                        # eof
                                        break
                                    if line.strip().startswith(END_OF_MOLECULE):
                                        fh_out.write("\n#ifdef POSRES\n")
                                        with open(args.res_file, 'r') as fh_res:
                                            for line2 in fh_res:
                                                fh_out.write(line2)
                                        fh_out.write("#endif\n\n")
                                        fh_out.write(line)
                                        not_found_molecule = False
                                        break
                                    fh_out.write(line)
                            else:
                                break


if __name__ == "__main__":
    __main__()
