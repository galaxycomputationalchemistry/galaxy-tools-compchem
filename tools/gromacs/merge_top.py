import argparse

import parmed as pmd


def merge_gro_files(prot_gro, lig_gro, cmplx_gro):
    prot = pmd.load_file(prot_gro)
    lig = pmd.load_file(lig_gro)
    cmplx = prot + lig
    cmplx.save(cmplx_gro)


def merge_top_files(prot_top, lig_top, cmplx_top):
    with open(lig_top, 'r') as f:
        lig_top_sections = f.read().split('\n[')

    # open ligand topology
    for n in range(len(lig_top_sections)):
        if 'atomtypes' in lig_top_sections[n][:10]:
            lig_atomtypes = lig_top_sections[n]
            del lig_top_sections[n]
            break
    else:
        lig_atomtypes = None
    lig_top_updated = '\n['.join(lig_top_sections)

    # open protein topology
    with open(prot_top, 'r') as f:
        prot_top_combined = f.read()
    if lig_atomtypes:
        prot_top_sections = prot_top_combined.split('[ moleculetype ]\n')
        prot_top_combined = (prot_top_sections[0] +
                             '; Include ligand atomtypes\n[' +
                             lig_atomtypes +
                             '\n[ moleculetype ]\n' +
                             prot_top_sections[1])
    prot_top_sections = prot_top_combined.split('; Include water topology')
    prot_top_combined = (prot_top_sections[0] +
                         '; Include ligand topology\n' +
                         lig_top_updated +
                         '\n; Include water topology' +
                         prot_top_sections[1])
    prot_top_combined += 'base     1\n'

    # save complex topology
    with open(cmplx_top, 'w') as f:
        f.write(prot_top_combined)


def main():
    parser = argparse.ArgumentParser(
        description='Perform SMD runs for dynamic undocking')
    parser.add_argument('--lig-top', help='Ligand TOP file.')
    parser.add_argument('--prot-top', help='Protein TOP file.')
    parser.add_argument('--lig-gro', help='Ligand GRO file.')
    parser.add_argument('--prot-gro', help='Protein GRO file.')
    parser.add_argument('--complex-top', help='Complex TOP file.')
    parser.add_argument('--complex-gro', help='Complex GRO file.')
    args = parser.parse_args()
    merge_gro_files(args.prot_gro, args.lig_gro, args.complex_gro)
    merge_top_files(args.prot_top, args.lig_top, args.complex_top)


if __name__ == "__main__":
    main()
