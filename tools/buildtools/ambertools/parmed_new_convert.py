import sys
import argparse
import io
import parmed
from parmed import gromacs, amber, unit as u
from parmed.tools.changeradii import ChRad
from contextlib import redirect_stdout


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--istr', help='input structure', required=True)
    parser.add_argument('--itop', help='input topology file', required=True)
    parser.add_argument('--istripmask', help='stripmask')
    parser.add_argument('--iradii', help='parmed radii are GB_RADII amber6, bondi, mbondi, mbondi2, mbondi3 - see https://parmed.github.io/ParmEd/html/_modules/parmed/tools/changeradii.html', required=True)
    parser.add_argument('--removedihe', action='store_true', default=False, help='remove dihedrals with zero periodicity')
    parser.add_argument('--removebox', action='store_true', default=False, help='remove periodic box info')
    parser.add_argument('--o_prmtop', help='AMBER output topology', required=True)
    return parser.parse_args()

def get_ids(dihedrals):
    """
    goes through dihedrals and looks for any with per=0. 
    returns a reverse sorted list of ids to be removed.
    """
    indices = []
    for k,v in enumerate(dihedrals):
        f = io.StringIO()
        with redirect_stdout(f):
            print(v)
        if f.getvalue().find("per=0") != -1:
            indices.append(k)
    indices.sort(reverse=True)
    return indices


args = parse_command_line(sys.argv)

gmx_top = gromacs.GromacsTopologyFile(args.itop)
gmx_gro = gromacs.GromacsGroFile.parse(args.istr)

if not args.removebox:
    # keep box info
    gmx_top.box = gmx_gro.box 
    gmx_top.positions = gmx_gro.positions


if args.removedihe:
    ids_to_remove = get_ids(gmx_top.dihedrals)
    print("Original number of dihedrals %i" % len(gmx_top.dihedrals))
    for i in ids_to_remove:
        gmx_top.dihedrals.pop(i)
    print("Update number of dihedrals %i" % len(gmx_top.dihedrals))

if args.istripmask is not None:
   if args.istripmask == "":
       pass
   else:
       gmx_top.strip(args.istripmask)

radii=str(args.iradii)
parmed.tools.changeRadii(gmx_top,radii)
amb_prm = amber.AmberParm.from_structure(gmx_top)
parmed.tools.changeRadii(amb_prm,radii)

if args.removebox:
    amb_prm.pointers['IFBOX'] = 0 


ChRad(amb_prm, radii)
for i, atom in enumerate(amb_prm.atoms):
    amb_prm.parm_data['RADII'][i] = atom.solvent_radius
    amb_prm.parm_data['SCREEN'][i] = atom.screen


amb_prm.write_parm(args.o_prmtop)
