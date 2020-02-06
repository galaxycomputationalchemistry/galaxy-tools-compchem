import argparse
import sys

import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis import distances


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--list1', help='list 2')
    parser.add_argument('--list2', help='list 2')
    parser.add_argument('--output', help='output')
    parser.add_argument('--header', dest='header', action='store_true')
    return parser.parse_args()


args = parse_command_line(sys.argv)

u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)

list1 = np.loadtxt(args.list1, dtype=str, delimiter="\t", ndmin=1)
list2 = np.loadtxt(args.list2, dtype=str, delimiter="\t", ndmin=1)

sel1 = [u.select_atoms(selection) for selection in list1]
sel2 = [u.select_atoms(selection) for selection in list2]

d = np.empty((u.trajectory.n_frames, list1.shape[0], list2.shape[0]),)

for ts in u.trajectory:
    c_o_m1 = np.array([selection.center_of_mass() for selection in sel1])
    c_o_m2 = np.array([selection.center_of_mass() for selection in sel2])
    distances.distance_array(c_o_m1, c_o_m2, result=d[ts.frame])

d = np.hstack((
    np.array(np.reshape(np.arange(
        0, d.shape[0]), (d.shape[0], 1)), dtype=int),  # add column w frame
    np.reshape(d, (d.shape[0], d.shape[1] * d.shape[2]))
))

if args.header:
    header = 'Frame\t' + '\t'.join(
        ['-'.join(pair) for pair in zip(
            sum([[n, ] * len(list2) for n in list1], []),
            list(list2) * len(list1),)]).replace(' ', '_')
else:
    header = ''

np.savetxt(args.output, d, header=header, comments='',
           fmt=['%d'] + ['%f'] * (d.shape[1] - 1), delimiter='\t')
