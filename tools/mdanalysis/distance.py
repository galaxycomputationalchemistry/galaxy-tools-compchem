#!/usr/bin/env python

import argparse
import sys

import MDAnalysis as mda

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import numpy as np


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--isegid1', help='segid 1')
    parser.add_argument('--iresid1', help='resid 1')
    parser.add_argument('--iname1', help='name 1')
    parser.add_argument('--isegid2', help='segid 2')
    parser.add_argument('--iresid2', help='resid 2')
    parser.add_argument('--iname2', help='name 2')
    parser.add_argument('--output', help='output')
    parser.add_argument('--odistance_plot', help='odistance plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)

atom1 = "(segid %s and resid %s and name %s)" % \
    (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % \
    (args.isegid2, args.iresid2, args.iname2)

u = mda.Universe(args.istr, args.itraj, 
        topology_format=args.istrext, format=args.itrajext)
x = u.select_atoms(atom1)
y = u.select_atoms(atom2)

with open(args.output, 'w') as f:
    for t in u.trajectory:
        r = x.positions - y.positions
        d = np.linalg.norm(r)
        f.write(str(t.frame) + '\t ')
        f.write(str(d) + '\n')

with open(args.output) as f:
    g = [xtmp.strip() for xtmp in f]
    data = [tuple(map(float, xtmp.split())) for xtmp in g[0:]]
    time = [xtmp[0] for xtmp in data]
    distance = [xtmp[1] for xtmp in data]
    plt.plot(time, distance)
    plt.xlabel('Frame No.')
    plt.ylabel(r'Distance ($\AA$)')
    plt.savefig(args.odistance_plot, format='png')
