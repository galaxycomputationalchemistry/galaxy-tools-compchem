#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import MDAnalysis
import numpy.linalg

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--idcd', help='input dcd')
    parser.add_argument('--ipdb', help='input pdb')
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

atom1 = "(segid %s and resid %s and name %s)" % (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % (args.isegid2, args.iresid2, args.iname2)

u = MDAnalysis.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")
x = u.select_atoms(atom1)
y = u.select_atoms(atom2)

with open(args.output, 'w') as f:
   for t in u.trajectory:
       r = x.positions - y.positions
       d = numpy.linalg.norm(r)
       f.write(str(t.frame)+'\t ')
       f.write(str(d)+'\n')

with open(args.output) as f:
     f=[x.strip() for x in f if x.strip()]
     data=[tuple(map(float,x.split())) for x in f[0:]]
     time=[x[0] for x in data]
     distance=[x[1] for x in data]
     plt.plot(time, distance)
     plt.xlabel('Frame No.')
     plt.ylabel('Distance ($\AA$)')
     plt.savefig(args.odistance_plot, format='png')


