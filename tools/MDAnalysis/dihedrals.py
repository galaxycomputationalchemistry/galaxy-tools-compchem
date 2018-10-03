#!/usr/bin/env python

import sys, os
import argparse
import MDAnalysis
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

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
    parser.add_argument('--isegid3', help='segid 3')
    parser.add_argument('--iresid3', help='resid 3')
    parser.add_argument('--iname3', help='name 3')
    parser.add_argument('--isegid4', help='segid 4')
    parser.add_argument('--iresid4', help='resid 4')
    parser.add_argument('--iname4', help='name 4')
    parser.add_argument('--output', help='output')
    parser.add_argument('--odihedral_plot', help='dihedral plot')

    return parser.parse_args()

args = parse_command_line(sys.argv)

atom1 = "(segid %s and resid %s and name %s)" % (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % (args.isegid2, args.iresid2, args.iname2)
atom3 = "(segid %s and resid %s and name %s)" % (args.isegid3, args.iresid3, args.iname3)
atom4 = "(segid %s and resid %s and name %s)" % (args.isegid4, args.iresid4, args.iname4)

def psi(u):
    A = u.select_atoms(atom1).positions
    B = u.select_atoms(atom2).positions
    C = u.select_atoms(atom3).positions
    D = u.select_atoms(atom4).positions

    psi = calc_dihedrals(A, B, C, D)
    return np.rad2deg(psi)

if __name__ == "__main__":
    import MDAnalysis
    from MDAnalysis.lib.distances import calc_dihedrals

    u = MDAnalysis.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")
    data = np.array([(u.trajectory.frame, psi(u)) for ts in u.trajectory])
    frame, psi = data.T

    zip(frame,psi)

    with open(args.output, 'w') as f:
       writer = csv.writer(f, delimiter='\t')
       writer.writerows(zip(frame,psi))


    with open(args.output) as f:
        f=[x.strip() for x in f if x.strip()]
        data=[tuple(map(float,x.split())) for x in f[0:]]
        time=[x[0] for x in data]
        dihedral=[x[1] for x in data]
        plt.plot(time, dihedral)
        plt.xlabel('Frame No.')
        plt.ylabel('Dihedral (degrees)')
        plt.savefig(args.odihedral_plot, format='png')

