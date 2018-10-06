#!/usr/bin/env python

import argparse
import csv
import sys

import MDAnalysis as mda

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
from numpy.linalg import norm

matplotlib.use('Agg')


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
    parser.add_argument('--output', help='output')
    parser.add_argument('--oangle_plot', help='angle plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)

atom1 = "(segid %s and resid %s and name %s)" % \
    (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % \
    (args.isegid2, args.iresid2, args.iname2)
atom3 = "(segid %s and resid %s and name %s)" % \
    (args.isegid3, args.iresid3, args.iname3)


def theta(u):
    A = u.select_atoms(atom1).center_of_geometry()
    B = u.select_atoms(atom2).center_of_geometry()
    C = u.select_atoms(atom3).center_of_geometry()
    BA = A - B
    BC = C - B
    theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
    return np.rad2deg(theta)


u = mda.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")
data = np.array([(u.trajectory.frame, theta(u)) for ts in u.trajectory])
frame, theta = data.T

with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(frame, theta))

with open(args.output) as f:
    g = [xtmp.strip() for xtmp in f]
    data = [tuple(map(float, xtmp.split())) for xtmp in g[0:]]
    time = [xtmp[0] for xtmp in data]
    angle = [xtmp[1] for xtmp in data]
    plt.plot(time, angle)
    plt.xlabel('Frame No.')
    plt.ylabel('Angle (degrees)')
    plt.savefig(args.oangle_plot, format='png')
