#!/usr/bin/env python

import argparse
import csv
import sys

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import pylab

from collections import namedtuple

import numpy as np

import seaborn as sns


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
    parser.add_argument('--isegid5', help='segid 1')
    parser.add_argument('--iresid5', help='resid 1')
    parser.add_argument('--iname5', help='name 1')
    parser.add_argument('--isegid6', help='segid 2')
    parser.add_argument('--iresid6', help='resid 2')
    parser.add_argument('--iname6', help='name 2')
    parser.add_argument('--isegid7', help='segid 3')
    parser.add_argument('--iresid7', help='resid 3')
    parser.add_argument('--iname7', help='name 3')
    parser.add_argument('--isegid8', help='segid 4')
    parser.add_argument('--iresid8', help='resid 4')
    parser.add_argument('--iname8', help='name 4')
    parser.add_argument('--output', help='output')
    parser.add_argument('--oramachandran_plot', help='dihedral plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)

Dihedral = namedtuple(
    'Dihedral', ['atom1', 'atom2', 'atom3', 'atom4'])

# order of dihedral atom is the crystallographic definition
# (see glycanstructure.org)

# phi
atom1 = "(segid %s and resid %s and name %s)" % \
    (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % \
    (args.isegid2, args.iresid2, args.iname2)
atom3 = "(segid %s and resid %s and name %s)" % \
    (args.isegid3, args.iresid3, args.iname3)
atom4 = "(segid %s and resid %s and name %s)" % \
    (args.isegid4, args.iresid4, args.iname4)

dihe_phi = Dihedral(atom1, atom2, atom3, atom4)

# psi
atom1 = "(segid %s and resid %s and name %s)" % \
    (args.isegid5, args.iresid5, args.iname5)
atom2 = "(segid %s and resid %s and name %s)" % \
    (args.isegid6, args.iresid6, args.iname6)
atom3 = "(segid %s and resid %s and name %s)" % \
    (args.isegid7, args.iresid7, args.iname7)
atom4 = "(segid %s and resid %s and name %s)" % \
    (args.isegid8, args.iresid8, args.iname8)

dihe_psi = Dihedral(atom1, atom2, atom3, atom4)


def calc_torsion(dihedral):
    """atom 1 -4 are valid atom selections. torsion in degrees is returned"""
    A = u.select_atoms(dihedral.atom1).positions
    B = u.select_atoms(dihedral.atom2).positions
    C = u.select_atoms(dihedral.atom3).positions
    D = u.select_atoms(dihedral.atom4).positions

    dihe = calc_dihedrals(A, B, C, D)
    return np.rad2deg(dihe)


u = mda.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")

phi_trajdata = np.array(
    [(u.trajectory.frame, calc_torsion(dihe_phi)) for ts in u.trajectory])
psi_trajdata = np.array(
    [(u.trajectory.frame, calc_torsion(dihe_psi)) for ts in u.trajectory])

phi_frame, phi_series = phi_trajdata.T
psi_frame, psi_series = psi_trajdata.T

phi_series = np.concatenate(phi_series, axis=0)
psi_series = np.concatenate(psi_series, axis=0)

zip(phi_frame, phi_series, psi_series)

with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(phi_frame, phi_series, psi_series))

with sns.axes_style("white"):
    figure()
    h = sns.jointplot(x=phi_series, y=psi_series, kind="kde", legend=True)
    h.set_axis_labels(r'$\Phi$ (degrees)', r'$\Psi$ (degrees)')
    h.ax_joint.set_xlim(-180, 180)
    h.ax_joint.set_ylim(-180, 180)
    plt.savefig(args.oramachandran_plot, format='png')
