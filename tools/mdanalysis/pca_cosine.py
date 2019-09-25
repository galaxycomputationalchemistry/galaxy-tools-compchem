#!/usr/bin/env python

import argparse
import csv
import sys

import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca

import numpy as np


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--icomponents', help='number of principle components')
    parser.add_argument('--iindex', help='index of the PC')
    parser.add_argument('--output', help='output')
    parser.add_argument('--cosout', help='cosine output')
    return parser.parse_args()


args = parse_command_line(sys.argv)

u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)

components = int(args.icomponents)
pca_index = int(args.iindex)

PSF_pca = pca.PCA(u, select='backbone')
PSF_pca.run()
n_pcs = np.where(PSF_pca.cumulated_variance > 0.95)[0][0]
atomgroup = u.select_atoms('backbone')

pca_space = PSF_pca.transform(atomgroup, n_components=components)
cosine = mda.analysis.pca.cosine_content(pca_space, pca_index)

PCA = list(pca_space)

with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(PCA)

with open(args.cosout, 'w') as f1:
    f1.write(str(cosine))
