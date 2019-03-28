#!/usr/bin/env python

import argparse
import sys

import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca

import numpy as np
import csv


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--idcd', help='input dcd')
    parser.add_argument('--ipdb', help='input pdb')
    parser.add_argument('--icomponents', help='number of principle components')
    parser.add_argument('--iindex', help='index of the PC')
    parser.add_argument('--output', help='output')
    parser.add_argument('--cosout', help='cosine output')
    return parser.parse_args()


args = parse_command_line(sys.argv)

u = mda.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")

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
