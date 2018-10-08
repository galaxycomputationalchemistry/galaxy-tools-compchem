#!/usr/bin/env python

import argparse
import csv
import sys

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

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
    parser.add_argument('--inbins', help='Number of bins in the histogram')
    parser.add_argument('--istart', help='Starting Point')
    parser.add_argument('--iend', help='End point')
    parser.add_argument('--output', help='output')
    parser.add_argument('--ordf_plot', help='RDF plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)

atom1 = "(segid %s and resid %s and name %s)" % \
    (args.isegid1, args.iresid1, args.iname1)
atom2 = "(segid %s and resid %s and name %s)" % \
    (args.isegid2, args.iresid2, args.iname2)
bins = int(args.inbins)
start = float(args.istart)
end = float(args.iend)

u = mda.Universe(args.ipdb, args.idcd, topology_format="PDB", format="DCD")
x = u.select_atoms(atom1)
y = u.select_atoms(atom2)

rdf = InterRDF(x, y, nbins=bins, range=(start, end))
rdf.run()
bins = rdf.bins
bins = np.around(bins, decimals=3)
RDF = rdf.rdf
zip(bins, RDF)

with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(bins, RDF))

with open(args.output) as f:
    g = [xtmp.strip() for xtmp in f]
    data = [tuple(map(float, xtmp.split())) for xtmp in g[0:]]
    time = [xtmp[0] for xtmp in data]
    rdf = [xtmp[1] for xtmp in data]
    plt.plot(time, rdf)
    plt.xlabel('r ($\AA$)')
    plt.ylabel('g(r)')
    plt.savefig(args.ordf_plot, format='png')
