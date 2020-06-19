#!/usr/bin/env python

import argparse
import csv
import sys

import itertools
import MDAnalysis as mda

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import numpy as np
import numpy.linalg

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--isegid1', help='segid 1')
    parser.add_argument('--ilabel', help='plot label')
    parser.add_argument('--ititle1', help='plot title')
    parser.add_argument('--output1', help='output1 - timeseries')
    parser.add_argument('--o_plot', help='End to End plot')
    return parser.parse_args()


args = parse_command_line(sys.argv)



u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)

ntermatoms = "(segid %s and name N)" % \
    (args.isegid1)
ctermatoms = "(segid %s and name C)" % \
    (args.isegid1)
# not sure how robust this selection really is 
nterm = u.select_atoms(ntermatoms)[0] # first atom named N
cterm = u.select_atoms(ctermatoms)[-1] # takes the last atom named 'C'
#print(nterm, cterm)




enddist = []

for ts in u.trajectory:  # iterate through all frames
    r = cterm.position - nterm.position  # end-to-end vector from atom positionitions
    d = numpy.linalg.norm(r)   # end-to-end distance
    enddist.append((ts.frame, d)) 

enddist = np.array(enddist)


color = itertools.cycle(['r', 'b', 'gold'])

fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, tight_layout=True)

params = {
   'axes.labelsize': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [4.5, 4.5],
   'figure.dpi':300
   }
plt.rcParams.update(params)

axs[0].plot(enddist[:,0], enddist[:,1], 'r-', lw=2, label=args.ilabel)
axs[0].set_xlabel("number of frames")
axs[0].set_ylabel(r"End to end distance  ($\AA$)")
axs[0].legend()

n, bins, patches = axs[1].hist(enddist[:,1], color=next(color), label=args.ilabel, alpha=0.5, density=True, stacked=True) #, bins=20)

axs[1].legend()
axs[1].set_ylabel('Density Normalised Frequency');
axs[1].set_xlabel(r'End to end distance ($\AA$)')
fig.suptitle(args.ititle1, fontsize=12, fontweight='bold')
fig.subplots_adjust(top=0.45)

print(" \n".join(['The End to End distance is measured between the following atoms:',str(nterm),str(cterm)]))

plt.savefig(args.o_plot, format='png') # svg is better but sticking with png for now


np.savetxt(args.output1, enddist, delimiter='\t')
