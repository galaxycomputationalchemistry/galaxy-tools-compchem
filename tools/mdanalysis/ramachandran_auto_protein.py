#!/usr/bin/env python

import argparse
import csv
import sys

import itertools
import MDAnalysis as mda

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
import numpy.linalg
from MDAnalysis.analysis.dihedrals import Ramachandran

import numpy.linalg
import numpy as np

import seaborn as sns
import importlib

import os
import sys
import json

import h5py
import base64
from jinja2 import Environment, FileSystemLoader


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--itraj', help='input traj')
    parser.add_argument('--istr', help='input str')
    parser.add_argument('--itrajext', help='input traj ext')
    parser.add_argument('--istrext', help='input str ext')
    parser.add_argument('--isegid1', help='segid 1')
    parser.add_argument('--iresid1', help='resid start')
    parser.add_argument('--iresid2', help='resid end')
    parser.add_argument('--iresname', help='resname e.g. ALA')
    parser.add_argument('--igroupby', help='groupby names or ids')
    parser.add_argument('--itemplatepath', help='template path')
    parser.add_argument('--o_plot1', help='MDA Ramachandran plot')
    parser.add_argument('--o_plot2', help='Seaborn Ramachandran plot')
    parser.add_argument('--o_data1', help='Timeseries in HDF5 format')
    parser.add_argument('--o_html1', help='Html overview output of all plots')
    return parser.parse_args()


args = parse_command_line(sys.argv)

currentpath = "."
if args.itemplatepath is not None:
    currentpath = args.itemplatepath


u = mda.Universe(args.istr, args.itraj,
                 topology_format=args.istrext, format=args.itrajext)
selection = "(segid %s)" % \
    (args.isegid1)

if args.iresname is not None:
    selection = "(segid %s and resname %s)" % \
        (args.isegid1, args.iresname)

if args.iresid1 is not None and args.iresid2 is not None:
    assert(int(args.iresid1) > 0), "ResID numbering starts at 1 for this tool."
    assert(int(args.iresid2) > 0), "ResID numbering starts at 1 for this tool."
    assert(int(args.iresid2) > int(args.iresid1)
           ), "ResID2 must be at least ResID1+1"
    selection = "(segid %s and resid %s-%s)" % \
        (args.isegid1, int(args.iresid1), int(args.iresid2))
    if args.iresname is not None:
        selection = "(segid %s and resid %s-%s and resname %s)" % \
            (args.isegid1, int(args.iresid1), int(args.iresid2), args.iresname)

r = u.select_atoms(selection)

assert(r != u.select_atoms('name thiscannotpossiblyexist')
       ), "The selection you specified returns an empty result. Check segment names and residue ID's. Also check the structure and trajectory file selected are the the correct ones"

if args.igroupby is not None:
    group_selections = {}  # dictionary of selections
    if args.igroupby == 'name':
        groupby = sorted(list(set(r.resnames)))
        for e in groupby:
            s = r & u.select_atoms("resname %s" % e)
            this_sel = "%s and resname %s" % (selection, e)
            group_selections[this_sel] = s
    elif args.igroupby == 'id':
        groupby = sorted(list(set(r.resids)))
        for e in groupby:
            s = r & u.select_atoms("resid %s" % e)
            this_sel = "%s and resid %s" % (selection, e)
            group_selections[this_sel] = s
    else:
        assert False, ("Invalid argument for igroupby. Only name and id are valid options.")


def ramachandran_plot(atomgroup, selection, outputfile1, outputfile2, image_format='png'):
    # plot standard mdanalysis and seaborn 2D with kde
    R = Ramachandran(atomgroup).run()
    fig, ax = plt.subplots(figsize=plt.figaspect(1))
    R.plot(ax=ax, color='k', marker='.', ref=True)

    a = R.angles.reshape(np.prod(R.angles.shape[:2]), 2)
    # open hdf file
    with h5py.File(args.o_data1, 'a') as f:
        setname = "%s" % (selection)
        f["/" + setname + "/ramachandran/phi"] = a[:, 0]
        f["/" + setname + "/ramachandran/psi"] = a[:, 1]
    plt.tight_layout()
    # svg is better but sticking with png for now
    plt.savefig(outputfile1, format=image_format)

    sns.reset_defaults()
    importlib.reload(plt)
    importlib.reload(sns)
    with sns.axes_style("white"):
        h = sns.jointplot(x=a[:, 0], y=a[:, 1],
                          kind="kde", space=0)
        h.set_axis_labels(r'$\phi$ (deg)', r'$\psi$ (deg)')
        h.ax_joint.set_xlim(-180, 180)
        h.ax_joint.set_ylim(-180, 180)
        h.ax_joint.xaxis.set_major_locator(ticker.MultipleLocator(60))
        h.ax_joint.yaxis.set_major_locator(ticker.MultipleLocator(60))
        # plt.tight_layout()  #<-- using this breaks the plot, but not using it causes a cut, so use bbox_inches
        plt.savefig(outputfile2, format=image_format, bbox_inches='tight')


def get_base64_encoded_image(image_path):
    """  encode image to string for use in html later"""
    with open(image_path, "rb") as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')


plots = []
if args.igroupby is not None:
    for k, v in group_selections.items():
        print(k, v)
        try:
            ramachandran_plot(v, str(k), "ramachandran1" +
                              str(k), "ramachandran2" + str(k))
            plots.append({'Name': "%s" % (k), 'plot1': get_base64_encoded_image(
                "ramachandran1" + str(k)), 'plot2': get_base64_encoded_image("ramachandran2" + str(k))})
        except Exception as einstance:
            print(type(einstance))
            print(einstance.args)
            print(einstance)

ramachandran_plot(r, selection, args.o_plot1, args.o_plot2)
plots.insert(0, {'Name': selection, 'plot1': get_base64_encoded_image(
    args.o_plot1), 'plot2': get_base64_encoded_image(args.o_plot2)})

template_environment = Environment(loader=FileSystemLoader(
    currentpath), lstrip_blocks=True, trim_blocks=True)
template = template_environment.get_template(
    'ramachandran_auto_protein_html.j2')
with open(args.o_html1, 'w+') as f:
    f.write(template.render(title="Ramachandran Plots", plots=plots))
