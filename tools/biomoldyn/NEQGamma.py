#!/usr/bin/env python

# coding: utf-8
# This script is a modified version of a script written by Steffen Wolf under the GPL v3.0.
# The original version can be accessed at https://github.com/moldyn/dcTMD/blob/master/NEQGamma.py
# # $\Gamma_{NEQ}(x)$ calculation routine
# * target: calculate $\Gamma_{NEQ}(t)=\frac{1}{k_BT}\int_0^tdt'\left<\delta f_c(t)\delta f_c(t')\right>$

import argparse
import json
import sys

import numpy as np

import pandas as pd

import scipy
import scipy.integrate
from scipy.ndimage.filters import gaussian_filter


def get_file_names(list_file):
    with open(list_file) as f:
        return [line for line in f.read().split('\n') if line]


def NEQGamma(file_names, output, output_frict, vel, T, av, sigma):
    N = len(file_names)
    RT = 0.0083144598 * T

    sys.stdout.write("reading data...\n")

    # read in initial data to get length of necessary array
    test_file = pd.read_csv(file_names[0], sep='\s+', header=None, skiprows=17, dtype=float)
    length_data = len(test_file[0].values)
    full_force_set = np.zeros((N, length_data))
    x = np.zeros(length_data)
    t = np.zeros(length_data)
    t = test_file[0].values
    x = test_file[0].values * vel

    # read in data
    for i in range(0, N):
        current_file_name = file_names[i]
        sys.stdout.write("reading file {}\n".format(current_file_name))
        input_file_data = pd.read_csv(current_file_name, sep='\s+', header=None, skiprows=17, dtype=float)
        full_force_set[i, :] = input_file_data[1].values

    # preprocessing
    # * force average: calculate $\left< f_c (t)\right>_N$. **Important:** this is an ensemble average over the trajectory ensemble $N$, not the time average over $t$
    av_force = np.zeros(length_data)
    av_forceintegral = np.zeros(length_data)
    for i in range(length_data):
        av_force[i] = np.mean(full_force_set[:, i])
    av_forceintegral[1:] = scipy.integrate.cumtrapz(av_force, x)

    # * calculate $\delta f_c(t) = f_c(t) - \left< f_c (t) \right>_N$ for all $t$
    sys.stdout.write("calculating fluctuations...\n")
    delta_force_set = np.zeros((N, length_data))
    for i in range(length_data):
        delta_force_set[:, i] = full_force_set[:, i] - av_force[i]

    # evaluation
    # * optimized algorithm for numerical evaluation:
    #     * integrate: $\int_0^t dt' \delta f_c(t')$ for all $t'$
    #     * multiply by $\delta f_c(t)$ to yield $\int_0^t dt'\delta f_c(t) \delta f_c(t')$ for $t$ with all $t' \leq t$ each
    #     * then calculate the ensemble average $\left< \int_0^t dt' \delta f_c(t) \delta f_c(t') \right>$
    int_delta_force_set = np.zeros((N, length_data))
    for n in range(N):
        int_delta_force_set[n, 1:] = scipy.integrate.cumtrapz(delta_force_set[n, :], t)

    sys.stdout.write("averaging and integrating...\n")
    intcorr = np.zeros((N, length_data))

    for n in range(N):
        for i in range(length_data):
            intcorr[n, i] = delta_force_set[n, i] * int_delta_force_set[n, i]
            if i % 1000 == 0:
                sys.stdout.write("Trajectory {:2d} {:3.1f} % done\r".format(n + 1, (i / length_data) * 100))

    # * shape of  $\int_0^t dt' \delta f_c(t) \delta f_c(t')$:
    sys.stdout.write("final average...\n")
    av_intcorr = np.zeros(length_data)
    for i in range(length_data):
        av_intcorr[i] = np.mean(intcorr[:, i]) / RT

    # * autocorrelation function evaluation:
    #     * calculate $\left< \delta f_c(t) \delta f_c(t') \right>$ for the last $t$

    corr_set = np.zeros((N, length_data))
    autocorr_set = np.zeros(length_data)

    sys.stdout.write("calculating and processing ACF...\n")
    for n in range(N):
        for i in range(length_data):
            corr_set[n, i] = delta_force_set[n, i] * delta_force_set[n, length_data - 1]

    for i in range(length_data):
        autocorr_set[i] = np.mean(corr_set[:, i])

    # * Gauss filter:
    sys.stdout.write("applying Gauss filter...\n")
    blurr = sigma
    blurred = np.zeros(length_data)
    blurred = gaussian_filter(av_intcorr, sigma=blurr)

    # * sliding window average:
    sys.stdout.write("applying sliding window average...\n")
    window = av
    runn_av = np.zeros(length_data)
    runn_av = np.convolve(av_intcorr, np.ones((window,)) / window, mode='same')

    #  runn_av[:-window+1] = np.convolve(av_intcorr, np.ones((window,))/window, mode='valid')
    # * $W_{diss}$ from integration:
    wdiss = np.zeros(length_data)
    wdiss[1:] = scipy.integrate.cumtrapz(av_intcorr, x) * vel

    sys.stdout.write("writing output...\n")
    distresult = open(output, "w")
    frictresult = open(output_frict, "w")

    distresult.write("#x  force_integral  frict_coeff   wdiss   corrected_force_integral\n")
    for i in range(length_data):
        distresult.write("{:15.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f}\n".format(x[i], av_forceintegral[i], av_intcorr[i], wdiss[i], av_forceintegral[i] - wdiss[i]))

    frictresult.write("#x   ACF   frict_coeff   gauss_filtered_frict_coeff   av_window_frict_coeff\n")
    for i in range(length_data):
        frictresult.write("{:15.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f}\n".format(x[i], autocorr_set[i], av_intcorr[i], blurred[i], runn_av[i]))

    distresult.close()
    frictresult.close()

    sys.stdout.write("Done!\n")


def main():
    parser = argparse.ArgumentParser(description='dcTMD friciton correction (please cite: Wolf, S., Stock, G. Targeted Molecular Dynamics Calculations of Free Energy Profiles Using a Nonequilibrium Friction Correction. J. Chem. Theory Comput. 2018, 14(12), 6175-6182, DOI: 10.1021/acs.jctc.8b00835). Integrates a constraint force file via trapezoid rule, calculates the NEQ memory friction kernel and friction factors, and performs a friction correction. First column: reaction coordinate in nm calculated via t * vel. Second column: force integral, i.e. the work profile. Third column: friction factors. Fourth column: trapezoid integral (final value) of friction work along reaction coordinate. Fourth column: friction corrected work profile. ATTENTION: Use with python3 or higher!')
    parser.add_argument('-i', metavar='<xvg force file>', type=str, help='List of xvg constraint force files prefix as given by Gromacs mdrun -pf option before running number.')
    parser.add_argument('-o', metavar='<combined results>', type=str, help='file to write x, dG(x), friction coefficeint by integration (time resolved), and the friction-corrected dG(x).')
    parser.add_argument('-ofrict', metavar='<combined friction results>', type=str, help='file to write x, ACF, friction coefficeint by integration (time resolved), gauss filtered friction coefficient, and slide window averaged friction.')
    parser.add_argument('-vel', metavar='<pull velocity>', type=float, help='pull velocity in nm/ps for converting simulation time t into distance x')
    parser.add_argument('-T', metavar='<temperature>', type=float, help='temperature in K')
    parser.add_argument('-av', metavar='<average window>', type=int, help='size of averaging window for displaying Gamma(x) (recommended: 4 to 20 per 100 data points)')
    parser.add_argument('-sigma', metavar='<gauss blurr>', type=int, help='sigma value for Gauss filter for displaying Gamma(x) (recommended: 4 per 100 data points)')
    parser.add_argument('-json', metavar='<json>', type=str, help='JSON file defining cluster membership')

    args = parser.parse_args()

    file_names = get_file_names(args.i)
    if args.json:
        with open(args.json) as f:
            j = json.load(f)
        file_names_dct = {n: [file_names[int(m)] for m in j[n]] for n in j}

        for cluster in file_names_dct:
            NEQGamma(file_names_dct[cluster], 'cluster{}_{}'.format(cluster, args.o), 'cluster{}_{}'.format(cluster, args.ofrict), args.vel, args.T, args.av, args.sigma)
    else:
        NEQGamma(file_names, args.o, args.ofrict, args.vel, args.T, args.av, args.sigma)


if __name__ == "__main__":
    main()
