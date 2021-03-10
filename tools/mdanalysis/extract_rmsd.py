import argparse
import json

import MDAnalysis as m
from MDAnalysis.analysis import rms

import numpy as np


def calc_rmsd(str_files, traj_files, str_format, traj_format, filepath_out,
              group, start, end, step):
    """
    the function will cycle through range 0 to no_t and load all files found.

    str_files: text file with filepaths for structures, one on each line
    traj_files: text file with filepaths for trajectories, one on each line
    filepath_in: directory where the files are located
    filepath_out: pickle file where results (3D matrix) should be saved to

    group: atoms for which RMSD should be calculated;
        use the MDAnalysis selection language

    start: first trajectory frame to calculate RMSD
    end: last trajectory frame to calculate RMSD
    step: how frequently frames are sampled between start and end; obviously,
        the larger the step, the quicker the script finishes
    """

    # open list of files
    with open(str_files) as f1, open(traj_files) as f2:
        str_file_list = f1.read().strip().split('\n')
        traj_file_list = f2.read().strip().split('\n')

        if sum(1 for line in f1) != sum(1 for line in f2):
            raise IOError('Number of structure and trajectory files unequal.')

    no_t = len(traj_file_list)

    # hard to find array size before loading files
    universe_coordinate_data = []

    for traj in range(no_t):
        # We no longer align here, users should do this themselves.
        u = m.Universe(str_file_list[traj], traj_file_list[traj],
                       format=traj_format, topology_format=str_format)
        u.transfer_to_memory()
        grp = u.select_atoms(group).universe
        coordinates = grp.trajectory.coordinate_array[start:end:step]
        universe_coordinate_data.append(coordinates)

    universe_coordinate_data = np.array(universe_coordinate_data)
    print("All trajs loaded by MDAnalysis")
    data = np.zeros((no_t, no_t, universe_coordinate_data.shape[1]))

    # calculate differences
    for traj1 in range(no_t):
        print("Calculating differences for traj {}".format(traj1))
        for traj2 in range(traj1):
            for frame in range(data.shape[2]):
                A = universe_coordinate_data[traj1, frame]
                B = universe_coordinate_data[traj2, frame]
                r = rms.rmsd(A, B)
                data[traj1, traj2, frame] = r
                data[traj2, traj1, frame] = r

    with open(filepath_out, 'w') as f:
        json.dump(data.tolist(), f, indent=4, sort_keys=True)

    print("Done!")
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trajs', required=True,
                        help='File containing trajectory filepaths.')
    parser.add_argument("--strs",
                        help='File containing structure filepaths.')
    parser.add_argument('--traj-format', required=True,
                        help='Trajectory format.')
    parser.add_argument("--str-format", help='Structure format.')
    parser.add_argument('-o', '--outfile',
                        help="Path to the output JSON file")
    parser.add_argument('--group', help="Atoms for which RMSD should be"
                        "calculated in MDAnalysis selection language")
    parser.add_argument('--start', type=int,
                        help="First trajectory frame to calculate RMSD")
    parser.add_argument('--end', type=int,
                        help="Last trajectory frame to calculate RMSD")
    parser.add_argument('--step', type=int,
                        help="Frame sampling frequency for RMSD calculation")
    args = parser.parse_args()
    calc_rmsd(args.strs, args.trajs, args.str_format,
              args.traj_format, args.outfile,
              args.group, args.start, args.end, args.step)


if __name__ == "__main__":
    main()
