import argparse
import json

import MDAnalysis as m
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

import numpy as np


def calc_rmsd(str_files, traj_files, ref_str, str_format, traj_format,
              ref_str_format, filepath_out, group, start, end, step,
              fitting_atoms):
    """
    this function does everything
    the function will cycle through range 0 to no_t and load all files found.

    str_files: text file with filepaths for structures, one on each line
    traj_files: text file with filepaths for trajectories, one on each line
    ref_str: reference structure for fitting
    filepath_in: directory where the files are located
    filepath_out: pickle file where results (3D matrix) should be saved to

    group: atoms for which RMSD should be calculated;
        use the MDAnalysis selection language
    fitting_atoms: atoms used for str alignment prior to RMSD calculation;
        use the MDAnalysis selection language

    start: first trajectory frame to calculate RMSD
    end: last trajectory frame to calculate RMSD
    step: how frequently frames are sampled between start and end; obviously,
        the larger the step, the quicker the script finishes
    """

    # open list of files
    with open(str_files) as f:
        str_file_list = f.read().strip().split('\n')

    with open(traj_files) as f:
        traj_file_list = f.read().strip().split('\n')

    if len(str_file_list) != len(traj_file_list):
        raise IOError('Number of structure and trajectory files not equal.')

    no_t = len(traj_file_list)

    data = np.zeros((no_t, no_t,
                    int((end - start)/step + ((end - start) % step > 0))))

    # load files
    universes = {}

    for traj in range(no_t):
        mobile = m.Universe(str_file_list[traj], traj_file_list[traj],
                            format=traj_format, topology_format=str_format)
        ref = m.Universe(ref_str, topology_format=ref_str_format)

        mobile.trajectory[-1]  # set mobile trajectory to last frame
        ref.trajectory[0]  # set reference trajectory to first frame

        # perform alignment
        align.AlignTraj(mobile, ref, select=fitting_atoms,
                        in_memory=True).run()

        grp = mobile.select_atoms(group)
        universes[traj] = m.core.universe.Merge(grp)  # create Universe w grp
        coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                           grp).run().results  # write to uv
        universes[traj].load_new(coordinates, format=MemoryReader)

    print("All trajs loaded by MDAnalysis")

    # calculate differences
    for traj1 in range(no_t):
        print("Calculating differences for traj {}".format(traj1))
        for traj2 in range(traj1):

            u1 = universes[traj1]
            u2 = universes[traj2]

            l1 = u1.select_atoms(group)
            l2 = u2.select_atoms(group)

            rmsd = rms.RMSD(l1, l2)

            rmsd.run()

            data[traj1, traj2] = rmsd.rmsd[:, 2]
            data[traj2, traj1] = rmsd.rmsd[:, 2]

    with open(filepath_out, 'w') as f:
        json.dump(data.tolist(), f)

    print("Done!")
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trajs', required=True,
                        help='File containing trajectory filepaths.')
    parser.add_argument("--strs",
                        help='File containing structure filepaths.')
    parser.add_argument("--ref-str",
                        help='File containing reference structure.')
    parser.add_argument('--traj-format', required=True,
                        help='Trajectory format.')
    parser.add_argument("--str-format", help='Structure format.')
    parser.add_argument("--ref-str-format",
                        help='Reference structure format.')
    parser.add_argument('-o', '--outfile',
                        help="Path to the output JSON file")
    parser.add_argument('--group', help="Atoms for which RMSD should be"
                        "calculated in MDAnalysis selection language")
    parser.add_argument('--fitting', help="Fitting atoms for alignment"
                        "prior to RMSD calculation")
    parser.add_argument('--start', type=int,
                        help="First trajectory frame to calculate RMSD")
    parser.add_argument('--end', type=int,
                        help="Last trajectory frame to calculate RMSD")
    parser.add_argument('--step', type=int,
                        help="Frame sampling frequency for RMSD calculation")
    args = parser.parse_args()

    calc_rmsd(args.strs, args.trajs, args.ref_str, args.str_format,
              args.traj_format, args.ref_str_format, args.outfile,
              args.group, args.start, args.end, args.step, args.fitting)


if __name__ == "__main__":
    main()
