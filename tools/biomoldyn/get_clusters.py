import argparse
import collections
import json

import numpy as np

from scipy.cluster.hierarchy import fcluster


def separate_clusters(Z_fpath, threshold, min_members, output):
    Z = np.loadtxt(Z_fpath)
    branch_assignments = fcluster(Z, threshold, criterion='distance')
    cluster_dict = collections.defaultdict(list)
    for n, val in enumerate(branch_assignments):
        cluster_dict[branch_assignments[n]].append(n)
    cluster_dict = {int(k): v for k, v in cluster_dict.items()
                    if len(v) >= min_members}
    with open(output, 'w') as f:
        json.dump(cluster_dict, f, indent=4, sort_keys=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--Z', required=True,
                        help='File for cluster linkage array.')
    parser.add_argument('--threshold', type=float, required=True,
                        help='Distance cutoff.')
    parser.add_argument('--min-members', type=int, required=True,
                        help='Minimum number of members of the cluster.')
    parser.add_argument('--output', required=True,
                        help='Output file.')
    args = parser.parse_args()

    separate_clusters(args.Z, args.threshold,
                      args.min_members, args.output)


if __name__ == "__main__":
    main()
