import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, cophenet
from scipy.spatial.distance import pdist

def json_to_np(fname, start=None, end=None):
    """
    Load json file and convert to numpy array
    """
    with open(fname) as f: 
        k = json.load(f)
    print(np.array(k)[:,:,start:end].shape)                                                                                                                                                                            
    return np.array(k)[:,:,start:end]

def flatten_tensor(tensor, normalize=True):
    """
    Flatten tensor to a 2D matrix along the time axis
    """
    av = np.mean(tensor, axis=(0,1)) if normalize else 1
    return np.mean(tensor/av, axis=2)

def get_cluster_linkage_array(mat, clustering_method='average'):
    Z = linkage(mat, clustering_method)
    c, coph_dists = cophenet(Z, pdist(mat))
    print('Cophenetic correlation coefficient: {}'.format(c))
    return Z

def plot_dist_mat(mat, output, cmap='plasma'):
    """
    Plot distance matrix as heatmap
    """
    fig, ax = plt.subplots(1)
    p = ax.pcolormesh(mat, cmap=cmap)
    plt.xlabel('Trajectory number')
    plt.ylabel('Trajectory number')
    plt.colorbar(p)
    plt.draw()
    plt.savefig(output, format='png')

def plot_dendrogram(Z, output):
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Trajectory index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.savefig(output, format='png')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', help='JSON input file (for 3D tensor).')
    parser.add_argument('--mat', help='Input tabular file (for 2D matrix).')
    parser.add_argument('--outp-mat', help='Tabular output file.')
    parser.add_argument('--Z', required=True, help='File for cluster linkage array.')
    parser.add_argument('--dendrogram', help="Path to the output dendrogram file")
    parser.add_argument('--heatmap', help="Path to the output distance matrix file")
    parser.add_argument('--clustering-method', default='average', choices=['single', 'complete', 'average', 'centroid', 'median', 'ward', 'weighted'], help="Method to use for clustering.")
    parser.add_argument('--cmap', type=str, default='plasma', help="Matplotlib colormap to use for plotting distance matrix.")
    parser.add_argument('--start', type=int, help="First trajectory frame to calculate distance matrix")
    parser.add_argument('--end', type=int, help="Last trajectory frame to calculate distance matrix")
    parser.add_argument('--normalize', action="store_true", help="Normalize the RMSD variation over the trajectories before averaging.")
    args = parser.parse_args()

    print(args)
    if args.json:
        tensor = json_to_np(args.json, args.start, args.end)
        mat = flatten_tensor(tensor, args.normalize)
        np.savetxt(args.outp_mat, mat)
    elif args.mat:
        mat = np.loadtxt(args.mat)
    else:
        print("Either --json or --mat must be specified.")
        exit(1)

    Z = get_cluster_linkage_array(mat, args.clustering_method)
    np.savetxt(args.Z, Z)

    if args.heatmap:
        plot_dist_mat(mat, args.heatmap, args.cmap)

    if args.dendrogram:
        plot_dendrogram(Z, args.dendrogram)

if __name__ == "__main__":
    main()