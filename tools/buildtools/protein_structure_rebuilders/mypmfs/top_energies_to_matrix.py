import argparse
import os

from MyPMFs_py.module import batch_download, training


def __main__():

    parser = argparse.ArgumentParser(description="Energy Matrix Generator")

    parser.add_argument(
        '--protein_file', default=None, help="Protein File Input")
    parser.add_argument(
        '--txtmatrix', default=None, help="MJ Style Energy Matrix")
    parser.add_argument(
        '--csvmatrix', default=None, help="CSV Energy Matrix")
    args = parser.parse_args()

    batch_download(args.protein_file, "pdb_files")

    directory = "./pdb_files"

    # create tsv input for MyPMFs training (excluding failed PDB downloads)
    new_file = ""
    for i in range(len(os.listdir(directory))):
        new_file += os.fsdecode(os.listdir(directory)[i])[:4] + "\n"

    new_file = new_file[:len(new_file)-1]

    with open("training_dataset.txt", "w") as file:
        file.write(new_file)

    training("training_dataset.txt", "pdb_files")

    aa = ["C", "M", "F", "I", "L", "V", "W", "Y", "A", "G",
          "T", "S", "N", "Q", "D", "E", "H", "R", "K", "P"]
    C = {}
    M = {}
    F = {}
    Iaa = {}
    L = {}
    V = {}
    W = {}
    Y = {}
    A = {}
    G = {}
    T = {}
    S = {}
    N = {}
    Q = {}
    D = {}
    E = {}
    H = {}
    R = {}
    K = {}
    P = {}
    dict = {"C": C, "M": M, "F": F, "I": Iaa, "L": L,
            "V": V, "W": W, "Y": Y, "A": A, "G": G,
            "T": T, "S": S, "N": N, "Q": Q, "D": D,
            "E": E, "H": H, "R": R, "K": K, "P": P}

    with open("data_potentials/top_energies.tsv") as f:
        for line in f:
            if line:
                line_data = line.strip().split('\t')
                (dict[line_data[0][0]])[line_data[0][3]] = line_data[1]
                (dict[line_data[0][3]])[line_data[0][0]] = line_data[1]

    # loop creates energy matrix in txt format
    txtmatrix = "    C       M       F       I       "
    txtmatrix += "L       V       W       Y       "
    txtmatrix += "A       G       T       S       "
    txtmatrix += "N       Q       D       E       H       R       K       P\n"
    for i in range(len(aa)):
        for k in range(i):
            txtmatrix = txtmatrix + "0e+00 "
        for j in range(i, len(aa)):
            txtmatrix = txtmatrix + str(dict[aa[i]][aa[j]]) + "e+00 "
        txtmatrix = txtmatrix + '\n'

    with open(args.txtmatrix, "w") as text_file:
        text_file.write(txtmatrix)

    # loop creates energy matrix in csv format
    csvmatrix = "C,M,F,I,L,V,W,Y,A,G,T,S,N,Q,D,E,H,R,K,P\n"
    for i in range(len(aa)):
        for j in range(len(aa)):
            csvmatrix = csvmatrix + str(dict[aa[i]][aa[j]])
            if j != len(aa) - 1:
                csvmatrix += ","
        csvmatrix = csvmatrix + '\n'

    with open(args.csvmatrix, "w") as file:
        file.write(csvmatrix)


if __name__ == "__main__":
    __main__()
