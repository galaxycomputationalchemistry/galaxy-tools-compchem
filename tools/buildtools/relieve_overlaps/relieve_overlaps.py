#!/usr/bin/env python
"""
Push apart near-overlapping atoms in a GROMACS structure before energy minimization.

Packed systems (membranes built with packmol-memgen, solvent boxes, inserted molecules)
sometimes contain atoms of different residues that (almost) coincide. The forces between them
are not finite in single precision, so minimization fails within a few steps. Every pair of
atoms from different residues closer than --min-distance is pushed apart symmetrically along
the line joining them to that distance; this is repeated until no such pair is left.
Minimization then removes the remaining strain.

Atoms of the same residue are never moved against each other, and the largest allowed
distance (1.0 Angstrom) is below the shortest bond between residues (about 1.33 Angstrom), so
bonded neighbours are not affected. Only the coordinates of moved atoms change; every other
character of the file is kept.
"""
import argparse
import collections
import math
import sys

CELL = 2.0  # Angstrom, grid cell for the neighbour search


def read_gro(path):
    with open(path) as fh:
        lines = fh.readlines()
    try:
        n = int(lines[1])
        body = lines[2:2 + n]
        box = [float(v) * 10 for v in lines[2 + n].split()[:3]]
        xyz = [[float(ln[20:28]) * 10, float(ln[28:36]) * 10, float(ln[36:44]) * 10] for ln in body]
    except (IndexError, ValueError) as exc:
        sys.exit(f"ERROR: {path} is not a valid GRO file ({exc}).")
    if len(body) != n:
        sys.exit(f"ERROR: {path} declares {n} atoms but contains {len(body)}.")
    return lines, n, xyz, box


def close_pairs(xyz, residue, box, limit):
    """Yield (a, b, d, r) for atoms of different residues closer than limit (minimum image)."""
    dims = [max(int(b // CELL), 1) for b in box]
    grid = collections.defaultdict(list)
    for k, x in enumerate(xyz):
        grid[tuple(int((x[i] % box[i]) // CELL) % dims[i] for i in range(3))].append(k)
    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    seen = set()
    for key, members in grid.items():
        for o in offsets:
            other = grid.get(tuple((key[i] + o[i]) % dims[i] for i in range(3)), ())
            for a in members:
                for b in other:
                    if b <= a or residue[a] == residue[b] or (a, b) in seen:
                        continue
                    d = [xyz[a][i] - xyz[b][i] for i in range(3)]
                    d = [d[i] - box[i] * round(d[i] / box[i]) for i in range(3)]
                    r = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                    if r < limit:
                        seen.add((a, b))
                        yield a, b, d, r


def describe(body, k):
    ln = body[k]
    return f"{ln[0:5].strip()}{ln[5:10].strip()}:{ln[10:15].strip()}({k + 1})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="GRO structure")
    ap.add_argument("--output", required=True, help="GRO structure with overlaps relieved")
    ap.add_argument("--report", required=True, help="text report")
    ap.add_argument("--min-distance", type=float, default=0.5,
                    help="push apart atoms of different residues closer than this (Angstrom, 0 < d <= 1.0)")
    ap.add_argument("--max-passes", type=int, default=20)
    args = ap.parse_args()
    if not 0 < args.min_distance <= 1.0:
        sys.exit("ERROR: --min-distance must be > 0 and <= 1.0 Angstrom (bonds between residues are longer).")

    lines, n, xyz, box = read_gro(args.input)
    body = lines[2:2 + n]
    residue = [ln[0:10] for ln in body]  # residue number + name
    start = [list(x) for x in xyz]
    initial = sorted(close_pairs(xyz, residue, box, args.min_distance - 1e-3), key=lambda p: p[3])

    pushed, moved, passes = 0, set(), 0
    for passes in range(1, args.max_passes + 1):
        pairs = list(close_pairs(xyz, residue, box, args.min_distance - 1e-3))
        if not pairs:
            passes -= 1
            break
        for a, b, d, r in pairs:
            u = [v / r for v in d] if r > 1e-6 else [1.0, 0.0, 0.0]
            shift = (args.min_distance - r) / 2 + 1e-3
            for i in range(3):
                xyz[a][i] += u[i] * shift
                xyz[b][i] -= u[i] * shift
            moved.update((a, b))
        pushed += len(pairs)
    else:
        if list(close_pairs(xyz, residue, box, args.min_distance - 1e-3)):
            sys.exit(f"ERROR: atoms closer than {args.min_distance} A remain after {args.max_passes} passes.")

    with open(args.output, "w") as fh:
        fh.writelines(lines[:2])
        for k, ln in enumerate(body):
            if k in moved:
                ln = f"{ln[:20]}{xyz[k][0] / 10:8.3f}{xyz[k][1] / 10:8.3f}{xyz[k][2] / 10:8.3f}{ln[44:]}"
            fh.write(ln)
        fh.writelines(lines[2 + n:])

    closest = min((r for *_, r in close_pairs(xyz, residue, box, 1.0)), default=None)
    largest = max((math.dist(xyz[k], start[k]) for k in moved), default=0.0)
    with open(args.report, "w") as fh:
        fh.write(f"Atoms: {n}\nMinimum distance between atoms of different residues: {args.min_distance} A\n\n")
        if initial:
            fh.write(f"{len(initial)} pairs closer than {args.min_distance} A in the input:\n")
            for a, b, _, r in initial[:50]:
                fh.write(f"  {r:6.3f} A  {describe(body, a)} -- {describe(body, b)}\n")
            if len(initial) > 50:
                fh.write(f"  ... {len(initial) - 50} more\n")
            fh.write(f"\n{pushed} pair moves in {passes} pass(es), {len(moved)} atoms moved, "
                     f"largest displacement {largest:.2f} A\n")
        else:
            fh.write("No overlapping atoms found; structure unchanged.\n")
        fh.write("Closest pair of atoms of different residues now: "
                 + (f"{closest:.2f} A\n" if closest is not None else "> 1.00 A\n"))
    sys.stdout.write(open(args.report).read())


if __name__ == "__main__":
    main()
