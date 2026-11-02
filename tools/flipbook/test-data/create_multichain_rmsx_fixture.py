#!/usr/bin/env python3
"""Create a compact two-chain RMSX regression fixture from the real two-chain protease trajectory."""

from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--output-topology", required=True)
    parser.add_argument("--output-trajectory", required=True)
    parser.add_argument("--frames", type=int, default=180)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.frames < 3:
        raise SystemExit("--frames must be at least 3")

    source = mda.Universe(args.topology, args.trajectory)
    frame_count = min(args.frames, len(source.trajectory))
    frame_indices = np.unique(np.linspace(0, len(source.trajectory) - 1, frame_count, dtype=int))

    if set(source.segments.segids) != {"A", "B"}:
        raise RuntimeError("Expected the source protease segments A and B.")
    if any(len(source.select_atoms(f"segid {chain} and name CA")) != 99 for chain in ["A", "B"]):
        raise RuntimeError("Expected 99 C-alpha residues in each protease chain.")

    output_topology = Path(args.output_topology)
    output_trajectory = Path(args.output_trajectory)
    output_topology.parent.mkdir(parents=True, exist_ok=True)
    output_trajectory.parent.mkdir(parents=True, exist_ok=True)
    source.trajectory[int(frame_indices[0])]
    source.atoms.write(str(output_topology))

    writer_kwargs = {"n_atoms": source.atoms.n_atoms}
    if output_trajectory.suffix.lower() == ".xtc":
        writer_kwargs["precision"] = 3
    with mda.Writer(str(output_trajectory), **writer_kwargs) as writer:
        for frame_index in frame_indices:
            source.trajectory[int(frame_index)]
            writer.write(source.atoms)

    print(
        f"Wrote real protease fixture with {source.atoms.n_atoms} atoms, "
        f"{len(frame_indices)} frames, and segment IDs A/B."
    )


if __name__ == "__main__":
    main()
