#!/usr/bin/env python3
"""Create a compact two-chain RMSX regression fixture from one trajectory."""

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
    parser.add_argument("--frames", type=int, default=36)
    parser.add_argument("--chain-offset", type=float, default=45.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.frames < 3:
        raise SystemExit("--frames must be at least 3")

    source = mda.Universe(args.topology, args.trajectory)
    frame_count = min(args.frames, len(source.trajectory))
    frame_indices = np.unique(np.linspace(0, len(source.trajectory) - 1, frame_count, dtype=int))

    combined = mda.Merge(source.atoms, source.atoms)
    if combined.atoms.n_segments != 2:
        raise RuntimeError(f"Expected two merged segments, found {combined.atoms.n_segments}.")
    combined.segments.segids = ["A", "B"]
    combined.atoms.chainIDs = np.concatenate(
        [
            np.full(source.atoms.n_atoms, "A", dtype=object),
            np.full(source.atoms.n_atoms, "B", dtype=object),
        ]
    )

    output_topology = Path(args.output_topology)
    output_trajectory = Path(args.output_trajectory)
    output_topology.parent.mkdir(parents=True, exist_ok=True)
    output_trajectory.parent.mkdir(parents=True, exist_ok=True)
    offset = np.array([args.chain_offset, 0.0, 0.0], dtype=np.float32)

    def set_combined_positions(frame_index: int) -> None:
        source.trajectory[frame_index]
        combined.trajectory.ts.time = source.trajectory.ts.time
        combined.trajectory.ts.frame = source.trajectory.ts.frame
        positions = source.atoms.positions.copy()
        combined.atoms.positions = np.concatenate([positions, positions + offset])
        if source.dimensions is not None:
            combined.dimensions = source.dimensions

    set_combined_positions(int(frame_indices[0]))
    combined.atoms.write(str(output_topology))

    writer_kwargs = {"n_atoms": combined.atoms.n_atoms}
    if output_trajectory.suffix.lower() == ".xtc":
        writer_kwargs["precision"] = 3
    with mda.Writer(str(output_trajectory), **writer_kwargs) as writer:
        for frame_index in frame_indices:
            set_combined_positions(int(frame_index))
            writer.write(combined.atoms)

    print(
        f"Wrote two-chain fixture with {combined.atoms.n_atoms} atoms, "
        f"{len(frame_indices)} frames, and segment IDs A/B."
    )


if __name__ == "__main__":
    main()
