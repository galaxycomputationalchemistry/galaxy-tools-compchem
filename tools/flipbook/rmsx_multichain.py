#!/usr/bin/env python3
"""Run RMSX for one or more chains and stage Galaxy collection outputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TypedDict


CHAIN_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


class ChainOutput(TypedDict):
    chain: str
    designation: str
    directory: Path
    rmsx: Path
    rmsd: Path
    rmsf: Path
    mask: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--galaxy-output-dir", required=True)
    parser.add_argument("--chain-mode", choices=("all", "selected"), required=True)
    parser.add_argument("--chains", default="", help="Comma-separated segment or chain IDs.")
    parser.add_argument("--num-slices", type=int, required=True)
    parser.add_argument("--start-frame", type=int, required=True)
    parser.add_argument("--end-frame", type=int)
    parser.add_argument("--analysis-type", choices=("protein", "dna", "rna", "generic"), required=True)
    parser.add_argument("--summary-n", type=int, required=True)
    parser.add_argument("--palette", required=True)
    parser.add_argument("--interpolate", action="store_true")
    parser.add_argument("--plot-helper", required=True)
    parser.add_argument("--rscript", default="Rscript")
    return parser.parse_args()


def unique_values(values) -> list[str]:
    seen = set()
    result = []
    for value in values:
        normalized = str(value).strip()
        if normalized and normalized not in seen:
            seen.add(normalized)
            result.append(normalized)
    return result


def parse_selected_chains(selector: str) -> list[str]:
    chains = unique_values(selector.split(","))
    if not chains:
        raise ValueError("Selected-chain mode requires at least one chain or segment ID.")
    invalid = [chain for chain in chains if not CHAIN_ID_RE.fullmatch(chain)]
    if invalid:
        raise ValueError(f"Invalid chain or segment ID(s): {', '.join(invalid)}")
    return chains


def chain_id_to_segids(universe, chain_id: str) -> list[str]:
    try:
        chain_ids = universe.atoms.chainIDs
        segids = universe.atoms.segids
    except Exception:
        return []
    return unique_values(
        segid
        for observed_chain, segid in zip(chain_ids, segids)
        if str(observed_chain).strip() == chain_id
    )


def resolve_chains(universe, mode: str, selector: str, analysis_type: str, selection_builder) -> list[str]:
    available_segids = unique_values(universe.atoms.segids)
    if not available_segids:
        raise ValueError("The topology does not contain any usable segment IDs.")

    if mode == "all":
        requested = available_segids
    else:
        requested = []
        missing = []
        for chain in parse_selected_chains(selector):
            if chain in available_segids:
                requested.append(chain)
                continue
            mapped = chain_id_to_segids(universe, chain)
            if len(mapped) == 1:
                requested.append(mapped[0])
            else:
                missing.append(chain)
        if missing:
            raise ValueError(
                f"Chain or segment ID(s) not found or ambiguous: {', '.join(missing)}. "
                f"Available segment IDs: {', '.join(available_segids)}"
            )

    valid = []
    skipped = []
    for chain in unique_values(requested):
        selection = selection_builder(analysis_type=analysis_type, chain_sele=chain)
        atoms = universe.select_atoms(selection).select_atoms("name CA")
        if len(atoms):
            valid.append(chain)
        else:
            skipped.append(chain)

    if mode == "selected" and skipped:
        raise ValueError(f"Selected segment(s) contain no analyzable CA atoms: {', '.join(skipped)}")
    if not valid:
        raise ValueError("No chains with analyzable CA atoms were found in the topology.")
    if skipped:
        print(f"Skipped segments without analyzable CA atoms: {', '.join(skipped)}")
    return valid


def one_matching_file(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {label} in {directory}, found {len(matches)}.")
    return matches[0]


def write_false_mask(rmsx_csv: Path, output: Path) -> None:
    with rmsx_csv.open(newline="", encoding="utf-8") as source, output.open(
        "w", newline="", encoding="utf-8"
    ) as destination:
        rows = csv.DictReader(source)
        writer = csv.DictWriter(
            destination,
            fieldnames=("ResidueID", "ChainID", "Masked"),
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "ResidueID": row.get("ResidueID", ""),
                    "ChainID": row.get("ChainID", ""),
                    "Masked": "False",
                }
            )


def merge_csv_files(inputs: list[Path], output: Path) -> None:
    fieldnames = None
    with output.open("w", newline="", encoding="utf-8") as destination:
        writer = None
        for path in inputs:
            with path.open(newline="", encoding="utf-8") as source:
                reader = csv.DictReader(source)
                current_fields = reader.fieldnames or []
                if fieldnames is None:
                    fieldnames = current_fields
                    writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
                    writer.writeheader()
                elif current_fields != fieldnames:
                    raise ValueError(f"CSV columns differ between chain outputs: {path}")
                writer.writerows(reader)
    if not fieldnames:
        raise ValueError("No CSV inputs were available to merge.")


def renumber_pdb_atom_serials(path: Path) -> None:
    serial = 0
    output_lines = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")):
            serial += 1
            if serial > 99999:
                raise ValueError(f"Combined PDB exceeds the 99,999 atom serial limit: {path}")
            line = f"{line[:6]}{serial:5d}{line[11:]}"
        output_lines.append(line)
    path.write_text("".join(output_lines), encoding="utf-8")


def pdb_slice_sort_key(path: Path) -> int:
    match = re.match(r"slice_(\d+)_first_frame\.pdb$", path.name)
    return int(match.group(1)) if match else 0


def stage_combined_pdb_slices(
    chain_outputs: list[ChainOutput],
    output_dir: Path,
    chain_index_path: Path,
    expected_slices: int,
) -> list[Path]:
    """Combine chain PDBs while recording logical chain atom serial ranges."""
    if not chain_outputs:
        raise ValueError("No chain outputs were available to stage PDB slices.")

    output_dir.mkdir(parents=True, exist_ok=True)
    source_paths = sorted(
        chain_outputs[0]["directory"].glob("slice_*_first_frame.pdb"),
        key=pdb_slice_sort_key,
    )
    if len(source_paths) != expected_slices:
        raise RuntimeError(f"Expected {expected_slices} PDB slices, found {len(source_paths)}.")

    slice_ranges = {}
    staged_paths = []
    for source_path in source_paths:
        header_lines = []
        first_lines = source_path.read_text(encoding="utf-8").splitlines()
        for line in first_lines:
            if line.startswith(("ATOM", "HETATM")):
                break
            if not line.startswith(("END", "MODEL", "ENDMDL", "TER")):
                header_lines.append(line)

        serial = 0
        atom_lines = []
        ranges = []
        for output in chain_outputs:
            chain_path = output["directory"] / source_path.name
            if not chain_path.is_file():
                raise FileNotFoundError(f"PDB slice missing for chain {output['chain']}: {chain_path}")
            source_atoms = [
                line
                for line in chain_path.read_text(encoding="utf-8").splitlines()
                if line.startswith(("ATOM", "HETATM"))
            ]
            if not source_atoms:
                raise ValueError(f"No atom records found for chain {output['chain']}: {chain_path}")

            start_serial = serial + 1
            for line in source_atoms:
                serial += 1
                if serial > 99999:
                    raise ValueError(f"Combined PDB exceeds the 99,999 atom serial limit: {source_path.name}")
                padded = line.ljust(11)
                atom_lines.append(f"{padded[:6]}{serial:5d}{padded[11:]}")
            ranges.append(
                {
                    "chain": output["chain"],
                    "atomSerialStart": start_serial,
                    "atomSerialEnd": serial,
                }
            )

        staged_path = output_dir / source_path.name
        staged_path.write_text(
            "\n".join([*header_lines, *atom_lines, "END"]) + "\n",
            encoding="utf-8",
        )
        staged_paths.append(staged_path)
        slice_ranges[source_path.name] = ranges

    chain_index = {
        "version": 1,
        "chains": [
            {"id": output["chain"], "designation": output["designation"]}
            for output in chain_outputs
        ],
        "slices": slice_ranges,
    }
    chain_index_path.write_text(f"{json.dumps(chain_index, indent=2)}\n", encoding="utf-8")
    return staged_paths


def canonicalize_pdb_as_single_structure(path: Path) -> None:
    """Remove concatenated-chain boundaries emitted by older RMSX releases."""
    header_lines = []
    atom_lines = []
    saw_atom = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            saw_atom = True
            atom_lines.append(line)
        elif not saw_atom and not line.startswith(("END", "MODEL")):
            header_lines.append(line)

    if not atom_lines:
        raise ValueError(f"No atom records found in staged PDB slice: {path}")
    path.write_text("\n".join([*header_lines, *atom_lines, "END"]) + "\n", encoding="utf-8")


def run_static_plots(
    plot_helper: Path,
    chain_outputs: list[ChainOutput],
    heatmap_dir: Path,
    triple_dir: Path,
    palette: str,
    interpolate: bool,
    minimum: float,
    maximum: float,
    rscript: str,
) -> None:
    for output in chain_outputs:
        designation = output["designation"]
        command = [
            sys.executable,
            str(plot_helper),
            "--rmsx-source",
            str(output["rmsx"]),
            "--rmsd-source",
            str(output["rmsd"]),
            "--rmsf-source",
            str(output["rmsf"]),
            "--palette",
            palette,
            "--heatmap-output",
            str(heatmap_dir / f"{designation}.png"),
            "--triple-output",
            str(triple_dir / f"{designation}.png"),
            "--min-value",
            str(minimum),
            "--max-value",
            str(maximum),
            "--rscript",
            rscript,
        ]
        if interpolate:
            command.append("--interpolate")
        subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    import MDAnalysis as mda
    from rmsx.core import compute_global_rmsx_min_max, get_selection_string, run_rmsx

    output_dir = Path(args.output_dir)
    galaxy_output_dir = Path(args.galaxy_output_dir)
    collection_dirs = {
        name: galaxy_output_dir / name
        for name in (
            "rmsx_tables",
            "rmsd_tables",
            "rmsf_tables",
            "mask_tables",
            "pdb_slices",
            "heatmap_plots",
            "triple_plots",
        )
    }
    for directory in collection_dirs.values():
        directory.mkdir(parents=True, exist_ok=True)

    topology = str(Path(args.topology))
    trajectory = str(Path(args.trajectory))
    universe = mda.Universe(topology)
    chains = resolve_chains(
        universe,
        args.chain_mode,
        args.chains,
        args.analysis_type,
        get_selection_string,
    )
    print(f"Resolved RMSX segment IDs: {', '.join(chains)}")

    chain_outputs: list[ChainOutput] = []
    for chain in chains:
        print(f"Running RMSX for segment {chain}")
        options = {
            "output_dir": str(output_dir),
            "num_slices": args.num_slices,
            "verbose": False,
            "interpolate": args.interpolate,
            "triple": False,
            "chain_sele": chain,
            "overwrite": True,
            "palette": args.palette,
            "start_frame": args.start_frame,
            "make_plot": False,
            "analysis_type": args.analysis_type,
            "summary_n": args.summary_n,
        }
        if args.end_frame is not None:
            options["end_frame"] = args.end_frame
        run_rmsx(topology, trajectory, **options)

        chain_dir = output_dir / f"chain_{chain}_rmsx"
        rmsx_csv = one_matching_file(chain_dir, "rmsx_*.csv", "RMSX CSV")
        rmsd_csv = chain_dir / "rmsd.csv"
        rmsf_csv = chain_dir / "rmsf.csv"
        for label, path in (("RMSD CSV", rmsd_csv), ("RMSF CSV", rmsf_csv)):
            if not path.is_file():
                raise FileNotFoundError(f"{label} not found: {path}")

        designation = f"chain_{chain}"
        staged_mask = collection_dirs["mask_tables"] / f"{designation}.csv"
        source_mask = chain_dir / "masked_residues.csv"
        if source_mask.is_file():
            shutil.copyfile(source_mask, staged_mask)
        else:
            write_false_mask(rmsx_csv, staged_mask)

        shutil.copyfile(rmsx_csv, collection_dirs["rmsx_tables"] / f"{designation}.csv")
        shutil.copyfile(rmsd_csv, collection_dirs["rmsd_tables"] / f"{designation}.csv")
        shutil.copyfile(rmsf_csv, collection_dirs["rmsf_tables"] / f"{designation}.csv")
        chain_outputs.append(
            {
                "chain": chain,
                "designation": designation,
                "directory": chain_dir,
                "rmsx": rmsx_csv,
                "rmsd": rmsd_csv,
                "rmsf": rmsf_csv,
                "mask": staged_mask,
            }
        )

    stage_combined_pdb_slices(
        chain_outputs,
        collection_dirs["pdb_slices"],
        galaxy_output_dir / "viewer_chain_index.json",
        args.num_slices,
    )

    rmsx_paths = [output["rmsx"] for output in chain_outputs]
    mask_paths = [output["mask"] for output in chain_outputs]
    merge_csv_files(rmsx_paths, galaxy_output_dir / "combined_rmsx.csv")
    merge_csv_files(mask_paths, galaxy_output_dir / "combined_mask.csv")
    global_min, global_max = compute_global_rmsx_min_max([str(path) for path in rmsx_paths])
    print(f"Shared RMSX range: {global_min:.6f} to {global_max:.6f}")

    plot_min, plot_max = global_min, global_max
    if plot_min == plot_max:
        plot_max += max(abs(plot_min) * 1e-6, 1e-6)

    run_static_plots(
        Path(args.plot_helper),
        chain_outputs,
        collection_dirs["heatmap_plots"],
        collection_dirs["triple_plots"],
        args.palette,
        args.interpolate,
        plot_min,
        plot_max,
        args.rscript,
    )
    print(f"Galaxy outputs staged for {len(chain_outputs)} chain(s).")


if __name__ == "__main__":
    main()
