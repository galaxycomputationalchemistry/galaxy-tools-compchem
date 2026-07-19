#!/usr/bin/env python3
"""Shared helpers for Flipbook Galaxy viewer manifest builders."""

import csv
import re
from pathlib import Path


SLICE_RE = re.compile(r"slice_(\d+)_first_frame\.pdb$")


def slice_sort_key(path):
    match = SLICE_RE.match(Path(path).name)
    return int(match.group(1)) if match else 0


def read_slices(pdb_dir):
    slices = []
    for path in sorted(Path(pdb_dir).glob("slice_*_first_frame.pdb"), key=slice_sort_key):
        match = SLICE_RE.match(path.name)
        if not match:
            continue
        slice_index = int(match.group(1))
        slices.append(
            {
                "index": slice_index,
                "id": f"slice_{slice_index}",
                "label": f"Slice {slice_index}",
                "filename": path.name,
                "rmsxColumn": f"slice_{slice_index}.dcd",
                "pdb": path.read_text(encoding="utf-8"),
            }
        )
    if not slices:
        raise ValueError(f"No slice_*_first_frame.pdb files found in {pdb_dir}")
    return slices


def read_rmsx_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"RMSX table is empty: {path}")
    slice_columns = [name for name in fieldnames if name.startswith("slice_") and name.endswith(".dcd")]
    if not slice_columns:
        raise ValueError(f"RMSX table has no slice columns: {path}")
    return rows, slice_columns


def summarize_slices(rmsx_rows, slice_columns):
    summaries = {}
    all_values = []
    for column in slice_columns:
        values = []
        max_residue = None
        max_value = None
        for row in rmsx_rows:
            try:
                value = float(row[column])
            except (KeyError, TypeError, ValueError):
                continue
            values.append(value)
            all_values.append(value)
            if max_value is None or value > max_value:
                max_value = value
                max_residue = row.get("ResidueID", "")
        if not values:
            continue
        summaries[column] = {
            "min": min(values),
            "max": max(values),
            "mean": sum(values) / len(values),
            "maxResidue": max_residue,
            "residueCount": len(values),
        }
    if not all_values:
        raise ValueError("No numeric RMSX values found in RMSX table.")
    return summaries, {"min": min(all_values), "max": max(all_values)}


def build_residue_payload(rmsx_rows, slice_columns):
    residues = []
    for row in rmsx_rows:
        residue_id = (row.get("ResidueID") or "").strip()
        if not residue_id:
            continue
        chain_id = (row.get("ChainID") or "").strip()
        values = {}
        for column in slice_columns:
            try:
                values[column] = float(row[column])
            except (KeyError, TypeError, ValueError):
                continue
        key = f"{chain_id}:{residue_id}" if chain_id else residue_id
        label = f"{residue_id} / chain {chain_id}" if chain_id else residue_id
        residues.append(
            {
                "id": residue_id,
                "chain": chain_id,
                "key": key,
                "label": label,
                "values": values,
            }
        )
    if not residues:
        raise ValueError("No residues found in RMSX table.")
    return residues
