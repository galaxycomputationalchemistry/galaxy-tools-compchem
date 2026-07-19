#!/usr/bin/env python3
"""Generate Galaxy-facing RMSX static plots from the packaged R plotting script."""

from __future__ import annotations

import argparse
import csv
import importlib.resources
import shutil
import subprocess
import sys
from pathlib import Path


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
R_PACKAGES = ("ggplot2", "viridis", "dplyr", "tidyr", "stringr", "readr", "gridExtra", "grid")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate RMSX heatmap and triple PNG plots.")
    parser.add_argument("--rmsx-source", required=True, help="Original RMSX CSV path with simulation time in its filename.")
    parser.add_argument("--rmsd-source", required=True, help="RMSD CSV path.")
    parser.add_argument("--rmsf-source", required=True, help="RMSF CSV path.")
    parser.add_argument("--palette", required=True, help="Viridis palette option to pass to RMSX R plotting.")
    parser.add_argument("--heatmap-output", required=True, help="Explicit Galaxy heatmap PNG output path.")
    parser.add_argument("--triple-output", required=True, help="Explicit Galaxy triple-plot PNG output path.")
    parser.add_argument("--interpolate", action="store_true", help="Interpolate the static RMSX heatmap raster.")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable path.")
    parser.add_argument("--plot-script", help="Override plot_rmsx.R path for tests/development.")
    return parser.parse_args()


def locate_plot_script(explicit: str | None = None) -> Path:
    if explicit:
        path = Path(explicit)
        if path.is_file():
            return path
        raise FileNotFoundError(f"Explicit RMSX plot script not found: {path}")

    try:
        candidate = importlib.resources.files("rmsx").joinpath("r_scripts", "plot_rmsx.R")
    except ModuleNotFoundError as exc:
        raise FileNotFoundError("Could not import the installed rmsx package to find plot_rmsx.R.") from exc

    with importlib.resources.as_file(candidate) as path:
        if path.is_file():
            return path
    raise FileNotFoundError("Installed rmsx package does not contain r_scripts/plot_rmsx.R.")


def read_single_chain_id(rmsx_csv: Path) -> str:
    with rmsx_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "ChainID" not in (reader.fieldnames or []):
            raise ValueError(f"RMSX CSV is missing required ChainID column: {rmsx_csv}")
        chain_ids = []
        for row in reader:
            chain_id = str(row.get("ChainID", "")).strip()
            if chain_id and chain_id not in chain_ids:
                chain_ids.append(chain_id)

    if not chain_ids:
        raise ValueError(f"RMSX CSV does not contain any ChainID values: {rmsx_csv}")
    if len(chain_ids) > 1:
        raise ValueError(
            "Static Galaxy plot outputs currently support one ChainID per run. "
            f"Found {', '.join(chain_ids)} in {rmsx_csv}; run chains separately for explicit PNG outputs."
        )
    return chain_ids[0]


def expected_r_plot_output(rmsx_csv: Path, chain_id: str) -> Path:
    return rmsx_csv.with_name(f"{rmsx_csv.stem}_rmsx_plot_chain_{chain_id}.png")


def verify_png(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} PNG was not created: {path}")
    if path.stat().st_size <= len(PNG_SIGNATURE):
        raise ValueError(f"{label} PNG is empty or truncated: {path}")
    with path.open("rb") as handle:
        signature = handle.read(len(PNG_SIGNATURE))
    if signature != PNG_SIGNATURE:
        raise ValueError(f"{label} output is not a PNG file: {path}")


def verify_r_packages(rscript: str) -> None:
    package_vector = ", ".join(f'"{package}"' for package in R_PACKAGES)
    command = [
        rscript,
        "-e",
        (
            f"needed <- c({package_vector}); "
            "missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]; "
            "if (length(missing)) stop(paste('Missing required R packages:', paste(missing, collapse=', ')), call. = FALSE)"
        ),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        details = "\n".join(part for part in (result.stdout.strip(), result.stderr.strip()) if part)
        if not details:
            details = f"Rscript package preflight exited with code {result.returncode}."
        raise RuntimeError(
            "Required R plotting packages are not available in the Galaxy runtime. "
            "Install them through Conda/container dependencies; runtime package installation is disabled.\n"
            f"{details}"
        )


def run_r_plot(
    *,
    rscript: str,
    plot_script: Path,
    rmsx_csv: Path,
    rmsd_csv: Path,
    rmsf_csv: Path,
    palette: str,
    interpolate: bool,
    triple: bool,
    generated_png: Path,
) -> None:
    if generated_png.exists():
        generated_png.unlink()
    command = [
        rscript,
        str(plot_script),
        str(rmsx_csv),
        str(rmsd_csv),
        str(rmsf_csv),
        "TRUE" if interpolate else "FALSE",
        "TRUE" if triple else "FALSE",
        palette,
        "",
        "",
        "FALSE",
        "",
        "FALSE",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        details = "\n".join(part for part in (result.stdout.strip(), result.stderr.strip()) if part)
        if not details:
            details = f"Rscript exited with code {result.returncode}."
        raise RuntimeError(f"RMSX static plot generation failed:\n{details}")
    verify_png(generated_png, "RMSX R plot")


def copy_verified_png(source: Path, destination: Path, label: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    verify_png(destination, label)


def main() -> None:
    args = parse_args()
    rmsx_csv = Path(args.rmsx_source)
    rmsd_csv = Path(args.rmsd_source)
    rmsf_csv = Path(args.rmsf_source)
    heatmap_output = Path(args.heatmap_output)
    triple_output = Path(args.triple_output)
    plot_script = locate_plot_script(args.plot_script)
    verify_r_packages(args.rscript)

    for label, path in (("RMSX CSV", rmsx_csv), ("RMSD CSV", rmsd_csv), ("RMSF CSV", rmsf_csv)):
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    chain_id = read_single_chain_id(rmsx_csv)
    generated_png = expected_r_plot_output(rmsx_csv, chain_id)

    run_r_plot(
        rscript=args.rscript,
        plot_script=plot_script,
        rmsx_csv=rmsx_csv,
        rmsd_csv=rmsd_csv,
        rmsf_csv=rmsf_csv,
        palette=args.palette,
        interpolate=args.interpolate,
        triple=True,
        generated_png=generated_png,
    )
    copy_verified_png(generated_png, triple_output, "RMSX triple plot")

    run_r_plot(
        rscript=args.rscript,
        plot_script=plot_script,
        rmsx_csv=rmsx_csv,
        rmsd_csv=rmsd_csv,
        rmsf_csv=rmsf_csv,
        palette=args.palette,
        interpolate=args.interpolate,
        triple=False,
        generated_png=generated_png,
    )
    copy_verified_png(generated_png, heatmap_output, "RMSX heatmap plot")

    print(f"RMSX static plots written: heatmap={heatmap_output}, triple={triple_output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
