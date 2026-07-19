#!/usr/bin/env python3
"""Validate Flipbook Galaxy wrapper inputs before launching RMSX."""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, help="Topology or structure file passed to RMSX.")
    parser.add_argument("--trajectory", required=True, help="Trajectory file passed to RMSX.")
    parser.add_argument("--selector", required=True, help="Comma-separated RMSX chain/segment selector.")
    parser.add_argument("--num-slices", required=True, type=int, help="Requested number of RMSX slices.")
    parser.add_argument("--start-frame", required=True, type=int, help="First trajectory frame to analyze.")
    parser.add_argument("--end-frame", type=int, help="Last trajectory frame to analyze.")
    parser.add_argument("--analysis-type", required=True, help="RMSX analysis type.")
    return parser.parse_args()


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def compact_list(values, limit=20):
    if not values:
        return "none"
    shown = values[:limit]
    suffix = "" if len(values) <= limit else f", ... ({len(values)} total)"
    return ", ".join(shown) + suffix


def unique_atom_values(atoms, attribute):
    try:
        raw_values = getattr(atoms, attribute)
    except Exception:
        return []
    values = sorted({str(value).strip() for value in raw_values if str(value).strip()})
    return values


def requested_selectors(selector):
    selectors = [part.strip() for part in selector.split(",") if part.strip()]
    return selectors


def main():
    args = parse_args()
    print("Flipbook Galaxy preflight")
    print(f"topology: {args.topology}")
    print(f"trajectory: {args.trajectory}")
    print(f"selector: {args.selector}")
    print(f"analysis type: {args.analysis_type}")
    print(f"requested slices: {args.num_slices}")
    print(f"start frame: {args.start_frame}")
    print(f"end frame: {args.end_frame if args.end_frame is not None else 'last'}")

    if args.num_slices < 1:
        return fail("Number of slices must be at least 1.")
    if args.start_frame < 0:
        return fail("Start frame must be 0 or greater.")
    if args.end_frame is not None and args.end_frame < args.start_frame:
        return fail("End frame must be greater than or equal to the start frame.")

    selectors = requested_selectors(args.selector)
    if not selectors:
        return fail("At least one chain/segment selector is required.")

    try:
        import MDAnalysis as mda
    except Exception as exc:
        return fail(f"MDAnalysis is not available for input validation: {exc}")

    try:
        universe = mda.Universe(args.topology, args.trajectory)
    except Exception as exc:
        return fail(f"MDAnalysis could not load the topology/trajectory pair: {exc}")

    atom_count = len(universe.atoms)
    try:
        frame_count = len(universe.trajectory)
    except Exception as exc:
        return fail(f"MDAnalysis loaded the files, but could not count trajectory frames: {exc}")

    print(f"atoms: {atom_count}")
    print(f"trajectory frames: {frame_count}")

    if atom_count < 1:
        return fail("The loaded topology contains no atoms.")
    if frame_count < 1:
        return fail("The loaded trajectory contains no frames.")
    if args.start_frame >= frame_count:
        return fail(f"Start frame {args.start_frame} is outside the trajectory frame range 0..{frame_count - 1}.")
    if args.end_frame is not None and args.end_frame >= frame_count:
        return fail(f"End frame {args.end_frame} is outside the trajectory frame range 0..{frame_count - 1}.")

    effective_end = args.end_frame if args.end_frame is not None else frame_count - 1
    frame_window = effective_end - args.start_frame + 1
    print(f"analysis frame window: {frame_window} frame(s)")
    if args.num_slices > frame_window:
        return fail(
            f"Requested {args.num_slices} slices, but only {frame_window} frame(s) are available "
            "in the selected frame window."
        )

    segids = unique_atom_values(universe.atoms, "segids")
    chain_ids = unique_atom_values(universe.atoms, "chainIDs")
    matched_segids = sorted(set(selectors).intersection(segids))
    matched_chain_ids = sorted(set(selectors).intersection(chain_ids))
    missing = [selector for selector in selectors if selector not in segids and selector not in chain_ids]

    print(f"available segment IDs: {compact_list(segids)}")
    print(f"available chain IDs: {compact_list(chain_ids)}")
    print(f"matched segment IDs: {compact_list(matched_segids)}")
    print(f"matched chain IDs: {compact_list(matched_chain_ids)}")

    if missing:
        return fail(
            "Selector value(s) not found in topology segment IDs or chain IDs: "
            f"{compact_list(missing)}. Available segment IDs: {compact_list(segids)}. "
            f"Available chain IDs: {compact_list(chain_ids)}."
        )
    if matched_chain_ids and not matched_segids:
        print(
            "WARNING: selector matched chain IDs but not segment IDs. RMSX receives the same value via "
            "--chain; if RMSX later reports an empty selection, try a segment ID from the list above."
        )

    print("preflight status: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
