#!/usr/bin/env python3
"""Focused tests for the optional RMSX Analysis manifest contract."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from flipbook_molstar_report import (  # noqa: E402
    annotate_slices_for_analysis,
    build_analysis_payload,
    downsample_rmsd_points,
    read_rmsd_points,
)
from rmsx_multichain import stage_combined_pdb_slices  # noqa: E402


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pdb_atom(serial, chain, residue, x):
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain:1s}{residue:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}          C\n"
    )


class AnalysisMetricTests(unittest.TestCase):
    def test_downsampling_preserves_order_endpoints_and_extrema(self):
        points = [(index, index * 0.01, float((index * 17) % 101)) for index in range(5000)]
        points[1777] = (1777, 17.77, -25.0)
        points[3888] = (3888, 38.88, 250.0)

        sampled, method = downsample_rmsd_points(points, maximum=64)

        self.assertEqual(method, "min-max-bin")
        self.assertLessEqual(len(sampled), 64)
        self.assertEqual(sampled[0], points[0])
        self.assertEqual(sampled[-1], points[-1])
        self.assertEqual([point[0] for point in sampled], sorted(point[0] for point in sampled))
        self.assertIn(points[1777], sampled)
        self.assertIn(points[3888], sampled)

        interior = points[1:-1]
        bin_count = (64 - 2) // 2
        sampled_set = set(sampled)
        for bin_index in range(bin_count):
            start = bin_index * len(interior) // bin_count
            end = (bin_index + 1) * len(interior) // bin_count
            bucket = interior[start:end]
            self.assertIn(min(bucket, key=lambda point: point[2]), sampled_set)
            self.assertIn(max(bucket, key=lambda point: point[2]), sampled_set)

    def test_payload_uses_nanoseconds_and_retains_rmsf_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rmsd_dir = root / "rmsd"
            rmsf_dir = root / "rmsf"
            rmsd_dir.mkdir()
            rmsf_dir.mkdir()
            write_csv(
                rmsd_dir / "chain_SYSTEM.csv",
                ["Frame", "Time", "RMSD"],
                [
                    {"Frame": 0, "Time": 0, "RMSD": 0.5},
                    {"Frame": 1, "Time": 1000, "RMSD": 1.5},
                    {"Frame": 2, "Time": 2500, "RMSD": 0.75},
                ],
            )
            write_csv(
                rmsf_dir / "chain_SYSTEM.csv",
                ["ResidueID", "RMSF"],
                [
                    {"ResidueID": "10", "RMSF": 0.2},
                    {"ResidueID": "2A", "RMSF": 0.4},
                    {"ResidueID": "11", "RMSF": 0.3},
                ],
            )
            chain_index = {
                "version": 1,
                "chains": [{"id": "SYSTEM", "designation": "chain_SYSTEM"}],
                "slices": {},
            }

            payload = build_analysis_payload(rmsd_dir, rmsf_dir, chain_index)

        self.assertEqual(payload["timeDomainNs"], [0.0, 2.5])
        self.assertEqual(payload["chains"][0]["rmsd"]["timeNs"], [0.0, 1.0, 2.5])
        self.assertEqual(payload["chains"][0]["rmsf"]["residueIds"], ["10", "2A", "11"])
        self.assertEqual(payload["chains"][0]["rmsf"]["values"], [0.2, 0.4, 0.3])

    def test_slice_annotations_share_time_domain_and_chain_ranges(self):
        slices = [
            {"filename": "slice_1_first_frame.pdb"},
            {"filename": "slice_2_first_frame.pdb"},
        ]
        analysis = {"timeDomainNs": [2.0, 6.0]}
        chain_index = {
            "slices": {
                "slice_1_first_frame.pdb": [
                    {"chain": "SYSTEM", "atomSerialStart": 1, "atomSerialEnd": 2}
                ],
                "slice_2_first_frame.pdb": [
                    {"chain": "SYSTEM", "atomSerialStart": 1, "atomSerialEnd": 2}
                ],
            }
        }

        annotate_slices_for_analysis(slices, analysis, chain_index)

        self.assertEqual(slices[0]["time"], {"startNs": 2.0, "endNs": 4.0, "centerNs": 3.0})
        self.assertEqual(slices[1]["time"], {"startNs": 4.0, "endNs": 6.0, "centerNs": 5.0})
        self.assertEqual(slices[0]["chainAtomRanges"][0]["chain"], "SYSTEM")


class MetricValidationTests(unittest.TestCase):
    def test_out_of_order_frames_fail_instead_of_drawing_misleading_time_axis(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rmsd.csv"
            path.write_text("Frame,Time,RMSD\n2,20,1\n1,10,2\n")
            with self.assertRaisesRegex(ValueError, "chronological"):
                read_rmsd_points(path)


class ChainIndexTests(unittest.TestCase):
    def test_combined_slices_record_logical_chain_atom_ranges(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outputs = []
            for chain, atom_count in (("SYSTEM", 2), ("LIGAND_SEGMENT", 3)):
                directory = root / chain
                directory.mkdir()
                for slice_index in (1, 2):
                    atoms = "".join(
                        pdb_atom(index, "A", index, float(index + slice_index))
                        for index in range(1, atom_count + 1)
                    )
                    (directory / f"slice_{slice_index}_first_frame.pdb").write_text(
                        f"HEADER TEST\n{atoms}END\n",
                        encoding="utf-8",
                    )
                outputs.append(
                    {
                        "chain": chain,
                        "designation": f"chain_{chain}",
                        "directory": directory,
                        "rmsx": root / "unused.csv",
                        "rmsd": root / "unused.csv",
                        "rmsf": root / "unused.csv",
                        "mask": root / "unused.csv",
                    }
                )

            (root / "unused.csv").write_text("Frame,Time,RMSD\n2,200,0.1\n3,300,0.2\n4,400,0.3\n5,500,0.4\n")
            output_dir = root / "combined"
            index_path = root / "viewer_chain_index.json"
            staged = stage_combined_pdb_slices(outputs, output_dir, index_path, expected_slices=2)
            index = json.loads(index_path.read_text(encoding="utf-8"))
            self.assertEqual(index["sliceTimes"][staged[0].name]["startFrame"], 2)
            self.assertEqual(index["sliceTimes"][staged[0].name]["endNs"], 0.3)
            annotated = [{"filename": path.name} for path in staged]
            annotate_slices_for_analysis(annotated, {"timeDomainNs": [0.2, 0.5]}, index)
            self.assertEqual(annotated[0]["time"]["endNs"], 0.3)
            atom_serials = [
                int(line[6:11])
                for line in staged[0].read_text(encoding="utf-8").splitlines()
                if line.startswith("ATOM")
            ]

        self.assertEqual(len(staged), 2)
        self.assertEqual([chain["id"] for chain in index["chains"]], ["SYSTEM", "LIGAND_SEGMENT"])
        self.assertEqual(
            index["slices"]["slice_1_first_frame.pdb"],
            [
                {"chain": "SYSTEM", "atomSerialStart": 1, "atomSerialEnd": 2},
                {"chain": "LIGAND_SEGMENT", "atomSerialStart": 3, "atomSerialEnd": 5},
            ],
        )
        self.assertEqual(atom_serials, [1, 2, 3, 4, 5])


if __name__ == "__main__":
    unittest.main()
