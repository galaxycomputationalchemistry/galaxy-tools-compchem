"""Regression tests for incomplete/duplicate CI reports."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "reports", Path(__file__).with_name("check_test_reports.py")
)
reports = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reports)


class CompletenessTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.first = self.tool("first", 5)
        self.second = self.tool("second", 2)
        self.first_cases = [("first", "1.0", index) for index in range(5)]
        self.second_cases = [("second", "1.0", index) for index in range(2)]

    def tool(self, name, count):
        path = self.root / f"{name}.xml"
        path.write_text(
            f'<tool id="{name}" version="@VERSION@"><macros>'
            '<token name="@VERSION@">1.0</token>'
            '<xml name="case"><test/></xml></macros><tests>'
            + '<expand macro="case"/>' * count
            + "</tests></tool>"
        )
        return path

    def write_report(self, path, cases, status="success"):
        path.write_text(
            json.dumps(
                {
                    "tests": [
                        dict(
                            data=dict(
                                tool_id=tool,
                                tool_version=version,
                                test_index=index,
                                status=status,
                            )
                        )
                        for tool, version, index in cases
                    ]
                }
            )
        )

    def artifact(self, name="Tool test output 0", groups=None):
        groups = groups or [[self.first]]
        directory = self.root / "artifacts" / name
        (directory / "raw-reports").mkdir(parents=True)
        (directory / "tool-groups.txt").write_text(
            "\n".join(" ".join(map(str, group)) for group in groups)
        )
        cases = []
        for index, group in enumerate(
            reports.expected_groups(directory / "tool-groups.txt")
        ):
            values = sorted(group["cases"])
            self.write_report(
                directory / "raw-reports" / f"{index}.json", values
            )
            cases.extend(values)
        self.write_report(directory / "tool_test_output.json", cases)
        return directory

    def test_complete_macro_expanded_group_passes(self):
        groups = reports.check_chunk(
            self.artifact(groups=[[self.first, self.second]])
        )
        self.assertEqual(sum(len(group["cases"]) for group in groups), 7)

    def test_one_or_four_of_five_successes_cannot_pass(self):
        directory = self.artifact()
        for count in (1, 4):
            with self.subTest(count=count):
                self.write_report(
                    directory / "raw-reports/0.json", self.first_cases[:count]
                )
                self.write_report(
                    directory / "tool_test_output.json",
                    self.first_cases[:count],
                )
                with self.assertRaises(ValueError):
                    reports.check_chunk(directory)

    def test_missing_group_report_cannot_disappear_in_merge(self):
        directory = self.artifact(groups=[[self.first], [self.second]])
        (directory / "raw-reports/1.json").unlink()
        self.write_report(
            directory / "tool_test_output.json", self.first_cases
        )
        with self.assertRaises(ValueError):
            reports.check_chunk(directory)

    def test_duplicate_group_cannot_replace_missing_group(self):
        directory = self.artifact(groups=[[self.first], [self.second]])
        self.write_report(directory / "raw-reports/1.json", self.first_cases)
        with self.assertRaises(ValueError):
            reports.check_chunk(directory)

    def test_duplicate_case_cannot_replace_missing_case(self):
        directory = self.artifact()
        self.write_report(
            directory / "raw-reports/0.json",
            self.first_cases[:4] + [self.first_cases[0]],
        )
        with self.assertRaises(ValueError):
            reports.check_chunk(directory)

    def test_merged_report_must_retain_all_cases_exactly_once(self):
        directory = self.artifact()
        for cases in (
            self.first_cases[:4],
            self.first_cases + [self.first_cases[0]],
        ):
            self.write_report(directory / "tool_test_output.json", cases)
            with self.assertRaises(ValueError):
                reports.check_chunk(directory)

    def test_wrong_tool_or_version_is_not_equivalent(self):
        directory = self.artifact()
        for tool, version in [("other", "1.0"), ("first", "2.0")]:
            self.write_report(
                directory / "raw-reports/0.json",
                [(tool, version, index) for index in range(5)],
            )
            with self.assertRaises(ValueError):
                reports.check_chunk(directory)

    def test_empty_skipped_and_malformed_reports_fail(self):
        directory = self.artifact()
        for contents in (
            {"tests": []},
            {"tests": [None]},
            {"tests": [{"data": {"status": "skip"}}]},
            [],
        ):
            (directory / "raw-reports/0.json").write_text(json.dumps(contents))
            with self.assertRaises(ValueError):
                reports.check_chunk(directory)

    def test_independent_repository_discovery_catches_omitted_tool(self):
        directory = self.artifact()
        selected = self.root / "selected.txt"
        selected.write_text(f"{self.first}\n{self.second}\n")
        expected = reports.expected_groups(selected)
        with patch.object(reports, "discover_groups", return_value=expected):
            with self.assertRaisesRegex(ValueError, "omit or duplicate"):
                reports.check_combined(
                    "repositories.txt",
                    directory.parent,
                    directory / "tool_test_output.json",
                    1,
                )

    def test_combined_requires_numbered_chunks_and_unique_cases(self):
        first = self.artifact()
        second = self.artifact("Tool test output 1", [[self.second]])
        selected = self.root / "selected.txt"
        selected.write_text(f"{self.first}\n{self.second}\n")
        expected = reports.expected_groups(selected)
        merged = self.root / "combined.json"
        self.write_report(merged, self.first_cases + self.second_cases)
        with patch.object(reports, "discover_groups", return_value=expected):
            self.assertEqual(
                reports.check_combined(
                    "repositories.txt", first.parent, merged, 2
                ),
                7,
            )
            with self.assertRaisesRegex(ValueError, "chunk artifacts"):
                reports.check_combined(
                    "repositories.txt", first.parent, merged, 3
                )
            self.write_report(
                merged,
                self.first_cases + self.second_cases + [self.first_cases[0]],
            )
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                reports.check_combined(
                    "repositories.txt", second.parent, merged, 2
                )

    def test_real_flipbook_xml_requires_all_five_cases(self):
        selected = self.root / "selected.txt"
        wrapper = (
            Path(__file__).resolve().parents[1] / "tools/flipbook/flipbook.xml"
        )
        selected.write_text(str(wrapper.relative_to(Path.cwd())))
        expected = reports.expected_groups(selected)
        self.assertEqual(
            expected[0]["cases"],
            frozenset(
                ("flipbook", "0.1.0+galaxy1", index) for index in range(5)
            ),
        )


if __name__ == "__main__":
    unittest.main()
