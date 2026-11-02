#!/usr/bin/env python3
"""Require every selected Galaxy test and every raw tool-group report."""

import argparse
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

from galaxy.tool_util.parser.factory import get_tool_source


def expected_groups(path):
    groups = []
    for line in Path(path).read_text().splitlines():
        paths = [str(Path(value).resolve()) for value in line.split()]
        if not paths:
            continue
        cases = []
        for tool_path in paths:
            source = get_tool_source(tool_path)
            tool_id, version = source.parse_id(), source.parse_version()
            tests = source.parse_tests_to_dict()["tests"]
            if not tool_id or not version or not tests:
                raise ValueError(
                    f"Tool has no identity/version/tests: {tool_path}"
                )
            cases.extend(
                (tool_id, version, index) for index in range(len(tests))
            )
        require_unique(cases, f"selected tool group {paths}")
        groups.append(dict(paths=paths, cases=frozenset(cases)))
    if not groups:
        raise ValueError(f"No tool groups selected in {path}")
    require_unique(
        [case for group in groups for case in group["cases"]], "selected tests"
    )
    return groups


def require_unique(values, label):
    duplicates = [
        value for value, count in Counter(values).items() if count > 1
    ]
    if duplicates:
        raise ValueError(f"Duplicate {label}: {duplicates}")


def report_cases(path):
    report = json.loads(Path(path).read_text())
    tests = report.get("tests") if isinstance(report, dict) else None
    if not isinstance(tests, list) or not tests:
        raise ValueError(f"No executed tests in {path}")
    cases = []
    for test in tests:
        data = test.get("data") if isinstance(test, dict) else None
        if not isinstance(data, dict) or data.get("status") != "success":
            raise ValueError(f"Failed, skipped or malformed test in {path}")
        tool_id, version, index = (
            data.get(key) for key in ("tool_id", "tool_version", "test_index")
        )
        if (
            not isinstance(tool_id, str)
            or not tool_id
            or not isinstance(version, str)
            or not version
            or type(index) is not int
            or index < 0
        ):
            raise ValueError(f"Missing exact test identity in {path}")
        cases.append((tool_id, version, index))
    require_unique(cases, f"test identities in {path}")
    return frozenset(cases)


def require_cases(actual, expected, label):
    if actual != expected:
        raise ValueError(
            f"{label}: missing={sorted(expected - actual)}, "
            f"unexpected={sorted(actual - expected)}"
        )


def check_chunk(artifact):
    artifact = Path(artifact)
    groups = expected_groups(artifact / "tool-groups.txt")
    raw_reports = sorted((artifact / "raw-reports").glob("*.json"))
    if len(raw_reports) != len(groups):
        raise ValueError(
            f"{artifact}: expected {len(groups)} group reports, "
            f"found {len(raw_reports)}"
        )
    actual_groups = [report_cases(report) for report in raw_reports]
    if Counter(actual_groups) != Counter(group["cases"] for group in groups):
        raise ValueError(
            f"{artifact}: missing, duplicate or incomplete tool-group report"
        )
    expected = frozenset(case for group in groups for case in group["cases"])
    require_cases(
        report_cases(artifact / "tool_test_output.json"),
        expected,
        str(artifact),
    )
    return groups


def discover_groups(repositories):
    paths = [
        line.strip()
        for line in Path(repositories).read_text().splitlines()
        if line.strip()
    ]
    if not paths:
        raise ValueError("No repositories selected")
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "all-tool-groups.txt"
        subprocess.run(
            [
                "planemo",
                "ci_find_tools",
                "--group_tools",
                "--output",
                str(output),
                *paths,
            ],
            check=True,
        )
        return expected_groups(output)


def check_combined(repositories, artifacts, report, chunks):
    if chunks < 1:
        raise ValueError("Expected a positive chunk count")
    expected = discover_groups(repositories)
    directories = sorted(
        path
        for path in Path(artifacts).glob("Tool test output *")
        if path.is_dir()
    )
    if {path.name for path in directories} != {
        f"Tool test output {index}" for index in range(chunks)
    }:
        raise ValueError("Missing or unexpected chunk artifacts")
    groups = [
        group for directory in directories for group in check_chunk(directory)
    ]
    actual_paths = Counter(path for group in groups for path in group["paths"])
    expected_paths = Counter(
        path for group in expected for path in group["paths"]
    )
    if actual_paths != expected_paths:
        raise ValueError("Chunk tool lists omit or duplicate selected tools")
    expected_cases = frozenset(
        case for group in expected for case in group["cases"]
    )
    require_cases(report_cases(report), expected_cases, "Combined report")
    return len(expected_cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    chunk = commands.add_parser("chunk")
    chunk.add_argument("artifact", type=Path)
    combined = commands.add_parser("combined")
    combined.add_argument("--repositories", required=True, type=Path)
    combined.add_argument("--artifacts", required=True, type=Path)
    combined.add_argument("--report", required=True, type=Path)
    combined.add_argument("--chunks", required=True, type=int)
    args = parser.parse_args()
    try:
        if args.mode == "chunk":
            groups = check_chunk(args.artifact)
            count = sum(len(group["cases"]) for group in groups)
        else:
            count = check_combined(
                args.repositories, args.artifacts, args.report, args.chunks
            )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.CalledProcessError,
    ) as exc:
        parser.exit(1, f"Test completeness check failed: {exc}\n")
    print(f"All {count} expected tool tests executed exactly once and passed")


if __name__ == "__main__":
    main()
