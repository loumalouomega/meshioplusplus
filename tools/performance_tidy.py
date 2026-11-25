"""Run the installed Clang-Tidy performance checks without modifying sources.

Thin wrapper around :mod:`tidy_report` defaulting to ``-*,performance-*``.
"""

from __future__ import annotations

import importlib.util
import pathlib
import subprocess
import sys

_TOOLS_DIR = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "tidy_report", _TOOLS_DIR / "tidy_report.py"
)
assert _spec is not None and _spec.loader is not None
_tidy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_tidy)

ROOT = _tidy.ROOT
CHECKS = _tidy.DEFAULT_CHECKS
DIAGNOSTIC = _tidy.DIAGNOSTIC


def _sync() -> None:
    _tidy.ROOT = ROOT
    _tidy.subprocess = subprocess


def tidy_config():
    _sync()
    return _tidy.tidy_config(CHECKS)


def tidy_arguments():
    _sync()
    return _tidy.tidy_arguments(CHECKS)


def first_party(path):
    _sync()
    return _tidy.first_party(path)


def source_path(entry):
    return _tidy.source_path(entry)


def configure(build_dir, backend):
    _sync()
    return _tidy.configure(build_dir, backend)


def configuration(build_dir, entries):
    return _tidy.configuration(build_dir, entries)


def parse_diagnostics(output, directory):
    _sync()
    return _tidy.parse_diagnostics(output, directory)


def scan(tool, build_dir, source, directory, log_dir):
    _sync()
    return _tidy.scan(tool, build_dir, source, directory, log_dir, checks=CHECKS)


def positive_int(value):
    return _tidy.positive_int(value)


def main(argv=None):
    _sync()
    return _tidy.main(
        argv,
        default_checks=CHECKS,
        default_out_name="performance-tidy.json",
        prefix="performance-tidy",
    )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        print(f"performance-tidy: {exc}", file=sys.stderr)
        sys.exit(1)
