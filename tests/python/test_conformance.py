"""The format conformance matrix (roadmap §3; declarations in conformance_spec.py).

Every writable format must have a declaration, and what a write and a read
actually keep must match it exactly -- in both directions, so a format that
starts keeping something is as visible as one that starts losing it. The
declarations also render doc/conformance.md and the "Round trip" column of
doc/formats.md; the last test fails when those pages are stale.
"""

import pathlib
import subprocess
import sys

import pytest

import meshioplusplus

from . import conformance_spec as cs

REPO = pathlib.Path(__file__).resolve().parents[2]
WRITABLE = sorted(meshioplusplus.formats()["writable"])


def test_every_writable_format_is_declared():
    missing = [f for f in WRITABLE if f not in cs.SPEC]
    assert not missing, (
        "no conformance declaration for: " + ", ".join(missing) + " -- run "
        "`python tools/gen_conformance_table.py --observe <fmt>` and review it"
    )


def test_no_declaration_for_a_format_that_is_gone():
    assert not sorted(set(cs.SPEC) - set(WRITABLE))


@pytest.mark.parametrize("fmt", WRITABLE)
def test_round_trip_matches_the_declaration(fmt):
    spec = cs.SPEC.get(fmt)
    if spec is None:
        pytest.skip("undeclared (reported by test_every_writable_format_is_declared)")
    declared = {k: v for k, v in spec.items() if k != "note"}
    cells = list(declared["cells"]) if "cells" in declared else None
    observed = cs.observe(fmt, cells, planar=declared.get("input") == "2d")
    if observed.get("error") == "ImportError":
        pytest.skip(f"{fmt}: an optional dependency is not installed")
    assert observed == declared


def test_the_generated_pages_are_current():
    run = subprocess.run(
        [sys.executable, str(REPO / "tools" / "gen_conformance_table.py"), "--check"],
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0, run.stdout + run.stderr
