"""The format conformance matrix (roadmap §2; declarations in conformance_spec.py).

Every writable format must have a declaration, and what a write and a read
actually keep must match it exactly -- in both directions, so a format that
starts keeping something is as visible as one that starts losing it. The
declarations also render doc/conformance.md and the "Round trip" column of
doc/formats.md; the last test fails when those pages are stale.
"""

import pathlib
import re
import subprocess
import sys

import pytest

import meshioplusplus

from . import conformance_spec as cs

REPO = pathlib.Path(__file__).resolve().parents[2]
WRITABLE = sorted(meshioplusplus.formats()["writable"])

# Writers (or the reader of a round trip) with no Python engine that exist only
# in some native builds, keyed to the core flag that says so. The declaration
# describes a build that has them; a build without one (the Windows CI job has
# no zlib, so no gidpost; the wheels have no OpenCASCADE, so STEP does not read
# back) skips.
_NEEDS_CORE = {"gid": "__has_gidpost__", "step": "__has_occt__"}


def test_every_writable_format_is_declared():
    missing = [f for f in WRITABLE if f not in cs.SPEC]
    assert not missing, (
        "no conformance declaration for: " + ", ".join(missing) + " -- run "
        "`python tools/gen_conformance_table.py --observe <fmt>` and review it"
    )


def test_no_declaration_for_a_format_that_is_gone():
    assert not sorted(set(cs.SPEC) - set(WRITABLE))


# Native registry keys that Python spells differently.
_NATIVE_ALIASES = {"ansysinp": "ansysInp", "dolfin": "dolfin-xml"}


def test_every_read_only_format_is_declared_with_its_reason():
    formats = meshioplusplus.formats()
    drift = cs.read_only_drift(formats["readable"], formats["writable"])
    assert drift == {"undeclared": [], "stale": []}, (
        "a format that reads without writing must name the roadmap item that "
        "closes the gap in conformance_spec.READ_ONLY (and one that writes must "
        "leave it)"
    )
    assert not set(cs.READ_ONLY) & set(cs.SPEC)
    assert all(reason.strip() for reason in cs.READ_ONLY.values())


def test_every_write_only_format_is_declared_with_its_item():
    formats = meshioplusplus.formats()
    drift = cs.write_only_drift(formats["readable"], formats["writable"])
    assert drift == {"undeclared": [], "stale": []}, (
        "a format that writes without reading must name the roadmap item that "
        "closes the gap in conformance_spec.WRITE_ONLY (and one that reads must "
        "leave it)"
    )
    assert all(reason.strip() for reason in cs.WRITE_ONLY.values())
    # The two lists and the round-trip matrix agree on what write-only means.
    assert {f for f, s in cs.SPEC.items() if s.get("error") == "write-only"} == set(
        cs.WRITE_ONLY
    )


def test_every_gap_names_a_roadmap_item():
    for fmt, reason in {**cs.READ_ONLY, **cs.WRITE_ONLY}.items():
        assert "roadmap" in reason, f"{fmt}: name the roadmap item in its reason"
        if reason.startswith("Exempt, not a gap"):
            assert "roadmap non-goals" in reason, f"{fmt}: an exemption is a decision"
        else:
            assert "roadmap §" in reason, f"{fmt}: a gap is closed by a roadmap item"


def test_the_native_registry_agrees():
    native = meshioplusplus._core.registry_formats()
    python = meshioplusplus.formats()

    def alias(names):
        return {_NATIVE_ALIASES.get(n, n) for n in names}

    readable, writable = alias(native["readable"]), alias(native["writable"])
    # Every native format is a Python one (Python adds its own formats and the
    # dependency-gated ones this build may lack).
    assert readable <= set(python["readable"])
    assert writable <= set(python["writable"])
    # Nothing declared read-only writes natively, and a format that reads
    # without writing natively is declared read-only or writes in Python.
    assert not writable & set(cs.READ_ONLY)
    assert readable - writable <= set(cs.READ_ONLY) | set(python["writable"])
    # The mirror: nothing declared write-only reads natively, and a format that
    # writes without reading natively is declared write-only or reads in Python.
    assert not readable & set(cs.WRITE_ONLY)
    read_anywhere = readable | set(python["readable"])
    assert {n for n in writable if cs._READ_AS.get(n, n) not in read_anywhere} <= set(
        cs.WRITE_ONLY
    )


def test_the_drift_check_fails_on_an_undeclared_write_only_format():
    # The probe: a format registered with a writer and no reader, and no
    # write-only reason, must be caught.
    from meshioplusplus._helpers import deregister_format, register_format

    register_format("probe-wo", [".probe-wo"], None, {"probe-wo": lambda *a: None})
    try:
        formats = meshioplusplus.formats()
        drift = cs.write_only_drift(formats["readable"], formats["writable"])
        assert drift["undeclared"] == ["probe-wo"]
        declared = {**cs.WRITE_ONLY, "probe-wo": "a test probe"}
        drift = cs.write_only_drift(formats["readable"], formats["writable"], declared)
        assert drift == {"undeclared": [], "stale": []}
    finally:
        deregister_format("probe-wo")
    assert "probe-wo" not in meshioplusplus.formats()["writable"]


def test_the_drift_check_fails_on_an_undeclared_read_only_format():
    # The probe: a format registered with a reader and no writer, and no
    # read-only reason, must be caught.
    from meshioplusplus._helpers import deregister_format, register_format

    register_format("probe-ro", [".probe-ro"], lambda filename: None, {})
    try:
        formats = meshioplusplus.formats()
        drift = cs.read_only_drift(formats["readable"], formats["writable"])
        assert drift["undeclared"] == ["probe-ro"]
        declared = {**cs.READ_ONLY, "probe-ro": "a test probe"}
        drift = cs.read_only_drift(formats["readable"], formats["writable"], declared)
        assert drift == {"undeclared": [], "stale": []}
    finally:
        deregister_format("probe-ro")
    assert "probe-ro" not in meshioplusplus.formats()["readable"]


@pytest.mark.parametrize("fmt", WRITABLE)
def test_round_trip_matches_the_declaration(fmt):
    spec = cs.SPEC.get(fmt)
    if spec is None:
        pytest.skip("undeclared (reported by test_every_writable_format_is_declared)")
    flag = _NEEDS_CORE.get(fmt)
    if flag and not getattr(meshioplusplus._core, flag, False):
        pytest.skip(f"{fmt}: this build is without {flag}")
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


# Read-only formats documented on another format's page.
_READ_ONLY_PAGE = {"ansys_rst_cyclic": "ansys_rst", "marc_t19": "marc"}


@pytest.mark.parametrize("fmt", sorted(cs.READ_ONLY))
def test_every_read_only_page_links_its_reason(fmt):
    page = REPO / "doc" / "formats" / f"{_READ_ONLY_PAGE.get(fmt, fmt)}.md"
    anchor = fmt.lower().replace("_", "-")
    assert f"conformance.md#{anchor}" in page.read_text(encoding="utf-8"), page.name


# Engine parity per direction (roadmap 3.4.13) ---------------------------------

_TWIN = re.compile(r"\b_(?:py|main)_(read|write)\b")
_PYTHON_ONLY_LAYERS = {"pmsh", "zarr", "cae", "usd"}
# Registered names that share a package with another format.
_ENGINE_ALIASES = {"ansys_rst_cyclic": "ansys_rst", "vtk42": "vtk", "vtk51": "vtk"}


def _package_of(fmt):
    table = (REPO / "src/python/meshioplusplus/_format_table.py").read_text()
    match = re.search(r'"%s": "([^"]+)"' % re.escape(fmt), table)
    return match.group(1) if match else fmt


def _twin_gaps():
    """``{(package, direction)}`` the Python twin of a core-served format lacks."""
    formats = meshioplusplus.formats()
    gaps = set()
    for fmt in sorted(set(formats["readable"]) | set(formats["writable"])):
        fmt = _ENGINE_ALIASES.get(fmt, fmt)
        if fmt in _PYTHON_ONLY_LAYERS:
            continue
        init = REPO / "src/python/meshioplusplus" / _package_of(fmt) / "__init__.py"
        found = set(_TWIN.findall(init.read_text())) if init.exists() else set()
        for direction in ("read", "write"):
            if direction not in found:
                gaps.add((fmt, direction))
    return gaps


def test_every_format_has_both_engines_or_a_recorded_exemption():
    declared = set(cs.ENGINE_EXEMPT) | {(f, "write") for f in cs.READ_ONLY}
    declared = {(_ENGINE_ALIASES.get(f, f), d) for f, d in declared}
    gaps = _twin_gaps()
    assert not gaps - declared, (
        "no Python twin for "
        + ", ".join(f"{f} ({d})" for f, d in sorted(gaps - declared))
        + " -- port it, or record the exemption in conformance_spec.ENGINE_EXEMPT"
    )
    assert (
        not set(cs.ENGINE_EXEMPT) - gaps
    ), "stale engine exemption (the twin now exists): " + ", ".join(
        f"{f} ({d})" for f, d in sorted(set(cs.ENGINE_EXEMPT) - gaps)
    )


def test_engine_exemptions_name_their_roadmap_item():
    for key, reason in cs.ENGINE_EXEMPT.items():
        assert "roadmap non-goals" in reason and reason.startswith("Exempt"), key
    for fmt, reason in cs.PARTIAL_ENGINE_EXEMPT.items():
        assert "roadmap non-goals" in reason, fmt
        assert (REPO / "doc/formats" / f"{fmt}.md").exists(), fmt
