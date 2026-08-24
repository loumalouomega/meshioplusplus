"""Provenance tag audit: doc/roadmap.md section 1's "audit and normalize" bullet.

Pins that the one-line provenance credit every writer emits
(``meshioplusplus._provenance.TAG`` on the Python side, ``detail::kProvenanceTag``
on the C++ side, see ``src/cpp/include/meshioplusplus/detail/provenance.hpp``) is
character-identical between the two engines, carries the release version, and
never regresses into the old per-writer drift this module replaced: a stale
``meshio`` name, a ``(C++ core)``-vs-``v{version}`` split, or a wall-clock
timestamp that made writing the same mesh twice produce different bytes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _provenance

_core = pytest.importorskip("meshioplusplus._core")

TRI = mp.Mesh(
    np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
    [("triangle", np.array([[0, 1, 2]]))],
)
TRI_2D = mp.Mesh(
    np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]),
    [("triangle", np.array([[0, 1, 2]]))],
)
TET = mp.Mesh(
    np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    [("tetra", np.array([[0, 1, 2, 3]]))],
)
HEX = mp.grid((1, 1, 1))


def _extract_tag(path: Path) -> bytes:
    """Pull the provenance tag out of a written file, whatever its format.

    Fixed-length extraction (not "scan to the next delimiter") is
    deliberate: the tag sits inside all sorts of surroundings across the
    ~25 formats that carry it -- an XML comment's closing ``-->``, a quoted
    ``"..."`` string, a parenthesised Ansys record, a NUL-padded binary
    slot -- and none of those characters may appear inside the tag itself,
    so reading exactly ``len(TAG)`` bytes from the match is the one rule
    that works uniformly everywhere.
    """
    data = path.read_bytes()
    idx = data.find(b"Written by meshio++")
    assert idx != -1, f"no provenance tag found in {path}"
    tag_len = len(_provenance.TAG.encode())
    return data[idx : idx + tag_len]


def _inspect_path(written_path: Path, name: str) -> Path:
    """The file to actually search for the tag in.

    Almost every format writes the tag into the file it was asked to write.
    EnSight is the one exception here: ``write(path.case, ...)`` writes the
    credit into a sibling ``.geo`` file, never into the ``.case`` file
    itself.
    """
    if name == "ensight":
        return written_path.with_suffix(".geo")
    return written_path


def _cases():
    """name -> (mesh, extension, cpp_write(path, mesh), py_write(path, mesh))."""
    from meshioplusplus.abaqus import _abaqus as py_abaqus
    from meshioplusplus.ansys import _ansys as py_ansys
    from meshioplusplus.avsucd import _avsucd as py_avsucd
    from meshioplusplus.ensight import _ensight as py_ensight
    from meshioplusplus.exodus import _exodus as py_exodus
    from meshioplusplus.flac3d import _flac3d as py_flac3d
    from meshioplusplus.flux import _flux as py_flux
    from meshioplusplus.mphtxt import _mphtxt as py_mphtxt
    from meshioplusplus.nastran import _nastran as py_nastran
    from meshioplusplus.netgen import _netgen as py_netgen
    from meshioplusplus.obj import _obj as py_obj
    from meshioplusplus.off import _off as py_off
    from meshioplusplus.permas import _permas as py_permas
    from meshioplusplus.ply import _ply as py_ply
    from meshioplusplus.stl import _stl as py_stl
    from meshioplusplus.tecplot import _tecplot as py_tecplot
    from meshioplusplus.tetgen import _tetgen as py_tetgen
    from meshioplusplus.triangle import _triangle as py_triangle
    from meshioplusplus.vti import _vti as py_vti
    from meshioplusplus.vtk import _vtk_42 as py_vtk42
    from meshioplusplus.vtk import _vtk_51 as py_vtk51
    from meshioplusplus.vtp import _vtp as py_vtp
    from meshioplusplus.vtu import _vtu as py_vtu

    return {
        "obj": (TRI, ".obj", _core.obj_write, py_obj.write),
        "off": (TRI, ".off", _core.off_write, py_off.write),
        "mphtxt": (TRI, ".mphtxt", _core.mphtxt_write, py_mphtxt.write),
        "avsucd": (TRI, ".avs", _core.avsucd_write, py_avsucd.write),
        "netgen": (
            TET,
            ".vol",
            lambda p, m: _core.netgen_write(p, m, ".16e"),
            py_netgen.write,
        ),
        "tetgen": (TET, ".node", _core.tetgen_write, py_tetgen.write),
        "triangle": (TRI_2D, ".node", _core.triangle_write, py_triangle.write),
        "abaqus": (TRI, ".inp", _core.abaqus_write, py_abaqus.write),
        "permas": (TRI, ".post", _core.permas_write, py_permas.write),
        "flac3d": (
            TET,
            ".f3grid",
            lambda p, m: _core.flac3d_write(p, m, "%.16e", False),
            lambda p, m: py_flac3d.write(p, m, binary=False),
        ),
        "flux": (TRI, ".pf3", _core.flux_write, py_flux.write),
        "tecplot": (TRI, ".dat", _core.tecplot_write, py_tecplot.write),
        "exodus": (TET, ".exo", _core.exodus_write, py_exodus.write),
        "ansys": (
            TRI,
            ".msh",
            lambda p, m: _core.ansys_write(p, m, False),
            lambda p, m: py_ansys.write(p, m, binary=False),
        ),
        "stl": (
            TRI,
            ".stl",
            lambda p, m: _core.stl_write(p, m, True, True),
            lambda p, m: py_stl.write(p, m, binary=True),
        ),
        "ensight": (
            TRI,
            ".case",
            lambda p, m: _core.ensight_write(p, m, False),
            lambda p, m: py_ensight.write(p, m, binary=False),
        ),
        "nastran": (TRI, ".bdf", _core.nastran_write, py_nastran.write),
        "ply": (
            TRI,
            ".ply",
            lambda p, m: _core.ply_write(p, m, False, True),
            lambda p, m: py_ply.write(p, m, binary=False),
        ),
        "vtk42": (
            TRI,
            ".vtk",
            lambda p, m: _core.vtk_write(p, m, False, False),
            py_vtk42.write,
        ),
        "vtk51": (
            TRI,
            ".vtk",
            lambda p, m: _core.vtk_write(p, m, False, True),
            py_vtk51.write,
        ),
        "vtu": (
            TRI,
            ".vtu",
            lambda p, m: _core.vtu_write(p, m, False, False),
            py_vtu.write,
        ),
        "vtp": (
            TRI,
            ".vtp",
            lambda p, m: _core.vtp_write(p, m, False, False),
            py_vtp.write,
        ),
        "vti": (
            HEX,
            ".vti",
            lambda p, m: _core.vti_write(p, m, False, False),
            py_vti.write,
        ),
    }


CASES = _cases()

# nastran's C++ reader is deliberately sentinel-gated (see doc/formats/nastran.md):
# only a file carrying the literal "meshioplusplus-cpp-nastran" comment -- which
# the Python writer never emits -- is accepted, so the two engines' *files*
# legitimately differ by that one extra line. The provenance tag itself is not
# an exception to anything: both engines still emit character-identical text.
NASTRAN_ONLY_CPP_HAS_SENTINEL = "nastran"


@pytest.mark.parametrize("name", sorted(CASES))
def test_engines_emit_character_identical_tag(tmp_path, name):
    mesh, ext, cpp_write, py_write = CASES[name]
    cpp_path = tmp_path / f"cpp_{name}{ext}"
    py_path = tmp_path / f"py_{name}{ext}"

    cpp_write(str(cpp_path), mesh)
    py_write(str(py_path), mesh)

    cpp_tag = _extract_tag(_inspect_path(cpp_path, name))
    py_tag = _extract_tag(_inspect_path(py_path, name))
    assert cpp_tag == py_tag == _provenance.TAG.encode()


def test_nastran_sentinel_is_cpp_only(tmp_path):
    """The one documented exception: the sentinel line, not the tag itself."""
    cpp_path = tmp_path / "cpp.bdf"
    py_path = tmp_path / "py.bdf"
    _core.nastran_write(str(cpp_path), TRI)
    from meshioplusplus.nastran import _nastran as py_nastran

    py_nastran.write(str(py_path), TRI)

    cpp_text = cpp_path.read_text()
    py_text = py_path.read_text()
    assert "meshioplusplus-cpp-nastran" in cpp_text
    assert "meshioplusplus-cpp-nastran" not in py_text
    assert _extract_tag(cpp_path) == _extract_tag(py_path) == _provenance.TAG.encode()


def test_openfoam_is_cpp_only(tmp_path):
    """OpenFOAM has no Python writer twin (documented in CLAUDE.md); the C++
    writer's banner still carries the tag. Unlike every other format, the
    credit cell is fixed-width (padded with trailing spaces up to the box's
    closing ``|``), so the match is a prefix check, not an exact one."""
    out = tmp_path / "case"
    _core.openfoam_write(str(out), HEX, {}, {})
    points_file = out / "constant" / "polyMesh" / "points"
    data = points_file.read_bytes()
    idx = data.find(b"Written by meshio++")
    assert idx != -1
    tag = _provenance.TAG.encode()
    assert data[idx : idx + len(tag)] == tag


def test_tag_contains_the_release_version():
    from meshioplusplus.__about__ import __version__

    assert _provenance.TAG == f"Written by meshio++ v{__version__}"
    assert _core.__version__ == __version__


@pytest.mark.parametrize("name", sorted(CASES))
def test_no_stale_drift_strings(tmp_path, name):
    """The regression guard for the drift this module fixed: no writer may
    reintroduce the old ``meshio`` (not ``meshio++``) name, the C++-only
    ``(C++ core)`` marker, or a version-less credit."""
    mesh, ext, cpp_write, py_write = CASES[name]
    cpp_path = tmp_path / f"cpp_{name}{ext}"
    py_path = tmp_path / f"py_{name}{ext}"
    cpp_write(str(cpp_path), mesh)
    py_write(str(py_path), mesh)

    for path in (_inspect_path(cpp_path, name), _inspect_path(py_path, name)):
        data = path.read_bytes()
        assert b"(C++ core)" not in data
        # "meshio v" without the "++" is the stale-fork spelling; a bare
        # "meshio++" with no following " v<digits>" is the version-less one.
        # The nastran sentinel legitimately contains "meshioplusplus" with no
        # "meshio v"/"meshio++ v" nearby, which is why it is not asserted on
        # beyond the two substring checks below.
        assert b"meshio v" not in data


def test_determinism_no_writer_embeds_a_timestamp(tmp_path):
    """The regression guard for the three deleted ``datetime.now()`` calls
    (obj, ply, exodus) plus flac3d's separate ``time.ctime()`` line: writing
    the same mesh twice must yield byte-identical output."""
    from meshioplusplus.exodus import _exodus as py_exodus
    from meshioplusplus.flac3d import _flac3d as py_flac3d
    from meshioplusplus.obj import _obj as py_obj
    from meshioplusplus.ply import _ply as py_ply

    checks = [
        ("obj", lambda p: py_obj.write(p, TRI), ".obj"),
        ("ply", lambda p: py_ply.write(p, TRI, binary=False), ".ply"),
        ("exodus", lambda p: py_exodus.write(p, TET), ".exo"),
        ("flac3d", lambda p: py_flac3d.write(p, TET, binary=False), ".f3grid"),
    ]
    for name, write, ext in checks:
        p1 = tmp_path / f"{name}_1{ext}"
        p2 = tmp_path / f"{name}_2{ext}"
        write(str(p1))
        write(str(p2))
        assert p1.read_bytes() == p2.read_bytes(), f"{name} is not deterministic"


# ---------------------------------------------------------------------------
# The opt-in record (v10.16.0): Mode/SlotTier semantics, the Python<->C++
# bridge, chain-sourcing, and conversion-assumption capture -- roadmap #1's
# bullets 2-7 (doc/provenance.md has the design; this pins the behaviour).
# ---------------------------------------------------------------------------


def test_off_mode_ignores_every_slot_tier():
    for tier in _provenance.SlotTier:
        result = _provenance.lines(tier)
        assert result == [_provenance.TAG]


def test_best_effort_renders_the_full_block_only_for_block_tier():
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        _provenance.set_source("in.vtu", "vtu")
        _provenance.add_operation("clean(weld=true)")
        _provenance.note("regions-dropped", "3 regions dropped")

        block = _provenance.lines(_provenance.SlotTier.BLOCK)
        assert block[0] == _provenance.TAG
        assert any("in.vtu" in line for line in block)
        assert any("clean(weld=true)" in line for line in block)
        assert any("3 regions dropped" in line for line in block)

        assert _provenance.lines(_provenance.SlotTier.SINGLE_LINE) == [_provenance.TAG]
        assert _provenance.lines(_provenance.SlotTier.BOUNDED) == [_provenance.TAG]
        assert s.get().notes == [
            _provenance.Note("regions-dropped", "3 regions dropped")
        ]


def test_required_throws_only_for_none_tier():
    with _provenance.scope(_provenance.Mode.REQUIRED):
        with pytest.raises(Exception):
            _provenance.lines(_provenance.SlotTier.NONE)
        # Degrading to the tag on a structurally smaller slot is not a failure.
        _provenance.lines(_provenance.SlotTier.SINGLE_LINE)
        _provenance.lines(_provenance.SlotTier.BOUNDED)
        _provenance.lines(_provenance.SlotTier.BLOCK)


def test_duplicate_notes_are_collapsed():
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        _provenance.note("dtype", "cast to int32")
        _provenance.note("dtype", "cast to int32")
        _provenance.note("dtype", "cast to int32")
        assert len(s.get().notes) == 1


def test_scopes_nest_and_restore():
    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("outer.vtu", "vtu")
        with _provenance.scope(_provenance.Mode.BEST_EFFORT):
            _provenance.set_source("inner.vtu", "vtu")
            assert _provenance.current_record().source_path == "inner.vtu"
        assert _provenance.current_record().source_path == "outer.vtu"
    assert _provenance.current_mode() is _provenance.Mode.OFF


def test_no_scope_means_off_and_note_is_a_no_op():
    assert _provenance.current_mode() is _provenance.Mode.OFF
    _provenance.note("x", "y")  # must not raise
    assert _provenance.current_record().notes == []


def test_timestamp_honours_source_date_epoch_and_the_off_switch(monkeypatch):
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "1000000000")
    assert _provenance.timestamp() == "2001-09-09T01:46:40Z"
    monkeypatch.delenv("SOURCE_DATE_EPOCH")

    monkeypatch.setenv("MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP", "off")
    assert _provenance.timestamp() == ""


def test_python_scope_drives_the_cpp_writer_via_the_bridge(tmp_path):
    """The bridge (bindings/python/_core.cpp's provenance_scope_push/_pop): a
    scope opened from Python must be honoured by a direct `_core.<fmt>_write`
    call, not just the pure-Python fallback writers."""
    out = tmp_path / "bridge.obj"
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        _provenance.set_source("in.vtu", "vtu")
        _provenance.note("regions-dropped", "Side regions have no OBJ equivalent")
        _core.obj_write(str(out), TRI)
        assert s.get().notes  # visible inside the scope too

    text = out.read_text()
    assert "in.vtu" in text
    assert "regions-dropped" in text

    # Scope closed: the C++ side's own thread-local state is back to Off.
    out2 = tmp_path / "bridge2.obj"
    _core.obj_write(str(out2), TRI)
    assert "in.vtu" not in out2.read_text()


def test_pipeline_records_source_target_and_operations(tmp_path):
    in_path = tmp_path / "in.obj"
    from meshioplusplus.obj import _obj as py_obj

    py_obj.write(str(in_path), TRI)
    out_path = tmp_path / "out.obj"

    settings = {
        "Version": 1,
        "Input": {"Path": str(in_path)},
        "Operations": [{"Op": "Clean", "Weld": True}],
        "Output": {"Path": str(out_path)},
    }
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        mp.run_pipeline(settings)
        rec = s.get()

    assert rec.source_path == str(in_path)
    assert rec.source_format == "obj"
    assert rec.target_format == "obj"
    assert len(rec.operations) == 1
    assert "Clean(" in rec.operations[0]
    assert "Weld=true" in rec.operations[0]


def test_cpp_and_python_pipeline_engines_agree_on_the_operation_chain(tmp_path):
    """Python's `_pipeline.py` is a separate, pure-Python engine
    (CLAUDE.md's own description) -- this pins that its rendering of a step
    still matches the C++ engine's `pipe_render_op`, for the common
    bool/int/string parameter shapes."""
    in_path = tmp_path / "in.obj"
    from meshioplusplus.obj import _obj as py_obj

    py_obj.write(str(in_path), TRI)

    settings = {
        "Version": 1,
        "Input": {"Path": str(in_path)},
        "Operations": [{"Op": "Clean", "Weld": True, "Atol": 1e-6}],
        "Output": {"Path": str(tmp_path / "py_out.obj")},
    }
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        mp.run_pipeline(settings)
        py_op = s.get().operations[0]

    cpp_out = tmp_path / "cpp_out.obj"
    from meshioplusplus._core import run_pipeline_file

    settings_path = tmp_path / "settings.json"
    import json

    settings["Output"]["Path"] = str(cpp_out)
    settings_path.write_text(json.dumps(settings))
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s2:
        run_pipeline_file(str(settings_path))
        cpp_op = s2.get().operations[0]

    assert py_op == cpp_op == "Clean(Atol=1e-06, Weld=true)"


def test_off_writer_records_dropped_cell_types(tmp_path):
    mesh = mp.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("triangle", np.array([[0, 1, 2]])), ("tetra", np.array([[0, 1, 2, 3]]))],
    )
    from meshioplusplus.off import _off as py_off

    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        py_off.write(str(tmp_path / "py.off"), mesh)
        py_notes = s.get().notes

    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s2:
        _core.off_write(str(tmp_path / "cpp.off"), mesh)
        cpp_notes = s2.get().notes

    assert any(n.category == "cells-dropped" and "tetra" in n.detail for n in py_notes)
    assert py_notes == cpp_notes


def test_warn_regions_dropped_records_a_note(tmp_path):
    """`detail::warn_regions_dropped` is the single choke point for the
    operations layer's own region-drop warning (9 call sites); pinned via
    `extract_surface`, which calls it whenever the input carries regions."""
    from meshioplusplus._surface import extract_surface

    mesh = mp.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
        point_sets={"corner": np.array([0])},
    )
    with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
        extract_surface(mesh)
        notes = s.get().notes
    assert any(n.category == "regions-dropped" for n in notes)


# ---------------------------------------------------------------------------
# Read-back (v10.17.0) -- roadmap #1's last bullet. See
# doc/provenance.md#reading-a-block-back.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label,text",
    [
        (
            "hash",
            "# Written by meshio++ v1.2.3\n# Converted from in.vtu (vtu)\nv 0 0 0\n",
        ),
        ("bang", "!PERMAS DataFile\n! Written by meshio++ v1.2.3\n$STRUCTURE\n"),
        ("xml-inline", '<?xml version="1.0"?>\n<!--Written by meshio++ v1.2.3-->\n'),
        ("xml-block", "<!--\nWritten by meshio++ v1.2.3\nTimestamp: X\n-->\n"),
        ("tecplot", 'TITLE = "Written by meshio++ v1.2.3"\nVARIABLES = "X"\n'),
        ("ansys", '(1 "Written by meshio++ v1.2.3")\n(2 3)\n'),
        ("openfoam", "|   \\\\  /    A nd  | Written by meshio++ v1.2.3        |\n"),
    ],
)
def test_scanner_handles_every_slot_shape(label, text):
    lines_found, recognised = _provenance.scan_provenance_text(text)
    assert recognised, label
    assert lines_found[0] == "Written by meshio++ v1.2.3", label


def test_scanner_keeps_parentheses_that_belong_to_the_content():
    """The bug this pins: an unconditional trailing-')' strip (needed for
    Ansys's `(1 "...")`) silently truncates every line legitimately ending in
    one."""
    lines_found, _ = _provenance.scan_provenance_text(
        "# Written by meshio++ v1.2.3\n"
        "# Converted from in.vtu (vtu)\n"
        "# Operation: Clean(Weld=true)\n"
    )
    assert lines_found[1] == "Converted from in.vtu (vtu)"
    assert lines_found[2] == "Operation: Clean(Weld=true)"


def test_scanner_is_honest_about_foreign_and_absent_blocks(tmp_path):
    assert _provenance.scan_provenance_text("v 0 0 0\n") == ([], False)
    assert _provenance.scan_provenance_text("# Created by SomeTool 3.2\n") == (
        [],
        False,
    )
    # A missing file is "nothing found", never an exception.
    assert _provenance.read_provenance_lines(tmp_path / "nope.obj") == ([], False)


def test_default_write_recovers_exactly_the_tag(tmp_path):
    out = tmp_path / "off.obj"
    mp.write(str(out), TRI)
    meta = mp.read_metadata(str(out))
    assert meta["provenance_recognised"] is True
    assert meta["provenance"] == [_provenance.TAG]


def test_scoped_write_round_trips_the_whole_block(tmp_path):
    out = tmp_path / "on.obj"
    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("in.msh", "gmsh")
        _provenance.note("regions-dropped", "Side regions dropped")
        mp.write(str(out), TRI)

    meta = mp.read_metadata(str(out))
    assert meta["provenance_recognised"] is True
    assert meta["provenance"][0] == _provenance.TAG
    assert any("in.msh" in line for line in meta["provenance"])
    assert any("regions-dropped" in line for line in meta["provenance"])


def test_a_read_block_is_never_re_emitted(tmp_path):
    """What makes "replace, never append" structural rather than a rule: a
    writer renders from the live record only, so converting N times leaves one
    block, not N."""
    first, second = tmp_path / "1.obj", tmp_path / "2.obj"
    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("original.msh", "gmsh")
        mp.write(str(first), TRI)

    mesh = mp.read(str(first))
    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("second.msh", "gmsh")
        mp.write(str(second), mesh)

    lines_found = mp.read_metadata(str(second))["provenance"]
    tags = [line for line in lines_found if line.startswith("Written by meshio++ v")]
    assert len(tags) == 1, "the block accumulated across a convert"
    assert not any("original.msh" in line for line in lines_found)


def test_cpp_and_python_scanners_agree(tmp_path):
    """`_provenance.scan_provenance_text` is a twin of the C++ scanner; a file
    scanned by either must yield the same lines."""
    out = tmp_path / "agree.obj"
    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("in.vtu", "vtu")
        _provenance.add_operation("Clean(Weld=true)")
        mp.write(str(out), TRI)

    cpp_lines, cpp_ok = _core.read_provenance_lines(str(out))
    py_lines, py_ok = _provenance.scan_provenance_text(
        out.read_bytes().replace(b"\0", b"\n").decode("utf-8", "replace")
    )
    assert list(cpp_lines) == py_lines
    assert bool(cpp_ok) == py_ok


# ---------------------------------------------------------------------------
# Cross-engine note parity (the guard for the assumption sweep). The character-
# identity guarantee means that wherever BOTH engines write a given mesh
# successfully, they must record the same conversion assumptions -- so widening
# `provenance_note` coverage one format at a time cannot silently desynchronise
# them.
# ---------------------------------------------------------------------------


# (format, _core writer, python writer, extension). Meshes are supplied per
# case so each can carry whatever the format's own lossy path needs.
def _note_parity_cases():
    from meshioplusplus.avsucd import _avsucd
    from meshioplusplus.cgns import _cgns
    from meshioplusplus.flac3d import _flac3d
    from meshioplusplus.mphtxt import _mphtxt
    from meshioplusplus.off import _off
    from meshioplusplus.ply import _ply
    from meshioplusplus.stl import _stl
    from meshioplusplus.tetgen import _tetgen
    from meshioplusplus.unv import _unv

    # A mesh mixing a surface block with a volume block: the shape that makes
    # surface-only formats drop something, and volume formats drop the rest.
    mixed = mp.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("triangle", np.array([[0, 1, 2]])), ("tetra", np.array([[0, 1, 2, 3]]))],
    )
    return [
        ("off", mixed, ".off", _core.off_write, _off.write),
        ("unv", mixed, ".unv", _core.unv_write, _unv.write),
        (
            "stl",
            mixed,
            ".stl",
            lambda p, m: _core.stl_write(p, m, False, True),
            lambda p, m: _stl.write(p, m, binary=False),
        ),
        (
            "ply",
            mixed,
            ".ply",
            lambda p, m: _core.ply_write(p, m, False, True),
            lambda p, m: _ply.write(p, m, binary=False),
        ),
        ("avsucd", mixed, ".avs", _core.avsucd_write, _avsucd.write),
        ("mphtxt", mixed, ".mphtxt", _core.mphtxt_write, _mphtxt.write),
        ("tetgen", mixed, ".node", _core.tetgen_write, _tetgen.write),
        (
            "flac3d",
            mixed,
            ".f3grid",
            lambda p, m: _core.flac3d_write(p, m, "%.16e", False),
            lambda p, m: _flac3d.write(p, m, binary=False),
        ),
        ("cgns", mixed, ".cgns", _core.cgns_write, _cgns.write),
    ]


@pytest.mark.parametrize("name", [c[0] for c in _note_parity_cases()])
def test_engines_record_the_same_notes(tmp_path, name):
    case = next(c for c in _note_parity_cases() if c[0] == name)
    _, mesh, ext, cpp_write, py_write = case

    def notes_from(write, path):
        with _provenance.scope(_provenance.Mode.BEST_EFFORT) as s:
            try:
                write(str(path), mesh)
            except Exception:
                return None  # this engine cannot write it; nothing to compare
            return [(n.category, n.detail) for n in s.get().notes]

    cpp_notes = notes_from(cpp_write, tmp_path / f"cpp{ext}")
    py_notes = notes_from(py_write, tmp_path / f"py{ext}")
    if cpp_notes is None or py_notes is None:
        pytest.skip(f"{name}: only one engine writes this mesh, nothing to compare")
    assert sorted(cpp_notes) == sorted(py_notes)
