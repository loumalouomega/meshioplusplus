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
