"""Tests for the C++ transient-XDMF writer, `_core.XdmfTimeSeriesWriter`.

This is the *explicitly reachable* C++ writer, not a replacement for the pure
Python `meshioplusplus.xdmf.TimeSeriesWriter` -- exactly the `.mdpa` precedent:
the two are not byte-for-byte interchangeable (different call signature, a
sibling rather than CWD-relative `.h5`, sorted-name array order), so the C++
one is exposed additionally rather than swapped in underneath a tested public
API. `test_python_writer_is_untouched` below pins that.

The read-back oracle is deliberately the *Python* `TimeSeriesReader` plus the
ordinary `meshioplusplus.read`, not the writer's own C++ reader -- a format's
own writer is not a sufficient test oracle for that format's reader (the
exodus lesson in CLAUDE.md).
"""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core

HAS_CPP_SERIES = hasattr(_core, "XdmfTimeSeriesWriter")
HAS_HDF5 = getattr(_core, "__has_hdf5__", False)

pytestmark = pytest.mark.skipif(not HAS_CPP_SERIES, reason="built without the C++ core")

POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)
CELLS = [("tetra", np.array([[0, 1, 2, 3]]))]
NUM_STEPS = 3


def _step_mesh(k: int) -> meshioplusplus.Mesh:
    """One step's mesh: same geometry, k-dependent data."""
    return meshioplusplus.Mesh(
        POINTS,
        CELLS,
        point_data={"phi": np.full(len(POINTS), float(k)) + np.arange(len(POINTS))},
        cell_data={"a": [np.array([10.0 * k])]},
    )


def _write_series(path, data_format: str, steps: int = NUM_STEPS) -> None:
    writer = _core.XdmfTimeSeriesWriter(str(path), data_format)
    writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
    for k in range(steps):
        writer.write_data(0.5 * k, _step_mesh(k))
    writer.finalize()


def _expected_phi(k: int) -> np.ndarray:
    return np.full(len(POINTS), float(k)) + np.arange(len(POINTS))


# ---------------------------------------------------------------------------
# Round-trips
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "data_format",
    [
        pytest.param(
            "HDF",
            marks=pytest.mark.skipif(not HAS_HDF5, reason="no HDF5 in this build"),
        ),
        "XML",
        "Binary",
    ],
)
def test_three_step_series_reads_back(tmp_path, data_format):
    """Write 3 steps, read every one back with the Python TimeSeriesReader."""
    path = tmp_path / "series.xdmf"
    _write_series(path, data_format)
    assert path.is_file()

    with meshioplusplus.xdmf.TimeSeriesReader(path) as reader:
        points, cells = reader.read_points_cells()
        assert reader.num_steps == NUM_STEPS
        assert np.allclose(points, POINTS)
        assert cells[0].type == "tetra"
        assert np.array_equal(cells[0].data, CELLS[0][1])

        for k in range(NUM_STEPS):
            t, point_data, cell_data = reader.read_data(k)
            assert t == pytest.approx(0.5 * k)
            assert np.allclose(point_data["phi"], _expected_phi(k))
            assert np.allclose(np.concatenate(cell_data["a"]), [10.0 * k])


@pytest.mark.skipif(not HAS_HDF5, reason="no HDF5 in this build")
def test_hdf_companion_is_a_sibling_of_the_xdmf(tmp_path):
    """The documented C++/Python difference: `<path minus ext>.h5`, not CWD.

    The Python writer resolves `filename.stem + ".h5"` against the process's
    working directory, which breaks as soon as the output path has a directory
    component. This is one of the reasons the C++ writer is exposed separately
    instead of replacing it.
    """
    subdir = tmp_path / "run" / "out"
    subdir.mkdir(parents=True)
    path = subdir / "series.xdmf"
    _write_series(path, "HDF")

    assert (subdir / "series.h5").is_file()
    # ... and the .xdmf references it by bare basename, so the pair is movable.
    text = path.read_text()
    assert "series.h5:/" in text
    assert str(subdir) not in text


def test_ordinary_read_sees_the_series(tmp_path):
    """`meshioplusplus.read` resolves the temporal collection structurally."""
    path = tmp_path / "series.xdmf"
    _write_series(path, "XML")

    mesh = meshioplusplus.read(path)
    assert np.allclose(mesh.points, POINTS)
    assert mesh.cells[0].type == "tetra"
    # No time_step given => the first step.
    assert np.allclose(mesh.point_data["phi"], _expected_phi(0))

    meta = meshioplusplus.read_metadata(path)
    assert np.allclose(meta["time_values"], [0.0, 0.5, 1.0])


# ---------------------------------------------------------------------------
# The object itself
# ---------------------------------------------------------------------------


def test_context_manager_finalizes(tmp_path):
    """`with` is how the Python writer is used; the C++ one matches."""
    path = tmp_path / "ctx.xdmf"
    with _core.XdmfTimeSeriesWriter(str(path), "XML") as writer:
        writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
        writer.write_data(0.0, _step_mesh(0))
        writer.write_data(1.0, _step_mesh(1))
        # Inherent to the format: the collection element has to enclose every
        # step, so nothing is on disk until the block ends.
        assert not path.exists()
        assert writer.num_steps == 2
        assert writer.finalized is False

    assert path.is_file()
    with meshioplusplus.xdmf.TimeSeriesReader(path) as reader:
        assert reader.num_steps == 2
        reader.read_points_cells()
        assert np.allclose(reader.read_data(1)[1]["phi"], _expected_phi(1))


def test_enter_returns_the_writer(tmp_path):
    writer = _core.XdmfTimeSeriesWriter(str(tmp_path / "e.xdmf"), "XML")
    with writer as entered:
        assert entered is writer


def test_exit_does_not_swallow_the_body_exception(tmp_path):
    path = tmp_path / "boom.xdmf"
    with pytest.raises(ZeroDivisionError):
        with _core.XdmfTimeSeriesWriter(str(path), "XML") as writer:
            writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
            raise ZeroDivisionError("from the body")
    # It still finalized, as the Python writer's __exit__ does.
    assert path.is_file()


def test_finalize_is_idempotent(tmp_path):
    path = tmp_path / "idem.xdmf"
    writer = _core.XdmfTimeSeriesWriter(str(path), "XML")
    writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
    writer.write_data(0.0, _step_mesh(0))
    writer.finalize()
    assert writer.finalized is True
    writer.finalize()  # no throw, no second write
    assert writer.num_steps == 1


@pytest.mark.skipif(not HAS_HDF5, reason="no HDF5 in this build")
def test_default_data_format_is_hdf(tmp_path):
    """No `data_format` => "HDF", matching the C++ and Python defaults."""
    path = tmp_path / "default.xdmf"
    writer = _core.XdmfTimeSeriesWriter(str(path))
    writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
    writer.write_data(0.0, _step_mesh(0))
    writer.finalize()
    assert (tmp_path / "default.h5").is_file()
    assert 'Format="HDF"' in path.read_text()


# ---------------------------------------------------------------------------
# Errors: every one must be a clean Python exception, never a crash
# ---------------------------------------------------------------------------


def test_unknown_data_format_raises(tmp_path):
    with pytest.raises(meshioplusplus.WriteError, match="unknown data format"):
        _core.XdmfTimeSeriesWriter(str(tmp_path / "x.xdmf"), "Parquet")


def test_write_data_before_points_cells_raises(tmp_path):
    writer = _core.XdmfTimeSeriesWriter(str(tmp_path / "x.xdmf"), "XML")
    with pytest.raises(meshioplusplus.WriteError, match="WritePointsCells"):
        writer.write_data(0.0, _step_mesh(0))


def test_write_points_cells_twice_raises(tmp_path):
    writer = _core.XdmfTimeSeriesWriter(str(tmp_path / "x.xdmf"), "XML")
    writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
    with pytest.raises(meshioplusplus.WriteError, match="more than once"):
        writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))


def test_write_after_finalize_raises(tmp_path):
    writer = _core.XdmfTimeSeriesWriter(str(tmp_path / "x.xdmf"), "XML")
    writer.write_points_cells(meshioplusplus.Mesh(POINTS, CELLS))
    writer.finalize()
    with pytest.raises(meshioplusplus.WriteError, match="finalized"):
        writer.write_data(0.0, _step_mesh(0))


@pytest.mark.skipif(HAS_HDF5, reason="this build has HDF5")
def test_hdf_without_hdf5_raises_by_name(tmp_path):
    """Compiled out => a named WriteError from the constructor, not a crash."""
    with pytest.raises(meshioplusplus.WriteError, match="MESHIOPLUSPLUS_WITH_HDF5"):
        _core.XdmfTimeSeriesWriter(str(tmp_path / "x.xdmf"), "HDF")


# ---------------------------------------------------------------------------
# The Python writer is NOT a shim over this one
# ---------------------------------------------------------------------------


def test_python_writer_is_untouched(tmp_path, monkeypatch):
    """`meshioplusplus.xdmf.TimeSeriesWriter` keeps its own API and behaviour.

    Guard against someone "unifying" the two: the Python writer takes raw
    `points, cells` / `point_data=, cell_data=` and writes its `.h5` relative
    to the CWD. Both are documented and tested elsewhere; changing them under
    users is the failure mode this test exists to catch.
    """
    assert meshioplusplus.xdmf.TimeSeriesWriter is not _core.XdmfTimeSeriesWriter
    assert meshioplusplus.xdmf.TimeSeriesWriter.__module__.startswith("meshioplusplus")

    monkeypatch.chdir(tmp_path)
    with meshioplusplus.xdmf.TimeSeriesWriter("py.xdmf", data_format="XML") as writer:
        # The raw-array signature, unchanged.
        writer.write_points_cells(POINTS, CELLS)
        writer.write_data(0.0, point_data={"phi": _expected_phi(0)})

    with meshioplusplus.xdmf.TimeSeriesReader(tmp_path / "py.xdmf") as reader:
        reader.read_points_cells()
        assert reader.num_steps == 1
        assert np.allclose(reader.read_data(0)[1]["phi"], _expected_phi(0))
