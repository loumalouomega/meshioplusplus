"""DOLFINx VTX output (ADIOS2 ``.bp`` directories): real DOLFINx 0.11 runs
checked against ParaView 6.1's ``ADIOS2VTXReader`` (and, for the steps ParaView
cannot read, against the raw ADIOS2 arrays); see ``meshes/vtx/README.md``.

The engine is the core when it was built with ADIOS2, else the ``adios2``
package: one process never holds both, since the package brings its own
ADIOS2 libraries."""

import pathlib
import shutil
import sys

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError
from meshioplusplus.vtx import _vtx

HAS_CORE = getattr(_core, "__has_adios2__", False)
if not HAS_CORE:
    pytest.importorskip("adios2")

VTX = pathlib.Path(__file__).parent / "meshes" / "vtx"
NAMES = ["heat", "elastic", "quads", "heat_np2", "complex"]

# VTK_LAGRANGE_* (type, nodes) -> meshio++ type, as the reader lowers them.
LOWERED = {
    (69, 3): "triangle",
    (70, 4): "quad",
    (71, 10): "tetra10",
    (5, 3): "triangle",
}


def _reference(name):
    return np.load(VTX / f"{name}.reference.npz")


def _raw(name):
    return np.load(VTX / f"{name}.raw.npz")


def _connectivity(mesh):
    return np.concatenate([c.data.ravel() for c in mesh.cells])


@pytest.mark.parametrize("name", NAMES)
def test_matches_paraview(name):
    ref = _reference(name)
    assert meshioplusplus.vtx.time_values(VTX / f"{name}.bp") == list(ref["times"])
    for k in ref["steps"]:
        mesh = meshioplusplus.read(VTX / f"{name}.bp", time_step=int(k))
        np.testing.assert_array_equal(mesh.points, ref[f"k{k}/points"])
        np.testing.assert_array_equal(_connectivity(mesh), ref[f"k{k}/connectivity"])
        sizes = np.diff(ref[f"k{k}/offsets"])
        types = [LOWERED[(t, n)] for t, n in zip(ref[f"k{k}/types"], sizes)]
        assert [c.type for c in mesh.cells] == sorted(set(types), key=types.index)
        assert mesh.field_data["meshio:time"][0] == ref["times"][k]
        for key in ref.files:
            if not key.startswith(f"k{k}/point/") and not key.startswith(f"k{k}/cell/"):
                continue
            kind, array = key.split("/")[1:]
            expected = ref[key]
            if kind == "point":
                got = mesh.point_data[array]
                # ParaView appends DOLFINx's stray one-value ghost blocks.
                expected = expected[: len(got)]
            else:
                got = np.concatenate(mesh.cell_data[array])
            np.testing.assert_array_equal(got, expected)


@pytest.mark.parametrize("name", ["heat", "heat_np2", "elastic", "complex"])
def test_every_step_matches_raw_arrays(name):
    """Every step, the ones ParaView aborts on included: the mesh carried from
    step 0 (`reuse`), the functions as ADIOS2 stored them."""
    raw = _raw(name)
    first = meshioplusplus.read(VTX / f"{name}.bp")
    steps = sorted({int(key.split("/")[0][1:]) for key in raw.files})
    for k in steps:
        mesh = meshioplusplus.read(VTX / f"{name}.bp", time_step=k)
        np.testing.assert_array_equal(mesh.points, first.points)
        np.testing.assert_array_equal(_connectivity(mesh), _connectivity(first))
        for key in raw.files:
            if not key.startswith(f"k{k}/"):
                continue
            array = key.split("/")[1]
            got = (
                mesh.point_data[array]
                if array in mesh.point_data
                else np.concatenate(mesh.cell_data[array])
            )
            np.testing.assert_array_equal(got, raw[key])


def test_ghost_arrays_skip_dolfinx_stray_blocks():
    """The first step of a `reuse` file holds malformed one-value ghost blocks;
    the arrays match the rank blocks and every step's are the same."""
    ids = [
        meshioplusplus.read(VTX / "heat_np2.bp", time_step=k).point_data[
            "vtkOriginalPointIds"
        ]
        for k in range(3)
    ]
    assert len(ids[0]) == 38
    for other in ids[1:]:
        np.testing.assert_array_equal(ids[0], other)
    assert set(ids[0]) == set(range(25))


def _by_position(mesh):
    return {tuple(np.round(p, 9)): i for i, p in enumerate(mesh.points)}


@pytest.mark.parametrize("step", [0, 2])
def test_drop_ghosts_gives_the_serial_mesh(step):
    parallel = meshioplusplus.vtx.read(
        VTX / "heat_np2.bp", time_step=step, ghosts="drop"
    )
    serial = meshioplusplus.read(VTX / "heat.bp", time_step=step)
    assert len(parallel.points) == len(serial.points)
    assert "vtkGhostType" not in parallel.point_data
    index = _by_position(parallel)
    order = np.array([index[tuple(np.round(p, 9))] for p in serial.points])
    np.testing.assert_allclose(
        parallel.point_data["u"][order], serial.point_data["u"], atol=1e-14
    )

    def cells(mesh):
        return {
            tuple(sorted(tuple(np.round(mesh.points[i], 9)) for i in c))
            for c in mesh.cells[0].data
        }

    assert cells(parallel) == cells(serial)


def test_points_only_and_array_filter():
    mesh = meshioplusplus.vtx.read(VTX / "elastic.bp", points_only=True)
    assert not mesh.point_data and not mesh.cell_data
    mesh = meshioplusplus.vtx.read(VTX / "elastic.bp", arrays=["density"])
    assert list(mesh.cell_data) == ["density"] and not mesh.point_data


def test_negative_and_out_of_range_steps():
    last = meshioplusplus.read(VTX / "heat.bp", time_step=-1)
    assert last.field_data["meshio:time"][0] == pytest.approx(0.02)
    with pytest.raises(ReadError, match="3"):
        meshioplusplus.read(VTX / "heat.bp", time_step=3)


def test_metadata_and_sequence():
    meta = meshioplusplus.read_metadata(VTX / "heat.bp")
    assert meta["time_values"] == pytest.approx([0.0, 0.01, 0.02])
    steps = list(meshioplusplus.read_sequence(VTX / "heat.bp"))
    assert [t for t, _ in steps] == pytest.approx([0.0, 0.01, 0.02])


def test_glob_keeps_bp_directories(tmp_path):
    for k in range(2):
        shutil.copytree(VTX / "quads.bp", tmp_path / f"run_{k}.bp")
    (tmp_path / "ignored").mkdir()
    entries = meshioplusplus.sequence_entries(str(tmp_path / "run_*.bp"))
    assert [pathlib.Path(e["path"]).name for e in entries] == ["run_0.bp", "run_1.bp"]


def test_sniffed_as_a_directory(tmp_path):
    target = tmp_path / "noext"
    shutil.copytree(VTX / "heat.bp", target)
    from meshioplusplus._sniff import sniff_format

    assert sniff_format(target) == "vtx"


def test_refuses_a_bp_without_schema():
    with pytest.raises(ReadError, match="adios4dolfinx"):
        meshioplusplus.read(VTX / "plain.bp")


def test_python_twin_names_its_dependency(monkeypatch):
    monkeypatch.setitem(sys.modules, "adios2", None)
    with pytest.raises(ImportError, match=r"meshioplusplus\[adios2\]"):
        _vtx.read(VTX / "heat.bp")
