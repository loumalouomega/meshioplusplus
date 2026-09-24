"""LS-DYNA binout (LSDA) files: both engines on a binout LS-DYNA wrote
(``glstat``) and on a ``nodout`` subset of another, every output checked against
lasso-python's ``Binout`` reading (frozen in ``binout_reference.npz`` by
``tools/gen_binout_reference.py``)."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.lsdyna_binout import _binout as py_binout

MESHES = pathlib.Path(__file__).parent / "meshes" / "lsdyna_binout"
GLSTAT = MESHES / "binout_glstat"
NODOUT = MESHES / "nodout" / "binout"
REFERENCE = MESHES / "binout_reference.npz"


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.lsdyna_binout
    return py_binout


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for mine, theirs in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(mine) == sorted(theirs)
        for name in mine:
            np.testing.assert_array_equal(mine[name], theirs[name])


@pytest.mark.parametrize("path", [GLSTAT, NODOUT], ids=["glstat", "nodout"])
@pytest.mark.parametrize("step", [0, 3, -1])
def test_engines_agree(path, step):
    _same(
        meshioplusplus.lsdyna_binout.read(path, time_step=step),
        py_binout.read(path, time_step=step),
    )


def test_nodout_steps_match_lasso(engine):
    """Each nodout output is a point cloud of its nodes at their coordinates,
    with displacement, rotation, velocity and acceleration as point data."""
    ref = np.load(REFERENCE)
    times = engine.time_values(NODOUT)
    np.testing.assert_array_equal(np.float32(times), ref["nodout|nodout|time"])
    for step in (0, 17, len(times) - 1):
        mesh = engine.read(NODOUT, time_step=step)
        assert [c.type for c in mesh.cells] == ["vertex"]
        np.testing.assert_array_equal(
            mesh.point_data["lsdyna:nid"], ref["nodout|nodout|ids"]
        )
        for axis, c in (("x", 0), ("y", 1), ("z", 2)):
            np.testing.assert_array_equal(
                np.float32(mesh.points[:, c]),
                ref[f"nodout|nodout|{axis}_coordinate"][step],
            )
            for name, var in (
                ("displacement", f"{axis}_displacement"),
                ("rotation", f"r{axis}_displacement"),
                ("velocity", f"{axis}_velocity"),
                ("rotational_velocity", f"r{axis}_velocity"),
                ("acceleration", f"{axis}_acceleration"),
                ("rotational_acceleration", f"r{axis}_acceleration"),
            ):
                np.testing.assert_array_equal(
                    np.float32(mesh.point_data[name][:, c]),
                    ref[f"nodout|nodout|{var}"][step],
                )
        assert mesh.field_data["meshio:time"][0] == pytest.approx(times[step])


def test_other_databases_attach_their_latest_output(engine):
    """glstat and matsum write every 1.0 s, nodout every 0.0979 s: a step
    carries each database's latest output at or before its time, with that
    output's own time; elout/beam writes with nodout."""
    ref = np.load(REFERENCE)
    glstat_times = ref["nodout|glstat|time"]
    for step in (0, 11, 59):
        mesh = engine.read(NODOUT, time_step=step)
        t = np.float32(mesh.field_data["meshio:time"][0])
        k = int(np.searchsorted(glstat_times, t, side="right")) - 1
        assert mesh.field_data["binout:glstat:time"][0] == glstat_times[k]
        assert (
            mesh.field_data["binout:glstat:kinetic_energy"][0]
            == ref["nodout|glstat|kinetic_energy"][k]
        )
        np.testing.assert_array_equal(
            np.float32(mesh.field_data["binout:matsum:internal_energy"]),
            ref["nodout|matsum|internal_energy"][k],
        )
        np.testing.assert_array_equal(
            mesh.field_data["binout:matsum:ids"], ref["nodout|matsum|ids"]
        )
        assert (
            mesh.field_data["binout:elout/beam:axial"][0]
            == ref["nodout|elout/beam|axial"][step]
        )


def test_without_nodout_the_steps_are_the_first_database(engine):
    """A binout LS-DYNA wrote with glstat and rwforc only: no points, the
    glstat outputs as steps."""
    ref = np.load(REFERENCE)
    times = engine.time_values(GLSTAT)
    np.testing.assert_array_equal(np.float32(times), ref["glstat|glstat|time"])
    mesh = engine.read(GLSTAT, time_step=4)
    assert len(mesh.points) == 0
    for var in ("kinetic_energy", "internal_energy", "total_energy", "time_step"):
        assert (
            mesh.field_data[f"binout:glstat:{var}"][0] == ref[f"glstat|glstat|{var}"][4]
        )
    assert mesh.field_data["binout:glstat:cycle"].dtype == np.int64


def test_selective_read(engine):
    mesh = engine.read(NODOUT, arrays=["displacement", "binout:glstat:kinetic_energy"])
    assert set(mesh.point_data) == {"lsdyna:nid", "displacement"}
    assert "binout:glstat:kinetic_energy" in mesh.field_data
    assert "binout:matsum:internal_energy" not in mesh.field_data
    bare = engine.read(NODOUT, points_only=True)
    assert set(bare.point_data) == {"lsdyna:nid"} and len(bare.points) == 32


def test_found_by_name_and_by_content(tmp_path):
    assert meshioplusplus.sniff_format(NODOUT) == "lsdyna_binout"
    renamed = tmp_path / "results.lsda"
    shutil.copy(GLSTAT, renamed)
    assert meshioplusplus.sniff_format(renamed) == "lsdyna_binout"
    assert _core.sniff_format(str(renamed)) == "lsdyna_binout"
    mesh = meshioplusplus.read(NODOUT, time_step=2)
    assert len(mesh.points) == 32
    assert meshioplusplus.read_metadata(NODOUT)["time_values"] == pytest.approx(
        meshioplusplus.lsdyna_binout.time_values(NODOUT)
    )


def test_refusals(engine, tmp_path):
    junk = tmp_path / "binout"
    junk.write_bytes(b"\x00" * 64)
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(junk)
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(NODOUT, time_step=60)
    truncated = tmp_path / "cut"
    truncated.write_bytes(NODOUT.read_bytes()[:5000])
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(truncated)
