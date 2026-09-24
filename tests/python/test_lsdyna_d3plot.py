"""LS-DYNA d3plot state databases: both engines, lasso-python's real test
families and the shell + solid families ``tools/gen_d3plot_reference.py`` writes
with lasso-python, every state checked against lasso-python's own reading
(frozen in ``lasso_reference.npz``)."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.lsdyna_d3plot import _d3plot as py_d3plot

MESHES = pathlib.Path(__file__).parent / "meshes" / "lsdyna_d3plot"
FAMILIES = sorted(MESHES.glob("*/*/d3plot"))
REFERENCE = MESHES / "lasso_reference.npz"
SHELL_SOLID = MESHES / "generated" / "shell_solid" / "d3plot"


def _id(path):
    return f"{path.parent.parent.name}/{path.parent.name}"


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.lsdyna_d3plot
    return py_d3plot


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for mine, theirs in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(mine) == sorted(theirs)
        for name in mine:
            assert np.asarray(mine[name]).dtype == np.asarray(theirs[name]).dtype
            np.testing.assert_array_equal(mine[name], theirs[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            assert x.dtype == y.dtype
            np.testing.assert_array_equal(x, y)
    assert sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in a.regions
    ) == sorted((r.kind, r.name, r.dim, r.tag, tuple(r.entries)) for r in b.regions)


@pytest.mark.parametrize("path", FAMILIES, ids=[_id(p) for p in FAMILIES])
@pytest.mark.parametrize("step", [0, 1, -1])
def test_engines_agree_on_every_family(path, step):
    n = len(py_d3plot.time_values(path))
    if step >= n:
        pytest.skip("fewer states")
    _same(
        meshioplusplus.lsdyna_d3plot.read(path, time_step=step),
        py_d3plot.read(path, time_step=step),
    )


@pytest.mark.parametrize("path", FAMILIES, ids=[_id(p) for p in FAMILIES])
def test_every_state_matches_lasso(engine, path):
    """The done-when: displacement, plastic strain, stress and deletion match
    lasso-python on every state (the generated shell + solid families delete
    shells and a solid along the way)."""
    ref = np.load(REFERENCE)
    key = _id(path)
    times = engine.time_values(path)
    np.testing.assert_allclose(times, ref[f"{key}:times"], rtol=1e-7)
    coords = ref[f"{key}:coords"].astype(np.float64)
    for step in range(len(times)):
        mesh = engine.read(path, time_step=step)
        np.testing.assert_array_equal(mesh.points, coords)
        if f"{key}:disp" in ref:
            np.testing.assert_array_equal(
                mesh.point_data["displacement"],
                ref[f"{key}:disp"][step].astype(np.float64) - coords,
            )
        for family, cell_types in (
            ("solid", {"hexahedron", "wedge", "tetra", "pyramid", "tetra10"}),
            ("shell", {"quad", "triangle", "quad8"}),
        ):
            if f"{key}:{family}_ids" not in ref or len(ref[f"{key}:{family}_ids"]) == 0:
                continue
            by_eid = {}
            for b, cells in enumerate(mesh.cells):
                if cells.type not in cell_types:
                    continue
                for i, e in enumerate(mesh.cell_data["lsdyna:eid"][b].tolist()):
                    by_eid[e] = (b, i)
            order = [by_eid[int(e)] for e in ref[f"{key}:{family}_ids"]]
            if f"{key}:{family}_eps" in ref:
                expected = ref[f"{key}:{family}_eps"][step]
                expected = expected.reshape(len(order), -1)
                got = np.array(
                    [mesh.cell_data["effective_plastic_strain"][b][i] for b, i in order]
                ).reshape(len(order), -1)
                np.testing.assert_array_equal(got[:, : expected.shape[1]], expected)
            if f"{key}:{family}_stress" in ref:
                expected = ref[f"{key}:{family}_stress"][step]
                n, p, w = expected.shape
                layout = mesh.field_data.get(
                    "lsdyna_d3plot:layout:stress", np.array([1, 6])
                )
                got = np.array([mesh.cell_data["stress"][b][i] for b, i in order])
                got = got.reshape(n, layout[0], layout[1])[:, :p, :w]
                np.testing.assert_array_equal(got, expected)
            if f"{key}:{family}_alive" in ref:
                expected = (ref[f"{key}:{family}_alive"][step] != 0).astype(np.int8)
                got = np.array([mesh.cell_data["lsdyna:alive"][b][i] for b, i in order])
                np.testing.assert_array_equal(got, expected)


def test_shell_solid_family_deletes_and_yields(engine):
    first = engine.read(SHELL_SOLID, time_step=0)
    last = engine.read(SHELL_SOLID, time_step=-1)
    assert [c.type for c in first.cells] == [
        "hexahedron",
        "wedge",
        "tetra",
        "line",
        "quad",
        "triangle",
    ]
    alive = np.concatenate(last.cell_data["lsdyna:alive"])
    assert alive.dtype == np.int8
    assert alive.sum() == len(alive) - 4  # three shells and the tetrahedron
    eps = np.concatenate([a[:, 0] for a in last.cell_data["effective_plastic_strain"]])
    assert np.nanmax(eps) > 0.0
    assert not np.any(np.concatenate(first.cell_data["lsdyna:alive"]) == 0)
    names = {r.name: r.tag for r in first.regions}
    assert names == {"block": 100, "plate": 200, "struts": 300}
    np.testing.assert_array_equal(
        last.field_data["lsdyna_d3plot:layout:stress"], [3, 6]
    )
    assert last.field_data["meshio:time"][0] == pytest.approx(3.0e-3)


def test_precision_and_splitting_do_not_change_the_reading(engine):
    single = engine.read(SHELL_SOLID, time_step=2)
    split = engine.read(
        MESHES / "generated" / "shell_solid_split" / "d3plot", time_step=2
    )
    double = engine.read(
        MESHES / "generated" / "shell_solid_double" / "d3plot", time_step=2
    )
    _same(single, split)
    # double precision carries more digits of the same (float32-born) values
    np.testing.assert_allclose(
        double.point_data["displacement"],
        single.point_data["displacement"],
        rtol=1e-6,
        atol=1e-7,
    )


def test_sequence_metadata_and_family_glob(tmp_path):
    folder = MESHES / "generated" / "shell_solid_split"
    times = meshioplusplus.read_metadata(folder / "d3plot")["time_values"]
    assert len(times) == 4
    assert _core.read_metadata(str(folder / "d3plot"))["time_values"] == times
    entries = meshioplusplus.sequence_entries(str(folder / "d3plot*"))
    assert [e["step"] for e in entries] == [0, 1, 2, 3]
    assert {pathlib.Path(e["path"]).name for e in entries} == {"d3plot"}


def test_selective_read(engine):
    mesh = engine.read(SHELL_SOLID, arrays=["displacement"], time_step=1)
    assert "displacement" in mesh.point_data and "stress" not in mesh.cell_data
    bare = engine.read(SHELL_SOLID, points_only=True)
    assert "displacement" not in bare.point_data
    assert "meshio:time" in bare.field_data


def test_found_by_name_and_by_content(tmp_path):
    assert meshioplusplus.read(SHELL_SOLID).cells
    path = tmp_path / "crash.bin"
    shutil.copy(SHELL_SOLID, path)
    assert meshioplusplus.sniff_format(path) == "lsdyna_d3plot"
    assert _core.sniff_format(str(path)) == "lsdyna_d3plot"
    big = MESHES / "generated" / "shell_solid_double" / "d3plot"
    assert _core.sniff_format(str(big)) == "lsdyna_d3plot"


def test_a_family_member_is_refused_by_name(engine):
    member = MESHES / "generated" / "shell_solid_split" / "d3plot02"
    with pytest.raises(meshioplusplus.ReadError, match="open the base file"):
        engine.read(member)
    with pytest.raises(meshioplusplus.ReadError, match="open the base file"):
        meshioplusplus.read(member)


def test_truncated_family_is_refused(engine, tmp_path):
    data = SHELL_SOLID.read_bytes()
    path = tmp_path / "d3plot"
    path.write_bytes(data[:300])
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(path)


def test_refused_features_are_named(engine, tmp_path):
    data = bytearray(SHELL_SOLID.read_bytes())
    path = tmp_path / "d3plot"
    femzip = bytearray(data)
    femzip[51 * 4 : 52 * 4] = (76_893_465).to_bytes(4, "little")
    path.write_bytes(bytes(femzip))
    with pytest.raises(meshioplusplus.ReadError, match="femzip"):
        engine.read(path)
    d3part = bytearray(data)
    d3part[11 * 4 : 12 * 4] = (5).to_bytes(4, "little")
    path.write_bytes(bytes(d3part))
    with pytest.raises(meshioplusplus.ReadError, match="d3part"):
        engine.read(path)
