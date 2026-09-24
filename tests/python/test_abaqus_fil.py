"""Abaqus results files (``.fil``): both engines, the ASCII and binary fixtures
written by ``tools/gen_abaqus_fil_fixtures.py`` (every value checked against its
closed form) and real Abaqus 2023 output from pybaqus (values checked against
pybaqus's own reading)."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.abaqus_fil import _abaqus_fil as py_fil

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "abaqus_fil"
MODELS = ["model.fil", "model_le.fil", "model_be.fil"]
FIXTURES = sorted(MESHES.glob("**/*.fil"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    if request.param == "core":
        return meshioplusplus.abaqus_fil
    return py_fil


def _same(a, b, rtol=0.0):
    """Equal meshes; ``rtol`` for ASCII (16 significant digits) against binary."""

    def eq(x, y):
        np.testing.assert_allclose(x, y, rtol=rtol, atol=0.0)

    eq(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for mine, theirs in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(mine) == sorted(theirs)
        for name in mine:
            eq(mine[name], theirs[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            eq(x, y)
    assert sorted(
        (r.kind, r.name, r.dim, tuple(r.entries)) for r in a.regions
    ) == sorted((r.kind, r.name, r.dim, tuple(r.entries)) for r in b.regions)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
@pytest.mark.parametrize("step", [0, -1])
def test_engines_agree_on_every_fixture(path, step):
    _same(
        meshioplusplus.abaqus_fil.read(path, time_step=step),
        py_fil.read(path, time_step=step),
    )


@pytest.mark.parametrize("name", MODELS[1:])
def test_binary_twins_read_as_the_ascii_file(engine, name):
    for step in (0, 1):
        _same(
            engine.read(MESHES / name, time_step=step),
            engine.read(MESHES / "model.fil", time_step=step),
            rtol=1e-15,
        )


# The closed forms tools/gen_abaqus_fil_fixtures.py writes.
def _u(node, comp, inc):
    return 0.001 * node + 0.1 * comp + inc


def _s(elem, point, sp, comp, inc):
    return 100.0 * elem + 10.0 * point + comp + 0.25 * sp + inc


@pytest.mark.parametrize("name", MODELS)
def test_model_and_increment_values(engine, name):
    mesh = engine.read(MESHES / name, time_step=1)
    assert [c.type for c in mesh.cells] == ["hexahedron20", "hexahedron", "quad"]
    for block in mesh.cells[:2]:
        for row in np.asarray(block.data):
            assert _is_valid(np.asarray(mesh.points)[row], block.type)
    assert [int(a[0]) for a in mesh.cell_data["abaqus:id"]] == [1, 2, 3]
    labels = mesh.point_data["abaqus:id"]
    assert float(mesh.field_data["meshio:time"]) == 1.0
    assert int(mesh.field_data["abaqus:increment"]) == 2
    assert int(mesh.field_data["abaqus:step"]) == 1
    u = mesh.point_data["U"]
    np.testing.assert_array_equal(u, [[_u(n, c, 2) for c in range(3)] for n in labels])
    rf = mesh.point_data["RF"]
    clamped = labels <= 9
    np.testing.assert_array_equal(rf[clamped][:, 0], -(labels[clamped]) * 2.0)
    assert np.isnan(rf[~clamped]).all()
    # Per-point data is one (cells, points * components) array in every block,
    # point-major, with its (points, components) layout in field_data.
    s = mesh.cell_data["S"]
    assert [a.shape for a in s] == [(1, 48)] * 3
    np.testing.assert_array_equal(mesh.field_data["abaqus:layout:S"], [8, 6])
    np.testing.assert_array_equal(
        s[0][0].reshape(8, 6),
        [[_s(1, p, 0, c, 2) for c in range(6)] for p in range(1, 9)],
    )
    np.testing.assert_array_equal(s[1][0, :6], [_s(2, 1, 0, c, 2) for c in range(6)])
    assert np.isnan(s[1][0, 6:]).all()
    np.testing.assert_array_equal(s[2][0, :3], [_s(3, 1, 1, c, 2) for c in range(3)])
    assert np.isnan(s[2][0, 3:]).all()
    sp5 = mesh.cell_data["S@sp5"]
    assert [a.shape for a in sp5] == [(1, 3)] * 3
    np.testing.assert_array_equal(sp5[2][0], [_s(3, 1, 5, c, 2) for c in range(3)])
    assert np.isnan(sp5[0]).all()
    sinv = mesh.cell_data["SINV"]
    assert sinv[1].shape == (1, 7)
    np.testing.assert_array_equal(sinv[1][0], [7.0 * c + 2 for c in range(7)])
    assert "abaqus:layout:SINV" not in mesh.field_data
    nforc = mesh.cell_data["NFORC"][1]
    assert nforc.shape == (1, 60)  # the widest cell, the hexahedron20, has 20 nodes
    np.testing.assert_array_equal(mesh.field_data["abaqus:layout:NFORC"], [20, 3])
    hex8 = np.asarray(mesh.cells[1].data)[0]
    per_node = nforc[0].reshape(20, 3)
    np.testing.assert_array_equal(per_node[:8, 0], labels[hex8] * 0.5 + 20.0)
    assert np.isnan(per_node[8:]).all()
    peeq = mesh.point_data["PEEQ"]
    np.testing.assert_allclose(peeq[:2], [0.02, 0.04])
    assert np.isnan(peeq[2:]).all()
    regions = {(r.kind, r.name): r for r in mesh.regions}
    clamped_set = regions[("point", "ASSEMBLY_PART-1-1_CLAMPED_NODES")]
    np.testing.assert_array_equal(sorted(labels[clamped_set.entries]), range(1, 10))
    assert regions[("cell", "SOLIDS")].dim == 3
    np.testing.assert_array_equal(regions[("cell", "SOLIDS")].entries, [0, 1])


def test_increments_are_steps(engine):
    assert list(engine.time_values(MESHES / "model_be.fil")) == [0.5, 1.0]
    first = engine.read(MESHES / "model.fil")
    assert float(first.field_data["meshio:time"]) == 0.5
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(MESHES / "model.fil", time_step=2)


def test_sequence_and_metadata():
    assert meshioplusplus.read_metadata(MESHES / "model_le.fil")["time_values"] == [
        0.5,
        1.0,
    ]
    assert _core.read_metadata(str(MESHES / "model.fil"))["time_values"] == [0.5, 1.0]


def test_selective_read(engine):
    mesh = engine.read(MESHES / "model.fil", arrays=["U"])
    assert (
        "U" in mesh.point_data
        and "S" not in mesh.cell_data
        and "RF" not in mesh.point_data
    )
    bare = engine.read(MESHES / "model.fil", points_only=True)
    assert "U" not in bare.point_data


# pybaqus's reading of its own hex_C3D8.fil (Abaqus 2023): node 2's U and the
# first integration point's S.
def test_real_abaqus_output_matches_pybaqus(engine):
    mesh = engine.read(MESHES / "pybaqus" / "hex_C3D8.fil")
    assert [c.type for c in mesh.cells] == ["hexahedron"]
    row = list(mesh.point_data["abaqus:id"]).index(2)
    np.testing.assert_array_equal(
        mesh.point_data["U"][row],
        [0.005484804966181764, 0.01164481342587608, 2.904946755494933e-33],
    )
    np.testing.assert_array_equal(
        mesh.cell_data["S"][0][0, :6],  # the first integration point
        [
            -1.781822547468652,
            6.695266022198746,
            3.419889858603343,
            23.52460259453869,
            3.390710085233756,
            52.63709925322325,
        ],
    )
    names = {r.name for r in mesh.regions}
    assert "ASSEMBLY_TEST_INSTANCE_SET-TEST_PART" in names
    gaps = engine.read(MESHES / "pybaqus" / "discontinuous_numbering_2D.fil")
    assert gaps.cell_data["S"][0].shape == (2, 12)
    np.testing.assert_array_equal(gaps.field_data["abaqus:layout:S"], [4, 3])


def test_truncated_binary_is_refused(engine, tmp_path):
    data = (MESHES / "model_le.fil").read_bytes()
    path = tmp_path / "cut.fil"
    path.write_bytes(data[:5000])
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(path)


@pytest.mark.parametrize("name", MODELS)
def test_sniffed_without_the_extension(tmp_path, name):
    path = tmp_path / "results.bin2"
    shutil.copy(MESHES / name, path)
    assert meshioplusplus.sniff_format(path) == "abaqus_fil"
    assert _core.sniff_format(str(path)) == "abaqus_fil"
