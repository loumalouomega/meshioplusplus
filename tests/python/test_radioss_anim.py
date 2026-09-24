"""OpenRadioss animation files (``<run>A001``...): both engines, the fixtures
written by ``tools/gen_radioss_anim_fixtures.py``, checked against what
OpenRadioss's own ``anim_to_vtk`` converter makes of them (``anim_to_vtk/``),
dispatch, sequences and refusals."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.radioss_anim import _anim as py_anim

MESHES = pathlib.Path(__file__).parent / "meshes" / "radioss_anim"
STATES = ["cubeA001", "cubeA002"]


@pytest.fixture(params=["core", "python"])
def read(request):
    if request.param == "core":
        return meshioplusplus.radioss_anim.read
    return py_anim.read


def _flat(arrays):
    return np.concatenate([np.asarray(a).reshape(len(a), -1) for a in arrays])


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


@pytest.mark.parametrize("name", STATES)
def test_engines_agree(name):
    a, b = _core.radioss_anim_read(str(MESHES / name)), py_anim.read(MESHES / name)
    np.testing.assert_array_equal(a.points, b.points)
    assert [(c.type, c.data.tolist()) for c in a.cells] == [
        (c.type, c.data.tolist()) for c in b.cells
    ]
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)
            assert np.asarray(x).dtype == np.asarray(y).dtype, k
    assert _regions(a) == _regions(b)


class _Ref:
    """What anim_to_vtk wrote: points, then point and cell data by name. It
    writes a triangle as a VTK triangle with four ids, which no VTK reader
    takes, so the sections are parsed here."""

    def __init__(self, path):
        tokens = path.read_text().split()
        self.point_data, self.cell_data = {}, {}
        k = 0
        target = None
        while k < len(tokens):
            t = tokens[k]
            if t == "POINTS":
                n = int(tokens[k + 1])
                self.points = np.array(tokens[k + 3 : k + 3 + 3 * n], float).reshape(
                    n, 3
                )
                self.n_points = n
                k += 3 + 3 * n
            elif t in ("POINT_DATA", "CELL_DATA"):
                count = int(tokens[k + 1])
                target = (
                    (self.point_data, count)
                    if t == "POINT_DATA"
                    else (self.cell_data, count)
                )
                k += 2
            elif t in ("SCALARS", "VECTORS", "TENSORS"):
                name = tokens[k + 1]
                width = {"SCALARS": 1, "VECTORS": 3, "TENSORS": 9}[t]
                k += 4 if t == "SCALARS" else 3
                if t == "SCALARS" and tokens[k - 1] != "1":
                    k -= 1
                if tokens[k] == "LOOKUP_TABLE":
                    k += 2
                store, count = target
                store[name] = [
                    np.array(tokens[k : k + width * count], float).reshape(count, width)
                ]
                k += width * count
            else:
                k += 1


def _anim_to_vtk(path):
    return _Ref(path)


@pytest.mark.parametrize("name", STATES)
def test_matches_anim_to_vtk(read, name):
    mesh = read(MESHES / name)
    ref = _anim_to_vtk(MESHES / "anim_to_vtk" / f"{name}.vtk")
    # anim_to_vtk prints six significant digits.
    np.testing.assert_allclose(mesh.points, ref.points, rtol=1e-5, atol=1e-6)
    for key in ("Velocity", "Displacement", "Temperature"):
        np.testing.assert_allclose(
            np.asarray(mesh.point_data[key]).reshape(len(mesh.points), -1),
            np.asarray(ref.point_data[key]).reshape(len(mesh.points), -1),
            rtol=1e-5,
            atol=1e-6,
        )
    np.testing.assert_array_equal(
        mesh.point_data["radioss:node_id"], np.ravel(ref.point_data["NODE_ID"])
    )
    # anim_to_vtk writes the 1-D, 2-D, 3-D and SPH cells in that order, which is
    # also the block order here (one type per family, but the 2-D and 3-D
    # families split into quad/triangle and hexahedron/tetra/wedge in order).
    np.testing.assert_array_equal(
        _flat(mesh.cell_data["radioss:element_id"]).ravel(),
        _flat(ref.cell_data["ELEMENT_ID"]).ravel(),
    )
    np.testing.assert_array_equal(
        _flat(mesh.cell_data["radioss:part"]).ravel(),
        _flat(ref.cell_data["PART_ID"]).ravel(),
    )
    types = np.concatenate([[c.type] * len(c.data) for c in mesh.cells])
    families = {
        "1DELEM_": ["line"],
        "2DELEM_": ["quad", "triangle"],
        "3DELEM_": ["hexahedron", "tetra", "wedge"],
        "SPHELEM_": ["vertex"],
    }
    for key in ref.cell_data:
        prefix = next((p for p in families if key.startswith(p)), None)
        if prefix is None or "Forces-Moments" in key:
            continue
        rows = np.isin(types, families[prefix])
        expected = _flat(ref.cell_data[key])[rows]
        # anim_to_vtk writes the names with underscores for spaces.
        (mine,) = [
            n for n in mesh.cell_data if n.replace(" ", "_") == key[len(prefix) :]
        ]
        got = _flat(mesh.cell_data[mine])[rows]
        if got.shape[1] == 6:
            # anim_to_vtk's 3x3 matrix holds the file's six values as
            # [[0, 3, 4], [3, 1, 5], [4, 5, 2]]; a 2-D tensor's (xx yy xy) as
            # [[0, 2, 0], [2, 1, 0], [0, 0, 0]].
            expected = expected[:, [0, 4, 8, 1, 2, 5]]
            got = np.nan_to_num(got)
        np.testing.assert_allclose(got, expected, rtol=1e-5, atol=1e-6, err_msg=key)
    forces = _flat(mesh.cell_data["Forces-Moments"])[types == "line"]
    for k, suffix in enumerate(["F1", "F2", "F3", "M1", "M2", "M3", "M4", "M5", "M6"]):
        np.testing.assert_allclose(
            forces[:, k],
            _flat(ref.cell_data[f"1DELEM_Forces-Moments{suffix}"])[:2, 0],
            rtol=1e-5,
        )


def test_cells_data_and_regions(read):
    mesh = read(MESHES / "cubeA002")
    assert [c.type for c in mesh.cells] == [
        "line",
        "quad",
        "triangle",
        "hexahedron",
        "tetra",
        "wedge",
        "vertex",
    ]
    # The triangle was a facet with a repeated node; the tetra and wedge
    # degenerate bricks.
    assert mesh.cells[2].data.tolist() == [[5, 10, 11]]
    assert mesh.cells[4].data.tolist() == [[1, 8, 9, 10]]
    assert float(mesh.field_data["meshio:time"][0]) == 0.5
    # One scalar name shared by the 2-D and 3-D families; NaN elsewhere.
    thickness = _flat(mesh.cell_data["Thickness"]).ravel()
    assert np.isnan(thickness[[0, 1, 7, 8]]).all()
    np.testing.assert_array_equal(thickness[2:7], [1.5, 2.5, 7.0, 8.0, 9.0])
    # A 2-D tensor keeps xx, yy, xy; its out-of-plane parts are NaN.
    stress = _flat(mesh.cell_data["Stress"])
    np.testing.assert_array_equal(stress[2][[0, 1, 3]], [10.0, 20.0, 5.0])
    assert np.isnan(stress[2][[2, 4, 5]]).all()
    np.testing.assert_array_equal(stress[4], [0, 1, 2, 3, 4, 5])
    assert _flat(mesh.cell_data["Forces-Moments"]).shape == (9, 9)
    # The second state deletes the tetra.
    np.testing.assert_array_equal(
        _flat(mesh.cell_data["radioss:alive"]).ravel(), [1, 1, 1, 1, 1, 0, 1, 1, 1]
    )
    np.testing.assert_array_equal(
        _flat(mesh.cell_data["radioss:material"]).ravel(),
        [7, 7, 7, 7, 8, 9, 9, 9, 9],
    )
    regions = {(r.name, r.tag): r for r in mesh.regions}
    assert sorted(regions) == [
        ("block", 21),
        ("frame", 31),
        ("side tria", 12),
        ("top shell", 11),
        ("water", 41),
        ("wedges", 22),
    ]
    assert regions[("wedges", 22)].dim == 3
    np.testing.assert_array_equal(regions[("wedges", 22)].entries, [5, 6])


def test_dispatch_sniffing_and_sequences(tmp_path):
    mesh = meshioplusplus.read(MESHES / "cubeA001")
    assert "radioss:alive" in mesh.cell_data
    other = tmp_path / "state.bin"
    shutil.copy(MESHES / "cubeA001", other)
    assert meshioplusplus.sniff_format(other) == "radioss_anim"
    assert _core.sniff_format(str(other)) == "radioss_anim"
    entries = meshioplusplus.sequence_entries(str(MESHES / "cubeA*"))
    assert [(e["time"], e["time_source"]) for e in entries] == [
        (0.0, "file"),
        (0.5, "file"),
    ]
    times = [t for t, _ in meshioplusplus.read_sequence(str(MESHES / "cubeA*"))]
    assert times == [0.0, 0.5]


def test_refusals(read, tmp_path):
    data = (MESHES / "cubeA001").read_bytes()
    bad = tmp_path / "badA001"
    bad.write_bytes(data[: len(data) // 2])
    with pytest.raises(meshioplusplus.ReadError, match="truncated"):
        read(bad)
    bad.write_bytes(b"\x00\x00\x54\x2a" + data[4:])
    with pytest.raises(meshioplusplus.ReadError, match="magic"):
        read(bad)
