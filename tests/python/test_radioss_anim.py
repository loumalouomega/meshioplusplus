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


def test_real_openradioss_run(read):
    """Animation files OpenRadioss wrote for a meshio++ deck (see
    radioss_th/column/README.md): displacement is the move from the first
    state, and the values are those anim_to_vtk prints."""
    column = MESHES.parent / "radioss_th" / "column"
    first = read(column / "columnA001")
    second = read(column / "columnA002")
    assert first.field_data["meshio:time"][0] == 0.0
    assert second.field_data["meshio:time"][0] == pytest.approx(0.00205397, rel=1e-5)
    assert [c.type for c in second.cells] == ["hexahedron"]
    np.testing.assert_allclose(
        second.point_data["Displacement"], second.points - first.points, atol=1e-6
    )
    np.testing.assert_allclose(
        second.cell_data["Von Mises"][0], [0.0110795, 0.0190114], rtol=1e-5
    )
    np.testing.assert_allclose(
        second.point_data["Velocity"][8],
        [-0.0497561, -0.0497561, -0.0165583],
        rtol=1e-5,
    )
    np.testing.assert_array_equal(second.cell_data["radioss:element_id"][0], [1, 2])


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


# ---------------------------------------------------------------------------
# The writer: both engines write the same bytes, and either reader gets the
# mesh back (to the 32-bit floats the file holds).
# ---------------------------------------------------------------------------

from meshioplusplus.radioss_anim import _write as py_write  # noqa: E402


@pytest.fixture(params=["core", "python"])
def write_engine(request):
    if request.param == "core":
        return lambda path, mesh: _core.radioss_anim_write(str(path), mesh)
    return py_write.write


def _census(mesh):
    counts = {}
    for block in mesh.cells:
        counts[block.type] = counts.get(block.type, 0) + len(block.data)
    return counts


def _all_families_mesh():
    points = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 0.5, 1],
            [3, 0, 0],
        ],
        dtype=float,
    )
    cells = [
        ("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]])),
        ("tetra", np.array([[1, 8, 9, 5]])),
        ("wedge", np.array([[1, 8, 9, 5, 10, 6]])),
        ("pyramid", np.array([[0, 1, 2, 3, 4]])),
        ("quad", np.array([[4, 5, 6, 7]])),
        ("triangle", np.array([[5, 10, 6]])),
        ("line", np.array([[0, 8], [8, 11]])),
        ("vertex", np.array([[11], [10]])),
    ]
    mesh = meshioplusplus.Mesh(points, cells)
    n = len(points)
    mesh.point_data["Temperature"] = np.arange(n, dtype=float) + 20.0
    mesh.point_data["Velocity"] = np.arange(3 * n, dtype=float).reshape(n, 3)
    ncell = sum(len(c.data) for c in mesh.cells)
    per = [len(c.data) for c in mesh.cells]
    split = np.cumsum(per)[:-1]
    mesh.cell_data["Thickness"] = list(
        np.split(np.arange(ncell, dtype=float) + 0.5, split)
    )
    mesh.cell_data["radioss:part"] = list(
        np.split(np.array([1, 1, 1, 1, 2, 2, 3, 3, 4, 4], dtype=np.int64), split)
    )
    mesh.field_data["meshio:time"] = np.array([0.25])
    return mesh


def test_the_engines_write_the_same_bytes(tmp_path):
    for name in STATES:
        mesh = meshioplusplus.read(MESHES / name)
        _core.radioss_anim_write(str(tmp_path / "cA001"), mesh)
        py_write.write(tmp_path / "pA001", mesh)
        assert (tmp_path / "cA001").read_bytes() == (tmp_path / "pA001").read_bytes()
    mesh = _all_families_mesh()
    _core.radioss_anim_write(str(tmp_path / "cA001"), mesh)
    py_write.write(tmp_path / "pA001", mesh)
    assert (tmp_path / "cA001").read_bytes() == (tmp_path / "pA001").read_bytes()


@pytest.mark.parametrize("name", STATES)
def test_a_state_survives_a_round_trip(name, write_engine, read, tmp_path):
    original = meshioplusplus.read(MESHES / name)
    out = tmp_path / "outA001"
    write_engine(out, original)
    again = read(out)
    np.testing.assert_allclose(again.points, original.points)
    assert _census(again) == _census(original)
    assert sorted(again.point_data) == sorted(original.point_data)
    for key, value in original.point_data.items():
        np.testing.assert_allclose(again.point_data[key], value, equal_nan=True)
    assert sorted(again.cell_data) == sorted(original.cell_data)
    for key in original.cell_data:
        np.testing.assert_allclose(
            _flat(again.cell_data[key]), _flat(original.cell_data[key]), equal_nan=True
        )
    assert _regions(again) == _regions(original)
    assert float(again.field_data["meshio:time"][0]) == pytest.approx(
        float(original.field_data["meshio:time"][0])
    )


def test_every_family_from_scratch(write_engine, read, tmp_path):
    mesh = _all_families_mesh()
    out = tmp_path / "scratchA001"
    write_engine(out, mesh)
    again = read(out)
    assert _census(again) == _census(mesh)
    np.testing.assert_allclose(again.point_data["Temperature"], 20 + np.arange(12))
    np.testing.assert_allclose(
        again.point_data["Velocity"], mesh.point_data["Velocity"]
    )
    assert float(again.field_data["meshio:time"][0]) == 0.25
    # Elements come back grouped by family, in the file's order; the part ids
    # and the thickness follow them.
    parts = {}
    for block, ids in zip(again.cells, again.cell_data["radioss:part"]):
        parts.setdefault(block.type, []).extend(ids.tolist())
    assert sorted(parts["hexahedron"]) == [1]
    assert sorted(parts["quad"] + parts["triangle"]) == [2, 2]
    assert sorted(parts["line"]) == [3, 3]
    assert sorted(parts["vertex"]) == [4, 4]


def test_cells_with_no_family_are_skipped_with_a_warning(write_engine, read, tmp_path):
    mesh = _all_families_mesh()
    mesh.cells.append(meshioplusplus.CellBlock("tetra10", np.arange(10)[None, :]))
    for name in list(mesh.cell_data):
        mesh.cell_data[name].append(np.zeros(1, dtype=mesh.cell_data[name][0].dtype))
    out = tmp_path / "skipA001"
    write_engine(out, mesh)
    assert "tetra10" not in _census(read(out))


def test_a_series_is_one_numbered_file_per_step(tmp_path):
    steps = []
    for k in range(3):
        mesh = _all_families_mesh()
        mesh.point_data["Temperature"] = mesh.point_data["Temperature"] + k
        steps.append((0.5 * k, mesh))
    written = meshioplusplus.write_sequence(
        tmp_path / "runA001", iter(steps), file_format="radioss_anim"
    )
    assert [pathlib.Path(p).name for p in written] == ["runA001", "runA002", "runA003"]
    back = list(meshioplusplus.read_sequence(str(tmp_path / "runA*")))
    assert [t for t, _ in back] == [0.0, 0.5, 1.0]
    for (_, got), (_, mesh) in zip(back, steps):
        np.testing.assert_allclose(
            got.point_data["Temperature"], mesh.point_data["Temperature"]
        )
    # The Python twin numbers the same files.
    with py_write.SeriesWriter(tmp_path / "pyA010") as w:
        w.write(0.0, steps[0][1])
        w.write(1.0, steps[1][1])
    assert (tmp_path / "pyA010").exists() and (tmp_path / "pyA011").exists()


def test_the_name_decides_the_format_and_a_plain_name_is_a_stem(tmp_path):
    formats = meshioplusplus.formats()
    assert "radioss_anim" in formats["writable"]
    meshioplusplus.write(tmp_path / "autoA001", _all_families_mesh())
    assert meshioplusplus.read(tmp_path / "autoA001").points.shape == (12, 3)
    # A name that is not an animation file's is a stem.
    written = meshioplusplus.write_sequence(
        tmp_path / "plain.out",
        [(0.0, _all_families_mesh())],
        file_format="radioss_anim",
    )
    assert [pathlib.Path(p).name for p in written] == ["plain.outA001"]
    with py_write.SeriesWriter(tmp_path / "stem") as w:
        w.write(0.0, _all_families_mesh())
    assert (tmp_path / "stemA001").exists()
