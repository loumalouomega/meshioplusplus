"""Ansys MAPDL results (``.rst``, ``.rth``): both engines against files MAPDL
wrote (and pymapdl-reader's frozen reading of them), plus a synthetic file that
pins what those files do not exercise -- a node rotated about all three axes, a
result set holding only some nodes, MAPDL's undefined value and the refusals."""

import collections
import math
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.ansys_rst import _ansys_rst as py_rst

RST = pathlib.Path(__file__).parent / "meshes" / "ansys" / "rst"
FIXTURES = sorted(RST.glob("*.rst")) + sorted(RST.glob("*.rth"))

_VTK = {
    1: "vertex",
    3: "line",
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
}
# pymapdl-reader's DOF labels as (meshio++ array, component).
_DOF = {
    "UX": ("U", 0),
    "UY": ("U", 1),
    "UZ": ("U", 2),
    "ROTX": ("ROT", 0),
    "ROTY": ("ROT", 1),
    "ROTZ": ("ROT", 2),
}


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    if request.param == "core":
        return meshioplusplus.ansys_rst
    return py_rst


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for data_a, data_b in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(data_a) == sorted(data_b)
        for name in data_a:
            np.testing.assert_array_equal(data_a[name], data_b[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert sorted((r.kind, r.name, tuple(r.entries)) for r in a.regions) == sorted(
        (r.kind, r.name, tuple(r.entries)) for r in b.regions
    )


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_set(path):
    times = meshioplusplus.ansys_rst.time_values(path)
    assert times == py_rst.time_values(path)
    for step in range(len(times)):
        _same(
            meshioplusplus.ansys_rst.read(path, time_step=step, lenient=True),
            py_rst.read(path, time_step=step, lenient=True),
        )


def _cells(points, cell_type, rows):
    out = collections.Counter()
    for row in rows:
        out[(cell_type, tuple(sorted(tuple(np.round(p, 7)) for p in points[row])))] += 1
    return out


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_matches_pymapdl_reader(engine, path):
    """The cells, set times and nodal solutions pymapdl-reader reads."""
    ref = np.load(RST / "pymapdl_reference.npz")
    key = path.name
    times = ref[f"{key}/times"]
    mesh = engine.read(path, lenient=True)

    expected = collections.Counter()
    points, offsets = ref[f"{key}/points"], ref[f"{key}/offsets"]
    conn = ref[f"{key}/connectivity"]
    for c, vtk in enumerate(ref[f"{key}/celltypes"]):
        if vtk:  # pymapdl-reader keeps elements it cannot map as empty cells
            row = conn[offsets[c] : offsets[c + 1]]
            expected += _cells(points, _VTK[int(vtk)], [row])
    got = collections.Counter()
    for block in mesh.cells:
        got += _cells(mesh.points, block.type, block.data)
    assert got == expected

    n_file = len(mesh.points)
    node_ids = np.zeros(n_file, dtype=np.int64)
    ids = py_rst._Results(str(path)).model(True)[1]
    for number, point in ids.items():
        node_ids[point] = number
    index = {int(n): k for k, n in enumerate(node_ids)}
    for step, time in enumerate(times):
        mesh = engine.read(path, time_step=step, lenient=True)
        assert mesh.field_data["meshio:time"][0] == time
        rows = [index[int(n)] for n in ref[f"{key}/{step}/nnum"]]
        values = ref[f"{key}/{step}/values"]
        for k, label in enumerate(ref[f"{key}/{step}/dofs"].tolist()):
            name, component = _DOF.get(label, (label, None))
            ours = mesh.point_data[name][rows]
            if component is not None:
                ours = ours[:, component]
            np.testing.assert_allclose(ours, values[:, k], rtol=1e-12, atol=0)


def test_modal_file(engine):
    """file.rst: a modal analysis, six modes of 40 SOLID186 bricks."""
    path = RST / "file.rst"
    np.testing.assert_allclose(
        engine.time_values(path),
        [7366.495, 7366.495, 11504.895, 17285.705, 17285.705, 20137.193],
        atol=1e-3,
    )
    mesh = engine.read(path, time_step=-1)
    assert [b.type for b in mesh.cells] == ["hexahedron20"]
    assert mesh.cells[0].data.shape == (40, 20)
    assert mesh.field_data["meshio:time"][0] == pytest.approx(20137.193, abs=1e-3)
    assert mesh.field_data["ansys:cumulative"].tolist() == [6]
    first = engine.read(path)
    np.testing.assert_allclose(
        first.point_data["U"][0], [28.962, -28.248, -0.309], atol=1e-3
    )
    np.testing.assert_array_equal(first.cell_data["ansys:element"][0], 186)


def test_thermal_file_has_temperatures(engine):
    mesh = engine.read(RST / "file.rth")
    assert sorted(mesh.point_data) == ["TEMP"]
    assert np.isfinite(mesh.point_data["TEMP"]).all()


def test_components_are_regions(engine):
    mesh = engine.read(RST / "hex_201.rst")
    names = {(r.kind, r.name): len(r.entries) for r in mesh.regions}
    assert names[("point", "NCOMP2")] == 98
    assert names[("cell", "ECOMP1")] == 22


def test_cyclic_model_reads_its_sector_with_a_warning(engine, capfd):
    mesh = engine.read(RST / "cyclic_v182.rst")
    assert "cyclic-symmetry" in capfd.readouterr().err
    assert len(mesh.points) == 460


def test_selective_reads(engine):
    path = RST / "shell181_2021R1.rst"
    assert not engine.read(path, points_only=True).point_data
    assert sorted(engine.read(path, arrays=["ROT"]).point_data) == ["ROT"]
    assert not engine.read(path, arrays=[]).point_data
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(path, time_step=1)


def test_sequences_and_metadata(tmp_path):
    path = RST / "file.rst"
    assert len(meshioplusplus.read_metadata(path)["time_values"]) == 6
    meshioplusplus.write_sequence(
        str(tmp_path / "mode_{step}.vtu"), meshioplusplus.read_sequence(path)
    )
    written = sorted(tmp_path.glob("mode_*.vtu"))
    assert len(written) == 6
    np.testing.assert_allclose(
        meshioplusplus.read(written[-1]).point_data["U"],
        meshioplusplus.read(path, time_step=5).point_data["U"],
    )


def test_registered_and_sniffed(tmp_path):
    out = meshioplusplus.formats()
    assert out["extensions"][".rst"] == ["ansys_rst"]
    assert out["extensions"][".rth"] == ["ansys_rst"]
    assert "ansys_rst" in out["readable"] and "ansys_rst" not in out["writable"]
    renamed = tmp_path / "results.bin"
    renamed.write_bytes((RST / "beam44.rst").read_bytes())
    assert meshioplusplus.sniff_format(renamed) == "ansys_rst"
    assert hasattr(_core, "ansys_rst_read")


# -- a synthetic results file ------------------------------------------------------


def _i32(value):
    return value - (1 << 32) if value >= 1 << 31 else value


class _Writer:
    """Records as MAPDL lays them out: length in words, flags, data, trailer."""

    def __init__(self):
        self.words = []

    def raw(self, data_words, flags):
        ptr = len(self.words)
        self.words += [len(data_words), _i32(flags << 24), *data_words, 0]
        return ptr

    def ints(self, values, size=None):
        values = list(values) + [0] * ((size or len(values)) - len(values))
        return self.raw(values, 0x80)

    def doubles(self, values):
        return self.raw(np.asarray(values, "<f8").view("<i4").tolist(), 0x00)

    def shorts(self, values):
        values = list(values) + [0] * (len(values) % 2)
        return self.raw(np.asarray(values, "<i2").view("<i4").tolist(), 0xC0)

    def bsparse_doubles(self, values):
        bits = sum(1 << k for k, v in enumerate(values) if v != 0)
        packed = np.asarray([v for v in values if v != 0], "<f8").view("<i4").tolist()
        return self.raw([len(values), _i32(bits)] + packed, 0x08)

    def patch(self, ptr, k, value):
        self.words[ptr + 2 + k] = _i32(value & 0xFFFFFFFF)

    def save(self, path):
        pathlib.Path(path).write_bytes(np.asarray(self.words, "<i4").tobytes())


_ANGLES = (30.0, 45.0, 60.0)  # node 3: THXY, THYZ, THZX
_NEQV = [8, 7, 6, 5, 4, 3, 2, 1]  # the solution's node order
_HEX = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]


def _write_synthetic(path, zlib=False, global_nnod=0, sectors=1):
    w = _Writer()
    w.ints([12], 100)  # the standard header: file 12
    header = w.ints([12, 8, 8, 10, 4, 1, 1, 0, 2], 80)
    w.patch(header, 20, sectors)
    w.patch(header, 48, global_nnod)
    w.patch(header, 14, w.ints(_NEQV))

    geometry = w.ints([0, 1, 0, 8, 1], 80)
    w.patch(header, 15, geometry)
    ety = w.ints([0])
    w.patch(geometry, 20, ety)
    solid185 = w.ints([1, 185, 0], 100)
    w.patch(solid185, 60, 8)  # nodelm
    w.patch(ety, 0, solid185 - ety)
    nodes = None
    for number, xyz in enumerate(_HEX, start=1):
        angles = _ANGLES if number == 3 else (0.0, 0.0, 0.0)
        values = [float(number), *map(float, xyz), *angles]
        ptr = w.bsparse_doubles(values) if number == 1 else w.doubles(values)
        nodes = ptr if nodes is None else nodes
    w.patch(geometry, 26, nodes)
    eid = w.ints([0, 0])
    w.patch(geometry, 28, eid)
    w.patch(eid, 0, w.shorts([2, 1, 3, 4, 0, 1, 0, 0, 11, 0, *range(1, 9)]) - eid)
    comp = w.ints([1] + [int.from_bytes(b"FACE", "big")] + [0x20202020] * 7 + [1, -4])
    w.patch(geometry, 50, comp)
    w.patch(geometry, 48, 1)

    sets = []
    for s in range(2):
        base = w.ints([0, 1, 8], 150)
        w.patch(base, 19, 4)
        for k, code in enumerate([1, 2, 3, 20]):
            w.patch(base, 20 + k, code)
        if s == 0:  # every node: U = (n, 0, 0), TEMP = 1000 + n
            rows = [[n, 0.0, 0.0, 1000.0 + n] for n in _NEQV]
            nsl = w.doubles(np.ravel(rows))
        else:  # nodes 5 and 2 only, node 5's TEMP undefined
            rows = [[0.0, 5.0, 0.0, 2.0**100], [0.0, 2.0, 0.0, 1002.0]]
            nsl = w.doubles(np.ravel(rows))
            w.ints([_NEQV.index(5) + 1, _NEQV.index(2) + 1])
        w.patch(base, 104, nsl - base)
        sets.append(base)
    dsi = w.ints(sets, 20)
    tim = w.doubles([0.5, 1.0] + [0.0] * 8)
    lsp = w.ints([1, 1, 1, 1, 2, 2], 30)
    for k, ptr in ((10, dsi), (11, tim), (12, lsp)):
        w.patch(header, k, ptr)
    if zlib:
        w.words[tim + 1] = _i32(0x20 << 24)
    w.save(path)


def _rotation(xy, yz, zx):
    def about(axis, degrees):
        c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
        i, j = [(1, 2), (2, 0), (0, 1)][axis]
        m = np.eye(3)
        m[i, i], m[i, j], m[j, i], m[j, j] = c, -s, s, c
        return m

    return about(2, xy) @ about(0, yz) @ about(1, zx)


def test_synthetic_file(engine, tmp_path):
    path = tmp_path / "synthetic.rst"
    _write_synthetic(path)
    assert engine.time_values(path) == [0.5, 1.0]

    mesh = engine.read(path)
    assert [b.type for b in mesh.cells] == ["hexahedron"]
    np.testing.assert_array_equal(mesh.cells[0].data, [list(range(8))])
    np.testing.assert_array_equal(mesh.points, np.array(_HEX, dtype=float))
    for name, value in (("ansys:mat", 2), ("ansys:real", 3), ("ansys:secnum", 4)):
        assert mesh.cell_data[name][0].tolist() == [value]
    assert [(r.kind, r.name, r.entries.tolist()) for r in mesh.regions] == [
        ("point", "FACE", [0, 1, 2, 3])
    ]
    assert mesh.field_data["meshio:time"].tolist() == [0.5]
    np.testing.assert_array_equal(mesh.point_data["TEMP"], 1000.0 + np.arange(1, 9))
    expected = np.array([[n, 0.0, 0.0] for n in range(1, 9)])
    # Node 3's axes are turned about Z, then the new X, then the newest Y.
    expected[2] = _rotation(*_ANGLES) @ expected[2]
    np.testing.assert_allclose(mesh.point_data["U"], expected, rtol=1e-14, atol=1e-14)

    partial = engine.read(path, time_step=1)
    assert partial.field_data["ansys:substep"].tolist() == [2]
    temp, u = partial.point_data["TEMP"], partial.point_data["U"]
    assert np.isnan(temp[[0, 2, 3, 4, 5, 6, 7]]).all()  # node 5's is undefined
    assert temp[1] == 1002.0
    np.testing.assert_array_equal(u[1], [0.0, 2.0, 0.0])
    np.testing.assert_array_equal(u[4], [0.0, 5.0, 0.0])
    assert np.isnan(u[[0, 2, 3, 5, 6, 7]]).all()


@pytest.mark.parametrize(
    "options, match",
    [
        ({"zlib": True}, "zlib"),
        ({"global_nnod": 16}, "distributed"),
    ],
)
def test_refusals(engine, tmp_path, options, match):
    path = tmp_path / "bad.rst"
    _write_synthetic(path, **options)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)
    junk = tmp_path / "junk.rst"
    junk.write_bytes(b"\x00" * 64)
    with pytest.raises(meshioplusplus.ReadError, match="not a MAPDL results file"):
        engine.read(junk)
    cut = tmp_path / "cut.rst"
    cut.write_bytes((RST / "beam44.rst").read_bytes()[:5000])
    with pytest.raises(meshioplusplus.ReadError, match="outside the file|past the end"):
        engine.read(cut)
