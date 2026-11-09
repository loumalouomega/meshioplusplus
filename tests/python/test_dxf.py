"""DXF R2000 drawings (roadmap 7.6.4): both engines, byte for byte."""

import math

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus.dxf import _dxf

needs_core = pytest.mark.skipif(
    not hasattr(_core, "dxf_write"), reason="needs the C++ core"
)


def _line_mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0.5]], float)
    mesh = mp.Mesh(
        pts, [("vertex", np.array([[0]])), ("line", np.array([[0, 1], [1, 2], [2, 3]]))]
    )
    mesh.field_data["meshio:length_unit"] = np.array([0.001])
    return mesh


def _slot_sketch():
    outline = []
    for t in np.linspace(-np.pi / 2, np.pi / 2, 25):
        outline.append((2 + 0.75 * np.cos(t), 0.75 * np.sin(t)))
    for t in np.linspace(np.pi / 2, 3 * np.pi / 2, 25):
        outline.append((-2 + 0.75 * np.cos(t), 0.75 * np.sin(t)))
    n = len(outline)
    pts = [(x, y, 0.0) for x, y in outline] + [(x, y, 1.0) for x, y in outline]
    tris = []
    for i in range(n):
        j = (i + 1) % n
        tris += [(i, j, j + n), (i, j + n, i + n)]
    prism = mp.Mesh(np.array(pts, float), [("triangle", np.array(tris))])
    section = mp.slice(prism, (0, 0, 0.5), (0, 0, 1), chain=True)
    return mp.fit_curves(section, tolerance=1e-6)


def _circle_sketch():
    t = np.linspace(0, 2 * math.pi, 41)[:-1]
    ring = np.c_[np.cos(t), np.sin(t), np.zeros_like(t)]
    segs = np.c_[np.arange(40), (np.arange(40) + 1) % 40]
    mesh = mp.Mesh(
        ring,
        [("line", segs)],
        cell_data={
            "section:chain": [np.zeros(40, np.int64)],
            "section:closed": [np.ones(40, np.int64)],
        },
    )
    return mp.fit_curves(mesh, tolerance=1e-3)


_ENTITIES = """0
SECTION
2
HEADER
9
$INSUNITS
70
4
0
ENDSEC
0
SECTION
2
ENTITIES
0
LWPOLYLINE
8
0
90
3
70
1
10
0.0
20
0.0
42
1.0
10
2.0
20
0.0
10
2.0
20
1.0
0
SPLINE
8
0
71
3
72
8
73
4
40
0
40
0
40
0
40
0
40
1
40
1
40
1
40
1
10
0
20
0
30
0
10
1
20
2
30
0
10
3
20
2
30
0
10
4
20
0
30
0
0
ARC
8
0
10
5
20
5
30
0
40
2
50
0
51
90
0
ENDSEC
0
EOF
"""


@needs_core
@pytest.mark.parametrize("make", [_line_mesh, _slot_sketch, _circle_sketch])
def test_writers_are_byte_identical(make, tmp_path):
    mesh = make()
    a, b = tmp_path / "py.dxf", tmp_path / "cpp.dxf"
    _dxf.write(a, mesh)
    _core.dxf_write(str(b), mesh)
    assert a.read_bytes() == b.read_bytes()


@needs_core
def test_readers_agree(tmp_path):
    path = tmp_path / "e.dxf"
    path.write_text(_ENTITIES)
    py = _dxf.read(path)
    cpp = mp.read(path)
    np.testing.assert_array_equal(py.points, cpp.points)
    np.testing.assert_array_equal(py.cells[0].data, cpp.cells[0].data)
    for name in ("dxf:entity", "dxf:type"):
        np.testing.assert_array_equal(py.cell_data[name][0], cpp.cell_data[name][0])
    assert sorted(set(cpp.cell_data["dxf:type"][0].tolist())) == [2, 4, 5]
    assert cpp.field_data["meshio:length_unit"][0] == 0.001
    # The bulge-1 span is a half circle of radius 1 about (1, 0).
    poly = cpp.cell_data["dxf:entity"][0] == 0
    ids = np.unique(cpp.cells[0].data[poly])
    p = cpp.points[ids]
    on_arc = p[:, 1] < -1e-12
    assert np.allclose(np.hypot(p[on_arc, 0] - 1, p[on_arc, 1]), 1.0)
    # The clamped cubic spline starts and ends on its end control points.
    spl = cpp.cells[0].data[cpp.cell_data["dxf:entity"][0] == 1]
    ends = cpp.points[[spl[0, 0], spl[-1, 1]]]
    np.testing.assert_allclose(ends, [[0, 0, 0], [4, 0, 0]], atol=1e-12)


@pytest.mark.parametrize("make", [_line_mesh, _slot_sketch])
def test_round_trip(make, tmp_path):
    mesh = make()
    path = tmp_path / "rt.dxf"
    mp.write(path, mesh)
    back = mp.read(path)
    lines = [c.data for c in mesh.cells if c.type == "line"][0]
    if "sketch:type" in mesh.field_data:
        # Lines come back exactly; arcs as 64 segments per full turn on the radius.
        types = back.cell_data["dxf:type"][0]
        assert sorted(set(types.tolist())) == [1, 2]
        arc_pts = back.points[np.unique(back.cells[0].data[types == 2])]
        r = np.minimum(
            np.hypot(arc_pts[:, 0] - 2, arc_pts[:, 1]),
            np.hypot(arc_pts[:, 0] + 2, arc_pts[:, 1]),
        )
        np.testing.assert_allclose(r, 0.75, atol=1e-9)
    else:
        assert len(back.cells[0].data) == len(lines)
        np.testing.assert_array_equal(back.points, mesh.points)
        assert back.field_data["meshio:length_unit"][0] == 0.001


def test_ezdxf_audit_is_clean(tmp_path):
    ezdxf = pytest.importorskip("ezdxf")
    for k, make in enumerate((_line_mesh, _slot_sketch, _circle_sketch)):
        path = tmp_path / f"a{k}.dxf"
        mp.write(path, make())
        doc = ezdxf.readfile(path)
        auditor = doc.audit()
        assert (len(auditor.errors), len(auditor.fixes)) == (0, 0)
        assert doc.dxfversion == "AC1015"


def test_refusals(tmp_path):
    bad = tmp_path / "t.dxf"
    bad.write_text("0\nSECTION\n2\nENTITIES\n0\nTEXT\n1\nhi\n0\nENDSEC\n0\nEOF\n")
    with pytest.raises(Exception, match="TEXT is not read"):
        mp.read(bad)
    with pytest.raises(Exception, match="TEXT is not read"):
        _dxf.read(bad)
    with pytest.raises(Exception):
        mp.write(tmp_path / "tri.dxf", mp.icosa_sphere(1.0, 0))
