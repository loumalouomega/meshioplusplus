import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _pvtk_structured as ps


def _grid(nx, ny, nz):
    x, y, z = np.meshgrid(
        np.arange(nx, dtype=float),
        np.arange(ny, dtype=float),
        np.arange(nz, dtype=float),
        indexing="ij",
    )
    pts = np.stack([x, y, z], -1).reshape(-1, 3)
    return pts


@pytest.mark.parametrize("kind", ["pvti", "pvtr", "pvts"])
def test_roundtrip_whole_lattice(tmp_path, kind):
    n = (4, 3, 2)
    pts = _grid(*n)
    if kind == "pvti":
        pytest.importorskip("numpy")
    # a lattice mesh in i-fastest order
    pts = np.stack([pts[:, 0], pts[:, 1], pts[:, 2]], -1)
    order = np.lexsort((pts[:, 0], pts[:, 1], pts[:, 2]))
    pts = pts[order]
    cell_type, conn = ps.lattice_cells(*n)
    mesh = mio.Mesh(pts, [(cell_type, conn)], point_data={"u": pts[:, 0].copy()})
    path = tmp_path / f"g.{kind}"
    ps.write(path, mesh, kind)
    out = mio.read(path)
    assert out.points.shape == pts.shape
    assert np.allclose(out.points, pts)
    assert np.allclose(out.point_data["u"], pts[:, 0])
    assert len(out.cells[0].data) == len(conn)


def test_two_pieces_share_a_plane(tmp_path):
    n = (5, 2, 2)
    pts = _grid(*n)
    order = np.lexsort((pts[:, 0], pts[:, 1], pts[:, 2]))
    pts = pts[order].reshape(2, 2, 5, 3)  # k, j, i
    extents = [(0, 2), (2, 4)]  # i ranges, sharing the i == 2 plane
    sources = []
    for number, (a, b) in enumerate(extents):
        sub = pts[:, :, a : b + 1].reshape(-1, 3)
        cell_type, conn = ps.lattice_cells(b - a + 1, 2, 2)
        mesh = mio.Mesh(sub, [(cell_type, conn)], point_data={"u": sub[:, 0].copy()})
        piece = tmp_path / f"p{number}.vts"
        ps.write(tmp_path / f"tmp{number}.pvts", mesh, "pvts")
        (tmp_path / f"tmp{number}" / f"tmp{number}_0000.vts").replace(piece)
        sources.append((f"p{number}.vts", f"{a} {b} 0 1 0 1"))
    body = "".join(f'<Piece Extent="{e}" Source="{s}"/>' for s, e in sources)
    (tmp_path / "w.pvts").write_text(
        '<?xml version="1.0"?><VTKFile type="PStructuredGrid" version="1.0">'
        f'<PStructuredGrid WholeExtent="0 4 0 1 0 1">{body}</PStructuredGrid></VTKFile>'
    )
    out = mio.read(tmp_path / "w.pvts")
    assert len(out.points) == 20
    assert np.allclose(np.sort(out.points[:, 0]), np.sort(_grid(*n)[:, 0]))
    assert len(out.cells[0].data) == 4
    assert {r.name for r in out.regions} == {"piece_0", "piece_1"}
    one = mio.pvts.read(tmp_path / "w.pvts", piece=1)
    assert len(one.points) == 12
