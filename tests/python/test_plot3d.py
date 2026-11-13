"""PLOT3D (roadmap 2.4.2): every variant generated here, read in both engines."""

import struct

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus._fortran_records import fortran_record
from meshioplusplus.plot3d import _plot3d


def _blocks(dim=3):
    """Two small blocks: 3x2x2 and 2x2x2 (3-D) or 3x2 and 2x2 (2-D)."""
    shapes = [(3, 2, 2), (2, 2, 2)] if dim == 3 else [(3, 2, 1), (2, 2, 1)]
    out = []
    for b, (ni, nj, nk) in enumerate(shapes):
        k, j, i = np.meshgrid(
            np.arange(nk), np.arange(nj), np.arange(ni), indexing="ij"
        )
        xyz = np.stack(
            [i.ravel() + 10.0 * b, j.ravel() * 0.5, k.ravel() * 0.25], axis=1
        ).astype(float)
        out.append(((ni, nj, nk), xyz))
    return out


def _encode(blocks, *, dim=3, w=8, iblank=False, multi=True, enc="fortran", big=False):
    o = ">" if big else "<"
    rf = o + ("f4" if w == 4 else "f8")
    dims = [d[:dim] for d, _ in blocks]
    coords = []
    for d, xyz in blocks:
        n = len(xyz)
        payload = xyz[:, :dim].T.astype(rf).tobytes()
        blank = (np.arange(n) % 2).astype(o + "i4")
        if iblank:
            payload += blank.tobytes()
        coords.append((xyz[:, :dim].T.ravel(), blank if iblank else None, payload))
    if enc == "ascii":
        words = ([str(len(blocks))] if multi else []) + [
            str(v) for d in dims for v in d
        ]
        for vals, blank, _ in coords:
            words += [repr(float(v)) for v in vals]
            if blank is not None:
                words += [str(int(v)) for v in blank]
        return (" ".join(words) + "\n").encode()
    header = [struct.pack(o + "i", len(blocks))] if multi else []
    header.append(np.asarray(dims, dtype=o + "i4").tobytes())
    if enc == "raw":
        return b"".join(header) + b"".join(c[2] for c in coords)
    layout = (4, o)
    return b"".join(fortran_record(h, layout) for h in header) + b"".join(
        fortran_record(c[2], layout) for c in coords
    )


@pytest.fixture(params=["core", "python"])
def engine(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
        return mio.plot3d.read
    return _plot3d.read


def _check(mesh, blocks, dim, iblank):
    expected = np.concatenate([xyz for _, xyz in blocks])
    if dim == 2:
        expected[:, 2] = 0.0
    assert np.allclose(mesh.points, expected)
    dims = np.asarray(mesh.field_data["plot3d:dims"])
    assert dims.tolist() == [list(d) for d, _ in blocks]
    want_type = "hexahedron" if dim == 3 else "quad"
    assert [c.type for c in mesh.cells] == [want_type] * len(blocks)
    assert [r.name for r in mesh.regions] == [f"block_{b}" for b in range(len(blocks))]
    if iblank:
        assert np.asarray(mesh.point_data["plot3d:iblank"]).tolist() == sum(
            [list(np.arange(len(x)) % 2) for _, x in blocks], []
        )
    else:
        assert "plot3d:iblank" not in mesh.point_data


@pytest.mark.parametrize("enc", ["fortran", "raw", "ascii"])
@pytest.mark.parametrize("dim", [3, 2])
@pytest.mark.parametrize("w", [4, 8])
@pytest.mark.parametrize("iblank", [False, True])
def test_variants(tmp_path, engine, enc, dim, w, iblank):
    if enc == "ascii" and w == 4:
        pytest.skip("ASCII has no precision")
    blocks = _blocks(dim)
    path = tmp_path / "g.x"
    path.write_bytes(_encode(blocks, dim=dim, w=w, iblank=iblank, enc=enc))
    _check(engine(path), blocks, dim, iblank)


@pytest.mark.parametrize("enc", ["fortran", "raw"])
def test_big_endian_and_single_block(tmp_path, engine, enc):
    blocks = _blocks(3)[:1]
    path = tmp_path / "g.p3d"
    path.write_bytes(_encode(blocks, multi=False, big=True, enc=enc))
    _check(engine(path), blocks, 3, False)


def test_options_pick_or_refuse_a_reading(tmp_path, engine):
    blocks = _blocks(3)[:1]
    path = tmp_path / "g.x"
    path.write_bytes(_encode(blocks, multi=False, enc="raw", w=8))
    _check(engine(path, precision="double", multiblock=False), blocks, 3, False)
    with pytest.raises(mio.ReadError, match="plot3d"):
        engine(path, precision="single", multiblock=False, iblank=False)


def test_garbage_is_refused(tmp_path, engine):
    path = tmp_path / "g.x"
    path.write_bytes(b"\x01\x02\x03")
    with pytest.raises(mio.ReadError, match="not a PLOT3D grid"):
        engine(path)


def test_solution_and_function_files(tmp_path, engine):
    blocks = _blocks(3)
    grid = tmp_path / "g.x"
    grid.write_bytes(_encode(blocks))
    layout = (4, "<")
    q = [fortran_record(struct.pack("<i", 2), layout)]
    q.append(
        fortran_record(np.asarray([d for d, _ in blocks], "<i4").tobytes(), layout)
    )
    f = [fortran_record(struct.pack("<i", 2), layout)]
    f.append(
        fortran_record(
            np.asarray([list(d) + [2] for d, _ in blocks], "<i4").tobytes(), layout
        )
    )
    rho, energy = [], []
    for b, (d, xyz) in enumerate(blocks):
        n = len(xyz)
        values = np.arange(5 * n, dtype=float).reshape(5, n) + 100 * b
        rho.append(values[0])
        energy.append(values[4])
        q.append(
            fortran_record(
                np.array([0.5, 2.0, 1e6, 3.0]).astype("<f8").tobytes(), layout
            )
        )
        q.append(fortran_record(values.astype("<f8").tobytes(), layout))
        f.append(
            fortran_record((np.ones((2, n)) * (b + 1)).astype("<f8").tobytes(), layout)
        )
    (tmp_path / "g.q").write_bytes(b"".join(q))
    (tmp_path / "g.f").write_bytes(b"".join(f))
    mesh = engine(grid, solution=tmp_path / "g.q", function=tmp_path / "g.f")
    assert np.allclose(mesh.point_data["density"], np.concatenate(rho))
    assert np.allclose(mesh.point_data["energy"], np.concatenate(energy))
    assert np.asarray(mesh.point_data["momentum"]).shape == (len(mesh.points), 3)
    assert np.allclose(
        mesh.field_data["plot3d:q_reference"], [[0.5, 2.0, 1e6, 3.0]] * 2
    )
    assert np.asarray(mesh.point_data["plot3d:f"]).shape == (len(mesh.points), 2)


def test_write_is_byte_equal_to_the_unformatted_form(tmp_path):
    blocks = _blocks(3)
    src = tmp_path / "g.x"
    src.write_bytes(_encode(blocks, iblank=True))
    mesh = mio.read(src)
    for writer, name in ((mio.plot3d.write, "core.x"), (_plot3d.write, "py.x")):
        out = tmp_path / name
        writer(out, mesh)
        assert out.read_bytes() == src.read_bytes()


def test_a_mesh_that_is_not_blocks_is_refused(tmp_path):
    mesh = mio.Mesh(np.zeros((3, 3)), [("triangle", np.array([[0, 1, 2]]))])
    for writer in (mio.plot3d.write, _plot3d.write):
        with pytest.raises(mio.WriteError, match="plot3d:dims"):
            writer(tmp_path / "x.x", mesh)
