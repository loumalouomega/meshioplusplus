"""Roadmap 3.3.7: field_data and regions in the structured-grid writers
(`.vti`, `.vtr`, `.vts`), both engines."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.vti import _vti
from meshioplusplus.vtr import _vtr
from meshioplusplus.vts import _vts

TWINS = {"vti": _vti, "vtr": _vtr, "vts": _vts}
HAS_ZLIB = getattr(_core, "__has_zlib__", False)
ZLIB_REASON = "the core writes binary VTK XML zlib-compressed; this build has no zlib"


def _lattice():
    n = 3
    points = np.array(
        [[i, j, k] for k in range(n) for j in range(n) for i in range(n)], float
    )

    def idx(i, j, k):
        return i + n * (j + n * k)

    cells = [
        [
            idx(i, j, k),
            idx(i + 1, j, k),
            idx(i + 1, j + 1, k),
            idx(i, j + 1, k),
            idx(i, j, k + 1),
            idx(i + 1, j, k + 1),
            idx(i + 1, j + 1, k + 1),
            idx(i, j + 1, k + 1),
        ]
        for k in range(2)
        for j in range(2)
        for i in range(2)
    ]
    return meshioplusplus.Mesh(
        points,
        [("hexahedron", cells)],
        field_data={"g": np.array([9.81]), "v": np.array([1, 2, 3])},
        regions=[
            meshioplusplus.Region("walls", "cell", [0, 3, 5], dim=3, tag=7),
            meshioplusplus.Region("inlet", "point", [1, 2]),
        ],
    )


def _check(back):
    assert sorted(back.field_data) == ["g", "v"]
    assert back.field_data["g"][0] == 9.81
    np.testing.assert_array_equal(back.field_data["v"], [1, 2, 3])
    regions = {r.name: r for r in back.regions}
    assert sorted(regions) == ["inlet", "walls"]
    assert list(regions["walls"].entries) == [0, 3, 5]
    assert (regions["walls"].dim, regions["walls"].tag) == (3, 7)
    assert list(regions["inlet"].entries) == [1, 2]


@pytest.mark.parametrize("fmt", ["vti", "vtr", "vts"])
@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("engine", ["core", "python"])
def test_field_data_and_regions_round_trip(fmt, binary, engine, tmp_path):
    if engine == "core" and binary and not HAS_ZLIB:
        pytest.skip(ZLIB_REASON)
    path = tmp_path / f"m.{fmt}"
    mesh = _lattice()
    if engine == "core":
        getattr(_core, f"{fmt}_write")(str(path), mesh, binary, True)
    else:
        TWINS[fmt].write(path, mesh, binary=binary)
    _check(meshioplusplus.read(path))
    _check(TWINS[fmt].read(path))  # the Python reader agrees with the core's


@pytest.mark.parametrize("fmt", ["vti", "vtr", "vts"])
def test_engines_write_the_same_binary_bytes(fmt, tmp_path):
    if not HAS_ZLIB:
        pytest.skip(ZLIB_REASON)
    a, b = tmp_path / f"core.{fmt}", tmp_path / f"twin.{fmt}"
    mesh = _lattice()
    getattr(_core, f"{fmt}_write")(str(a), mesh, True, True)
    TWINS[fmt].write(b, mesh, binary=True)
    assert a.read_bytes() == b.read_bytes()


# --------------------------------------------------------------------------- #
# Roadmap 3.3.7: pvtu / vtm piece regions survive a write                      #
# --------------------------------------------------------------------------- #


def _pieces_mesh():
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [2, 0, 0], [2, 1, 0.0]]
    )
    return meshioplusplus.Mesh(
        points,
        [("triangle", [[0, 1, 2], [1, 3, 2]]), ("triangle", [[1, 4, 3], [4, 5, 3]])],
        regions=[
            meshioplusplus.Region("left", "cell", [0, 1]),
            meshioplusplus.Region("right", "cell", [2, 3]),
        ],
    )


def test_vtm_pieces_are_named_after_their_regions(tmp_path):
    from meshioplusplus.vtm import _vtm

    for name, write in (
        ("core", lambda p, m: meshioplusplus.vtm.write(p, m)),
        ("twin", lambda p, m: _vtm.write(p, m)),
    ):
        path = tmp_path / f"{name}.vtm"
        write(path, _pieces_mesh())
        back = meshioplusplus.read(path)
        assert sorted(r.name for r in back.regions) == ["left", "right"], name


def test_vtm_region_names_are_escaped(tmp_path):
    mesh = _pieces_mesh()
    mesh.regions = [
        meshioplusplus.Region('a"<b>&', "cell", [0, 1]),
        meshioplusplus.Region("right", "cell", [2, 3]),
    ]
    path = tmp_path / "q.vtm"
    meshioplusplus.vtm.write(path, mesh)
    back = meshioplusplus.read(path)
    assert 'a"<b>&' in {r.name for r in back.regions}


@pytest.mark.parametrize("engine", ["core", "python"])
def test_pvtu_write_keeps_the_pieces_of_a_read_file(engine, tmp_path):
    from meshioplusplus.pvtu import _pvtu

    if engine == "core" and not HAS_ZLIB:
        pytest.skip(ZLIB_REASON)
    mesh = _pieces_mesh()
    mesh.regions = [
        meshioplusplus.Region("piece_0", "cell", [0, 1]),
        meshioplusplus.Region("piece_1", "cell", [2, 3]),
    ]
    mesh.cell_data["c"] = [np.arange(2.0), np.arange(2.0, 4.0)]
    path = tmp_path / f"{engine}.pvtu"
    if engine == "core":
        _core.pvtu_write_codec(str(path), mesh, True, "zlib", "partition:part")
    else:
        _pvtu.write(path, mesh)
    pieces = sorted((tmp_path / engine).glob("*.vtu"))
    assert len(pieces) == 2
    back = meshioplusplus.read(path)
    assert [(r.name, list(r.entries)) for r in back.regions] == [
        ("piece_0", [0, 1]),
        ("piece_1", [2, 3]),
    ]
    assert list(back.cell_data) == ["c"]


def test_the_python_vtu_writer_takes_a_piece_that_opens_with_an_empty_block(tmp_path):
    # carve_by_part keeps the parts with no cells, so a piece can open with an empty
    # cell block; the offsets of the blocks after it must still start from zero.
    from meshioplusplus.vtu import _vtu

    points = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [2, 1, 0.0], [1, 1, 0]])
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("triangle", np.empty((0, 3), dtype=np.int64)),
            ("triangle", [[0, 1, 4], [1, 2, 3]]),
        ],
    )
    path = tmp_path / "empty_first.vtu"
    _vtu.write(path, mesh)
    back = meshioplusplus.read(path)
    cells = np.concatenate([np.asarray(c.data) for c in back.cells if len(c.data)])
    np.testing.assert_array_equal(cells, [[0, 1, 4], [1, 2, 3]])
