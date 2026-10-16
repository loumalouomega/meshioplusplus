import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.add_cell_data(helpers.tri_mesh, [("medit:ref", (), int)]),
    ],
)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.medit.write, meshioplusplus.medit.read, mesh, 1.0e-15
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.mesh")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.mesh")
    # same for binary files
    helpers.generic_io(tmp_path / "test.meshb")
    helpers.generic_io(tmp_path / "test.0.meshb")


# same tests with ugrid format files converted with UGC from
# https://www.simcenter.msstate.edu/


@pytest.mark.parametrize(
    "filename, ref_num_points, ref_num_triangle, ref_num_quad, ref_num_wedge, ref_num_tet, ref_num_hex, ref_tag_counts",
    [
        (
            "sphere_mixed.1.meshb",
            3270,
            864,
            0,
            3024,
            9072,
            0,
            {1: 432, 2: 216, 3: 216},
        ),
        ("hch_strct.4.meshb", 306, 12, 178, 96, 0, 144, {1: 15, 2: 15, 3: 160}),
        ("hch_strct.4.be.meshb", 306, 12, 178, 96, 0, 144, {1: 15, 2: 15, 3: 160}),
        ("cube86.mesh", 39, 72, 0, 0, 86, 0, {1: 14, 2: 14, 3: 14, 4: 8, 5: 14, 6: 8}),
    ],
)
def test_reference_file(
    filename,
    ref_num_points,
    ref_num_triangle,
    ref_num_quad,
    ref_num_wedge,
    ref_num_tet,
    ref_num_hex,
    ref_tag_counts,
):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "medit" / filename

    mesh = meshioplusplus.read(filename)
    assert mesh.points.shape[0] == ref_num_points
    assert mesh.points.shape[1] == 3

    medit_meshio_id = {
        "triangle": None,
        "quad": None,
        "tetra": None,
        "pyramid": None,
        "wedge": None,
        "hexahedron": None,
    }

    for i, cell_block in enumerate(mesh.cells):
        if cell_block.type in medit_meshio_id:
            medit_meshio_id[cell_block.type] = i

    # validate element counts
    if ref_num_triangle > 0:
        c = mesh.cells[medit_meshio_id["triangle"]]
        assert c.data.shape == (ref_num_triangle, 3)
    else:
        assert medit_meshio_id["triangle"] is None

    if ref_num_quad > 0:
        c = mesh.cells[medit_meshio_id["quad"]]
        assert c.data.shape == (ref_num_quad, 4)
    else:
        assert medit_meshio_id["quad"] is None

    if ref_num_tet > 0:
        c = mesh.cells[medit_meshio_id["tetra"]]
        assert c.data.shape == (ref_num_tet, 4)
    else:
        assert medit_meshio_id["tetra"] is None

    if ref_num_wedge > 0:
        c = mesh.cells[medit_meshio_id["wedge"]]
        assert c.data.shape == (ref_num_wedge, 6)
    else:
        assert medit_meshio_id["wedge"] is None

    if ref_num_hex > 0:
        c = mesh.cells[medit_meshio_id["hexahedron"]]
        assert c.data.shape == (ref_num_hex, 8)
    else:
        assert medit_meshio_id["hexahedron"] is None

    # validate boundary tags

    # gather tags
    all_tags = []
    for k, c in enumerate(mesh.cells):
        if c.type not in ["triangle", "quad"]:
            continue
        all_tags.append(mesh.cell_data["medit:ref"][k])

    all_tags = np.concatenate(all_tags)

    # validate against known values
    unique, counts = np.unique(all_tags, return_counts=True)
    tags = dict(zip(unique, counts))
    assert tags.keys() == ref_tag_counts.keys()
    for key in tags.keys():
        assert tags[key] == ref_tag_counts[key]


# --- malformed-input / error-path coverage (Python reference reader) ---
# These drive the reader's explicit `raise ReadError` branches directly, since
# the C++ core (tried first by the shim) may accept some malformed files.
from meshioplusplus.medit._medit import read as _medit_py_read  # noqa: E402


def test_medit_vertices_before_version_raises(tmp_path):
    p = tmp_path / "bad.mesh"
    p.write_text("Dimension 3\nVertices\n1\n0 0 0 0\n")
    with pytest.raises(meshioplusplus.ReadError):
        _medit_py_read(str(p))


def test_medit_non_alpha_keyword_raises(tmp_path):
    p = tmp_path / "bad.mesh"
    p.write_text("123 456\n")
    with pytest.raises(meshioplusplus.ReadError):
        _medit_py_read(str(p))


@pytest.mark.parametrize("engine", ["core", "python"])
def test_missing_dimension_is_inferred(engine):
    # FEconv writes no `Dimension` keyword (tools/gen_feconv_quirk_fixtures.py).
    from meshioplusplus import _core

    path = (
        pathlib.Path(__file__).resolve().parent
        / "meshes"
        / "medit"
        / "no_dimension.mesh"
    )
    mesh = (
        _core.medit_read_ascii(str(path)) if engine == "core" else _medit_py_read(path)
    )
    assert mesh.points.shape == (5, 3)
    np.testing.assert_array_equal(mesh.points[4], [1, 1, 1])
    tetra = [c for c in mesh.cells if c.type == "tetra"]
    np.testing.assert_array_equal(tetra[0].data, [[0, 1, 2, 3], [1, 4, 2, 3]])
    np.testing.assert_array_equal(mesh.cell_data["medit:ref"][0], [1, 2])


def test_medit_file_ending_inside_vertices_raises(tmp_path):
    # No `Dimension` and nothing after the vertex count: the blank-line skip
    # used to spin for ever on readline()'s "" at EOF.
    p = tmp_path / "truncated.mesh"
    p.write_text("MeshVersionFormatted 2\nVertices\n3\n\n")
    with pytest.raises(meshioplusplus.ReadError, match="ends inside"):
        _medit_py_read(str(p))


_KEYWORD_FILE = """MeshVersionFormatted 2
Dimension 3

Vertices
4
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 1

Edges
3
1 2 1
2 3 1
3 1 1

Triangles
1
1 2 3 1

Corners
2
1
3

RequiredVertices
1
2

Ridges
2
1
3

RequiredEdges
1
2

RequiredTriangles
1
1

End
"""


@pytest.mark.parametrize("engine", ["core", "python"])
def test_keyword_id_lists_are_regions_and_come_back(engine, tmp_path):
    """Roadmap 3.3.9: Corners, Ridges and the Required* lists are kept."""
    from meshioplusplus import _core
    from meshioplusplus.medit import _medit

    path = tmp_path / "keywords.mesh"
    path.write_text(_KEYWORD_FILE)
    read = _core.medit_read_ascii if engine == "core" else _medit.read
    mesh = read(str(path))
    regions = {r.name: (r.kind, list(r.entries)) for r in mesh.regions}
    assert regions == {
        "medit:Corners": ("point", [0, 2]),
        "medit:RequiredVertices": ("point", [1]),
        "medit:Ridges": ("cell", [0, 2]),
        "medit:RequiredEdges": ("cell", [1]),
        "medit:RequiredTriangles": ("cell", [3]),
    }

    # Written back (by either engine), they are read back the same by both.
    out = tmp_path / "again.mesh"
    write = _core.medit_write_ascii if engine == "core" else _medit.write
    write(str(out), mesh)
    for reader in (_core.medit_read_ascii, _medit.read):
        again = reader(str(out))
        assert {r.name: list(r.entries) for r in again.regions} == {
            n: e for n, (_, e) in regions.items()
        }


# -- binary, both engines ------------------------------------------------------

import glob  # noqa: E402

from meshioplusplus import _core  # noqa: E402
from meshioplusplus.medit import _medit as py_medit  # noqa: E402

BINARY = sorted(
    glob.glob(str(pathlib.Path(__file__).parent / "meshes" / "medit" / "*.meshb"))
)


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    np.testing.assert_array_equal(a.point_data["medit:ref"], b.point_data["medit:ref"])
    for x, y in zip(a.cell_data["medit:ref"], b.cell_data["medit:ref"]):
        np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("path", BINARY)
def test_binary_engines_read_and_write_the_same(tmp_path, path):
    core, twin = _core.medit_read_binary(path), py_medit.read(path)
    _same(core, twin)
    _core.medit_write_binary(str(tmp_path / "c.meshb"), core)
    py_medit.write_binary_file(str(tmp_path / "p.meshb"), core)
    assert (tmp_path / "c.meshb").read_bytes() == (tmp_path / "p.meshb").read_bytes()
    _same(meshioplusplus.read(tmp_path / "c.meshb"), core)


def _lists_mesh(refs):
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0.0]]),
        [
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("line", np.array([[0, 1], [1, 2]])),
        ],
    )
    mesh.point_data["medit:ref"] = np.array([1, 2, 3, 4])
    mesh.cell_data["medit:ref"] = [np.array([5, 6]), np.array(refs)]
    mesh.regions = [
        meshioplusplus.Region("medit:Corners", "point", np.array([0, 2])),
        meshioplusplus.Region("medit:Ridges", "cell", np.array([2, 3])),
        meshioplusplus.Region("medit:RequiredTriangles", "cell", np.array([1])),
    ]
    return mesh


@pytest.mark.parametrize("refs", [[7, 8], [7, 2**40]])
def test_binary_id_lists_and_wide_integers(tmp_path, refs):
    mesh = _lists_mesh(refs)
    _core.medit_write_binary(str(tmp_path / "c.meshb"), mesh)
    py_medit.write_binary_file(str(tmp_path / "p.meshb"), mesh)
    data = (tmp_path / "c.meshb").read_bytes()
    assert data == (tmp_path / "p.meshb").read_bytes()
    assert int.from_bytes(data[4:8], "little") == (4 if refs[1] > 2**31 else 3)
    for read in (_core.medit_read_binary, py_medit.read):
        back = read(str(tmp_path / "c.meshb"))
        assert {r.name: r.entries.tolist() for r in back.regions} == {
            "medit:Corners": [0, 2],
            "medit:Ridges": [2, 3],
            "medit:RequiredTriangles": [1],
        }
        np.testing.assert_array_equal(back.cell_data["medit:ref"][1], refs)


def test_binary_refusals(tmp_path):
    raw = tmp_path / "ok.meshb"
    _core.medit_write_binary(str(raw), _lists_mesh([7, 8]))
    data = raw.read_bytes()
    for cut in (10, 40, 100):
        (tmp_path / "t.meshb").write_bytes(data[:cut])
        for read in (_core.medit_read_binary, py_medit.read):
            with pytest.raises(meshioplusplus.ReadError):
                read(str(tmp_path / "t.meshb"))
    (tmp_path / "bad.meshb").write_bytes(b"\x07\x00\x00\x00" + data[4:])
    for read in (_core.medit_read_binary, py_medit.read):
        with pytest.raises(meshioplusplus.ReadError, match="nvalid code"):
            read(str(tmp_path / "bad.meshb"))
