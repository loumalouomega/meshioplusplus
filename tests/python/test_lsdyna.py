import io
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus._fallback import set_strict_core
from meshioplusplus.lsdyna import _cards
from meshioplusplus.lsdyna import _lsdyna as py_lsdyna

from . import helpers


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.lsdyna
    return py_lsdyna


def fmt(values, widths):
    return "".join(f"{v:>{w}}" for v, w in zip(values, widths))


def std_node(nid, x, y, z):
    return fmt([nid, f"{x:.8e}", f"{y:.8e}", f"{z:.8e}"], [8, 16, 16, 16])


def long_node(nid, x, y, z):
    return fmt([nid, f"{x:.12e}", f"{y:.12e}", f"{z:.12e}"], [20, 20, 20, 20])


CUBE = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
    (0.5, 0.5, 2),
]


def cube_nodes(node_fn=std_node, ids=range(1, 10)):
    return [node_fn(i, *CUBE[i - 1]) for i in ids]


def write_deck(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def canon(mesh):
    """Everything two engines must agree on, region order aside."""
    return (
        mesh.points.dtype,
        mesh.points.tolist(),
        [(c.type, np.asarray(c.data).tolist()) for c in mesh.cells],
        sorted(
            (r.kind, r.name, r.dim, r.tag, r.entries.tolist()) for r in mesh.regions
        ),
    )


def region(mesh, name, kind=None):
    hits = [r for r in mesh.regions if r.name == name and kind in (None, r.kind)]
    assert len(hits) == 1, f"{name!r}: {[r.name for r in mesh.regions]}"
    return hits[0]


# A deck exercising most of the reader: every element family, both solid layouts, the
# degenerate solids, parts, and every set flavour.
MIXED = [
    "*KEYWORD",
    "$ a comment",
    "*TITLE",
    "mixed deck",
    "*NODE",
    *cube_nodes(),
    std_node(10, 2, 0, 0),
    std_node(11, 2, 1, 0),
    std_node(12, 0, 0, 2),
    "*ELEMENT_SOLID",
    fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [8] * 10),
    fmt([2, 1, 1, 2, 3, 4, 9, 9, 9, 9], [8] * 10),
    fmt([3, 2, 1, 2, 3, 3, 5, 6, 7, 7], [8] * 10),
    fmt([4, 2, 1, 2, 3, 9, 9, 9, 9, 9], [8] * 10),
    "*ELEMENT_SHELL",
    fmt([1, 3, 1, 2, 3, 3], [8] * 6),
    fmt([2, 3, 1, 2, 3, 4], [8] * 6),
    "*ELEMENT_BEAM",
    fmt([1, 4, 1, 2, 9], [8] * 5),
    "*ELEMENT_DISCRETE",
    fmt([2, 4, 2, 3], [8] * 4),
    "*ELEMENT_MASS",
    fmt([1, 5, "1.5", 6], [8, 8, 16, 8]),
    "*PART",
    "top",
    fmt([1, 1, 1], [10] * 3),
    "*PART",
    "",
    fmt([2, 1, 1], [10] * 3),
    "*PART",
    "skin",
    fmt([3, 1, 1], [10] * 3),
    "*PART",
    "rod",
    fmt([4, 1, 1], [10] * 3),
    "*PART",
    "lonely",
    fmt([5, 1, 1], [10] * 3),
    "*SET_NODE_LIST_TITLE",
    "fixed",
    fmt([10], [10]),
    fmt([1, 2, 3, 4], [10] * 4),
    "*SET_NODE_LIST_GENERATE",
    fmt([20, 0, 0, 0, 0, 0], [10] * 6),
    fmt([1, 3, 6, 8], [10] * 4),
    "*SET_SOLID_LIST",
    fmt([30], [10]),
    fmt([1, 3], [10] * 2),
    "*SET_SHELL_LIST_TITLE",
    "shells",
    fmt([30], [10]),
    fmt([2], [10]),
    "*SET_PART_LIST",
    fmt([40], [10]),
    fmt([1, 3], [10] * 2),
    "*SET_SEGMENT_TITLE",
    "faces",
    fmt([50], [10]),
    fmt([1, 2, 3, 4], [10] * 4),
    fmt([1, 2, 9, 9], [10] * 4),
    "*END",
    "*NODE",
    "999,ignored after END",
]


def test_reads_the_mixed_deck(engine, tmp_path):
    mesh = engine.read(write_deck(tmp_path / "mixed.k", MIXED))
    types = {c.type: np.asarray(c.data).tolist() for c in mesh.cells}
    # Cells are grouped by type in order of first appearance, file order within.
    assert [c.type for c in mesh.cells] == [
        "hexahedron",
        "pyramid",
        "wedge",
        "tetra",
        "triangle",
        "quad",
        "line",
        "vertex",
    ]
    assert types["hexahedron"] == [[0, 1, 2, 3, 4, 5, 6, 7]]
    assert types["pyramid"] == [[0, 1, 2, 3, 8]]
    assert types["wedge"] == [[0, 1, 2, 4, 5, 6]]
    assert types["tetra"] == [[0, 1, 2, 8]]
    assert types["triangle"] == [[0, 1, 2]]
    assert types["quad"] == [[0, 1, 2, 3]]
    assert types["line"] == [[0, 1], [1, 2]]  # a beam and a discrete element
    assert types["vertex"] == [[4]]
    assert len(mesh.points) == 12  # nothing after *END is read

    top = region(mesh, "top", "cell")
    assert (top.tag, top.dim, top.entries.tolist()) == (1, 3, [0, 1])
    # A blank title falls back to "Part <pid>".
    assert region(mesh, "Part 2", "cell").entries.tolist() == [2, 3]
    assert region(mesh, "skin", "cell").dim == 2
    assert region(mesh, "rod", "cell").entries.tolist() == [6, 7]
    lonely = region(mesh, "lonely", "cell")
    assert (lonely.tag, len(lonely.entries)) == (5, 0)  # an empty part is kept

    assert region(mesh, "fixed", "point").entries.tolist() == [0, 1, 2, 3]
    assert region(mesh, "fixed", "point").tag == 10
    assert region(mesh, "NODE set 20", "point").entries.tolist() == [0, 1, 2, 5, 6, 7]
    # Solid 1 is the hexahedron, solid 3 the wedge; a shell set of the same id 30
    # is a different region.
    assert region(mesh, "SOLID set 30", "cell").entries.tolist() == [0, 2]
    assert region(mesh, "shells", "cell").entries.tolist() == [5]
    assert region(mesh, "shells", "cell").tag == 30
    # A part set is the union of its parts' elements.
    assert region(mesh, "PART set 40", "cell").entries.tolist() == [0, 1, 4, 5]
    faces = region(mesh, "faces", "side")
    assert faces.tag == 50
    # Segment 1 is the hexahedron's bottom face (facet 4); segment 2 is a
    # triangular face of the pyramid (cell 1), its facet 1.
    assert faces.entries.tolist() == [[0, 4], [1, 1]]


def test_engines_agree(tmp_path):
    path = write_deck(tmp_path / "mixed.k", MIXED)
    assert canon(meshioplusplus.lsdyna.read(path)) == canon(py_lsdyna.read(path))


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.quad_mesh,
        helpers.tri_quad_mesh,
        helpers.line_tri_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
    ],
)
def test_round_trip(engine, tmp_path, mesh):
    helpers.write_read(tmp_path, engine.write, engine.read, mesh, 1e-8, ".k")


@pytest.mark.parametrize("ext", [".k", ".key", ".dyn"])
def test_extensions_are_registered(tmp_path, ext):
    path = tmp_path / f"cube{ext}"
    meshioplusplus.write(path, helpers.hex_mesh)
    assert meshioplusplus.read(path).cells[0].type == "hexahedron"
    assert "lsdyna" in meshioplusplus.extension_to_filetypes[ext]
    assert "lsdyna" in meshioplusplus.formats()["readable"]
    assert "lsdyna" in meshioplusplus.formats()["writable"]


def test_engines_write_the_same_bytes(tmp_path):
    mesh = py_lsdyna.read(write_deck(tmp_path / "mixed.k", MIXED))
    a, b = tmp_path / "a.k", tmp_path / "b.k"
    meshioplusplus.lsdyna.write(a, mesh)
    py_lsdyna.write(b, mesh)
    assert a.read_bytes() == b.read_bytes()
    assert b"Written by meshio++" in a.read_bytes()


def test_round_trip_keeps_regions(engine, tmp_path):
    mesh = py_lsdyna.read(write_deck(tmp_path / "mixed.k", MIXED))
    out = tmp_path / "out.k"
    engine.write(out, mesh)
    back = engine.read(out)
    for r in mesh.regions:
        if r.entries.size == 0:
            continue  # an empty part comes back as an empty set
        twins = [
            q
            for q in back.regions
            if q.kind == r.kind and q.tag == r.tag and q.name.startswith(r.name)
        ]
        assert twins, r
        # A region spanning several element families comes back as one set each.
        assert sum(len(q.entries) for q in twins) == len(r.entries), r.name

    # Cells are regrouped by type, so compare the geometry cell by cell.
    def cells_of(m):
        return sorted(
            (c.type, tuple(sorted(row))) for c in m.cells for row in c.data.tolist()
        )

    assert cells_of(back) == cells_of(mesh)


def test_coordinates_survive(engine, tmp_path):
    rng = np.random.default_rng(3)
    pts = np.vstack(
        [
            rng.uniform(-1e3, 1e3, (20, 3)),
            [[0.1, 0.25, -0.5], [1e-30, -1e30, 3.0], [1e100, -1e-100, 0.0]],
        ]
    )
    mesh = meshioplusplus.Mesh(pts, [("vertex", np.arange(len(pts))[:, None])])
    engine.write(tmp_path / "p.k", mesh)
    back = engine.read(tmp_path / "p.k")
    assert np.allclose(back.points, pts, rtol=1e-9, atol=0.0)
    # Short decimals are exact, not rounded to nine digits.
    assert back.points[-3].tolist() == [0.1, 0.25, -0.5]


class TestWedgeOrientation:
    """A wedge collapsed from a hexahedron must have a positive volume."""

    def wedge_volume(self, p):
        p = np.asarray(p, dtype=float)
        n = np.cross(p[1] - p[0], p[2] - p[0])
        return n @ (p[3:].mean(0) - p[:3].mean(0)) / 2

    def read_cell(self, engine, tmp_path, coords, ids, kind="wedge"):
        lines = ["*NODE"] + [std_node(i + 1, *c) for i, c in enumerate(coords)]
        lines += ["*ELEMENT_SOLID", fmt([1, 1, *ids], [8] * 10)]
        mesh = engine.read(write_deck(tmp_path / f"{kind}.k", lines))
        assert [c.type for c in mesh.cells] == [kind]
        return mesh, np.asarray(mesh.cells[0].data[0])

    def test_n3_n3_form(self, engine, tmp_path):
        # n1 n2 n3 n3 n5 n6 n7 n7
        coords = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
        mesh, cell = self.read_cell(engine, tmp_path, coords, [1, 2, 3, 3, 4, 5, 6, 6])
        assert self.wedge_volume(mesh.points[cell]) == pytest.approx(0.5)

    def test_n5_n5_form(self, engine, tmp_path):
        # n1 n2 n3 n4 n5 n5 n6 n6: triangles (n1 n2 n5) and (n4 n3 n6)
        coords = [
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0.5, 0, 1),
            (0.5, 1, 1),
        ]
        mesh, cell = self.read_cell(engine, tmp_path, coords, [1, 2, 3, 4, 5, 5, 6, 6])
        assert self.wedge_volume(mesh.points[cell]) == pytest.approx(0.5)

    def test_written_wedge_reads_back_positive(self, engine, tmp_path):
        engine.write(tmp_path / "w.k", helpers.wedge_mesh)
        back = engine.read(tmp_path / "w.k")
        for row in back.cells[0].data:
            assert self.wedge_volume(back.points[row]) > 0


class TestCollapse:
    def cell_of(self, engine, tmp_path, ids):
        lines = ["*NODE"] + [std_node(i, i, i * i, 0) for i in range(1, 9)]
        lines += ["*ELEMENT_SOLID", fmt([1, 1, *ids], [8] * 10)]
        mesh = engine.read(write_deck(tmp_path / "c.k", lines))
        return mesh.cells[0].type, mesh.cells[0].data[0].tolist()

    def test_tetra(self, engine, tmp_path):
        assert self.cell_of(engine, tmp_path, [1, 2, 3, 4, 4, 4, 4, 4]) == (
            "tetra",
            [0, 1, 2, 3],
        )

    def test_tetra_with_repeated_third_node(self, engine, tmp_path):
        assert self.cell_of(engine, tmp_path, [1, 2, 3, 3, 4, 4, 4, 4]) == (
            "tetra",
            [0, 1, 2, 3],
        )

    def test_pyramid(self, engine, tmp_path):
        assert self.cell_of(engine, tmp_path, [1, 2, 3, 4, 5, 5, 5, 5]) == (
            "pyramid",
            [0, 1, 2, 3, 4],
        )

    def test_hexahedron_is_untouched(self, engine, tmp_path):
        assert self.cell_of(engine, tmp_path, [1, 2, 3, 4, 5, 6, 7, 8]) == (
            "hexahedron",
            [0, 1, 2, 3, 4, 5, 6, 7],
        )

    def test_triangle_shell_with_blank_fourth_node(self, engine, tmp_path):
        lines = ["*NODE"] + [std_node(i, i, 0, 0) for i in range(1, 5)]
        lines += ["*ELEMENT_SHELL", fmt([1, 1, 1, 2, 3], [8] * 5)]
        mesh = engine.read(write_deck(tmp_path / "s.k", lines))
        assert mesh.cells[0].type == "triangle"

    def test_written_solids_are_degenerate_hexahedra(self, tmp_path):
        for mesh, repeat in [(helpers.tet_mesh, 3), (helpers.pyramid_mesh, 4)]:
            py_lsdyna.write(tmp_path / "d.k", mesh)
            text = (tmp_path / "d.k").read_text().splitlines()
            first = text[text.index("*ELEMENT_SOLID") + 1]
            row = [int(v) + 1 for v in mesh.cells[0].data[0]]
            apex = row[repeat]
            assert first.split()[2:] == [
                str(v) for v in row[:repeat] + [apex] * (8 - repeat)
            ]


class TestCardFormats:
    def test_mixed_formats_in_one_deck(self, engine, tmp_path):
        lines = [
            "*KEYWORD",
            "*NODE",
            *cube_nodes(std_node, range(1, 5)),
            "*NODE+",  # long: every field 20 columns
            *cube_nodes(long_node, range(5, 9)),
            "*NODE",
            "9, 0.5, 0.5, 2.0",  # free format, per line
            "*ELEMENT_SOLID",
            fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [8] * 10),
            "*ELEMENT_SOLID%",  # i10: the 8-column integers widen to 10
            fmt([2, 1, 1, 2, 3, 4, 9, 9, 9, 9], [10] * 10),
            "*ELEMENT_SOLID",
            "3,2,1,2,3,3,5,6,7,7",
        ]
        mesh = engine.read(write_deck(tmp_path / "m.k", lines))
        assert len(mesh.points) == 9
        assert mesh.points[8].tolist() == [0.5, 0.5, 2.0]
        assert mesh.points[4].tolist() == [0.0, 0.0, 1.0]
        assert [c.type for c in mesh.cells] == ["hexahedron", "pyramid", "wedge"]

    def test_deck_wide_long(self, engine, tmp_path):
        lines = [
            "*KEYWORD LONG=Y",
            "*NODE",
            *cube_nodes(long_node),
            "*ELEMENT_SOLID",
            fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [20] * 10),
            "*PART",
            "wide",
            fmt([7, 1, 1], [20] * 3),
            "*SET_NODE_LIST",
            fmt([3], [20]),
            fmt([1, 2], [20] * 2),
        ]
        mesh = engine.read(write_deck(tmp_path / "long.k", lines))
        assert region(mesh, "wide").tag == 7
        assert region(mesh, "NODE set 3").entries.tolist() == [0, 1]
        assert mesh.points[8].tolist() == [0.5, 0.5, 2.0]

    def test_deck_wide_i10(self, engine, tmp_path):
        node = lambda i, x, y, z: fmt(  # noqa: E731
            [i, f"{x:.8e}", f"{y:.8e}", f"{z:.8e}"], [10, 16, 16, 16]
        )
        lines = [
            "*KEYWORD I10=Y",
            "*NODE",
            *cube_nodes(node),
            "*ELEMENT_SOLID",
            fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [10] * 10),
        ]
        mesh = engine.read(write_deck(tmp_path / "i10.k", lines))
        assert mesh.cells[0].data.tolist() == [[0, 1, 2, 3, 4, 5, 6, 7]]
        assert mesh.points[8].tolist() == [0.5, 0.5, 2.0]

    def test_long_suffix_is_undone_by_minus(self, engine, tmp_path):
        lines = [
            "*KEYWORD LONG=Y",
            "*NODE-",
            *cube_nodes(std_node),
        ]
        mesh = engine.read(write_deck(tmp_path / "minus.k", lines))
        assert mesh.points[8].tolist() == [0.5, 0.5, 2.0]

    def test_fortran_real_spellings(self, engine, tmp_path):
        lines = [
            "*NODE",
            fmt([1, "1.5D+01", "2.5-1", "-3.0E0"], [8, 16, 16, 16]),
        ]
        lines += ["*ELEMENT_MASS", fmt([1, 1, "0.0", 1], [8, 8, 16, 8])]
        mesh = engine.read(write_deck(tmp_path / "f.k", lines))
        assert mesh.points[0].tolist() == [15.0, 0.25, -3.0]

    def test_lowercase_keywords(self, engine, tmp_path):
        lines = [
            "*node",
            *cube_nodes()[:4],
            "*element_shell",
            fmt([1, 1, 1, 2, 3, 4], [8] * 6),
        ]
        mesh = engine.read(write_deck(tmp_path / "lc.k", lines))
        assert mesh.cells[0].type == "quad"

    def test_two_line_tetra10(self, engine, tmp_path):
        lines = ["*NODE"] + [std_node(i, i, i * i % 7, 0) for i in range(1, 11)]
        lines += [
            "*ELEMENT_SOLID",
            fmt([1, 1], [8, 8]),
            fmt(range(1, 11), [8] * 10),
        ]
        mesh = engine.read(write_deck(tmp_path / "t10.k", lines))
        assert mesh.cells[0].type == "tetra10"
        assert mesh.cells[0].data.tolist() == [list(range(10))]

    def test_variant_cards_with_extra_lines(self, engine, tmp_path):
        lines = [
            "*NODE",
            *cube_nodes(),
            "*ELEMENT_SHELL_THICKNESS",
            fmt([1, 1, 1, 2, 3, 4], [8] * 6),
            fmt(["1.0", "1.0", "1.0", "1.0"], [16] * 4),
            fmt([2, 1, 5, 6, 7, 8], [8] * 6),
            fmt(["2.0", "2.0", "2.0", "2.0"], [16] * 4),
            "*ELEMENT_BEAM",
            fmt([3, 1, 1, 5, 9], [8] * 5),
            fmt(["1.0", "1.0"], [16, 16]),  # an optional second card
            fmt([4, 1, 2, 6, 9], [8] * 5),
        ]
        mesh = engine.read(write_deck(tmp_path / "v.k", lines))
        counts = {c.type: len(c.data) for c in mesh.cells}
        assert counts == {"quad": 2, "line": 2}


class TestIncludes:
    def deck(self, tmp_path):
        write_deck(
            tmp_path / "main.k",
            [
                "*KEYWORD",
                "*INCLUDE",
                "sub/nodes.k",
                "*INCLUDE",
                "sub/elements.k",
                "*SET_SOLID_LIST_TITLE",
                "everything",
                fmt([1], [10]),
                fmt([1, 2], [10] * 2),
                "*SET_NODE_LIST",
                fmt([2], [10]),
                fmt([9], [10]),
                "*END",
            ],
        )
        write_deck(tmp_path / "sub" / "nodes.k", ["*NODE", *cube_nodes()])
        write_deck(
            tmp_path / "sub" / "elements.k",
            [
                "*ELEMENT_SOLID",
                fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [8] * 10),
                fmt([2, 1, 1, 2, 3, 4, 9, 9, 9, 9], [8] * 10),
                "*PART",
                "block",
                fmt([1, 1, 1], [10] * 3),
            ],
        )
        return tmp_path / "main.k"

    def test_sets_may_refer_to_included_elements(self, engine, tmp_path):
        mesh = engine.read(self.deck(tmp_path))
        assert len(mesh.points) == 9
        assert region(mesh, "everything").entries.tolist() == [0, 1]
        assert region(mesh, "block").entries.tolist() == [0, 1]
        assert region(mesh, "NODE set 2").entries.tolist() == [8]

    def test_include_path(self, engine, tmp_path):
        write_deck(tmp_path / "lib" / "nodes.k", ["*NODE", *cube_nodes()])
        main = write_deck(
            tmp_path / "main.k",
            ["*INCLUDE_PATH", "lib", "*INCLUDE", "nodes.k"],
        )
        assert len(engine.read(main).points) == 9

    def test_include_path_relative(self, engine, tmp_path):
        write_deck(tmp_path / "a" / "lib" / "nodes.k", ["*NODE", *cube_nodes()])
        write_deck(
            tmp_path / "a" / "main.k",
            ["*INCLUDE_PATH_RELATIVE", "lib", "*INCLUDE", "nodes.k"],
        )
        assert len(engine.read(tmp_path / "a" / "main.k").points) == 9

    def test_included_file_resolves_relative_to_its_parent(self, engine, tmp_path):
        write_deck(tmp_path / "d" / "nodes.k", ["*NODE", *cube_nodes()])
        write_deck(tmp_path / "d" / "mid.k", ["*INCLUDE", "nodes.k"])
        write_deck(tmp_path / "top.k", ["*INCLUDE", "d/mid.k"])
        assert len(engine.read(tmp_path / "top.k").points) == 9

    def test_windows_separators(self, engine, tmp_path):
        write_deck(tmp_path / "d" / "nodes.k", ["*NODE", *cube_nodes()])
        write_deck(tmp_path / "top.k", ["*INCLUDE", "d\\nodes.k"])
        assert len(engine.read(tmp_path / "top.k").points) == 9

    def test_continued_file_name(self, engine, tmp_path):
        write_deck(tmp_path / "longname.k", ["*NODE", *cube_nodes()])
        write_deck(tmp_path / "top.k", ["*INCLUDE", "long+", "name.k"])
        assert len(engine.read(tmp_path / "top.k").points) == 9

    def test_recursion_is_bounded(self, engine, tmp_path):
        write_deck(tmp_path / "loop.k", ["*INCLUDE", "loop.k"])
        with pytest.raises(ReadError, match="deeper"):
            engine.read(tmp_path / "loop.k")

    def test_missing_include_is_a_warning(self, engine, tmp_path, capfd):
        write_deck(
            tmp_path / "main.k",
            ["*INCLUDE", "materials.k", "*NODE", *cube_nodes()],
        )
        assert len(engine.read(tmp_path / "main.k").points) == 9
        assert "materials.k" in capfd.readouterr().err

    def test_transform_is_ignored_with_a_warning(self, engine, tmp_path, capfd):
        write_deck(tmp_path / "n.k", ["*NODE", *cube_nodes()])
        write_deck(
            tmp_path / "main.k",
            ["*INCLUDE_TRANSFORM", "n.k", fmt(["1.0", "1"], [10, 10])],
        )
        mesh = engine.read(tmp_path / "main.k")
        assert mesh.points[8].tolist() == [0.5, 0.5, 2.0]
        assert "INCLUDE_TRANSFORM" in capfd.readouterr().err


class TestSets:
    def base(self):
        return ["*NODE", *cube_nodes(), "*ELEMENT_SOLID"] + [
            fmt([1, 1, 1, 2, 3, 4, 5, 6, 7, 8], [8] * 10)
        ]

    def test_generate_ranges(self, engine, tmp_path):
        lines = self.base() + [
            "*SET_NODE_LIST_GENERATE",
            fmt([1], [10]),
            fmt([1, 3, 7, 9], [10] * 4),
        ]
        mesh = engine.read(write_deck(tmp_path / "g.k", lines))
        assert region(mesh, "NODE set 1").entries.tolist() == [0, 1, 2, 6, 7, 8]

    def test_free_format_set(self, engine, tmp_path):
        lines = self.base() + ["*SET_NODE_LIST", "4", "1,2,3,4,5,6,7,8,9,1"]
        mesh = engine.read(write_deck(tmp_path / "g.k", lines))
        assert len(region(mesh, "NODE set 4").entries) == 9

    def test_same_title_and_id_in_two_families_do_not_collide(self, engine, tmp_path):
        lines = self.base() + [
            "*ELEMENT_SHELL",
            fmt([1, 1, 1, 2, 3, 4], [8] * 6),
            "*SET_SOLID_LIST_TITLE",
            "same",
            fmt([1], [10]),
            fmt([1], [10]),
            "*SET_SHELL_LIST_TITLE",
            "same",
            fmt([1], [10]),
            fmt([1], [10]),
        ]
        mesh = engine.read(write_deck(tmp_path / "d.k", lines))
        names = sorted(r.name for r in mesh.regions)
        assert names == ["same", "same [shell]"]

    def test_segment_sets_match_faces(self, engine, tmp_path):
        lines = self.base() + [
            "*SET_SEGMENT",
            fmt([1], [10]),
            fmt([1, 2, 3, 4], [10] * 4),  # bottom quad -> facet 4
            fmt([1, 2, 6, 5], [10] * 4),  # front quad -> facet 2
            fmt([1, 2, 3, 3], [10] * 4),  # not a face of the hexahedron
        ]
        mesh = engine.read(write_deck(tmp_path / "s.k", lines))
        assert region(mesh, "SEGMENT set 1").entries.tolist() == [[0, 2], [0, 4]]

    def test_dangling_ids_are_dropped_with_a_warning(self, engine, tmp_path, capfd):
        lines = self.base() + ["*SET_NODE_LIST", fmt([1], [10]), fmt([1, 99], [10] * 2)]
        mesh = engine.read(write_deck(tmp_path / "x.k", lines))
        assert region(mesh, "NODE set 1").entries.tolist() == [0]
        assert "undefined ids" in capfd.readouterr().err

    def test_unsupported_set_variant_is_a_warning(self, engine, tmp_path, capfd):
        lines = self.base() + ["*SET_NODE_COLUMN", fmt([1], [10])]
        mesh = engine.read(write_deck(tmp_path / "x.k", lines))
        assert not [r for r in mesh.regions if r.kind == "point"]
        assert "SET_NODE_COLUMN" in capfd.readouterr().err


class TestSkippedContent:
    def test_parameter_references_are_skipped_with_a_warning(
        self, engine, tmp_path, capfd
    ):
        lines = ["*NODE", *cube_nodes(), "&loadnode,1.0,2.0,3.0"]
        mesh = engine.read(write_deck(tmp_path / "p.k", lines))
        assert len(mesh.points) == 9
        assert "PARAMETER" in capfd.readouterr().err

    def test_pgp_blocks_are_skipped(self, engine, tmp_path, capfd):
        lines = [
            "*NODE",
            *cube_nodes(),
            "-----BEGIN PGP MESSAGE-----",
            "garbage that is not a card",
            "-----END PGP MESSAGE-----",
            "*ELEMENT_SHELL",
            fmt([1, 1, 1, 2, 3, 4], [8] * 6),
        ]
        mesh = engine.read(write_deck(tmp_path / "pgp.k", lines))
        assert mesh.cells[0].type == "quad"
        assert "PGP" in capfd.readouterr().err

    def test_conversion_directives_warn(self, engine, tmp_path, capfd):
        lines = ["*NODE", *cube_nodes(), "*ELEMENT_SOLID_TET4TOTET10", "1,1,1,2,3,4"]
        mesh = engine.read(write_deck(tmp_path / "c.k", lines))
        assert not mesh.cells
        assert "TET4TOTET10" in capfd.readouterr().err

    def test_unrelated_keywords_are_ignored(self, engine, tmp_path):
        lines = [
            "*CONTROL_TERMINATION",
            fmt(["1.0"], [10]),
            "*MAT_ELASTIC",
            fmt([1, "7.8e-6", "210.0", "0.3"], [10] * 4),
            "*NODE",
            *cube_nodes(),
            "*DEFINE_CURVE",
            "1",
            "0.0,1.0",
        ]
        mesh = engine.read(write_deck(tmp_path / "u.k", lines))
        assert len(mesh.points) == 9


class TestErrors:
    def test_undefined_node(self, engine, tmp_path):
        lines = [
            "*NODE",
            *cube_nodes(),
            "*ELEMENT_SHELL",
            fmt([1, 1, 1, 2, 3, 77], [8] * 6),
        ]
        with pytest.raises(ReadError, match="undefined node 77"):
            engine.read(write_deck(tmp_path / "e.k", lines))

    def test_duplicate_node(self, engine, tmp_path):
        lines = ["*NODE", *cube_nodes(), std_node(3, 9, 9, 9)]
        with pytest.raises(ReadError, match="duplicate node id 3"):
            engine.read(write_deck(tmp_path / "e.k", lines))

    def test_duplicate_element_id_within_a_family(self, engine, tmp_path):
        lines = ["*NODE", *cube_nodes(), "*ELEMENT_SHELL"]
        lines += [fmt([1, 1, 1, 2, 3, 4], [8] * 6)] * 2
        with pytest.raises(ReadError, match="duplicate shell element id 1"):
            engine.read(write_deck(tmp_path / "e.k", lines))

    def test_the_same_id_in_two_families_is_fine(self, engine, tmp_path):
        lines = ["*NODE", *cube_nodes(), "*ELEMENT_SHELL"]
        lines += [
            fmt([1, 1, 1, 2, 3, 4], [8] * 6),
            "*ELEMENT_BEAM",
            fmt([1, 1, 1, 2], [8] * 4),
        ]
        assert len(engine.read(write_deck(tmp_path / "e.k", lines)).cells) == 2

    def test_garbage_in_a_numeric_field(self, engine, tmp_path):
        lines = ["*NODE", fmt([1, "abc", "0.0", "0.0"], [8, 16, 16, 16])]
        with pytest.raises(ReadError, match="invalid real field"):
            engine.read(write_deck(tmp_path / "e.k", lines))

    def test_missing_file(self, engine, tmp_path):
        with pytest.raises(ReadError):
            engine.read(tmp_path / "nope.k")

    def test_unsupported_write_type(self, engine, tmp_path):
        with pytest.raises(WriteError, match="triangle6"):
            engine.write(tmp_path / "t6.k", helpers.triangle6_mesh)


class TestWriter:
    def test_parts_from_dimensioned_regions(self, engine, tmp_path):
        mesh = meshioplusplus.Mesh(
            helpers.tri_quad_mesh.points,
            helpers.tri_quad_mesh.cells,
            regions=[
                meshioplusplus._regions.Region("tris", "cell", [0, 1, 3], dim=2, tag=7),
                meshioplusplus._regions.Region("quad", "cell", [2], dim=2, tag=7),
            ],
        )
        engine.write(tmp_path / "p.k", mesh)
        back = py_lsdyna.read(tmp_path / "p.k")
        tags = {r.name: r.tag for r in back.regions if r.kind == "cell"}
        # A repeated tag cannot be a second pid: the later part is renumbered.
        assert sorted(tags) == ["quad", "tris"]
        assert len(set(tags.values())) == 2

    def test_cells_no_part_claims_get_one_part_per_block(self, engine, tmp_path):
        engine.write(tmp_path / "b.k", helpers.tri_quad_mesh)
        back = py_lsdyna.read(tmp_path / "b.k")
        names = sorted(r.name for r in back.regions)
        assert names == sorted(c.type for c in helpers.tri_quad_mesh.cells)

    def test_overlapping_regions_are_sets_not_lost(self, engine, tmp_path):
        R = meshioplusplus._regions.Region
        mesh = meshioplusplus.Mesh(
            helpers.hex_mesh.points,
            helpers.hex_mesh.cells,
            regions=[
                R("a", "cell", [0], dim=3, tag=1),
                R("b", "cell", [0], dim=3, tag=2),
            ],
        )
        engine.write(tmp_path / "o.k", mesh)
        back = py_lsdyna.read(tmp_path / "o.k")
        assert sorted(r.name for r in back.regions) == ["a", "b"]
        assert {r.name: r.dim for r in back.regions} == {"a": 3, "b": -1}

    def test_titles_that_look_like_keywords(self, engine, tmp_path):
        R = meshioplusplus._regions.Region
        mesh = meshioplusplus.Mesh(
            helpers.hex_mesh.points,
            helpers.hex_mesh.cells,
            regions=[R("*NODE", "cell", [0], dim=3, tag=1)],
        )
        engine.write(tmp_path / "k.k", mesh)
        back = py_lsdyna.read(tmp_path / "k.k")
        assert [r.name for r in back.regions] == ["*NODE"]
        assert len(back.points) == len(helpers.hex_mesh.points)

    def test_side_regions_use_face_nodes(self, engine, tmp_path):
        R = meshioplusplus._regions.Region
        mesh = meshioplusplus.Mesh(
            helpers.hex_mesh.points,
            helpers.hex_mesh.cells,
            regions=[R("top", "side", [[0, 5]], tag=3)],
        )
        engine.write(tmp_path / "s.k", mesh)
        back = py_lsdyna.read(tmp_path / "s.k")
        side = [r for r in back.regions if r.kind == "side"][0]
        assert (side.name, side.tag, side.entries.tolist()) == ("top", 3, [[0, 5]])


class TestTokenizer:
    def test_widths_per_mode(self):
        i8 = ("i", 8)
        assert _cards.field_width(*i8, _cards.STD) == 8
        assert _cards.field_width(*i8, _cards.I10) == 10
        assert _cards.field_width(*i8, _cards.LONG) == 20
        assert _cards.field_width("r", 16, _cards.I10) == 16
        assert _cards.field_width("i", 10, _cards.I10) == 10

    def test_split_and_free_format(self):
        layout = _cards.ELEMENT
        assert _cards.split_card("       1       2", layout, _cards.STD)[:3] == [
            "1",
            "2",
            "",
        ]
        assert _cards.split_card("1, 2,,4", layout, _cards.STD)[:4] == [
            "1",
            "2",
            "",
            "4",
        ]

    @pytest.mark.parametrize(
        "x", [0.0, 1.0, -2.5, 0.1, 1 / 3, 123456.789, 1e-5, -1e-30, 1e100, 6.02e23]
    )
    def test_format_real16_fits_and_round_trips_when_it_can(self, x):
        s = _cards.format_real16(x)
        assert len(s) <= 16
        assert float(s) == pytest.approx(x, rel=2e-9)
        if len(repr(x)) <= 8:
            assert float(s) == x


class TestEngineWiring:
    def test_native_path_is_used(self, tmp_path):
        path = tmp_path / "cube.k"
        py_lsdyna.write(path, helpers.hex_mesh)
        set_strict_core(True)
        try:
            mesh = meshioplusplus.lsdyna.read(path)
            meshioplusplus.lsdyna.write(tmp_path / "again.k", mesh)
        finally:
            set_strict_core(None)
        assert mesh.cells[0].type == "hexahedron"

    def test_core_has_the_functions(self):
        for name in ("lsdyna_read", "lsdyna_write"):
            assert hasattr(_core, name)

    def test_buffers_use_the_python_reader(self):
        buf = io.StringIO()
        py_lsdyna.write(buf, helpers.tet_mesh)
        buf.seek(0)
        mesh = meshioplusplus.lsdyna.read(buf)
        assert mesh.cells[0].type == "tetra"

    def test_core_and_python_read_pathlib_and_str(self, tmp_path):
        py_lsdyna.write(tmp_path / "cube.k", helpers.hex_mesh)
        p = pathlib.Path(tmp_path / "cube.k")
        assert canon(meshioplusplus.lsdyna.read(p)) == canon(
            meshioplusplus.lsdyna.read(str(p))
        )


FIXTURES = pathlib.Path(__file__).parent / "meshes" / "lsdyna"


@pytest.mark.skipif(not FIXTURES.exists(), reason="fixtures not present")
class TestFixtureDecks:
    def test_multi_file_deck(self, engine, capfd):
        mesh = engine.read(FIXTURES / "crash_main.k")
        # Long-format grid nodes, free-format extras, I10 solids, free-format
        # shells, all pulled in through *INCLUDE / *INCLUDE_PATH.
        assert len(mesh.points) == 44
        counts = {c.type: len(c.data) for c in mesh.cells}
        assert counts == {
            "hexahedron": 8,
            "wedge": 4,
            "tetra": 2,
            "quad": 12,
            "line": 3,
        }
        parts = {r.name: r for r in mesh.regions if r.kind == "cell" and r.dim >= 0}
        assert {n: p.tag for n, p in parts.items()} == {
            "matrix": 1,
            "skin": 2,
            "plate": 3,
            "loose tets": 4,
            "rods": 5,
        }
        assert {n: len(p.entries) for n, p in parts.items()} == {
            "matrix": 8,
            "skin": 12,
            "plate": 4,
            "loose tets": 2,
            "rods": 3,
        }
        assert sum(len(p.entries) for p in parts.values()) == sum(counts.values())
        assert len(region(mesh, "NODE set 1").entries) == 9
        assert len(region(mesh, "bottom skin").entries) == 4
        assert len(region(mesh, "solid parts").entries) == 8 + 4 + 2
        assert region(mesh, "loaded face").entries.tolist() == [[0, 4]]
        # The one include the fixture deliberately does not ship is a warning.
        assert "materials_not_shipped.k" in capfd.readouterr().err

    def test_rewritten_fixture_reads_back(self, engine, tmp_path):
        mesh = engine.read(FIXTURES / "crash_main.k")
        engine.write(tmp_path / "out.k", mesh)
        back = engine.read(tmp_path / "out.k")
        assert {c.type: len(c.data) for c in back.cells} == {
            c.type: len(c.data) for c in mesh.cells
        }
        assert sorted(
            (r.name, r.kind, r.tag, len(r.entries)) for r in mesh.regions
        ) == sorted((r.name, r.kind, r.tag, len(r.entries)) for r in back.regions)

    def test_engines_agree_on_the_fixture(self):
        assert canon(meshioplusplus.lsdyna.read(FIXTURES / "crash_main.k")) == canon(
            py_lsdyna.read(FIXTURES / "crash_main.k")
        )
