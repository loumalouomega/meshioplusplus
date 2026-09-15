"""Gates for the generated documentation diagrams (``doc/diagrams/``).

The committed SVGs under ``doc/public/diagrams/`` are generated, not drawn,
so this pins three things: the generator still reproduces them byte for byte
(staleness), every figure has a PNG twin and a page that embeds it, and the
cell-type figures place every node where the code's own topology tables say
it belongs. Nothing here writes into ``doc/``.
"""

import importlib.util
import pathlib
import sys
import xml.etree.ElementTree as ET

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO / "doc" / "diagrams" / "gen_diagrams.py"
OUT = REPO / "doc" / "public" / "diagrams"


@pytest.fixture(scope="module")
def gen():
    spec = importlib.util.spec_from_file_location("gen_diagrams", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tables(gen):
    sys.path.insert(0, str(SCRIPT.parent))
    from diaglib.tables import Tables

    return Tables()


def _close(a, b, tol=1e-12):
    return all(abs(x - y) <= tol for x, y in zip(a, b))


def _mid(a, b):
    return tuple((x + y) / 2.0 for x, y in zip(a, b))


def _centroid(pts):
    n = float(len(pts))
    return tuple(sum(p[i] for p in pts) / n for i in range(3))


# --------------------------------------------------------------------------- #
# the committed outputs                                                        #
# --------------------------------------------------------------------------- #
def test_lists_every_committed_diagram(gen):
    committed = sorted(p.stem for p in OUT.glob("*.svg"))
    assert committed == sorted(gen.REGISTRY), "run doc/diagrams/gen_diagrams.py"


def test_generated_svgs_are_current(gen, tmp_path):
    gen.write_svgs(tmp_path)
    for name in gen.REGISTRY:
        fresh = (tmp_path / f"{name}.svg").read_text(encoding="utf-8")
        committed = (OUT / f"{name}.svg").read_text(encoding="utf-8")
        assert (
            fresh == committed
        ), f"{name}.svg is stale: run doc/diagrams/gen_diagrams.py"


def test_png_twins_exist_at_2x(gen):
    for name in gen.REGISTRY:
        svg = (OUT / f"{name}.svg").read_text(encoding="utf-8")
        png = OUT / f"{name}.png"
        assert png.exists(), png
        w, h = gen.viewbox_size(svg)
        assert gen.png_size(png) == (round(w * gen.SCALE), round(h * gen.SCALE)), png


def test_every_diagram_is_embedded(gen):
    used = gen.referenced_names()
    unused = sorted(set(gen.REGISTRY) - used)
    assert not unused, f"embedded by no page under doc/ nor README.md: {unused}"
    dangling = sorted(used - set(gen.REGISTRY))
    assert not dangling, f"referenced but produced by no figure: {dangling}"


def test_svgs_are_well_formed_and_on_paper(gen):
    from diaglib import palette

    for name in gen.REGISTRY:
        path = OUT / f"{name}.svg"
        root = ET.parse(path).getroot()
        assert root.get("viewBox", "").startswith("0 0 "), name
        children = [c for c in root if c.tag.split("}")[-1] not in ("title", "defs")]
        first = children[0]
        assert first.tag.endswith("rect") and first.get("fill") == palette.PAPER, name
        text = path.read_text(encoding="utf-8")
        import re

        colours = {c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}\b", text)}
        assert colours <= palette.ALL, f"{name}: {sorted(colours - palette.ALL)}"


def test_check_passes_on_the_committed_tree(gen):
    assert gen.check(OUT) == []


# --------------------------------------------------------------------------- #
# the cell-type figures derive from the code's tables                          #
# --------------------------------------------------------------------------- #
DRAWN_TYPES = [
    "vertex",
    "line",
    "line3",
    "triangle",
    "triangle6",
    "quad",
    "quad8",
    "quad9",
    "tetra",
    "tetra10",
    "hexahedron",
    "hexahedron20",
    "hexahedron24",
    "hexahedron27",
    "wedge",
    "wedge15",
    "wedge18",
    "pyramid",
    "pyramid13",
    "pyramid14",
]


@pytest.mark.parametrize("cell_type", DRAWN_TYPES)
def test_node_counts_match_num_nodes_per_cell(tables, cell_type):
    assert len(tables.node_coordinates(cell_type)) == tables.num_nodes[cell_type]


@pytest.mark.parametrize(
    "base", ["line", "triangle", "quad", "tetra", "hexahedron", "wedge", "pyramid"]
)
def test_elevate_mid_edge_nodes_are_edge_midpoints(tables, base):
    target, ncorners, edges = tables.elevate[base]
    nodes = tables.node_coordinates(target)
    for k, (a, b) in enumerate(edges):
        assert _close(nodes[ncorners + k], _mid(nodes[a], nodes[b])), (target, k)


@pytest.mark.parametrize(
    "cell_type",
    [
        t
        for t in DRAWN_TYPES
        if t
        in (
            "tetra10",
            "hexahedron20",
            "hexahedron27",
            "wedge15",
            "wedge18",
            "pyramid13",
            "pyramid14",
        )
    ],
)
def test_cell_faces_rows_place_their_quadratic_nodes_at_midpoints(tables, cell_type):
    nodes = tables.node_coordinates(cell_type)
    for out_type, nc, row in tables.cell_faces[cell_type]:
        corners = row[:nc]
        for k in range(nc):
            mid = row[nc + k]
            a, b = corners[k], corners[(k + 1) % nc]
            assert _close(nodes[mid], _mid(nodes[a], nodes[b])), (
                cell_type,
                out_type,
                row,
            )
        if out_type == "quad9":
            assert _close(nodes[row[8]], _centroid([nodes[i] for i in corners])), (
                cell_type,
                row,
            )


@pytest.mark.parametrize("cell_type", ["triangle6", "quad8", "quad9"])
def test_cell_edges_rows_place_their_mid_nodes_at_midpoints(tables, cell_type):
    nodes = tables.node_coordinates(cell_type)
    for _, _, (a, b, m) in tables.cell_edges[cell_type]:
        assert _close(nodes[m], _mid(nodes[a], nodes[b])), (cell_type, a, b, m)


def test_face_and_body_centres(tables):
    hexa = tables.node_coordinates("hexahedron27")
    for k, row in enumerate(tables.refine.QUAD_FACES["hexahedron"]):
        assert _close(hexa[20 + k], _centroid([hexa[i] for i in row])), (k, row)
    assert _close(hexa[26], _centroid(hexa[:8]))
    hexa24 = tables.node_coordinates("hexahedron24")
    assert hexa24[:24] == hexa[:24]
    wedge = tables.node_coordinates("wedge18")
    for k, row in enumerate(tables.refine.QUAD_FACES["wedge"]):
        assert _close(wedge[15 + k], _centroid([wedge[i] for i in row])), (k, row)
    quad = tables.node_coordinates("quad9")
    assert _close(quad[8], _centroid(quad[:4]))
    pyr = tables.node_coordinates("pyramid14")
    assert _close(pyr[13], _centroid(pyr[:4]))
    assert pyr[:13] == tables.node_coordinates("pyramid13")


@pytest.mark.parametrize("cell_type", ["triangle", "quad"])
def test_refine_children_tile_the_parent(tables, cell_type):
    layout = tables.node_coordinates(
        {"triangle": "triangle6", "quad": "quad9"}[cell_type]
    )

    def area(ring):
        s = 0.0
        for k in range(len(ring)):
            (x0, y0, _), (x1, y1, _) = (
                layout[ring[k]],
                layout[ring[(k + 1) % len(ring)]],
            )
            s += x0 * y1 - x1 * y0
        return s / 2.0

    parent = area(list(range(tables.num_nodes[cell_type])))
    for mask, (children, alt, _, _) in tables.refine.TABLES[cell_type].items():
        for variant in (children, alt):
            if not variant:
                continue
            areas = [area(child) for child in variant]
            assert all(a > 0 for a in areas), (cell_type, mask, variant)
            assert abs(sum(areas) - parent) < 1e-12, (cell_type, mask)


def test_tables_match_live_package(tables):
    skin = pytest.importorskip("meshioplusplus._skin")
    surface = pytest.importorskip("meshioplusplus._surface")
    convert_cells = pytest.importorskip("meshioplusplus._convert_cells")
    common = pytest.importorskip("meshioplusplus._common")
    assert tables.cell_faces == skin._CELL_FACES
    assert tables.cell_edges == surface._CELL_EDGES
    assert tables.elevate == convert_cells._ELEVATE
    assert tables.num_nodes == common.num_nodes_per_cell
