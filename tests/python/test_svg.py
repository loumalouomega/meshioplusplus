from xml.etree import ElementTree as ET

import pytest

import meshioplusplus
from meshioplusplus import _colormap

from . import helpers, helpers_coloring

SVG_NS = "{http://www.w3.org/2000/svg}"

test_set = [
    helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
]


def _drawable_cell_count(mesh):
    return sum(
        len(cb.data) for cb in mesh.cells if cb.type in ("line", "triangle", "quad")
    )


@pytest.mark.parametrize("mesh", test_set)
def test(mesh, tmp_path):
    filepath = tmp_path / "out.svg"
    meshioplusplus.write_points_cells(filepath, mesh.points, mesh.cells)
    # Output is valid SVG with one <path> per drawable cell.
    paths = ET.parse(filepath).getroot().findall(f"{SVG_NS}path")
    assert len(paths) == _drawable_cell_count(mesh)


@pytest.mark.parametrize("mesh", test_set)
def test_cpp_matches_python(mesh, tmp_path):
    # SVG is plain text: the C++ core writer must be byte-identical to the
    # pure-Python reference, as TikZ already is.
    cpp = tmp_path / "cpp.svg"
    py = tmp_path / "py.svg"
    meshioplusplus._core.svg_write(str(cpp), mesh)
    meshioplusplus.svg._svg.write(str(py), mesh)
    assert cpp.read_bytes() == py.read_bytes()

    n_cpp = len(ET.parse(cpp).getroot().findall(f"{SVG_NS}path"))
    assert n_cpp == _drawable_cell_count(mesh)


def _count_paths(filepath):
    return len(ET.parse(filepath).getroot().findall(f"{SVG_NS}path"))


def test_non_flat_shell_is_projected(tmp_path):
    # A non-flat 3D surface mesh is projected (isometric camera) and drawn.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.svg"
    meshioplusplus.svg.write(filepath, mesh)
    assert _count_paths(filepath) == 1


@pytest.mark.parametrize(
    "mesh, num_faces",
    [
        (helpers.tet_mesh, 6),
        (helpers.hex_mesh, 6),
        (helpers.wedge_mesh, 5),
        (helpers.pyramid_mesh, 5),
    ],
)
def test_volume_mesh_draws_skin(mesh, num_faces, tmp_path):
    import copy

    filepath = tmp_path / "skin.svg"
    meshioplusplus.svg.write(filepath, copy.deepcopy(mesh))
    assert _count_paths(filepath) == num_faces


def test_camera_angles_change_projection(tmp_path):
    import copy

    p1 = tmp_path / "iso.svg"
    p2 = tmp_path / "other.svg"
    meshioplusplus.svg.write(p1, copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.svg.write(
        p2, copy.deepcopy(helpers.hex_mesh), azimuth=10.0, elevation=60.0
    )
    assert _count_paths(p1) == _count_paths(p2) == 6
    assert p1.read_bytes() != p2.read_bytes()


def test_cpp_matches_python_3d(tmp_path):
    import copy

    # The byte-identity guarantee extends to the 3D projected path: skin
    # extraction, camera arithmetic, depth sort, and formatting must all
    # agree between the C++ core and the Python reference.
    cpp = tmp_path / "cpp.svg"
    py = tmp_path / "py.svg"
    meshioplusplus._core.svg_write(str(cpp), copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.svg._svg.write(str(py), copy.deepcopy(helpers.hex_mesh))
    assert cpp.read_bytes() == py.read_bytes()


# --- data-driven colouring ---------------------------------------------------


@pytest.mark.parametrize(
    "make_mesh, kwargs", helpers_coloring.CASE_ARGS, ids=helpers_coloring.CASE_IDS
)
def test_colored_cpp_matches_python(make_mesh, kwargs, tmp_path):
    # The byte-identity guarantee extends to the coloured path: array lookup,
    # scalarisation, the per-face mean, the auto range, the colormap index and
    # the colorbar geometry must all agree with the pure-Python reference.
    cpp = tmp_path / "cpp.svg"
    py = tmp_path / "py.svg"
    meshioplusplus._core.svg_write(str(cpp), make_mesh(), **kwargs)
    meshioplusplus.svg._svg.write(str(py), make_mesh(), **kwargs)
    assert cpp.read_bytes() == py.read_bytes()


@pytest.mark.parametrize(
    "make_mesh, kwargs", helpers_coloring.CASE_ARGS, ids=helpers_coloring.CASE_IDS
)
def test_colored_shim_matches_python(make_mesh, kwargs, tmp_path):
    # The shim passes everything POSITIONALLY into _core, so a mis-ordered
    # argument would silently colour by the wrong thing rather than fail.
    shim = tmp_path / "shim.svg"
    py = tmp_path / "py.svg"
    meshioplusplus.svg.write(shim, make_mesh(), **kwargs)
    meshioplusplus.svg._svg.write(str(py), make_mesh(), **kwargs)
    assert shim.read_bytes() == py.read_bytes()


def test_color_by_unset_is_unchanged(tmp_path):
    # The regression guard: passing the colouring parameters explicitly at
    # their defaults must reproduce the legacy output exactly.
    mesh = helpers_coloring.volume_mesh()
    plain = tmp_path / "plain.svg"
    explicit = tmp_path / "explicit.svg"
    meshioplusplus._core.svg_write(str(plain), mesh)
    meshioplusplus._core.svg_write(
        str(explicit),
        mesh,
        color_by="",
        component=None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        nan_color="#808080",
        colorbar=True,  # ignored: no colouring means no bar
    )
    assert plain.read_bytes() == explicit.read_bytes()


def _fills(path):
    return [p.get("fill") for p in ET.parse(path).getroot().findall(f"{SVG_NS}path")]


def test_colored_paths_carry_a_winning_inline_style(tmp_path):
    # Regression guard: a bare `fill="..."` attribute on a <path> has ZERO
    # CSS specificity and loses to ANY stylesheet rule, including the plain
    # "path {fill: ...}" type selector this writer always emits - every
    # cascade-honouring renderer (browsers, cairosvg) then paints every
    # coloured face in the flat fallback colour instead. The `fill`
    # attribute alone is therefore not enough; an inline `style` (which
    # always wins the cascade) must carry the same colour.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
        cell_data={"tag": [[1.0]]},
    )
    out = tmp_path / "styled.svg"
    meshioplusplus.svg.write(out, mesh, color_by="tag")
    paths = ET.parse(out).getroot().findall(f"{SVG_NS}path")
    assert len(paths) == 1
    fill = paths[0].get("fill")
    assert fill is not None
    assert paths[0].get("style") == f"fill:{fill}"


def test_cell_data_lands_on_the_right_facet(tmp_path):
    # A cube of two hexes tagged 0 and 1: every drawn skin facet must carry its
    # OWNING hex's colour, which is what "surface:parent_cell" provenance buys.
    import numpy as np

    points = [
        [x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0, 2.0)
    ]
    hexes = [
        [0, 1, 4, 3, 6, 7, 10, 9],
        [1, 2, 5, 4, 7, 8, 11, 10],
    ]
    mesh = meshioplusplus.Mesh(
        points, [("hexahedron", hexes)], cell_data={"tag": [np.array([0.0, 1.0])]}
    )
    out = tmp_path / "cube.svg"
    meshioplusplus.svg.write(out, mesh, color_by="tag")

    lo = _colormap.colormap_lookup(_colormap.colormap_table("viridis"), 0.0)
    hi = _colormap.colormap_lookup(_colormap.colormap_table("viridis"), 1.0)
    lo_hex = "#%02x%02x%02x" % lo
    hi_hex = "#%02x%02x%02x" % hi
    fills = _fills(out)
    # Only the two endpoint colours appear, and each hex contributes 5 of the
    # 10 boundary faces (the shared interior face is not on the boundary).
    assert set(fills) == {lo_hex, hi_hex}
    assert fills.count(lo_hex) == 5
    assert fills.count(hi_hex) == 5


def test_point_data_averages_over_corners(tmp_path):
    # One triangle whose corner values are 0, 3, 6: the face value is the mean,
    # 3, which sits exactly at the middle of an explicit 0..6 range.
    import numpy as np

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
        point_data={"T": np.array([0.0, 3.0, 6.0])},
    )
    out = tmp_path / "tri.svg"
    meshioplusplus.svg.write(out, mesh, color_by="T", vmin=0.0, vmax=6.0)
    expected = _colormap.colormap_lookup(_colormap.colormap_table("viridis"), 0.5)
    assert _fills(out) == ["#%02x%02x%02x" % expected]


def test_auto_range_is_the_drawn_finite_range(tmp_path):
    # Auto vmin/vmax span the drawn faces' finite values, so the extremes land
    # on the colormap endpoints and the NaN face takes nan_color.
    import numpy as np

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2], [1, 4, 3]])],
        cell_data={"tag": [np.array([-7.0, np.nan, 11.0])]},
    )
    out = tmp_path / "range.svg"
    meshioplusplus.svg.write(out, mesh, color_by="tag", nan_color="#123456")
    table = _colormap.colormap_table("viridis")
    lo = "#%02x%02x%02x" % _colormap.colormap_lookup(table, 0.0)
    hi = "#%02x%02x%02x" % _colormap.colormap_lookup(table, 1.0)
    assert _fills(out) == [lo, "#123456", hi]


def test_clamping_outside_an_explicit_range(tmp_path):
    import numpy as np

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2]])],
        cell_data={"tag": [np.array([-100.0, 100.0])]},
    )
    out = tmp_path / "clamp.svg"
    meshioplusplus.svg.write(out, mesh, color_by="tag", vmin=0.0, vmax=1.0)
    table = _colormap.colormap_table("viridis")
    assert _fills(out) == [
        "#%02x%02x%02x" % _colormap.colormap_lookup(table, 0.0),
        "#%02x%02x%02x" % _colormap.colormap_lookup(table, 1.0),
    ]


def test_colorbar_widens_only_the_viewbox(tmp_path):
    mesh = helpers_coloring.flat_mesh()
    without = tmp_path / "without.svg"
    with_bar = tmp_path / "with.svg"
    meshioplusplus.svg.write(without, mesh, color_by="tag")
    meshioplusplus.svg.write(with_bar, mesh, color_by="tag", colorbar=True)

    root_a = ET.parse(without).getroot()
    root_b = ET.parse(with_bar).getroot()
    box_a = root_a.get("viewBox").split()
    box_b = root_b.get("viewBox").split()
    # x, y and height are untouched; only the width grows.
    assert (box_a[0], box_a[1], box_a[3]) == (box_b[0], box_b[1], box_b[3])
    assert float(box_b[2]) > float(box_a[2])
    # Every mesh path is byte-identical; the bar is pure append.
    assert [p.get("d") for p in root_a.findall(f"{SVG_NS}path")] == [
        p.get("d") for p in root_b.findall(f"{SVG_NS}path")
    ]
    assert len(root_b.findall(f"{SVG_NS}rect")) == 32
    assert len(root_b.findall(f"{SVG_NS}text")) == 2
    assert root_a.findall(f"{SVG_NS}rect") == []


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"color_by": "nope"}, "no point_data or cell_data array"),
        ({"color_by": "tag", "cmap": "nope"}, "unknown colormap"),
        ({"color_by": "T", "component": 9}, "out of range"),
        ({"color_by": "tag", "vmin": 5.0, "vmax": 1.0}, "vmin must not exceed vmax"),
    ],
)
def test_invalid_options_raise(kwargs, match, tmp_path):
    # The shim falls back to Python on a C++ exception, and the twin raises the
    # same error -- so a genuine user error still surfaces, from either side.
    mesh = helpers_coloring.flat_mesh()
    with pytest.raises(ValueError, match=match):
        meshioplusplus.svg.write(tmp_path / "bad.svg", mesh, **kwargs)
    with pytest.raises(Exception, match=match):
        meshioplusplus._core.svg_write(str(tmp_path / "bad.svg"), mesh, **kwargs)
