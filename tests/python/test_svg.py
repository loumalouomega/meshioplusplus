import copy
from xml.etree import ElementTree as ET

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _colormap

from . import helpers, helpers_coloring, helpers_drawing

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


# --- roadmap §1 "writers that drop data without a warning" ---


def _mesh_with_unsupported_cell():
    # Flat (z=0) so the 3-D projection path isn't taken; a bare tetra block
    # is not one of SVG's line/triangle/quad primitives.
    return meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("tetra", [[0, 1, 2, 3]])],
    )


def test_warns_when_skipping_an_unsupported_cell_type(tmp_path, capfd):
    mesh = _mesh_with_unsupported_cell()
    meshioplusplus.svg.write(tmp_path / "skip.svg", mesh)
    err = capfd.readouterr().err
    assert "tetra" in err


def test_warns_when_skipping_python_engine(tmp_path, capsys):
    from meshioplusplus.svg._svg import write as _py_write

    mesh = _mesh_with_unsupported_cell()
    _py_write(str(tmp_path / "skip_py.svg"), mesh)
    err = capsys.readouterr().err
    assert "tetra" in err


# --- reading ------------------------------------------------------------------

_PLANAR_MESHES = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
    helpers.line_tri_mesh,
]


def _written(mesh, tmp_path, name="out.svg", **kwargs):
    path = tmp_path / name
    meshioplusplus.svg.write(path, copy.deepcopy(mesh), **kwargs)
    return path


@pytest.mark.parametrize("mesh", _PLANAR_MESHES)
def test_planar_round_trip(mesh, tmp_path):
    # Written, read through both engines, the drawing comes back in model units:
    # same cells, same corners (to the 3 decimals the file keeps), z = 0.
    path = _written(mesh, tmp_path)
    back = helpers_drawing.read_both("svg", path)
    assert back.points.shape[1] == 3
    assert not back.points[:, 2].any()
    expected = meshioplusplus.Mesh(
        copy.deepcopy(mesh.points), [c for c in mesh.cells if len(c.data)]
    )
    assert helpers_drawing.cell_coordinates(back) == helpers_drawing.cell_coordinates(
        expected
    )


def test_dyadic_planar_input_is_exact(tmp_path):
    # The scale (50 here) and the y-flip are restored from the marker at full
    # precision, so coordinates that survive the 3 decimals come back bit-exact.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [2.0, 0.0], [1.0, 0.5], [2.0, 1.0], [0.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2]]), ("quad", [[0, 2, 3, 4]])],
    )
    back = helpers_drawing.read_both("svg", _written(mesh, tmp_path))
    assert np.array_equal(back.points[:, :2], mesh.points)
    assert [c.type for c in back.cells] == ["triangle", "quad"]
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
    assert np.array_equal(back.cells[1].data, mesh.cells[1].data)


def test_shared_corners_come_back_shared(tmp_path):
    back = helpers_drawing.read_both("svg", _written(helpers.tri_mesh, tmp_path))
    assert len(back.points) == 4
    assert np.array_equal(back.cells[0].data, [[0, 1, 2], [0, 2, 3]])


def test_welding_is_on_the_exact_coordinate_text(tmp_path):
    path = tmp_path / "text.svg"
    path.write_text(
        helpers_drawing.marked_svg(
            '<path d="M 0 0L 1 0L 1 1Z" /><path d="M 0.0 0L 1 0L 1 1Z" />'
        )
    )
    back = helpers_drawing.read_both("svg", path)
    assert len(back.points) == 4  # "0" and "0.0" are two points
    assert [c.type for c in back.cells] == ["triangle"]


def test_a_polygon_path_is_a_polygon_cell(tmp_path):
    path = tmp_path / "poly.svg"
    path.write_text(
        helpers_drawing.marked_svg('<path d="M 0 0L 4 0L 5 2L 2 4L -1 2Z" />')
    )
    back = helpers_drawing.read_both("svg", path)
    assert [c.type for c in back.cells] == ["polygon"]
    assert back.cells[0].data.shape == (1, 5)


@pytest.mark.parametrize("image_width", [None, 40, 250.5])
def test_the_scale_is_restored_whatever_image_width_was(image_width, tmp_path):
    mesh = meshioplusplus.Mesh(
        [[10.0, -4.0], [14.0, -4.0], [14.0, -2.0], [10.0, -2.0]],
        [("quad", [[0, 1, 2, 3]])],
    )
    path = _written(mesh, tmp_path, image_width=image_width)
    back = helpers_drawing.read_both("svg", path)
    assert np.allclose(back.points[:, :2], mesh.points, atol=1e-3)


def test_a_projected_drawing_comes_back_flat(tmp_path, capfd):
    path = _written(helpers.hex_mesh, tmp_path)
    capfd.readouterr()
    core, _ = helpers_drawing.readers("svg")
    back = core(str(path))
    assert back.points.shape[1] == 3 and not back.points[:, 2].any()
    assert len(back.cells) >= 1
    err = capfd.readouterr().err
    assert err.count("projected from 3-D") == 1
    assert "depth lost" in err


def test_a_planar_drawing_warns_once(tmp_path, capfd):
    path = _written(helpers.tri_mesh, tmp_path)
    capfd.readouterr()
    core, py = helpers_drawing.readers("svg")
    core(str(path))
    assert capfd.readouterr().err.count("planar drawing, z = 0 assumed") == 1
    py(str(path))
    assert capfd.readouterr().err.count("planar drawing, z = 0 assumed") == 1


@pytest.mark.parametrize(
    "make_mesh, kwargs", helpers_coloring.CASE_ARGS, ids=helpers_coloring.CASE_IDS
)
def test_colored_drawings_read_the_same_in_both_engines(make_mesh, kwargs, tmp_path):
    path = _written(make_mesh(), tmp_path, **kwargs)
    back = helpers_drawing.read_both("svg", path)
    assert len(back.cells) >= 1
    # Fills come back as RGB, one row per cell, never as the data they came from.
    assert "svg:fill" in back.cell_data
    for block, fill in zip(back.cells, back.cell_data["svg:fill"]):
        assert fill.dtype == np.int32
        assert fill.shape == (len(block.data), 3)
        assert ((fill >= 0) & (fill <= 255)).all()
    assert not back.point_data and "tag" not in back.cell_data


def test_fills_match_the_file(tmp_path):
    path = _written(
        helpers_coloring.flat_mesh(), tmp_path, color_by="tag", colorbar=True
    )
    paths = ET.parse(path).getroot().findall(f"{SVG_NS}path")
    expected = [[int(p.get("fill")[i : i + 2], 16) for i in (1, 3, 5)] for p in paths]
    back = helpers_drawing.read_both("svg", path)
    assert len(back.points) == 4  # the 32 colorbar swatches added no points
    assert back.cell_data["svg:fill"][0].tolist() == expected


def test_a_cell_with_no_fill_is_minus_one(tmp_path):
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [("triangle", [[0, 1, 2]]), ("line", [[1, 3]])],
        cell_data={"tag": [np.array([0.0]), np.array([1.0])]},
    )
    back = helpers_drawing.read_both("svg", _written(mesh, tmp_path, color_by="tag"))
    triangle_fill, line_fill = back.cell_data["svg:fill"]
    assert (triangle_fill >= 0).all()
    assert line_fill.tolist() == [[-1, -1, -1]]


def test_a_flat_fill_gives_no_array(tmp_path):
    back = helpers_drawing.read_both("svg", _written(helpers.tri_mesh, tmp_path))
    assert not back.cell_data


def test_points_only_and_arrays_narrow_the_fill_away(tmp_path):
    path = _written(helpers_coloring.flat_mesh(), tmp_path, color_by="tag")
    assert "svg:fill" in helpers_drawing.read_both("svg", path).cell_data
    assert not helpers_drawing.read_both("svg", path, points_only=True).cell_data
    assert not helpers_drawing.read_both("svg", path, arrays=[]).cell_data
    assert not helpers_drawing.read_both("svg", path, arrays=["x"]).cell_data
    got = helpers_drawing.read_both("svg", path, arrays=["svg:fill"])
    assert "svg:fill" in got.cell_data


def test_the_dispatch_reads_it(tmp_path):
    path = _written(helpers.quad_mesh, tmp_path)
    assert meshioplusplus.read(path).cells[0].type == "quad"
    assert meshioplusplus.svg.read(path).cells[0].type == "quad"
    with open(path, "rb") as f:
        assert meshioplusplus.svg.read(f).cells[0].type == "quad"
    with open(path) as f:
        assert meshioplusplus.svg.read(f).cells[0].type == "quad"


def test_a_file_without_the_marker_needs_strict_false(tmp_path):
    path = tmp_path / "plain.svg"
    path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10">'
        '<path d="M 0 10L 10 10L 10 0Z" /></svg>'
    )
    core, py = helpers_drawing.readers("svg")
    for read in (core, py, meshioplusplus.svg.read, meshioplusplus.read):
        with pytest.raises(meshioplusplus.ReadError, match="no meshio\\+\\+ drawing"):
            read(str(path))
    # strict=False reads it as a plain drawing, the y-flip undone from the viewBox.
    for read in (core, py, meshioplusplus.svg.read):
        back = read(str(path), strict=False)
        assert back.points.tolist() == [[0, 0, 0], [10, 0, 0], [10, 10, 0]]
        assert back.cells[0].data.tolist() == [[0, 1, 2]]


def test_a_file_from_a_release_before_the_marker_is_unmarked(tmp_path):
    path = tmp_path / "old.svg"
    path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1">'
        '<!--Written by meshio++ v10.0.0--><path d="M 0 1L 1 1L 1 0Z" /></svg>'
    )
    core, py = helpers_drawing.readers("svg")
    for read in (core, py):
        with pytest.raises(meshioplusplus.ReadError, match="no meshio\\+\\+ drawing"):
            read(str(path))


_REFUSED = [
    ('<polygon points="0,0 1,0 1,1" />', "", "<polygon>"),
    ('<path d="M 0 0L 1 0L 1 1Z" transform="translate(1)" />', "", "transform"),
    ('<path d="M 0 0C 1 0 1 1 0 1Z" />', "", "curve command 'C'"),
    ('<path d="M 0 0Q 1 0 1 1Z" />', "", "curve command 'Q'"),
    ('<path d="M 0 0A 1 1 0 0 1 1 1Z" />', "", "curve command 'A'"),
    ('<path d="M 0 0l 1 0l 0 1z" />', "", "relative path command 'l'"),
    ('<path d="M 0 0H 1V 1Z" />', "", "path command 'H'"),
    ('<path d="M 0 0L 1 0Z M 2 2L 3 3" />', "", "subpaths"),
    ('<path d="M 0 0 1 0 1 1Z" />', "", "repeated coordinates"),
    ('<path d="M 0 0L 1" />', "", "missing its coordinates"),
    ('<path d="M 0 0L 1 0L 1 1" />', "", "open path of 3"),
    ('<path d="M 0 0L 1 0Z" />', "", "closed path of 2"),
    ('<path d="L 0 0" />', "", "must start with M"),
    ('<path d="" />', "", "empty path"),
    ('<path d="M 1e999 0L 1 0L 1 1Z" />', "", "bad number"),
    ("<path />", "", "no d attribute"),
    ('<g><path d="M 0 0L 1 0L 1 1Z" /></g>', "", "<g>"),
    ('<circle cx="1" cy="1" r="1" />', "", "<circle>"),
    ('<path d="M 0 0L 1 0L 1 1Z" />', ' transform="scale(2)"', "transform"),
    (
        '<path d="M 0 0L 1 0L 1 1Z" />',
        ' xmlns:xlink="http://www.w3.org/1999/xlink"',
        "namespace prefixes",
    ),
]


@pytest.mark.parametrize("body, root_extra, match", _REFUSED)
def test_unsupported_constructs_are_refused_by_name(body, root_extra, match, tmp_path):
    path = tmp_path / "bad.svg"
    path.write_text(helpers_drawing.marked_svg(body, root_extra))
    core, py = helpers_drawing.readers("svg")
    for read in (core, py):
        with pytest.raises(meshioplusplus.ReadError, match=match):
            read(str(path))


@pytest.mark.parametrize(
    "text",
    [
        "",
        "not xml",
        "<svg",
        "<html></html>",
        '<svg xmlns="http://example.org/" viewBox="0 0 1 1"></svg>',
        '<svg xmlns="http://www.w3.org/2000/svg"><!--\nWritten by meshio++ v1\n'
        "Drawing: planar scale=1 offset=0,0\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: planar scale=1 offset=0,0\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: planar scale=0 offset=0,0\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: planar\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: projected azimuth=1\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: sideways\n--></svg>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1 1"><!--\nWritten by '
        "meshio++ v1\nDrawing: planar scale=nan offset=0,0\n--></svg>",
    ],
)
def test_malformed_input_is_always_a_read_error(text, tmp_path):
    path = tmp_path / "bad.svg"
    path.write_text(text)
    core, py = helpers_drawing.readers("svg")
    for read in (core, py):
        with pytest.raises(meshioplusplus.ReadError):
            read(str(path), strict=False)


def test_an_empty_drawing_reads_as_an_empty_mesh(tmp_path):
    back = helpers_drawing.read_both("svg", _written(helpers.empty_mesh, tmp_path))
    assert back.points.shape == (0, 3) and not back.cells


def test_the_marker_lines_follow_the_tag(tmp_path):
    text = _written(helpers.tri_mesh, tmp_path).read_text()
    head = text.split("-->")[0].splitlines()
    assert head[0].endswith("><!--")
    assert head[1].startswith("Written by meshio++ v")
    assert head[2] == "Drawing: planar scale=100 offset=0,1"
    text = _written(helpers.hex_mesh, tmp_path, name="proj.svg").read_text()
    assert (
        "\nDrawing: projected azimuth=45 elevation=35.264389682754654 roll=0 scale="
        in text
    )


def test_the_marker_survives_a_provenance_block(tmp_path):
    # Notes, the source and a timestamp join the comment after the marker line;
    # a "--" in one is broken apart so the document stays well formed.
    from meshioplusplus import _provenance

    with _provenance.scope(_provenance.Mode.BEST_EFFORT):
        _provenance.set_source("in.vtu", "vtu")
        _provenance.note("regions-dropped", "side regions -- no equivalent")
        path = _written(helpers.tri_mesh, tmp_path)
    text = path.read_text()
    assert "Converted from in.vtu (vtu)" in text
    assert "side regions - - no equivalent" in text
    ET.parse(path)  # well formed
    back = helpers_drawing.read_both("svg", path)
    assert len(back.points) == 4


def test_the_marker_is_there_with_provenance_off(tmp_path):
    from meshioplusplus import _provenance

    with _provenance.scope(_provenance.Mode.OFF):
        path = _written(helpers.tri_mesh, tmp_path)
    assert "Drawing: planar" in path.read_text()
    assert len(helpers_drawing.read_both("svg", path).points) == 4


def test_cpp_and_python_writers_agree_with_the_marker(tmp_path):
    # The byte-identity guarantee covers the new comment, planar and projected.
    for mesh in (helpers.tri_mesh_2d, helpers.hex_mesh, helpers.empty_mesh):
        cpp = tmp_path / "cpp.svg"
        py = tmp_path / "py.svg"
        meshioplusplus._core.svg_write(str(cpp), copy.deepcopy(mesh), image_width=None)
        meshioplusplus.svg._svg.write(str(py), copy.deepcopy(mesh), image_width=None)
        assert cpp.read_bytes() == py.read_bytes()
        assert b"Drawing: " in cpp.read_bytes()
