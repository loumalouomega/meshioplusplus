import pytest

import meshioplusplus
from meshioplusplus import _colormap

from . import helpers, helpers_coloring

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
    filepath = tmp_path / "out.tikz"
    meshioplusplus.write_points_cells(filepath, mesh.points, mesh.cells)

    content = filepath.read_text()
    assert "\\documentclass{standalone}" in content
    assert "\\begin{tikzpicture}" in content
    assert "\\end{tikzpicture}" in content
    assert content.count("\\draw") == _drawable_cell_count(mesh)


@pytest.mark.parametrize("mesh", test_set)
def test_cpp_matches_python(mesh, tmp_path):
    # TikZ is plain text: the C++ core writer must be byte-identical to the
    # pure-Python reference.
    cpp = tmp_path / "cpp.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus._core.tikz_write(str(cpp), mesh)
    meshioplusplus.tikz._tikz.write(str(py), mesh)
    assert cpp.read_text() == py.read_text()


def test_standalone_false(tmp_path):
    mesh = helpers.tri_mesh_2d
    filepath = tmp_path / "snippet.tikz"
    meshioplusplus.tikz.write(filepath, mesh, standalone=False)

    content = filepath.read_text()
    assert "\\documentclass" not in content
    assert "\\begin{tikzpicture}" in content
    assert "\\draw" in content


def test_non_flat_shell_is_projected(tmp_path):
    # A non-flat 3D surface mesh is projected (isometric camera) and drawn.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2]])],
    )
    filepath = tmp_path / "out.tikz"
    meshioplusplus.tikz.write(filepath, mesh)
    content = filepath.read_text()
    assert "\\documentclass{standalone}" in content
    assert content.count("\\draw") == 1


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

    filepath = tmp_path / "skin.tikz"
    meshioplusplus.tikz.write(filepath, copy.deepcopy(mesh))
    content = filepath.read_text()
    assert content.count("\\draw") == num_faces
    assert content.count("-- cycle;") == num_faces


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
    ],
)
def test_cpp_matches_python_3d(mesh, tmp_path):
    import copy

    # The byte-identity guarantee extends to the 3D projected path: skin
    # extraction, camera arithmetic, depth sort, and formatting must all
    # agree between the C++ core and the Python reference.
    cpp = tmp_path / "cpp.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus._core.tikz_write(str(cpp), copy.deepcopy(mesh))
    meshioplusplus.tikz._tikz.write(str(py), copy.deepcopy(mesh))
    assert cpp.read_text() == py.read_text()


def test_camera_angles_change_projection(tmp_path):
    import copy

    p1 = tmp_path / "iso.tikz"
    p2 = tmp_path / "other.tikz"
    meshioplusplus.tikz.write(p1, copy.deepcopy(helpers.hex_mesh))
    meshioplusplus.tikz.write(
        p2, copy.deepcopy(helpers.hex_mesh), azimuth=10.0, elevation=60.0
    )
    assert p1.read_text() != p2.read_text()
    assert p1.read_text().count("\\draw") == p2.read_text().count("\\draw") == 6


# --- data-driven colouring ---------------------------------------------------


@pytest.mark.parametrize(
    "make_mesh, kwargs", helpers_coloring.CASE_ARGS, ids=helpers_coloring.CASE_IDS
)
def test_colored_cpp_matches_python(make_mesh, kwargs, tmp_path):
    # The byte-identity guarantee extends to the coloured path: array lookup,
    # scalarisation, the per-face mean, the auto range, the colormap index and
    # the colorbar geometry must all agree with the pure-Python reference.
    cpp = tmp_path / "cpp.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus._core.tikz_write(str(cpp), make_mesh(), **kwargs)
    meshioplusplus.tikz._tikz.write(str(py), make_mesh(), **kwargs)
    assert cpp.read_bytes() == py.read_bytes()


@pytest.mark.parametrize(
    "make_mesh, kwargs", helpers_coloring.CASE_ARGS, ids=helpers_coloring.CASE_IDS
)
def test_colored_shim_matches_python(make_mesh, kwargs, tmp_path):
    # The shim passes everything POSITIONALLY into _core, so a mis-ordered
    # argument would silently colour by the wrong thing rather than fail.
    shim = tmp_path / "shim.tikz"
    py = tmp_path / "py.tikz"
    meshioplusplus.tikz.write(shim, make_mesh(), **kwargs)
    meshioplusplus.tikz._tikz.write(str(py), make_mesh(), **kwargs)
    assert shim.read_bytes() == py.read_bytes()


def test_color_by_unset_is_unchanged(tmp_path):
    # The regression guard: passing the colouring parameters explicitly at
    # their defaults must reproduce the legacy output exactly.
    mesh = helpers_coloring.volume_mesh()
    plain = tmp_path / "plain.tikz"
    explicit = tmp_path / "explicit.tikz"
    meshioplusplus._core.tikz_write(str(plain), mesh)
    meshioplusplus._core.tikz_write(
        str(explicit),
        mesh,
        color_by="",
        component=None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        nan_color="gray",
        colorbar=True,  # ignored: no colouring means no bar
    )
    assert plain.read_bytes() == explicit.read_bytes()


def _rgb(t):
    r, g, b = _colormap.colormap_lookup(_colormap.colormap_table("viridis"), t)
    return f"{{rgb,255:red,{r};green,{g};blue,{b}}}"


def test_golden_fill_vocabulary(tmp_path):
    # Pins the exact TikZ colour spelling: the braces matter, since without
    # them the inner commas would split the surrounding option list.
    import numpy as np

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2]])],
        cell_data={"tag": [np.array([0.0, 1.0])]},
    )
    out = tmp_path / "golden.tikz"
    meshioplusplus.tikz.write(out, mesh, color_by="tag")
    text = out.read_text()
    assert f"\\draw[fill={_rgb(0.0)}, draw=black]" in text
    assert f"\\draw[fill={_rgb(1.0)}, draw=black]" in text
    assert "fill={rgb,255:red,68;green,1;blue,84}" in text  # viridis low = #440154


def test_nan_color_and_clamping(tmp_path):
    import numpy as np

    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2], [1, 4, 3]])],
        cell_data={"tag": [np.array([-100.0, np.nan, 100.0])]},
    )
    out = tmp_path / "nan.tikz"
    meshioplusplus.tikz.write(
        out, mesh, color_by="tag", vmin=0.0, vmax=1.0, nan_color="red!50"
    )
    text = out.read_text()
    assert "fill=red!50, draw=black" in text
    assert f"fill={_rgb(0.0)}, draw=black" in text  # clamped from -100
    assert f"fill={_rgb(1.0)}, draw=black" in text  # clamped from +100


def test_colorbar_appends_only(tmp_path):
    mesh = helpers_coloring.flat_mesh()
    without = tmp_path / "without.tikz"
    with_bar = tmp_path / "with.tikz"
    meshioplusplus.tikz.write(without, mesh, color_by="tag")
    meshioplusplus.tikz.write(with_bar, mesh, color_by="tag", colorbar=True)

    a = without.read_text().splitlines()
    b = with_bar.read_text().splitlines()
    # TikZ has no viewBox, so the mesh lines are untouched and the bar is a
    # pure insertion before \end{tikzpicture}.
    draws_a = [line for line in a if "\\draw[" in line]
    draws_b = [line for line in b if "\\draw[" in line]
    assert draws_a == draws_b
    assert len([line for line in b if "\\fill[" in line]) == 32
    assert len([line for line in b if "\\node[" in line]) == 2
    assert [line for line in a if "\\fill[" in line] == []


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
        meshioplusplus.tikz.write(tmp_path / "bad.tikz", mesh, **kwargs)
    with pytest.raises(Exception, match=match):
        meshioplusplus._core.tikz_write(str(tmp_path / "bad.tikz"), mesh, **kwargs)
