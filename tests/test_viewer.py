"""Tests for the mesh -> renderer mapping.

Everything here except the ``viewer``-marked tests is pure: it exercises
``_to_polyscope_payload`` with no polyscope, no display and no browser, which
is where essentially all of the logic lives.
"""

import base64
import copy
import io
import json
import re
import struct
from pathlib import Path

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._viewer import (
    _flatten_cell_data,
    _quantities_from_array,
    _to_polyscope_payload,
    has_viewer,
)
from meshioplusplus._viewer_browser import _renderable_surface

from . import helpers


def _q(payload, name):
    """The one quantity named ``name``, or None."""
    for q in payload.quantities:
        if q.name == name:
            return q
    return None


def _names(payload):
    return sorted(q.name for q in payload.quantities)


# --- kind routing --------------------------------------------------------- #


def test_surface_mesh_routes_to_surface():
    p = _to_polyscope_payload(helpers.tri_mesh)
    assert p.kind == "surface"
    assert p.faces == [[0, 1, 2], [0, 2, 3]]
    assert p.tets is None and p.hexes is None and p.mixed_cells is None


def test_volume_mesh_routes_to_volume():
    p = _to_polyscope_payload(helpers.tet_mesh)
    assert p.kind == "volume"
    assert p.tets.shape == (2, 4)
    assert p.mixed_cells is None
    assert p.faces is None


def test_hex_mesh_uses_the_hex_fast_path_and_is_not_simplexified():
    p = _to_polyscope_payload(helpers.hex_mesh)
    assert p.kind == "volume"
    assert p.hexes.shape == (1, 8)
    assert p.tets is None
    # Decomposing a hex mesh would destroy exactly the structure a user opens a
    # viewer to look at.
    assert p.notes == []


def test_line_mesh_routes_to_curve():
    p = _to_polyscope_payload(helpers.line_mesh)
    assert p.kind == "curve"
    assert p.edges.shape[1] == 2


def test_points_kind_ignores_cells():
    p = _to_polyscope_payload(helpers.tet_mesh, kind="points")
    assert p.kind == "points"
    assert p.num_primitives == 0
    assert len(p.vertices) == len(helpers.tet_mesh.points)


def test_unknown_kind_raises():
    with pytest.raises(ValueError, match="unknown kind"):
        _to_polyscope_payload(helpers.tri_mesh, kind="wireframe")


def test_volume_kind_on_a_surface_mesh_raises():
    with pytest.raises(ValueError, match="needs 3D cells"):
        _to_polyscope_payload(helpers.tri_mesh, kind="volume")


def test_polyhedron_volume_raises_by_name():
    with pytest.raises(ValueError, match="polyhedron"):
        _to_polyscope_payload(helpers.polyhedron_mesh, kind="volume")


# --- geometry ------------------------------------------------------------- #


def test_2d_points_are_padded_to_three_columns():
    p = _to_polyscope_payload(helpers.tri_mesh_2d)
    assert p.vertices.shape[1] == 3
    assert np.all(p.vertices[:, 2] == 0.0)
    assert any("padded to 3D" in n for n in p.notes)


def test_mixed_triangle_quad_faces_stay_ragged():
    # Three blocks: triangle, quad, triangle. Faces follow block order, so the
    # two triangle blocks are NOT grouped together.
    p = _to_polyscope_payload(helpers.tri_quad_mesh)
    assert [len(f) for f in p.faces] == [3, 3, 4, 3]


def test_ragged_polygon_block_becomes_a_ragged_face_list():
    p = _to_polyscope_payload(helpers.polygon_mesh)
    assert p.kind == "surface"
    assert len({len(f) for f in p.faces}) > 1
    assert all(isinstance(i, int) for f in p.faces for i in f)


def test_tet_and_hex_together_become_padded_mixed_cells():
    mesh = meshioplusplus.Mesh(
        np.vstack([helpers.hex_mesh.points, [[0.5, 0.5, 2.0]]]),
        [
            ("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7]]),
            ("tetra", [[4, 5, 6, 8]]),
        ],
    )
    p = _to_polyscope_payload(mesh)
    assert p.mixed_cells.shape == (2, 8)
    assert p.tets is None and p.hexes is None
    # Row order is ours (block order), not polyscope's merge order.
    assert np.all(p.mixed_cells[0] >= 0)  # the hexahedron
    assert np.all(p.mixed_cells[1, 4:] < 0)  # the tetrahedron, padded
    assert p.mixed_cells[1, :4].tolist() == [4, 5, 6, 8]
    # No conversion happened: both types are natively renderable.
    assert p.notes == []


def test_wedge_is_simplexified_with_a_note():
    p = _to_polyscope_payload(helpers.wedge_mesh)
    assert p.kind == "volume"
    assert p.tets.shape == (3, 4)
    assert any("split into tetrahedra" in n for n in p.notes)


def test_pyramid_is_simplexified():
    p = _to_polyscope_payload(helpers.pyramid_mesh)
    assert p.tets.shape == (2, 4)


def test_simplexifying_a_mixed_mesh_says_hexes_were_collateral():
    mesh = meshioplusplus.Mesh(
        np.vstack([helpers.hex_mesh.points, [[2.0, 0.0, 0.0], [2.0, 1.0, 0.0]]]),
        [
            ("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7]]),
            ("wedge", [[1, 8, 2, 5, 9, 6]]),
        ],
    )
    p = _to_polyscope_payload(mesh)
    assert p.tets is not None
    assert any("hexahedra in the same mesh were split too" in n for n in p.notes)


def test_higher_order_volume_is_linearized_then_simplexified():
    p = _to_polyscope_payload(helpers.tet10_mesh)
    assert p.kind == "volume"
    assert p.tets.shape[1] == 4
    # The mid-side nodes are orphaned by linearization and pruned.
    assert len(p.vertices) < len(helpers.tet10_mesh.points)


def test_higher_order_surface_is_linearized():
    p = _to_polyscope_payload(helpers.triangle6_mesh)
    assert p.kind == "surface"
    assert all(len(f) == 3 for f in p.faces)
    assert any("linear bases" in n for n in p.notes)


def test_surface_kind_extracts_the_boundary_of_a_volume_mesh():
    p = _to_polyscope_payload(helpers.hex_mesh, kind="surface")
    assert p.kind == "surface"
    assert len(p.faces) == 6
    assert all(len(f) == 4 for f in p.faces)


def test_input_mesh_is_never_modified():
    before = copy.deepcopy(helpers.tet10_mesh)
    _to_polyscope_payload(helpers.tet10_mesh)
    assert np.allclose(helpers.tet10_mesh.points, before.points)
    assert len(helpers.tet10_mesh.cells) == len(before.cells)
    assert helpers.tet10_mesh.cells[0].type == before.cells[0].type


# --- cell_data alignment (the reason cell_source exists) ------------------- #


def _with_cell_values(mesh, values_per_block, name="tag"):
    m = copy.deepcopy(mesh)
    m.cell_data[name] = [np.asarray(v) for v in values_per_block]
    return m


def test_cell_data_follows_block_order_both_ways_round():
    """The regression test for the whole cell_source design.

    ``tri_quad_mesh`` (triangle, quad, triangle) and ``quad_tri_mesh`` (quad,
    triangle) put the same cell types in different block orders. A per-face
    quantity must follow the faces, not be regrouped by type.
    """
    tq = _with_cell_values(helpers.tri_quad_mesh, [[10, 11], [20], [30]])
    qt = _with_cell_values(helpers.quad_tri_mesh, [[20], [10, 11]])

    p_tq = _to_polyscope_payload(tq)
    p_qt = _to_polyscope_payload(qt)

    assert [len(f) for f in p_tq.faces] == [3, 3, 4, 3]
    assert [len(f) for f in p_qt.faces] == [4, 3, 3]
    assert _q(p_tq, "tag").values.tolist() == [10.0, 11.0, 20.0, 30.0]
    assert _q(p_qt, "tag").values.tolist() == [20.0, 10.0, 11.0]


def test_a_skipped_lower_dimensional_block_does_not_offset_cell_data():
    """A gmsh-style boundary-marker `line` block sits before the triangles.

    A naive per-dimension concatenation would read the line block's values as
    the triangles' and colour the whole mesh wrong.
    """
    mesh = meshioplusplus.Mesh(
        helpers.tri_mesh.points,
        [("line", [[0, 1], [1, 2]]), ("triangle", [[0, 1, 2], [0, 2, 3]])],
    )
    mesh.cell_data["tag"] = [np.array([-1, -2]), np.array([7, 8])]
    p = _to_polyscope_payload(mesh)
    assert p.kind == "surface"
    assert len(p.faces) == 2
    assert _q(p, "tag").values.tolist() == [7.0, 8.0]


def test_cell_data_is_replicated_to_simplexified_children():
    mesh = _with_cell_values(helpers.wedge_mesh, [[5]])
    p = _to_polyscope_payload(mesh)
    assert p.tets.shape[0] == 3
    assert _q(p, "tag").values.tolist() == [5.0, 5.0, 5.0]


def test_surface_cell_data_is_sampled_from_the_owning_volume_cell():
    mesh = _with_cell_values(helpers.hex_mesh, [[42]])
    p = _to_polyscope_payload(mesh, kind="surface")
    assert len(p.faces) == 6
    assert _q(p, "tag").values.tolist() == [42.0] * 6
    assert any("owning volume cell" in n for n in p.notes)
    # The provenance array is plumbing, not something to colour by.
    assert "surface:parent_cell" not in _names(p)


def test_point_data_on_an_extracted_surface_is_the_subset_not_the_original():
    mesh = copy.deepcopy(helpers.tet_mesh)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    p = _to_polyscope_payload(mesh, kind="surface")
    assert len(_q(p, "t").values) == len(p.vertices)


def test_malformed_cell_data_is_dropped_with_a_note_not_guessed():
    mesh = copy.deepcopy(helpers.tri_quad_mesh)
    mesh.cell_data["bad"] = [np.array([1, 2])]  # one block for two
    p = _to_polyscope_payload(mesh)
    assert _q(p, "bad") is None
    assert any("bad" in n and "dropped" in n for n in p.notes)


def test_flatten_cell_data_rejects_inconsistent_widths():
    mesh = copy.deepcopy(helpers.tri_quad_mesh)
    mesh.cell_data["w"] = [np.zeros((2, 3)), np.zeros((1, 2)), np.zeros((1, 3))]
    notes = []
    assert _flatten_cell_data(mesh, "w", notes) is None
    assert any("disagree on component count" in n for n in notes)


# --- data array -> quantity rules ----------------------------------------- #


def _map(name, arr):
    notes = []
    return (
        _quantities_from_array(name, np.asarray(arr), "vertices", "point_data", notes),
        notes,
    )


def test_scalar_array_becomes_one_scalar_quantity():
    qs, _ = _map("t", [1.0, 2.0, 3.0])
    assert len(qs) == 1
    assert qs[0].kind == "scalar" and qs[0].values.shape == (3,)
    assert qs[0].values.dtype == np.float64


def test_two_wide_array_is_padded_to_a_vector():
    qs, notes = _map("v", [[1.0, 2.0]])
    assert len(qs) == 1 and qs[0].kind == "vector"
    assert qs[0].values.tolist() == [[1.0, 2.0, 0.0]]
    assert any("zero-padded" in n for n in notes)


def test_three_wide_array_becomes_a_vector_plus_a_magnitude():
    qs, _ = _map("disp", [[3.0, 4.0, 0.0]])
    assert [q.kind for q in qs] == ["vector", "scalar"]
    assert qs[1].name == "disp:magnitude"
    assert qs[1].values.tolist() == [5.0]


def test_color_is_detected_by_name_not_by_range():
    """A normalized field in [0,1] must not silently render as RGB."""
    qs, _ = _map("normalized_displacement", [[0.1, 0.2, 0.3]])
    assert qs[0].kind == "vector"

    qs, _ = _map("color", [[0.1, 0.2, 0.3]])
    assert len(qs) == 1 and qs[0].kind == "color"

    qs, _ = _map("cell_rgb", [[0.1, 0.2, 0.3]])
    assert qs[0].kind == "color"


def test_uint8_colors_are_scaled_to_unit_range():
    qs, _ = _map("colors", np.array([[255, 128, 0]], dtype=np.uint8))
    assert qs[0].kind == "color"
    assert qs[0].values[0, 0] == pytest.approx(1.0)
    assert qs[0].values[0, 2] == pytest.approx(0.0)


def test_a_color_named_array_out_of_range_falls_back_to_a_vector():
    qs, _ = _map("color", [[5.0, -3.0, 12.0]])
    assert qs[0].kind == "vector"


def test_tensor_becomes_components_plus_a_norm():
    qs, notes = _map("stress", np.arange(9, dtype=float).reshape(1, 9))
    names = [q.name for q in qs]
    assert names == [f"stress:{i}" for i in range(9)] + ["stress:norm"]
    assert any("no tensor quantity" in n for n in notes)


def test_very_wide_array_keeps_only_a_magnitude():
    qs, notes = _map("modes", np.ones((2, 40)))
    assert [q.name for q in qs] == ["modes:magnitude"]
    assert any("only its magnitude" in n for n in notes)


def test_bool_array_is_cast_not_dropped():
    qs, _ = _map("flag", np.array([True, False]))
    assert qs[0].kind == "scalar"
    assert qs[0].values.tolist() == [1.0, 0.0]


def test_non_numeric_array_is_dropped_with_a_note():
    qs, notes = _map("label", np.array(["a", "b"]))
    assert qs == []
    assert any("non-numeric" in n for n in notes)


def test_nan_forces_an_explicit_colormap_range():
    qs, _ = _map("t", [1.0, np.nan, 3.0])
    assert qs[0].vminmax == (1.0, 3.0)


def test_all_finite_array_leaves_the_range_to_polyscope():
    qs, _ = _map("t", [1.0, 2.0])
    assert qs[0].vminmax is None


def test_all_nan_array_is_dropped():
    qs, notes = _map("t", [np.nan, np.nan])
    assert qs == []
    assert any("no finite values" in n for n in notes)


def test_point_and_cell_data_name_collision_is_disambiguated():
    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["t"] = np.zeros(len(mesh.points))
    mesh.cell_data["t"] = [np.array([1, 2])]
    p = _to_polyscope_payload(mesh)
    assert _q(p, "t").defined_on == "vertices"
    assert _q(p, "t (cells)").defined_on == "faces"
    assert any("both point and cell data" in n for n in p.notes)


def test_notes_are_empty_when_nothing_lossy_happened():
    p = _to_polyscope_payload(helpers.tri_mesh)
    assert p.notes == []


def test_quantity_values_are_contiguous_float64():
    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=np.float32)
    p = _to_polyscope_payload(mesh)
    q = _q(p, "t")
    assert q.values.dtype == np.float64
    assert q.values.flags["C_CONTIGUOUS"]


# --- the browser backend (no browser needed) ------------------------------ #


def _payload_of(page):
    """The base64 VTP embedded in a generated page, decoded."""
    match = re.search(r'id="vtp-payload">([^<]*)</script>', page)
    assert match, "the generated page has no payload"
    return base64.b64decode(match.group(1))


def test_browser_page_embeds_the_mesh_and_the_bundle():
    from meshioplusplus._viewer_browser import build_page

    page = build_page(helpers.tri_mesh, title="a mesh")
    assert "<title>a mesh</title>" in page
    # Self-contained: no external script, style or wasm to fetch. That is what
    # makes it work over file:// with no local HTTP server.
    assert "<script" in page and "src=http" not in page.replace('"', "")
    assert ".wasm" not in page
    assert len(page) > 100_000, "the vtk.js bundle does not look inlined"

    vtp = _payload_of(page)
    assert vtp.startswith(b"<?xml") or vtp.startswith(b"<VTKFile")


def test_browser_payload_round_trips_through_meshioplusplus():
    """The embedded bytes must be a mesh meshio++ itself can read back."""
    from meshioplusplus._viewer_browser import build_page

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    vtp = _payload_of(build_page(mesh))

    back = meshioplusplus.read(io.BytesIO(vtp), file_format="vtp")
    assert len(back.points) == len(mesh.points)
    assert sum(len(b) for b in back.cells) == 2
    assert np.allclose(back.point_data["t"], mesh.point_data["t"])


def test_browser_payload_keeps_multi_component_data():
    """A vector array must survive to the renderer.

    This is the property the whole VTP-in-C++ path exists to preserve; the flat
    WASM mesh representation cannot carry it.
    """
    from meshioplusplus._viewer_browser import build_page

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["disp"] = np.arange(len(mesh.points) * 3, dtype=float).reshape(
        -1, 3
    )
    back = meshioplusplus.read(
        io.BytesIO(_payload_of(build_page(mesh))), file_format="vtp"
    )
    assert back.point_data["disp"].shape == (len(mesh.points), 3)


def test_browser_renders_a_volume_mesh_as_its_boundary():
    """vtk.js has no volume mesh type, so a solid is shown by its skin."""
    from meshioplusplus._viewer_browser import build_page

    vtp = _payload_of(build_page(helpers.hex_mesh))
    back = meshioplusplus.read(io.BytesIO(vtp), file_format="vtp")
    assert sum(len(b) for b in back.cells) == 6
    assert all(b.type == "quad" for b in back.cells)


def test_browser_linearizes_higher_order_surface_cells():
    from meshioplusplus._viewer_browser import build_page

    vtp = _payload_of(build_page(helpers.triangle6_mesh))
    back = meshioplusplus.read(io.BytesIO(vtp), file_format="vtp")
    assert all(b.type == "triangle" for b in back.cells)


def test_browser_rejects_kind_volume_rather_than_ignoring_it():
    from meshioplusplus._viewer_browser import build_page

    with pytest.raises(ValueError, match="no volume mesh type"):
        build_page(helpers.hex_mesh, kind="volume")


def test_view_browser_writes_a_page_without_opening_one(tmp_path, monkeypatch):
    import webbrowser

    opened = []
    monkeypatch.setattr(webbrowser, "open", lambda url: opened.append(url) or True)

    path = meshioplusplus.view(helpers.tri_mesh, backend="browser")
    assert len(opened) == 1
    assert opened[0].startswith("file://")
    assert _payload_of(Path(path).read_text(encoding="utf-8"))


def test_bundled_viewer_asset_is_present_and_current():
    """The committed asset must exist and carry the slot Python fills.

    Deliberately structural rather than a byte comparison against a fresh
    build: Rollup output is not reproducible across npm patch releases, so a
    `git diff --exit-code` here would be a recurring red build unrelated to any
    change of ours.
    """
    from meshioplusplus._viewer_browser import _ASSET, _PAYLOAD_SLOT

    assert _ASSET.exists(), f"{_ASSET} is missing; see src/viewer/README.md"
    page = _ASSET.read_text(encoding="utf-8")
    assert _PAYLOAD_SLOT in page
    assert ".wasm" not in page, "the embedded viewer must not carry WebAssembly"


# --- the polyscope backend, when installed -------------------------------- #


def test_has_viewer_is_a_bool():
    assert isinstance(has_viewer(), bool)


def test_color_by_an_unknown_array_lists_what_is_available():
    from meshioplusplus._viewer import _resolve_color_by

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["t"] = np.zeros(len(mesh.points))
    p = _to_polyscope_payload(mesh)
    with pytest.raises(ValueError, match=r"no data array named 'nope'.*available: t"):
        _resolve_color_by(p, "nope")


def test_view_rejects_an_unknown_backend():
    from meshioplusplus import view

    with pytest.raises(ValueError, match="unknown backend"):
        view(helpers.tri_mesh, backend="opengl")


def test_polyscope_backend_error_names_the_install_command(monkeypatch):
    """The missing-dependency path must say how to fix itself."""
    import builtins

    from meshioplusplus import _viewer

    real = builtins.__import__

    def blocked(name, *a, **k):
        if name == "polyscope":
            raise ImportError("no polyscope here")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", blocked)
    assert _viewer.has_viewer() is False
    with pytest.raises(ImportError, match=r"pip install meshioplusplus\[viewer\]"):
        _viewer._import_polyscope()


@pytest.fixture
def ps_headless():
    """An initialized headless polyscope, or a skip."""
    ps = pytest.importorskip("polyscope")
    from meshioplusplus._viewer import _init_polyscope

    try:
        _init_polyscope(ps, headless=True)
    except Exception as e:  # pragma: no cover - depends on the GL runtime
        pytest.skip(f"no headless polyscope backend: {e}")
    ps.remove_all_structures()
    yield ps
    ps.remove_all_structures()


#: Each polyscope structure names its counts differently.
_POINT_COUNT = {
    "surface": "n_vertices",
    "volume": "n_vertices",
    "curve": "n_nodes",
    "points": "n_points",
}
_PRIMITIVE_COUNT = {"surface": "n_faces", "volume": "n_cells", "curve": "n_edges"}


@pytest.mark.viewer
@pytest.mark.parametrize(
    "fixture", ["tet_mesh", "hex_mesh", "tri_quad_mesh", "wedge_mesh", "line_mesh"]
)
def test_register_matches_the_payload_counts(ps_headless, fixture):
    from meshioplusplus._viewer import _register

    payload = _to_polyscope_payload(getattr(helpers, fixture))
    struct = _register(ps_headless, payload, "test")
    assert getattr(struct, _POINT_COUNT[payload.kind])() == len(payload.vertices)
    if payload.kind in _PRIMITIVE_COUNT:
        got = getattr(struct, _PRIMITIVE_COUNT[payload.kind])()
        assert got == payload.num_primitives


@pytest.mark.viewer
def test_every_quantity_kind_registers(ps_headless):
    """Scalars, vectors and colours must all survive the real polyscope call.

    Each structure type has its own `defined_on` vocabulary -- point clouds
    take none at all -- so this is the only place a wrong one shows up.
    """
    from meshioplusplus._viewer import _register

    mesh = copy.deepcopy(helpers.tet_mesh)
    n = len(mesh.points)
    mesh.point_data["t"] = np.arange(n, dtype=float)
    mesh.point_data["disp"] = np.zeros((n, 3))
    mesh.point_data["color"] = np.zeros((n, 3))
    mesh.cell_data["tag"] = [np.array([1, 2])]
    payload = _to_polyscope_payload(mesh)
    _register(ps_headless, payload, "test", color_by="t")

    for kind in ("points", "curve"):
        p = _to_polyscope_payload(helpers.line_mesh, kind=kind)
        _register(ps_headless, p, f"as-{kind}")


@pytest.mark.viewer
def test_screenshot_writes_a_png_of_the_requested_size(ps_headless, tmp_path):
    from meshioplusplus import screenshot

    mesh = copy.deepcopy(helpers.hex_mesh)
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    out = tmp_path / "shot.png"
    screenshot(mesh, out, color_by="t", size=(320, 240))

    data = out.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert struct.unpack(">II", data[16:24]) == (320, 240)


@pytest.mark.viewer
def test_screenshot_actually_renders_the_mesh_it_was_given(ps_headless, tmp_path):
    """Two different meshes must not produce the same image.

    Catches both failure modes at once: a blank canvas, and a pipeline that
    runs end to end while rendering something unrelated to its input.
    """
    from meshioplusplus import screenshot

    a, b = tmp_path / "a.png", tmp_path / "b.png"
    screenshot(helpers.tet_mesh, a, size=(160, 120))
    screenshot(helpers.tri_mesh, b, size=(160, 120))
    assert a.read_bytes() != b.read_bytes()

    # And each has more than a background's worth of colour in it.
    buf = np.asarray(ps_headless.screenshot_to_buffer())
    assert buf.shape[:2] == (120, 160)


@pytest.mark.viewer
def test_consecutive_screenshots_do_not_accumulate_structures(ps_headless, tmp_path):
    """The second image shows only the second mesh."""
    from meshioplusplus import screenshot

    screenshot(helpers.tet_mesh, tmp_path / "1.png", size=(160, 120))
    screenshot(helpers.tri_mesh, tmp_path / "2.png", size=(160, 120))
    solo = (tmp_path / "solo.png").as_posix()
    ps_headless.remove_all_structures()
    screenshot(helpers.tri_mesh, solo, size=(160, 120))
    assert (tmp_path / "2.png").read_bytes() == (tmp_path / "solo.png").read_bytes()


@pytest.mark.viewer
def test_view_with_show_false_registers_without_blocking(ps_headless):
    from meshioplusplus import view

    struct = view(helpers.tet_mesh, backend="polyscope", show=False, name="m")
    assert struct.n_vertices() == len(helpers.tet_mesh.points)


def test_browser_carries_cell_data_onto_the_boundary():
    """A solid's per-cell tag must reach the renderer.

    `extract_surface` drops cell data, so without the parent-id gather this
    array would simply vanish -- and colouring a solid by its material is the
    common case. The WASM `convertSurface` binding does the same thing, and the
    two must agree.
    """
    from meshioplusplus._viewer_browser import build_page

    mesh = copy.deepcopy(helpers.hex_mesh)
    mesh.cell_data["material"] = [np.array([7])]
    back = meshioplusplus.read(
        io.BytesIO(_payload_of(build_page(mesh))), file_format="vtp"
    )
    assert np.array_equal(np.concatenate(back.cell_data["material"]), [7] * 6)
    # The provenance array is plumbing; it must not clutter a colour-by menu.
    assert "surface:parent_cell" not in back.cell_data


# --- what Python bakes into the offline page ------------------------------ #


def _slot(page, tag_id):
    """The content of one payload slot."""
    match = re.search(rf'id="{tag_id}">([^<]*)</script>', page)
    assert match, f"the generated page has no {tag_id} slot"
    return match.group(1)


def test_offline_page_can_bake_in_quality_metrics():
    """The page has no WASM, so anything it should show is computed here.

    Decoded and read back rather than grepped: a substring check would pass on
    a corrupt payload.
    """
    from meshioplusplus._viewer_browser import build_page

    page = build_page(helpers.hex_mesh, quality=True)
    back = meshioplusplus.read(io.BytesIO(_payload_of(page)), file_format="vtp")
    assert any(name.startswith("quality:") for name in back.cell_data)
    # The metrics ride the existing parent-id gather onto the boundary, so
    # every boundary face must carry one.
    metric = next(k for k in back.cell_data if k.startswith("quality:"))
    assert sum(len(b) for b in back.cell_data[metric]) == 6


def test_offline_page_bakes_stats_the_browser_could_not_derive():
    """The page renders only the boundary; volume and area come from here."""
    from meshioplusplus._viewer_browser import build_page

    stats = json.loads(_slot(build_page(helpers.hex_mesh), "stats-payload"))
    assert stats["Volume"] == "1"
    assert stats["Surface area"] == "6"


def test_offline_page_can_preselect_a_colour_array():
    from meshioplusplus._viewer_browser import build_page

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["t"] = np.zeros(len(mesh.points))
    config = json.loads(_slot(build_page(mesh, color_by="t"), "view-config"))
    assert config["colorBy"] == "point:t:-1"


def test_offline_page_raises_on_an_unknown_colour_array():
    """Raise rather than warn: the mesh is right here to check against."""
    from meshioplusplus._viewer_browser import build_page

    with pytest.raises(ValueError, match=r"no data array named 'nope'"):
        build_page(helpers.tri_mesh, color_by="nope")


def test_python_and_wasm_surfaces_agree():
    """`_renderable_surface` is the twin of the WASM `convertSurfaceOps([])`.

    Only their agreement makes a mesh look the same whether it reached the
    viewer through Python or through the hosted demo. Until now that was
    implied by two Playwright suites happening to assert the same numbers.
    """
    surface = _renderable_surface(helpers.hex_mesh, "auto")
    assert sum(len(b) for b in surface.cells) == 6
    assert all(b.type == "quad" for b in surface.cells)
    assert len(surface.points) == 8
