"""The software rasterizer through Python (meshioplusplus._render) and both CLIs.

The raster's own probes (top-left rule, depth, encodings) are C++ tests in
tests/cpp/test_render.cpp; this file covers what the Python surface adds: the
option plumbing, the file forms, the screenshot fallback, the core-only
contract and byte-identity between the Python and native ``snapshot`` verbs.
"""

import glob
import json
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

import meshioplusplus as mio

from . import helpers

pytest.importorskip("meshioplusplus._core")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _native_cli():
    env = os.environ.get("MESHIOPLUSPLUS_NATIVE_CLI")
    if env and os.path.isfile(env):
        return env
    found = [
        p
        for p in glob.glob(os.path.join(ROOT, "build", "*", "meshioplusplus"))
        if os.path.isfile(p) and os.access(p, os.X_OK)
    ]
    return max(found, key=os.path.getmtime) if found else None


NATIVE = _native_cli()


def _cube():
    pts = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ]
    quads = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    mesh = mio.Mesh(np.array(pts, dtype=float), [("quad", np.array(quads))])
    mesh.point_data["x"] = np.array(pts, dtype=float)[:, 0]
    mesh.cell_data["id"] = [np.arange(6, dtype=float)]
    return mesh


def test_render_image_is_an_rgba_array_with_a_cell_id_buffer():
    image, info = mio.render_image(_cube(), 64, 48, return_info=True)
    assert image.shape == (48, 64, 4)
    assert image.dtype == np.uint8
    assert info["cell_ids"].shape == (48, 64)
    # The isometric camera sees the faces x = 1, y = 1 and z = 1 only.
    assert set(np.unique(info["cell_ids"])) == {-1, 1, 3, 5}
    assert info["colored"] is False
    assert info["vmin"] is None


def test_render_image_is_deterministic_and_options_change_it():
    a = mio.render_image(_cube(), 80, 60, supersample=2, edges="all")
    b = mio.render_image(_cube(), 80, 60, supersample=2, edges="all")
    assert np.array_equal(a, b)
    c = mio.render_image(_cube(), 80, 60, supersample=2, edges="all", view="+x")
    assert not np.array_equal(a, c)


def test_colouring_reports_the_range():
    image, info = mio.render_image(_cube(), 40, 40, color_by="id", return_info=True)
    assert info["colored"] is True
    assert (info["vmin"], info["vmax"]) == (0.0, 5.0)
    assert info["notes"] == ["id: 0 .. 5 (viridis)"]
    with pytest.raises(ValueError, match="nope"):
        mio.render_image(_cube(), 40, 40, color_by="nope")


def test_unknown_options_are_named():
    with pytest.raises(ValueError, match="unknown render option 'bogus'"):
        mio.render_image(_cube(), bogus=1)
    with pytest.raises(TypeError, match="encoding"):
        mio.render_image(_cube(), encoding="braille")
    with pytest.raises(ValueError, match="supersample"):
        mio.render_image(_cube(), supersample=3)


def test_render_text_fills_the_requested_cells():
    text = mio.render_text(_cube(), cols=30, rows=10, color_depth="mono")
    lines = text.splitlines()
    assert len(lines) == 10
    assert all(len(line) == 30 for line in lines)
    assert "\x1b" not in text
    coloured = mio.render_text(
        _cube(), cols=30, rows=10, color_depth="256", color_by="x"
    )
    assert "\x1b[38;5;" in coloured or "\x1b[48;5;" in coloured
    assert coloured.rstrip("\n").endswith("x: 0 .. 1 (viridis)")


@pytest.mark.parametrize(
    "encoding", ["halfblock", "quadrant", "sextant", "braille", "ascii"]
)
def test_every_cell_encoding_renders(encoding):
    text = mio.render_text(
        _cube(), cols=12, rows=6, encoding=encoding, color_depth="truecolor"
    )
    assert len(text.splitlines()) == 6


def test_snapshot_writes_each_form(tmp_path):
    mesh = _cube()
    png = mio.snapshot(mesh, tmp_path / "a.png", width=50, height=40)
    data = open(png, "rb").read()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert int.from_bytes(data[16:20], "big") == 50
    assert int.from_bytes(data[20:24], "big") == 40
    mio.snapshot(mesh, tmp_path / "a.txt", cols=20, rows=8)
    txt = (tmp_path / "a.txt").read_text(encoding="utf-8")
    assert len(txt.splitlines()) == 8 and "\x1b" not in txt
    mio.snapshot(mesh, tmp_path / "a.ansi", cols=20, rows=8)
    assert "\x1b[" in (tmp_path / "a.ansi").read_text(encoding="utf-8")
    mio.snapshot(mesh, tmp_path / "a.html", cols=20, rows=8)
    assert (
        (tmp_path / "a.html").read_text(encoding="utf-8").startswith("<!DOCTYPE html>")
    )
    mio.snapshot(
        mesh, tmp_path / "a.cast", cols=20, rows=8, cast_frames=3, cast_fps=3.0
    )
    lines = (tmp_path / "a.cast").read_text(encoding="utf-8").splitlines()
    assert json.loads(lines[0]) == {"height": 8, "version": 2, "width": 20}
    assert [json.loads(line)[0] for line in lines[1:]] == [
        0.0,
        pytest.approx(1 / 3),
        pytest.approx(2 / 3),
    ]
    with pytest.raises(ValueError, match="expected .png"):
        mio.snapshot(mesh, tmp_path / "a.vtu")


def test_png_matches_the_image(tmp_path):
    zlib = pytest.importorskip("zlib")
    mesh = _cube()
    image = mio.render_image(mesh, 30, 20)
    mio.snapshot(mesh, tmp_path / "a.png", width=30, height=20)
    data = open(tmp_path / "a.png", "rb").read()
    length = int.from_bytes(data[33:37], "big")
    assert data[37:41] == b"IDAT"
    raw = zlib.decompress(data[41 : 41 + length])
    rows = np.frombuffer(raw, dtype=np.uint8).reshape(20, 30 * 4 + 1)
    assert (rows[:, 0] == 0).all()
    assert np.array_equal(rows[:, 1:].reshape(20, 30, 4), image)


def test_screenshot_falls_back_without_polyscope(tmp_path, monkeypatch):
    from meshioplusplus import _viewer

    monkeypatch.setattr(_viewer, "has_viewer", lambda: False)
    out = tmp_path / "s.png"
    with pytest.warns(RuntimeWarning, match="software rasterizer"):
        mio.screenshot(_cube(), out, size=(64, 48), camera=((3, 3, 3), (0.5, 0.5, 0.5)))
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mio.screenshot(_cube(), tmp_path / "t.png", size=(64, 48), transparent=True)


def test_without_the_core_it_raises_by_name(monkeypatch):
    from meshioplusplus import _render

    monkeypatch.setitem(sys.modules, "meshioplusplus._core", None)
    monkeypatch.delattr(mio, "_core", raising=False)
    with pytest.raises(NotImplementedError, match="C\\+\\+-core only"):
        _render.render_image(_cube())


def test_mesh_kinds_render(tmp_path):
    for mesh in (
        helpers.tri_mesh,
        helpers.tet_mesh,
        helpers.line_mesh,
        helpers.hex_mesh,
    ):
        image = mio.render_image(mesh, 32, 32)
        assert image[..., 3].any()
    cloud = mio.Mesh(np.random.default_rng(0).random((20, 3)), [])
    assert mio.render_image(cloud, 32, 32)[..., 3].any()


# --- both CLIs ---------------------------------------------------------------


def _python_cli(args):
    return mio._cli.main(args) or 0


@pytest.mark.skipif(NATIVE is None, reason="no native CLI build found")
@pytest.mark.parametrize(
    "name, extra",
    [
        (
            "a.png",
            ["--width", "80", "--height", "60", "--supersample", "2", "--edges", "all"],
        ),
        ("a.ansi", ["--cols", "30", "--rows", "12", "--color-by", "x", "--colorbar"]),
        (
            "a.txt",
            ["--cols", "30", "--rows", "12", "--encoding", "braille", "--view", "+z"],
        ),
        (
            "a.html",
            [
                "--cols",
                "20",
                "--rows",
                "8",
                "--color-depth",
                "256",
                "--shading",
                "smooth",
            ],
        ),
    ],
)
def test_both_clis_write_the_same_bytes(tmp_path, name, extra):
    infile = tmp_path / "in.vtu"
    mio.write(infile, _cube())
    py_out = tmp_path / ("py_" + name)
    native_out = tmp_path / ("native_" + name)
    assert _python_cli(["snapshot", str(infile), str(py_out), *extra]) == 0
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), str(native_out), *extra],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert py_out.read_bytes() == native_out.read_bytes()


@pytest.mark.skipif(NATIVE is None, reason="no native CLI build found")
def test_both_clis_print_the_same_text(tmp_path, capsys):
    infile = tmp_path / "in.vtu"
    mio.write(infile, _cube())
    args = ["--cols", "24", "--rows", "8", "--color-depth", "16"]
    assert _python_cli(["snapshot", str(infile), "-", *args]) == 0
    py_text = capsys.readouterr().out
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), "-", *args], capture_output=True
    )
    assert r.returncode == 0
    assert py_text.encode("utf-8") == r.stdout


def test_cli_rejects_modifiers_without_color_by(tmp_path):
    infile = tmp_path / "in.vtu"
    mio.write(infile, _cube())
    with pytest.raises(SystemExit, match="--cmap requires --color-by"):
        _python_cli(["snapshot", str(infile), "-", "--cmap", "turbo"])


# --- field rendering (v16.34.0) ------------------------------------------------


def _grid(n=4):
    """An n x n grid of quads in z = 0 with point arrays x, y, v (a vector) and
    cell arrays id and mat, and a cell region for each half."""
    xs, ys = np.meshgrid(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1))
    pts = np.stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)], axis=1)
    quads = []
    for j in range(n):
        for i in range(n):
            p = j * (n + 1) + i
            quads.append([p, p + 1, p + n + 2, p + n + 1])
    mesh = mio.Mesh(pts, [("quad", np.array(quads))])
    mesh.point_data["x"] = pts[:, 0].copy()
    mesh.point_data["v"] = np.tile([0.0, 1.0, 0.0], (len(pts), 1))
    mesh.point_data["u"] = np.tile([0.0, 0.0, 0.25], (len(pts), 1))
    mesh.cell_data["mat"] = [np.arange(len(quads), dtype=float) % 2]
    return mesh


def _count(image, color):
    return int((image == np.array(color, dtype=np.uint8)).all(axis=-1).sum())


def test_clip_symmetric_and_log_report_what_they_did():
    mesh = _grid()
    _, plain = mio.render_image(mesh, 40, 40, color_by="x", return_info=True, view="+z")
    _, clipped = mio.render_image(
        mesh, 40, 40, color_by="x", clip=(30, 70), return_info=True, view="+z"
    )
    assert plain["vmin"] < clipped["vmin"] < clipped["vmax"] < plain["vmax"]
    _, sym = mio.render_image(
        mesh, 40, 40, color_by="x", symmetric=True, return_info=True, view="+z"
    )
    assert sym["vmin"] == -sym["vmax"]
    _, logged = mio.render_image(
        mesh, 40, 40, color_by="mat", scale="log", return_info=True, view="+z"
    )
    assert any("log" in n for n in logged["notes"])
    with pytest.raises(ValueError, match="percentiles"):
        mio.render_image(mesh, color_by="x", clip=(90, 10))
    with pytest.raises(ValueError, match="scale must be"):
        mio.render_image(mesh, color_by="x", scale="sqrt")


def test_isolines_vectors_and_warp_draw_in_their_colors():
    mesh = _grid(8)
    img = mio.render_image(
        mesh,
        80,
        80,
        color_by="x",
        isolines=3,
        iso_color=(255, 0, 255),
        cmap="grey",
        shading="none",
        view="+z",
    )
    assert _count(img, (255, 0, 255, 255)) > 80
    with pytest.raises(ValueError, match="point"):
        mio.render_image(mesh, color_by="mat", isolines=2)
    arrows = mio.render_image(
        mesh, 100, 100, vectors="v", vector_count=9, vector_color=(255, 0, 0), view="+z"
    )
    assert _count(arrows, (255, 0, 0, 255)) > 30
    mesh.point_data["swirl"] = np.stack(
        [
            -(mesh.points[:, 1] - 0.5),
            mesh.points[:, 0] - 0.5,
            np.zeros(len(mesh.points)),
        ],
        axis=1,
    )
    flow = mio.render_image(
        mesh,
        100,
        100,
        streamlines="swirl",
        stream_color=(255, 0, 0),
        view="+z",
        shading="none",
    )
    assert _count(flow, (255, 0, 0, 255)) > 100
    again = mio.render_image(
        mesh,
        100,
        100,
        streamlines="swirl",
        stream_color=(255, 0, 0),
        view="+z",
        shading="none",
    )
    assert (flow == again).all()
    few = mio.render_image(
        mesh,
        100,
        100,
        streamlines="swirl",
        stream_seeds=3,
        stream_color=(255, 0, 0),
        view="+z",
        shading="none",
    )
    assert 0 < _count(few, (255, 0, 0, 255)) < _count(flow, (255, 0, 0, 255))
    with pytest.raises(ValueError, match="streamline array 'nope'"):
        mio.render_image(mesh, streamlines="nope")
    with pytest.raises(ValueError, match="stream seeds"):
        mio.render_image(mesh, streamlines="swirl", stream_seeds=0)
    warped = mio.render_image(
        mesh, 60, 60, warp="u", warp_outline=True, outline_color=(0, 255, 0), view="+x"
    )
    assert _count(warped, (0, 255, 0, 255)) > 0
    with pytest.raises(ValueError, match="warp array 'nope'"):
        mio.render_image(mesh, warp="nope")


def test_categories_regions_and_diagnostics_carry_keys():
    mesh = _grid()
    mesh.regions = []
    _, info = mio.render_image(
        mesh, 40, 40, color_by="mat", categorical=True, return_info=True, view="+z"
    )
    assert any(n.startswith("mat = 0: #0072b2") for n in info["notes"])
    _, orient = mio.render_image(
        mesh, 40, 40, diagnostic="orientation", return_info=True, view="+z"
    )
    assert any(n.startswith("orientation:") for n in orient["notes"])
    _, qual = mio.render_image(
        mesh,
        40,
        40,
        diagnostic="quality",
        quality_metric="scaled_jacobian",
        return_info=True,
    )
    assert qual["colored"] is True
    with pytest.raises(ValueError, match="diagnostic must be"):
        mio.render_image(mesh, diagnostic="plastic")
    with pytest.raises(ValueError, match="excludes"):
        mio.render_image(mesh, diagnostic="orientation", color_by="x")


def test_expr_and_reduce_colour_like_arrays():
    mesh = _grid()
    _, doubled = mio.render_image(
        mesh, 30, 30, expr="x * 2", return_info=True, view="+z"
    )
    _, plain = mio.render_image(mesh, 30, 30, color_by="x", return_info=True, view="+z")
    assert doubled["vmax"] == pytest.approx(2 * plain["vmax"])
    assert doubled["notes"][0].startswith("x * 2:")
    mesh.point_data["stress"] = np.tile([5.0, 0, 0, 0, 0, 0], (len(mesh.points), 1))
    _, mises = mio.render_image(
        mesh, 30, 30, color_by="stress", reduce="mises", return_info=True, view="+z"
    )
    assert mises["vmax"] == pytest.approx(5.0)
    with pytest.raises(ValueError, match="either color_by or expr"):
        mio.render_image(mesh, color_by="x", expr="x * 2")
    with pytest.raises(ValueError, match="reduce must be"):
        mio.render_image(mesh, color_by="stress", reduce="trace")


def test_field_flags_match_between_the_two_clis(tmp_path):
    if NATIVE is None:
        pytest.skip("no native CLI build found")
    infile = tmp_path / "in.vtu"
    mio.write(infile, _grid(6))
    flags = [
        "--color-by",
        "x",
        "--clip",
        "5,95",
        "--symmetric",
        "--colorbar",
        "--isolines",
        "2",
        "--vectors",
        "v",
        "--vector-count",
        "9",
        "--warp",
        "u",
        "--warp-outline",
        "--view",
        "+z",
        "--width",
        "64",
        "--height",
        "48",
    ]
    py_out = tmp_path / "py.png"
    native_out = tmp_path / "native.png"
    assert _python_cli(["snapshot", str(infile), str(py_out), *flags]) == 0
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), str(native_out), *flags],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert py_out.read_bytes() == native_out.read_bytes()
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), "-", "--clip", "5,95"],
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0 and "requires --color-by" in r.stderr


# --- cut-aways (v16.36.0, ABI 25) ----------------------------------------------


def _centre(image, info):
    h, w = image.shape[:2]
    return int(info["cell_ids"][h // 2, w // 2]), tuple(
        int(c) for c in image[h // 2, w // 2]
    )


def test_a_cutaway_shows_the_tinted_inside_of_a_cube():
    cube = _cube()
    kw = dict(view="-x", shading="none", background=(0, 0, 0, 255), return_info=True)
    image, info = mio.render_image(cube, 64, 64, **kw)
    assert _centre(image, info)[0] == 4  # the x = 0 face
    for plane in (
        "+x:0.5",
        (0.5, 0, 0, 1, 0, 0),
        [(0.5, 0, 0, 1, 0, 0)],
        "0.5,0,0,1,0,0",
    ):
        image, info = mio.render_image(
            cube, 64, 64, cutaway=plane, cutaway_tint=(10, 200, 30), **kw
        )
        assert _centre(image, info) == (5, (10, 200, 30, 255)), plane
    assert any(n.startswith("cutaway: 1 plane") for n in info["notes"])


def test_cutaway_forms_are_validated():
    cube = _cube()
    with pytest.raises(ValueError, match="AXIS:OFFSET"):
        mio.render_image(cube, cutaway="x:1")
    with pytest.raises(ValueError, match="six numbers"):
        mio.render_image(cube, cutaway=(1, 2, 3))
    with pytest.raises(ValueError, match="at most two"):
        mio.render_image(cube, cutaway=["+x:0.1", "+y:0.1", "+z:0.1"])
    with pytest.raises(ValueError, match="non-zero normal"):
        mio.render_image(cube, cutaway=(0, 0, 0, 0, 0, 0))
    a = mio.render_image(cube, 32, 32, cutaway=["+x:0.5", "+y:0.5"], view="+z")
    b = mio.render_image(
        cube, 32, 32, cutaway=[(0.5, 0, 0, 1, 0, 0), "+y:0.5"], view="+z"
    )
    assert (a == b).all()


def test_a_cut_keeps_the_framing_of_the_whole_model():
    cube = _cube()
    whole, wi = mio.render_image(
        cube, 64, 64, view="+z", shading="none", return_info=True
    )
    cut, ci = mio.render_image(
        cube, 64, 64, view="+z", shading="none", cutaway="+x:0.5", return_info=True
    )
    ids_whole = wi["cell_ids"]
    ids_cut = ci["cell_ids"]
    # Where the cut picture draws the top face, the whole picture drew it too.
    drawn = ids_cut == 1
    assert drawn.any() and (ids_whole[drawn] == 1).all()
    assert 0.3 < drawn.sum() / (ids_whole == 1).sum() < 0.7


def test_cutaway_flags_match_between_the_two_clis(tmp_path):
    if NATIVE is None:
        pytest.skip("no native CLI build found")
    infile = tmp_path / "cube.vtu"
    mio.write(infile, _cube())
    flags = [
        "--cutaway",
        "+x:0.4",
        "--cutaway",
        "0,0.3,0,0,-1,0",
        "--cutaway-tint",
        "#00ff80",
        "--view",
        "iso",
        "--shading",
        "smooth",
        "--width",
        "64",
        "--height",
        "48",
    ]
    py_out = tmp_path / "py.png"
    native_out = tmp_path / "native.png"
    assert _python_cli(["snapshot", str(infile), str(py_out), *flags]) == 0
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), str(native_out), *flags],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert py_out.read_bytes() == native_out.read_bytes()
    plain = tmp_path / "plain.png"
    assert _python_cli(["snapshot", str(infile), str(plain), "--width", "64"]) == 0
    assert plain.read_bytes() != py_out.read_bytes()
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), "-", "--cutaway-tint", "#ff0000"],
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0 and "requires --cutaway" in r.stderr
    r = subprocess.run(
        [NATIVE, "snapshot", str(infile), "-", "--cutaway", "sideways:1"],
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0 and "AXIS:OFFSET" in r.stderr


# --- the synthwave theme (v16.38.0) ----------------------------------------------


def test_the_synthwave_theme_paints_a_banded_sunset_and_defaults_its_colormap():
    mesh = _grid()
    img = mio.render_image(mesh, 64, 96, theme="synthwave", view="+z", zoom=0.4)
    assert tuple(img[0, 0]) == (20, 8, 60, 255)  # the top of the sky, opaque
    assert len({tuple(img[y, 0, :3]) for y in range(96)}) <= 12  # banded
    # color_by with no cmap takes the theme's colormap; naming one wins.
    a = mio.render_image(mesh, 48, 48, color_by="x", theme="synthwave", view="+z")
    b = mio.render_image(
        mesh, 48, 48, color_by="x", theme="synthwave", cmap="synthwave", view="+z"
    )
    c = mio.render_image(
        mesh, 48, 48, color_by="x", theme="synthwave", cmap="viridis", view="+z"
    )
    assert (a == b).all() and not (a == c).all()
    opaque = mio.render_image(
        mesh, 64, 96, theme="synthwave", background=(1, 2, 3, 255), zoom=0.4
    )
    assert tuple(opaque[0, 0]) == (1, 2, 3, 255)


def test_the_theme_options_are_validated_and_the_post_processes_deterministic():
    mesh = _grid()
    with pytest.raises(ValueError, match="theme"):
        mio.render_image(mesh, grid_floor=True)
    with pytest.raises(ValueError, match="theme must be"):
        mio.render_image(mesh, theme="vaporwave")
    base = dict(theme="synthwave", color_by="x", view="+z")
    plain = mio.render_image(mesh, 96, 64, **base)
    for option in ("bloom", "fringe", "scanlines", "grid_floor"):
        out = mio.render_image(mesh, 96, 64, **{**base, option: True})
        assert not (out == plain).all(), option
        assert (out == mio.render_image(mesh, 96, 64, **{**base, option: True})).all()
    # The post-processes need no theme.
    assert not (
        mio.render_image(mesh, 96, 64, scanlines=True, view="+z")
        == mio.render_image(mesh, 96, 64, view="+z")
    ).all()
    phase0 = mio.render_image(mesh, 96, 64, grid_floor=True, **base)
    phase1 = mio.render_image(mesh, 96, 64, grid_floor=True, theme_phase=1, **base)
    assert not (phase0 == phase1).all()
