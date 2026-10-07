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
