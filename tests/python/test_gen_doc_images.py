"""Smoke test for ``tools/gen_doc_images.py``.

The committed doc figures are generated, not hand-drawn, so the generator has
to keep working. This runs it against the bundled sample into a temp directory
and asserts the expected files come out non-empty -- it never writes into
``doc/``.
"""

import importlib.util
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO / "tools" / "gen_doc_images.py"


@pytest.fixture(scope="module")
def gen():
    spec = importlib.util.spec_from_file_location("gen_doc_images", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_lists_its_outputs(gen):
    # --list mode touches nothing and still names every figure.
    names = [p.name for p in gen.vector_figures(dry_run=True)]
    assert names == [
        "color_by_quality.svg",
        "color_by_partition.svg",
        "color_by_diff_error.svg",
    ]


def test_vector_figures_are_written_and_nonempty(gen, tmp_path, monkeypatch):
    monkeypatch.setattr(gen, "IMAGES", tmp_path)
    written = gen.vector_figures(dry_run=False)
    assert written
    for path in written:
        assert path.exists(), path
        assert path.stat().st_size > 0, path
        text = path.read_text()
        # Each figure is a real coloured SVG with a colorbar, not a blank frame.
        assert text.startswith("<svg ")
        assert 'fill="#' in text
        assert "<rect " in text


def test_raster_figures_need_polyscope(gen, tmp_path, monkeypatch):
    pytest.importorskip("polyscope")
    monkeypatch.setattr(gen, "VIEWER", tmp_path)
    written = gen.raster_figures(dry_run=False)
    assert written
    for path in written:
        assert path.exists() and path.stat().st_size > 0, path


TERMINAL = [
    "tui_loop.svg",
    "tui_help.svg",
    "tui_braille.svg",
    "tui_probe.svg",
    "tui_command.svg",
    "tui_cutaway.svg",
    "tui_compare.svg",
    "tui_series.svg",
    "tui_streamlines.svg",
]


def test_terminal_figures_are_listed(gen):
    assert [p.name for p in gen.terminal_figures(dry_run=True)] == TERMINAL


def test_terminal_figures_are_deterministic_screens(gen, tmp_path, monkeypatch):
    pytest.importorskip("meshioplusplus._core")
    monkeypatch.setattr(gen, "IMAGES", tmp_path)
    first = {p.name: p.read_bytes() for p in gen.terminal_figures(dry_run=False)}
    assert list(first) == TERMINAL
    for name, data in first.items():
        text = data.decode()
        assert text.startswith("<svg ")
        assert "<rect " in text
        assert "q quit" in text or "this help" in text or ":" in text  # a status line
    assert b"<circle" in first["tui_braille.svg"]  # Braille dots are drawn as shapes
    assert "cutaway: 1 plane" in first["tui_cutaway.svg"].decode()
    again = {p.name: p.read_bytes() for p in gen.terminal_figures(dry_run=False)}
    assert again == first  # the same bytes twice: a replay has no clock in it
