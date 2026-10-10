"""The interactive terminal viewer (roadmap 7.2, 7.3): ``tui``, ``view(backend=...)``.

Two kinds of test. *Replays* feed a recorded input stream to the loop on a
screen of a stated size (no terminal needed) and compare what it wrote, byte for
byte between the C++ CLI, the Python CLI and ``_core.tui`` and through a small
terminal emulator (``tools/ansi_screen.py``). *Pty tests* run the real binary on
a pseudo-terminal and check what a user would see: the terminal comes back in
cooked mode with the cursor visible however the viewer ends, and a resize never
leaves a picture bigger than the screen.
"""

import fcntl
import importlib.util
import os
import pathlib
import select
import signal
import subprocess
import sys
import termios
import time

import numpy as np
import pytest

import meshioplusplus

pytest.importorskip("meshioplusplus._core")

REPO = pathlib.Path(__file__).resolve().parents[2]


def _load_screen():
    spec = importlib.util.spec_from_file_location(
        "ansi_screen", REPO / "tools" / "ansi_screen.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["ansi_screen"] = (
        module  # dataclasses resolve annotations via sys.modules
    )
    spec.loader.exec_module(module)
    return module


ansi_screen = _load_screen()


def _native_cli():
    env = os.environ.get("MESHIOPLUSPLUS_NATIVE_CLI")
    if env and os.path.exists(env):
        return env
    found = sorted(
        (REPO / "build").glob("*/meshioplusplus"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return str(found[0]) if found else None


NATIVE = _native_cli()
posix_only = pytest.mark.skipif(
    sys.platform == "win32", reason="pseudo-terminals are POSIX"
)
needs_native = pytest.mark.skipif(NATIVE is None, reason="native CLI is not built")

CLEAR = b"\x1b[2J"
LEAVE = b"\x1b[?1049l"
SHOW_CURSOR = b"\x1b[?25h"
ENTER = b"\x1b[?1049h"


@pytest.fixture(scope="module")
def mesh_file(tmp_path_factory):
    """A small coloured cube of quads."""
    path = tmp_path_factory.mktemp("tui") / "cube.vtu"
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        float,
    )
    quads = np.array(
        [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [3, 7, 6, 2],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
        ]
    )
    mesh = meshioplusplus.Mesh(
        pts, [("quad", quads)], point_data={"u": np.arange(8, dtype=float)}
    )
    meshioplusplus.write(str(path), mesh)
    return str(path)


def _mouse(code, x, y, release=False):
    return f"\x1b[<{code};{x};{y}{'m' if release else 'M'}".encode()


_CLI = "import sys, meshioplusplus._cli as c; sys.exit(c.main(sys.argv[1:]) or 0)"


def _python_cli(*args):
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    return subprocess.run(
        [sys.executable, "-c", _CLI, *args],
        capture_output=True,
        env=env,
        timeout=60,
    )


# --------------------------------------------------------------------------
# Replays
# --------------------------------------------------------------------------

DRAG = (
    _mouse(0, 20, 10)
    + _mouse(32, 30, 10)
    + _mouse(32, 40, 14)
    + _mouse(0, 40, 14, release=True)
    + b"+q"
)


def _replay_core(mesh_file, script, cols=80, rows=24, **options):
    mesh = meshioplusplus.read(mesh_file)
    return meshioplusplus.tui(
        mesh,
        replay=script,
        cols=cols,
        rows=rows,
        color_depth="truecolor",
        title="cube.vtu",
        **options,
    )


def test_a_replay_draws_a_picture_and_a_status_line(mesh_file):
    result = _replay_core(mesh_file, b"q", color_by="u")
    screen = ansi_screen.replay(result["output"], 80, 24)
    text = screen.text()
    assert result["exit"] == 0
    assert sum(1 for line in text[:-1] if line) > 5  # the picture
    assert "cube.vtu" in text[-1] and "az 45" in text[-1] and "q quit" in text[-1]
    # The status line is reverse video across the whole width.
    assert screen.cell(23, 79).reverse


def test_a_replay_is_deterministic(mesh_file):
    first = _replay_core(mesh_file, DRAG, color_by="u")
    second = _replay_core(mesh_file, DRAG, color_by="u")
    assert first["output"] == second["output"]
    assert first["frames"] >= 2


def test_dragging_and_zooming_change_the_camera_and_the_picture(mesh_file):
    still = _replay_core(mesh_file, b"q")
    moved = _replay_core(mesh_file, DRAG)
    assert moved["azimuth"] != still["azimuth"]
    assert moved["zoom"] > 1.0
    final = ansi_screen.replay(moved["output"], 80, 24).text()
    initial = ansi_screen.replay(still["output"], 80, 24).text()
    assert final != initial


def test_the_view_keys_set_the_named_views(mesh_file):
    for key, azimuth, elevation in (
        (b"2", 0.0, 0.0),
        (b"3", 180.0, 0.0),
        (b"6", -90.0, 90.0),
    ):
        result = _replay_core(mesh_file, key + b"q")
        assert (result["azimuth"], result["elevation"]) == (azimuth, elevation)


def test_reset_returns_to_the_start_and_help_is_shown(mesh_file):
    result = _replay_core(mesh_file, b"2+r?q")
    assert (result["azimuth"], result["zoom"]) == (45.0, 1.0)
    screen = ansi_screen.replay(result["output"], 80, 24)
    assert any("this help" in line for line in screen.text())


def test_the_viewer_never_draws_outside_a_small_screen(mesh_file):
    result = _replay_core(mesh_file, b"q", cols=30, rows=10)
    screen = ansi_screen.replay(result["output"], 30, 10)
    assert screen.non_blank_rows() <= 10
    result = _replay_core(mesh_file, b"q", cols=5, rows=3)
    assert b"terminal too small" in result["output"]


def test_field_options_reach_the_viewer(mesh_file):
    plain = _replay_core(mesh_file, b"q")
    colored = _replay_core(mesh_file, b"q", color_by="u", cmap="magma", colorbar=True)
    assert plain["output"] != colored["output"]
    screen = ansi_screen.replay(colored["output"], 80, 24)
    assert any(line.startswith("u: ") for line in screen.text())  # the colour range


def test_a_bad_option_is_a_value_error(mesh_file):
    with pytest.raises(ValueError):
        _replay_core(mesh_file, b"q", color_by="no_such_array")
    with pytest.raises(ValueError):
        _replay_core(mesh_file, b"q", color_by="u", cmap="no_such_map")


def test_tui_without_a_terminal_refuses_by_name(mesh_file):
    with pytest.raises(RuntimeError, match="snapshot"):
        meshioplusplus.tui(meshioplusplus.read(mesh_file))


def test_tui_reads_a_path_and_titles_it(mesh_file):
    result = meshioplusplus.tui(mesh_file, replay=b"q", cols=60, rows=14)
    assert b"cube.vtu" in result["output"]


def test_the_cli_replays_match_the_core_byte_for_byte(mesh_file, tmp_path):
    script = tmp_path / "script.bin"
    script.write_bytes(DRAG)
    expected = _replay_core(mesh_file, DRAG, color_by="u")["output"]
    py = _python_cli(
        "tui",
        mesh_file,
        "--replay",
        str(script),
        "--cols",
        "80",
        "--rows",
        "24",
        "--color-depth",
        "truecolor",
        "--color-by",
        "u",
    )
    assert py.returncode == 0, py.stderr.decode()
    assert py.stdout == expected
    if NATIVE:
        native = subprocess.run(
            [
                NATIVE,
                "tui",
                mesh_file,
                "--replay",
                str(script),
                "--cols",
                "80",
                "--rows",
                "24",
                "--color-depth",
                "truecolor",
                "--color-by",
                "u",
            ],
            capture_output=True,
            timeout=60,
        )
        assert native.returncode == 0, native.stderr.decode()
        assert native.stdout == expected


# --------------------------------------------------------------------------
# view(backend=...)
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "display, viewer, terminal, expected",
    [
        (True, True, True, "polyscope"),
        (True, False, True, "browser"),  # a window can open: nothing changes
        (True, True, False, "polyscope"),
        (False, True, True, "terminal"),  # no display, a terminal: the new case
        (False, False, True, "terminal"),
        (False, True, False, "browser"),  # a notebook or a pipe
        (False, False, False, "browser"),
    ],
)
def test_auto_backend_table(monkeypatch, display, viewer, terminal, expected):
    from meshioplusplus import _viewer

    monkeypatch.setattr(_viewer, "_has_display", lambda: display)
    monkeypatch.setattr(_viewer, "has_viewer", lambda: viewer)
    monkeypatch.setattr(_viewer, "_terminal_available", lambda: terminal)
    assert _viewer._auto_backend() == expected


def test_view_rejects_an_unknown_backend_and_accepts_terminal(mesh_file):
    from meshioplusplus import _viewer

    with pytest.raises(ValueError, match="unknown backend"):
        _viewer.view(meshioplusplus.read(mesh_file), backend="opengl")
    # terminal is a backend; without a tty it refuses by name rather than ValueError
    with pytest.raises(RuntimeError, match="snapshot"):
        _viewer.view(meshioplusplus.read(mesh_file), backend="terminal")


def test_the_cli_accepts_terminal_as_a_backend(mesh_file):
    result = _python_cli("view", mesh_file, "--backend", "terminal")
    assert result.returncode != 0
    assert b"snapshot" in result.stderr or b"terminal" in result.stderr
    result = _python_cli("view", mesh_file, "--backend", "opengl")
    assert result.returncode == 2  # argparse still rejects what it does not know


# --------------------------------------------------------------------------
# Pseudo-terminals
# --------------------------------------------------------------------------


class Pty:
    """A process on a pseudo-terminal of a stated size."""

    def __init__(self, argv, rows=24, cols=80, env=None):
        self.master, slave = os.openpty()
        self.set_size(rows, cols)
        full_env = dict(os.environ if env is None else env)
        full_env.setdefault("TERM", "xterm-256color")
        full_env.setdefault("COLORTERM", "truecolor")
        self.proc = subprocess.Popen(
            argv,
            stdin=slave,
            stdout=slave,
            stderr=slave,
            start_new_session=True,
            preexec_fn=lambda: fcntl.ioctl(0, termios.TIOCSCTTY, 0),
            env=full_env,
        )
        os.close(slave)
        self.output = b""

    def set_size(self, rows, cols):
        fcntl.ioctl(
            self.master,
            termios.TIOCSWINSZ,
            rows.to_bytes(2, "little")
            + cols.to_bytes(2, "little")
            + (0).to_bytes(4, "little"),
        )

    def pump(self, seconds=0.2):
        end = time.time() + seconds
        while time.time() < end:
            ready, _, _ = select.select([self.master], [], [], 0.05)
            if ready:
                try:
                    data = os.read(self.master, 65536)
                except OSError:
                    return
                if not data:
                    return
                self.output += data

    def wait_for(self, marker, timeout=15):
        end = time.time() + timeout
        while time.time() < end:
            if marker in self.output:
                return True
            self.pump(0.1)
            if self.proc.poll() is not None:
                self.pump(0.2)
                return marker in self.output
        return False

    def send(self, data):
        os.write(self.master, data)

    def finish(self, timeout=15):
        end = time.time() + timeout
        while self.proc.poll() is None and time.time() < end:
            self.pump(0.1)
        self.pump(0.2)
        return self.proc.poll()

    def cooked(self):
        lflag = termios.tcgetattr(self.master)[3]
        return bool(lflag & termios.ICANON) and bool(lflag & termios.ECHO)

    def close(self):
        if self.proc.poll() is None:
            self.proc.kill()
            self.proc.wait()
        try:
            os.close(self.master)
        except OSError:
            pass


@pytest.fixture
def pty_run():
    started = []

    def run(*args, **kwargs):
        p = Pty(*args, **kwargs)
        started.append(p)
        return p

    yield run
    for p in started:
        p.close()


def _tui_argv(mesh_file, *extra):
    return [NATIVE, "tui", mesh_file, *extra]


def _python_tui_argv(mesh_file, *extra):
    return [sys.executable, "-c", _CLI, "tui", mesh_file, *extra]


def _both():
    params = [pytest.param("python", id="python")]
    if NATIVE:
        params.append(pytest.param("native", id="native"))
    return params


def _argv_for(kind, mesh_file, *extra):
    return (
        _tui_argv(mesh_file, *extra)
        if kind == "native"
        else _python_tui_argv(mesh_file, *extra)
    )


@posix_only
@pytest.mark.parametrize("kind", _both())
def test_quit_restores_the_terminal(pty_run, mesh_file, kind):
    p = pty_run(_argv_for(kind, mesh_file))
    assert p.wait_for(ENTER)
    assert p.wait_for(b"q quit")
    assert not p.cooked()  # raw while the viewer runs
    p.send(b"q")
    assert p.finish() == 0
    assert p.cooked()
    assert LEAVE in p.output and SHOW_CURSOR in p.output
    # The alternate screen was entered before it was left.
    assert p.output.index(ENTER) < p.output.index(LEAVE)


@posix_only
@pytest.mark.parametrize("kind", _both())
@pytest.mark.parametrize(
    "sig, code",
    [(signal.SIGTERM, 128 + 15), (signal.SIGHUP, 128 + 1), (signal.SIGINT, 128 + 2)],
)
def test_a_signal_restores_the_terminal(pty_run, mesh_file, kind, sig, code):
    p = pty_run(_argv_for(kind, mesh_file))
    assert p.wait_for(b"q quit")
    p.proc.send_signal(sig)
    status = p.finish()
    assert p.cooked()
    assert LEAVE in p.output and SHOW_CURSOR in p.output
    if kind == "native":
        assert status == code
    else:
        # Python reports SIGINT as KeyboardInterrupt (exit 1 with a traceback or
        # the interpreter's own status) and the others as the loop's exit code.
        assert status != 0 or sig == signal.SIGHUP


@posix_only
@needs_native
def test_end_of_input_ends_the_session_cleanly(pty_run, mesh_file):
    p = pty_run(_tui_argv(mesh_file))
    assert p.wait_for(b"q quit")
    p.send(b"\x04")  # Ctrl-D
    assert p.finish() == 0
    assert p.cooked()


@posix_only
@needs_native
def test_a_resize_redraws_within_the_new_size(pty_run, mesh_file):
    p = pty_run(_tui_argv(mesh_file, "--color-depth", "truecolor"), rows=30, cols=100)
    assert p.wait_for(b"q quit")
    p.pump(0.3)
    p.output = b""
    p.set_size(10, 40)  # the kernel sends SIGWINCH to the foreground group
    assert p.wait_for(b"cube.vtu", timeout=10)  # the status line, cut to 40 columns
    p.pump(0.3)
    p.send(b"q")
    assert p.finish() == 0
    screen = ansi_screen.replay(p.output, 40, 10)
    assert screen.non_blank_rows() <= 10
    assert "cube.vtu" in screen.text()[-1]


@posix_only
@needs_native
def test_keys_and_mouse_reach_the_loop(pty_run, mesh_file):
    p = pty_run(_tui_argv(mesh_file, "--color-depth", "truecolor"))
    assert p.wait_for(b"q quit")
    p.pump(0.2)
    p.output = b""
    p.send(b"3")  # the -x view
    assert p.wait_for(b"view -x")
    p.send(_mouse(65, 10, 10))  # wheel down: zoom out
    p.pump(0.3)
    p.send(b"q")
    assert p.finish() == 0
    assert b"view -x" in p.output


@posix_only
@needs_native
def test_without_a_terminal_the_native_viewer_names_snapshot(mesh_file):
    result = subprocess.run(
        [NATIVE, "tui", mesh_file], capture_output=True, stdin=subprocess.DEVNULL
    )
    assert result.returncode != 0
    assert b"snapshot" in result.stderr
    assert ENTER not in result.stdout  # nothing about the terminal changed


@posix_only
def test_auto_picks_the_terminal_on_a_tty_without_a_display_and_the_browser_on_a_pipe(
    pty_run,
):
    code = "from meshioplusplus import _viewer; print('BACKEND=' + _viewer._auto_backend())"
    env = {
        k: v for k, v in os.environ.items() if k not in ("DISPLAY", "WAYLAND_DISPLAY")
    }
    if sys.platform == "darwin":
        pytest.skip("macOS always has a window system")
    p = pty_run([sys.executable, "-c", code], env=env)
    assert p.finish() == 0
    assert b"BACKEND=terminal" in p.output
    piped = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, env=env, timeout=60
    )
    assert b"BACKEND=browser" in piped.stdout


@posix_only
def test_view_on_a_tty_takes_over_the_terminal_and_quits(pty_run, mesh_file):
    code = (
        "import meshioplusplus as m; "
        f"m.view(m.read({mesh_file!r}), backend='terminal'); print('RETURNED')"
    )
    p = pty_run([sys.executable, "-c", code])
    assert p.wait_for(b"q quit")
    p.send(b"q")
    assert p.finish() == 0
    assert b"RETURNED" in p.output
    assert p.cooked()
    assert p.output.index(LEAVE) < p.output.index(b"RETURNED")


# --------------------------------------------------------------------------
# Probing, commands, sessions, comparison, series (roadmap 7.2)
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def series_files(tmp_path_factory):
    """Four steps of a cube whose field grows with the step."""
    folder = tmp_path_factory.mktemp("series")
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        float,
    )
    quads = np.array(
        [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [3, 7, 6, 2],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
        ]
    )
    paths = []
    for k in range(4):
        mesh = meshioplusplus.Mesh(
            pts,
            [("quad", quads)],
            point_data={"u": np.arange(8, dtype=float) * (k + 1)},
        )
        path = folder / f"out_{k + 1}.vtu"
        meshioplusplus.write(str(path), mesh)
        paths.append(str(path))
    return folder, paths


def test_a_click_probes_a_cell_and_a_pin_keeps_it(mesh_file):
    click = _mouse(0, 40, 10) + _mouse(0, 40, 10, release=True)
    result = _replay_core(mesh_file, click + b"iq", view="+z", color_by="u")
    assert result["probe"][0].startswith("cell ")
    assert any(line.startswith("pin A") for line in result["probe"])
    assert any("u " in line for line in result["probe"])
    screen = ansi_screen.replay(result["output"], 80, 24)
    assert any(line.startswith("cell ") for line in screen.text())
    cleared = _replay_core(mesh_file, click + b"i0q", view="+z", color_by="u")
    assert cleared["probe"] == []


def test_the_command_line_speaks_the_snapshot_flags(mesh_file):
    result = _replay_core(mesh_file, b":bogus\r", color_by="u")
    assert "unknown command" in result["status"]
    result = _replay_core(mesh_file, b":cmap magma\r", color_by="u")
    plain = _replay_core(mesh_file, b"q", color_by="u")
    assert result["output"] != plain["output"]
    result = _replay_core(mesh_file, b":zoom 2\r", color_by="u")
    assert result["zoom"] == 2.0
    result = _replay_core(mesh_file, b":view +x\r", color_by="u")
    assert (result["azimuth"], result["elevation"]) == (0.0, 0.0)
    result = _replay_core(mesh_file, b":cmap magma\x1b", color_by="u")  # Escape cancels
    assert result["output"] == plain["output"] or "preparing" not in result["status"]


def test_cutaway_keys_and_options_clip_the_picture(mesh_file):
    plain = _replay_core(mesh_file, b"q")
    keyed = _replay_core(mesh_file, b"xq")
    assert keyed["output"] != plain["output"]
    assert "cut 1" in keyed["status"]
    optioned = _replay_core(mesh_file, b"q", cutaway="+x:0.5")
    assert "cut 1" in optioned["status"]
    commanded = _replay_core(mesh_file, b":clip +x 0.5\r")
    assert "cut 1" in commanded["status"]
    assert (
        _replay_core(mesh_file, b":clip +x 0.5\r:clip off\r")["status"].count("cut")
        == 0
    )


def test_a_session_file_is_written_and_read_back(mesh_file, tmp_path):
    session = tmp_path / "view.json"
    first = _replay_core(mesh_file, b":zoom 2.5\r", color_by="u", session=str(session))
    assert first["zoom"] == 2.5
    text = session.read_text()
    assert '"version": 1' in text and "--zoom=2.5" in text
    again = _replay_core(mesh_file, b"q", color_by="u", session=str(session))
    assert again["zoom"] == 2.5
    session.write_text('{"version": 1, "mystery": 3}')
    with pytest.raises(ValueError, match="mystery"):
        _replay_core(mesh_file, b"q", session=str(session))


def test_a_series_is_stepped_and_followed_by_name(series_files):
    folder, paths = series_files
    glob = str(folder / "out_*.vtu")
    result = meshioplusplus.tui(
        series=glob,
        replay=b"]]q",
        cols=100,
        rows=24,
        color_depth="truecolor",
        color_by="u",
    )
    assert result["step"] == 2
    assert "step 3/4" in result["status"]
    listed = meshioplusplus.tui(
        series=paths, replay=b"}q", cols=100, rows=24, color_depth="truecolor"
    )
    assert listed["step"] == 3
    # The colour range is the first step's: the last step's own range never shows.
    first = meshioplusplus.tui(
        series=glob,
        replay=b"q",
        cols=100,
        rows=24,
        color_depth="truecolor",
        color_by="u",
    )
    later = meshioplusplus.tui(
        series=glob,
        replay=b"}q",
        cols=100,
        rows=24,
        color_depth="truecolor",
        color_by="u",
    )
    s_first = "\n".join(ansi_screen.replay(first["output"], 100, 24).text())
    s_later = "\n".join(ansi_screen.replay(later["output"], 100, 24).text())
    note = [line for line in s_first.splitlines() if line.startswith("u: ")][0]
    assert note in s_later
    from meshioplusplus._exceptions import ReadError

    with pytest.raises(ReadError, match="matched no files"):
        meshioplusplus.tui(series=str(folder / "nothing_*.vtu"), replay=b"q")
    with pytest.raises(TypeError):
        meshioplusplus.tui(replay=b"q")


def test_two_meshes_side_by_side_and_their_difference(mesh_file, series_files):
    _, paths = series_files
    both = meshioplusplus.tui(
        paths[0],
        compare=paths[3],
        replay=b"q",
        cols=101,
        rows=24,
        color_depth="truecolor",
        color_by="u",
    )
    screen = ansi_screen.replay(both["output"], 101, 24)
    assert all(screen.cell(r, 50).char == "\u2502" for r in range(20))  # (101 - 1) // 2
    note = "\n".join(screen.text())
    assert "out_1.vtu" in note and "out_4.vtu" in note
    diff = meshioplusplus.tui(
        paths[0],
        compare=paths[3],
        diff=True,
        replay=b"q",
        cols=101,
        rows=24,
        color_depth="truecolor",
        color_by="u",
    )
    assert diff["output"] != both["output"]
    with pytest.raises(ValueError, match="difference"):
        meshioplusplus.tui(
            paths[0], compare=paths[3], diff=True, replay=b"q", color_by="missing"
        )
    with pytest.raises(ValueError, match="cell encoding"):
        meshioplusplus.tui(
            paths[0], compare=paths[3], replay=b"q", encoding="kitty", color_by="u"
        )


def test_the_cli_replays_of_a_series_and_a_comparison_match(series_files, tmp_path):
    if NATIVE is None:
        pytest.skip("native CLI is not built")
    folder, paths = series_files
    script = tmp_path / "keys.bin"
    script.write_bytes(b"]x:cmap magma\r\x1b")
    common = [
        "--replay",
        str(script),
        "--cols",
        "101",
        "--rows",
        "30",
        "--color-depth",
        "truecolor",
        "--color-by",
        "u",
    ]
    for args in (
        [str(folder / "out_*.vtu")],
        paths[:3],
        [paths[0], "--compare", paths[3]],
        [paths[0], "--compare", paths[3], "--diff"],
    ):
        py = _python_cli("tui", *args, *common)
        native = subprocess.run(
            [NATIVE, "tui", *args, *common], capture_output=True, timeout=60
        )
        assert py.returncode == native.returncode == 0, (py.stderr, native.stderr)
        assert py.stdout == native.stdout, args


@posix_only
@pytest.mark.parametrize("kind", _both())
def test_the_command_line_works_on_a_real_terminal(pty_run, mesh_file, kind):
    p = pty_run(_argv_for(kind, mesh_file, "--color-depth", "truecolor"))
    assert p.wait_for(b"q quit")
    p.send(b":zoom 2\r")
    assert p.wait_for(b"zoom 2.00")
    p.send(b":nonsense\r")
    assert p.wait_for(b"unknown command")
    p.send(b"q")
    assert p.finish() == 0
    assert p.cooked()


@posix_only
@pytest.mark.parametrize("kind", _both())
def test_following_a_series_on_a_real_terminal(pty_run, series_files, tmp_path, kind):
    _, paths = series_files
    live = tmp_path / "live"
    live.mkdir()
    for k in (0, 1):
        meshioplusplus.write(
            str(live / f"run_{k + 1}.vtu"), meshioplusplus.read(paths[k])
        )
    p = pty_run(
        _argv_for(
            kind,
            str(live / "run_*.vtu"),
            "--follow",
            "--follow-interval",
            "100",
            "--settle",
            "100",
            "--color-by",
            "u",
        )
    )
    assert p.wait_for(b"step 2/2")  # a followed series starts on the newest step
    meshioplusplus.write(str(live / "run_3.vtu"), meshioplusplus.read(paths[2]))
    assert p.wait_for(b"step 3/3", timeout=20)
    assert p.wait_for(b"followed")
    p.send(b"q")
    assert p.finish() == 0
    assert p.cooked()
