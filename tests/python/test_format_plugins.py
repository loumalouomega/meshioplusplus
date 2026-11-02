"""Third-party formats by entry point (roadmap 7.1.3): a separate package that
declares a ``meshioplusplus.formats`` entry point adds its format by being
installed -- here, a throwaway distribution on ``sys.path``."""

import importlib
import sys
import textwrap

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _helpers, _plugins

FOO_MODULE = """
import numpy as np

import meshioplusplus

REGISTERED = []


def _read(filename):
    with open(filename) as fh:
        rows = [line.split() for line in fh if line.strip()]
    pts = np.array([[float(v) for v in r[1:]] for r in rows if r[0] == "v"])
    tris = np.array([[int(v) for v in r[1:]] for r in rows if r[0] == "t"])
    return meshioplusplus.Mesh(pts, [("triangle", tris)])


def _write(filename, mesh, **kwargs):
    with open(filename, "w") as fh:
        for p in mesh.points:
            fh.write("v " + " ".join(repr(float(x)) for x in p) + "\\n")
        for c in mesh.cells[0].data:
            fh.write("t " + " ".join(str(int(i)) for i in c) + "\\n")


def register():
    REGISTERED.append(True)
    meshioplusplus.register_format("foo", [".foo"], _read, {"foo": _write})
"""

BROKEN_MODULE = """
def register():
    raise RuntimeError("this plugin is broken on purpose")
"""


def _make_dist(root, name, module_name, source, entry):
    (root / f"{module_name}.py").write_text(textwrap.dedent(source))
    info = root / f"{name.replace('-', '_')}-1.0.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(
        f"Metadata-Version: 2.1\nName: {name}\nVersion: 1.0\n"
    )
    (info / "entry_points.txt").write_text(f"[meshioplusplus.formats]\n{entry}\n")
    importlib.invalidate_caches()


@pytest.fixture
def plugin_path(tmp_path, monkeypatch):
    """A directory on sys.path holding the plugin distributions, and a fresh
    plugin state before and after."""
    monkeypatch.delenv("MESHIOPLUSPLUS_NO_PLUGINS", raising=False)
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.invalidate_caches()

    def reset():
        # Deregister first: it reads the registry, which discovers plugins.
        mp.deregister_format("foo")
        for mod in ("mio_foo_plugin", "mio_broken_plugin"):
            sys.modules.pop(mod, None)
        _plugins._format_state.update(done=False, loading=False, loaded={}, errors={})
        importlib.invalidate_caches()

    reset()
    yield tmp_path
    reset()


def test_an_installed_package_adds_its_format(plugin_path):
    _make_dist(
        plugin_path,
        "mio-foo",
        "mio_foo_plugin",
        FOO_MODULE,
        "foo = mio_foo_plugin:register",
    )
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    mesh = mp.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    # The probe: nothing but the installed entry point makes `.foo` work.
    path = plugin_path / "x.foo"
    mp.write(str(path), mesh)
    back = mp.read(str(path))
    np.testing.assert_allclose(back.points, pts)
    assert "foo" in mp.formats()["readable"]
    assert sys.modules["mio_foo_plugin"].REGISTERED == [True]
    # Loaded once, however often the registry is used.
    mp.read(str(path))
    assert sys.modules["mio_foo_plugin"].REGISTERED == [True]
    assert list(_plugins.load_format_plugins()) == ["foo"]


def test_a_broken_plugin_warns_and_is_skipped(plugin_path):
    _make_dist(
        plugin_path,
        "mio-broken",
        "mio_broken_plugin",
        BROKEN_MODULE,
        "broken = mio_broken_plugin:register",
    )
    _make_dist(
        plugin_path,
        "mio-foo",
        "mio_foo_plugin",
        FOO_MODULE,
        "foo = mio_foo_plugin:register",
    )
    with pytest.warns(RuntimeWarning, match="'broken'.*broken on purpose"):
        _helpers._ensure_plugins()
    assert "broken" in _plugins.format_plugin_errors()
    assert "foo" in _plugins.load_format_plugins()  # the others still load


def test_no_plugins_environment_skips_discovery(plugin_path, monkeypatch):
    _make_dist(
        plugin_path,
        "mio-foo",
        "mio_foo_plugin",
        FOO_MODULE,
        "foo = mio_foo_plugin:register",
    )
    monkeypatch.setenv("MESHIOPLUSPLUS_NO_PLUGINS", "1")
    assert _plugins.load_format_plugins() == {}
    with pytest.raises(mp.ReadError):
        mp.read(str(plugin_path / "nothing.foo"))
