"""Add-on discovery (``meshioplusplus._plugins``)."""

import sys
import types

import pytest

from meshioplusplus import _plugins
from meshioplusplus.__about__ import __version__


class _EntryPoint:
    def __init__(self, name, value):
        self.name, self.value = name, value


@pytest.fixture
def plugins(monkeypatch):
    saved = dict(_plugins._state)
    monkeypatch.delenv("MESHIOPLUSPLUS_NO_PLUGINS", raising=False)

    def install(*eps):
        monkeypatch.setattr(_plugins, "_entry_points", lambda: list(eps))
        return _plugins.load_plugins(force=True)

    yield install
    _plugins._state.clear()
    _plugins._state.update(saved)


def _module(monkeypatch, name, version):
    mod = types.ModuleType(name)
    mod.__version__ = version
    monkeypatch.setitem(sys.modules, name, mod)
    return mod


def test_matching_addon_is_loaded(plugins, monkeypatch):
    mod = _module(monkeypatch, "fake_addon_ok", __version__)
    assert plugins(_EntryPoint("occt", "fake_addon_ok")) == {"occt": mod}


def test_version_mismatch_is_refused_by_name(plugins, monkeypatch):
    _module(monkeypatch, "fake_addon_old", "0.0.1")
    with pytest.warns(UserWarning, match="fake_addon_old|occt"):
        assert plugins(_EntryPoint("occt", "fake_addon_old")) == {}
    assert "0.0.1" in _plugins.plugin_errors()["occt"]


def test_broken_addon_does_not_break_discovery(plugins):
    with pytest.warns(UserWarning, match="not loaded"):
        assert plugins(_EntryPoint("occt", "no_such_addon_module")) == {}


def test_opt_out(plugins, monkeypatch):
    _module(monkeypatch, "fake_addon_ok", __version__)
    monkeypatch.setenv("MESHIOPLUSPLUS_NO_PLUGINS", "1")
    assert plugins(_EntryPoint("occt", "fake_addon_ok")) == {}
