"""Every per-format shim honours ``core_declined``, not just the three tested by hand.

``test_no_broad_core_swallow.py`` proves each shim *asks* ``core_declined``; that
is a statement about syntax. This drives all of them (90 ``read``/``write``
functions across 47 format packages) with a stand-in for the compiled ``_core``
whose every entry point raises, and pins the outcome of each:

* a ``ReadError``/``WriteError`` falls back to the Python twin;
* a ``TypeError`` (never a decline) propagates and the twin is not called;
* under ``MESHIOPLUSPLUS_STRICT_CORE=1`` the core's own error is re-raised.

Each case also asserts the stand-in was actually reached. Without that, a shim
whose pre-check skips the core would "fall back" trivially and prove nothing.
"""

import ast
import importlib
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import ReadError, WriteError
from meshioplusplus._fallback import set_strict_core

PKG = pathlib.Path(meshioplusplus.__file__).resolve().parent
_SKIP_DIRS = {"_cli", "mcp", "_viewer_assets", "_cxml", "__pycache__"}

# The one shim whose Python twin is not spelled `_py_<read|write>`.
_TWIN = {("vtk", "write"): "_main_write"}

# Entry points that take more than the defaults to reach every `_core` call.
_VARIANTS = {("gmsh", "write"): [{"fmt_version": "2.2"}, {"fmt_version": "4.1"}]}


def _calls_core(node):
    return any(
        isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and isinstance(n.func.value, ast.Name)
        and n.func.value.id == "_core"
        for n in ast.walk(node)
    )


def _discover():
    """(package, "read"|"write", extra kwargs) for every shim that wraps a `_core` call."""
    found = []
    for init in sorted(PKG.glob("*/__init__.py")):
        if init.parent.name in _SKIP_DIRS:
            continue
        for fn in ast.parse(init.read_text(encoding="utf-8")).body:
            if not (isinstance(fn, ast.FunctionDef) and fn.name in ("read", "write")):
                continue
            if any(
                isinstance(t, ast.Try) and any(_calls_core(s) for s in t.body)
                for t in ast.walk(fn)
            ):
                for kwargs in _VARIANTS.get((init.parent.name, fn.name), [{}]):
                    found.append((init.parent.name, fn.name, kwargs))
    return found


SHIMS = _discover()
_IDS = [
    f"{pkg}.{fn}" + (f"[{kw['fmt_version']}]" if kw else "") for pkg, fn, kw in SHIMS
]

_MESH = meshioplusplus.Mesh(
    np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    [("triangle", np.array([[0, 1, 2]]))],
)


class _RaisingCore:
    """A stand-in for ``meshioplusplus._core`` whose entry points all raise."""

    def __init__(self, exc):
        self._exc = exc
        self.calls = []

    def __getattr__(self, name):
        if name.startswith("__has_"):
            return True  # build-capability probes: say every codec/library exists
        if name.startswith("__"):
            raise AttributeError(name)

        def entry(*args, **kwargs):
            self.calls.append(name)
            raise self._exc

        return entry


@pytest.fixture(autouse=True)
def _reset_strict(monkeypatch):
    monkeypatch.delenv("MESHIOPLUSPLUS_STRICT_CORE", raising=False)
    set_strict_core(None)
    yield
    set_strict_core(None)


def _arm(monkeypatch, pkg, fn, exc):
    """Swap in the raising core and a recording twin; return (module, core, twin_calls)."""
    mod = importlib.import_module(f"meshioplusplus.{pkg}")
    core = _RaisingCore(exc)
    twin_calls = []
    monkeypatch.setattr(mod, "_core", core)
    monkeypatch.setattr(
        mod, _TWIN.get((pkg, fn), f"_py_{fn}"), lambda *a, **k: twin_calls.append(1)
    )
    for name in dir(mod):
        if name.startswith("_HAS_"):  # optional-library gates decided at import time
            monkeypatch.setattr(mod, name, True)
    if (
        pkg == "cgns"
    ):  # an ADF-container file must not fall back; pretend this one is HDF5
        monkeypatch.setattr(mod, "_is_hdf5", lambda filename: True)
    return mod, core, twin_calls


def _call(mod, fn, **kwargs):
    if fn == "read":
        return mod.read("model.dat", **kwargs)
    return mod.write("model.out", _MESH, **kwargs)


def _decline(fn):
    return ReadError("declined") if fn == "read" else WriteError("declined")


def test_the_discovery_finds_every_shim():
    # 44 read + 46 write functions, and gmsh.write is run down both of its paths.
    assert len(SHIMS) >= 90, len(SHIMS)
    assert {"gmsh", "vtu", "xdmf", "cgns", "exodus"} <= {pkg for pkg, _, _ in SHIMS}


@pytest.mark.parametrize("pkg, fn, kwargs", SHIMS, ids=_IDS)
def test_a_recognised_decline_falls_back_to_the_python_twin(
    monkeypatch, pkg, fn, kwargs
):
    mod, core, twin_calls = _arm(monkeypatch, pkg, fn, _decline(fn))
    _call(mod, fn, **kwargs)
    assert core.calls, f"{pkg}.{fn} never reached the core; the stand-in proved nothing"
    assert twin_calls == [1]


@pytest.mark.parametrize("pkg, fn, kwargs", SHIMS, ids=_IDS)
def test_an_exception_that_is_never_a_decline_propagates(monkeypatch, pkg, fn, kwargs):
    mod, core, twin_calls = _arm(monkeypatch, pkg, fn, TypeError("stale build"))
    with pytest.raises(TypeError, match="stale build"):
        _call(mod, fn, **kwargs)
    assert core.calls
    assert twin_calls == []


@pytest.mark.parametrize("pkg, fn, kwargs", SHIMS, ids=_IDS)
def test_strict_mode_reraises_the_cores_own_error(monkeypatch, pkg, fn, kwargs):
    set_strict_core(True)
    mod, core, twin_calls = _arm(monkeypatch, pkg, fn, _decline(fn))
    with pytest.raises((ReadError, WriteError), match="declined"):
        _call(mod, fn, **kwargs)
    assert core.calls
    assert twin_calls == []


# Formats whose Python twin cannot answer the same question for these arguments,
# so falling back would quietly return step 0 (or a single region) under the
# name of what was asked for: a wrong answer, not a slower one.
_NO_TWIN_FOR = [
    ("cgns", {"time_step": 1}),
    ("ensight", {"time_step": 1}),
    ("gmsh", {"time_step": 1}),
    ("med", {"time_step": 1}),
    ("xdmf", {"time_step": 1}),
    ("openfoam", {"region": "fluid"}),
    ("openfoam", {"time_step": 1}),
]


@pytest.mark.parametrize(
    "pkg, kwargs", _NO_TWIN_FOR, ids=[f"{p}:{list(k)[0]}" for p, k in _NO_TWIN_FOR]
)
def test_a_request_the_twin_cannot_answer_reraises_instead_of_falling_back(
    monkeypatch, pkg, kwargs
):
    mod, core, twin_calls = _arm(monkeypatch, pkg, "read", ReadError("no such step"))
    with pytest.raises(ReadError, match="no such step"):
        mod.read("model.dat", **kwargs)
    assert core.calls
    assert twin_calls == []


def test_openfoam_multi_region_case_reraises_instead_of_falling_back(monkeypatch):
    mod, core, twin_calls = _arm(
        monkeypatch, "openfoam", "read", ReadError("found multi-region case")
    )
    with pytest.raises(ReadError, match="multi-region"):
        mod.read("case")
    assert twin_calls == []
