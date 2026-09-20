"""The C++ fast path's fall-back decision (roadmap "reader fallback swallows the reason").

Every format shim used to wrap ``_core.<fmt>_read`` / ``_write`` in a bare
``except Exception: pass``. ``meshioplusplus._fallback.core_declined`` now decides:
a recognised decline falls back quietly, an unexpected failure falls back
*loudly*, and an exception that is never a decline propagates. See that module's
docstring for the table these tests pin.
"""

import importlib
import logging
from pathlib import Path

import pytest

import meshioplusplus
from meshioplusplus import ReadError, WriteError, _core
from meshioplusplus._fallback import core_declined, set_strict_core, strict_core

MSH_PATH = Path(__file__).resolve().parent / "meshes" / "msh" / "insulated-2.2.msh"


@pytest.fixture(autouse=True)
def _reset_strict(monkeypatch):
    monkeypatch.delenv("MESHIOPLUSPLUS_STRICT_CORE", raising=False)
    set_strict_core(None)
    yield
    set_strict_core(None)


def _records(caplog):
    return [r for r in caplog.records if r.name == "meshioplusplus"]


# --- the decision itself -------------------------------------------------


@pytest.mark.parametrize("exc", [ReadError("x"), WriteError("x")])
def test_recognised_decline_falls_back_quietly(exc, caplog):
    with caplog.at_level(logging.DEBUG, logger="meshioplusplus"):
        assert core_declined(exc, "fmt", "read", "a.file") is True
    recs = _records(caplog)
    assert [r.levelno for r in recs] == [logging.DEBUG]
    assert "declined" in recs[0].getMessage()


@pytest.mark.parametrize(
    "exc",
    [
        ValueError("x"),
        IndexError("x"),
        KeyError("x"),
        RuntimeError("x"),
        AttributeError("x"),
    ],
)
def test_unexpected_failure_falls_back_but_warns(exc, caplog):
    with caplog.at_level(logging.DEBUG, logger="meshioplusplus"):
        assert core_declined(exc, "fmt", "write", "a.file") is True
    recs = _records(caplog)
    assert [r.levelno for r in recs] == [logging.WARNING]
    msg = recs[0].getMessage()
    assert "failed unexpectedly" in msg and type(exc).__name__ in msg
    assert "a.file" in msg and "MESHIOPLUSPLUS_STRICT_CORE" in msg


@pytest.mark.parametrize("exc", [TypeError("x"), MemoryError(), RecursionError()])
def test_never_a_decline_propagates_and_is_not_logged(exc, caplog):
    # RecursionError subclasses RuntimeError: it must not land in the warn tier.
    with caplog.at_level(logging.DEBUG, logger="meshioplusplus"):
        assert core_declined(exc, "fmt", "read", "a.file") is False
    assert _records(caplog) == []


@pytest.mark.parametrize("exc", [ReadError("x"), ValueError("x")])
def test_strict_mode_reraises_everything(exc, caplog):
    set_strict_core(True)
    with caplog.at_level(logging.DEBUG, logger="meshioplusplus"):
        assert core_declined(exc, "fmt", "read", "a.file") is False
    assert [r.levelno for r in _records(caplog)] == [logging.WARNING]


@pytest.mark.parametrize(
    "value, expected",
    [
        ("1", True),
        ("on", True),
        ("TRUE", True),
        (" yes ", True),
        ("0", False),
        ("off", False),
        ("", False),
    ],
)
def test_env_var_parsing(monkeypatch, value, expected):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", value)
    set_strict_core(None)
    assert strict_core() is expected


def test_env_var_unset_is_off():
    assert strict_core() is False


def test_env_var_is_read_once(monkeypatch):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    set_strict_core(None)
    assert strict_core() is True
    monkeypatch.delenv("MESHIOPLUSPLUS_STRICT_CORE")
    assert strict_core() is True  # cached, like MESHIOPLUSPLUS_PROVENANCE


# --- through a real shim -------------------------------------------------

SENTINEL = object()


@pytest.fixture
def ansys(monkeypatch):
    mod = importlib.import_module("meshioplusplus.ansys")
    monkeypatch.setattr(mod, "_py_read", lambda filename: SENTINEL)
    monkeypatch.setattr(mod, "_py_write", lambda *a, **k: SENTINEL)
    return mod


def _raising(exc):
    def fn(*args, **kwargs):
        raise exc

    return fn


def test_shim_read_falls_back_on_read_error(ansys, monkeypatch):
    monkeypatch.setattr(_core, "ansys_read", _raising(ReadError("declined")))
    assert ansys.read("x.msh") is SENTINEL


def test_shim_read_falls_back_on_unexpected_failure(ansys, monkeypatch, caplog):
    monkeypatch.setattr(_core, "ansys_read", _raising(ValueError("stoll")))
    with caplog.at_level(logging.WARNING, logger="meshioplusplus"):
        assert ansys.read("x.msh") is SENTINEL
    assert any("failed unexpectedly" in r.getMessage() for r in _records(caplog))


@pytest.mark.parametrize("exc", [TypeError("bad arg"), MemoryError()])
def test_shim_read_propagates_never_a_decline(ansys, monkeypatch, exc):
    monkeypatch.setattr(_core, "ansys_read", _raising(exc))
    with pytest.raises(type(exc)):
        ansys.read("x.msh")


def test_shim_write_falls_back_on_write_error(ansys, monkeypatch):
    monkeypatch.setattr(_core, "ansys_write", _raising(WriteError("declined")))
    assert ansys.write("x.msh", object()) is SENTINEL


def test_shim_strict_mode_reraises(ansys, monkeypatch):
    monkeypatch.setattr(_core, "ansys_read", _raising(ReadError("declined")))
    set_strict_core(True)
    with pytest.raises(ReadError, match="declined"):
        ansys.read("x.msh")


def test_shim_buffer_skips_the_core(ansys, monkeypatch):
    import io

    monkeypatch.setattr(
        _core, "ansys_read", _raising(AssertionError("core must not run"))
    )
    assert ansys.read(io.StringIO("")) is SENTINEL


def test_time_step_reraise_is_preserved(monkeypatch):
    # gmsh/xdmf/med/... re-raise instead of falling back when a step was
    # asked for: the Python twin would hand back step 0 under step N's name.
    mod = importlib.import_module("meshioplusplus.gmsh")
    monkeypatch.setattr(mod, "_py_read", lambda filename: SENTINEL)
    monkeypatch.setattr(_core, "gmsh_read", _raising(ReadError("no such step")))
    assert mod.read("x.msh") is SENTINEL
    with pytest.raises(ReadError, match="no such step"):
        mod.read("x.msh", time_step=1)


def test_openfoam_region_reraise_is_preserved(monkeypatch):
    mod = importlib.import_module("meshioplusplus.openfoam")
    monkeypatch.setattr(mod, "_py_read", lambda filename: SENTINEL)
    monkeypatch.setattr(
        _core, "openfoam_read", _raising(ReadError("multi-region case"))
    )
    with pytest.raises(ReadError):
        mod.read("case")  # message names a multi-region case
    monkeypatch.setattr(_core, "openfoam_read", _raising(ReadError("declined")))
    assert mod.read("case") is SENTINEL
    with pytest.raises(ReadError):
        mod.read("case", region="fluid")


# --- the regression the tiering exists to prevent -----------------------


def test_ambiguous_msh_still_resolves_without_warnings(caplog):
    # `.msh` tries ansys, gmsh, freefem in turn; the non-gmsh candidates decline
    # with ReadError, which must stay on the quiet (DEBUG) tier.
    with caplog.at_level(logging.DEBUG, logger="meshioplusplus"):
        mesh = meshioplusplus.read(MSH_PATH)
    assert len(mesh.points) > 0
    assert [r for r in _records(caplog) if r.levelno >= logging.WARNING] == []


def test_ambiguous_msh_still_resolves_in_strict_mode():
    # Strict mode re-raises out of each shim, but `read()`'s candidate loop
    # catches ReadError and moves on, so a real gmsh file still reads.
    set_strict_core(True)
    assert len(meshioplusplus.read(MSH_PATH).points) > 0
