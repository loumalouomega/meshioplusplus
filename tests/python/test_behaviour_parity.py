"""Pin behaviour parity decisions from roadmap 3.4.12."""

from __future__ import annotations

import ast
import importlib.util
import io
import pathlib
import sys
import types

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
PKG = REPO / "src" / "python" / "meshioplusplus"


def _load_shims():
    """Load ``_shims.py`` in isolation without requiring ``numpy`` or ``_core``."""
    pkg = types.ModuleType("meshioplusplus")
    pkg.__path__ = [str(PKG)]
    sys.modules.setdefault("meshioplusplus", pkg)
    # Provide minimal _core stub only if not already loaded.
    if "meshioplusplus._core" not in sys.modules:
        sys.modules["meshioplusplus._core"] = types.ModuleType("meshioplusplus._core")
    spec = importlib.util.spec_from_file_location(
        "meshioplusplus._shims", PKG / "_shims.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["meshioplusplus._shims"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_time_values_shim_guards_file_like_buffers():
    shims = _load_shims()
    fake_pkg = types.ModuleType("meshioplusplus._fake_tv_pkg")
    core_calls: list[str] = []
    py_calls: list[object] = []

    class _Core:
        @staticmethod
        def fake_time_values(path: str) -> list[float]:
            core_calls.append(path)
            return [1.0, 2.0]

    fake_pkg._core = _Core  # type: ignore[attr-defined]
    fake_pkg._py_time_values = lambda src: (py_calls.append(src), [3.0])[1]  # type: ignore[attr-defined]
    sys.modules[fake_pkg.__name__] = fake_pkg
    try:
        fn = shims.time_values_shim(fake_pkg.__name__, "fake", "doc")
        buf = io.BytesIO(b"data")
        assert fn(buf) == [3.0]
        assert core_calls == []
        assert py_calls == [buf]

        assert fn("run.dat") == [1.0, 2.0]
        assert core_calls == ["run.dat"]
    finally:
        sys.modules.pop(fake_pkg.__name__, None)


@pytest.mark.parametrize(
    "fmt",
    [
        "nastran_op2",
        "radioss_th",
        "marc",
        "ansys_rst",
        "lsdyna_d3plot",
        "lsdyna_binout",
        "femap",
        "abaqus_fil",
        "selafin",
        "ugrid_nc",
    ],
)
def test_every_time_values_format_uses_guarded_shim(fmt: str):
    src = (PKG / fmt / "__init__.py").read_text(encoding="utf-8")
    uses_tv_shim = "time_values = time_values_shim(" in src
    uses_series_shims = "series_shims(" in src and "time_values" in src
    assert uses_tv_shim or uses_series_shims


def test_medmulti_writes_med_num():
    src = (PKG / "med" / "_medmulti.py").read_text(encoding="utf-8")
    tree = ast.parse(src)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "_write_mesh_step"
    ]
    assert len(calls) == 1
    kw = {k.arg: k.value for k in calls[0].keywords}
    assert isinstance(kw["write_num"], ast.Constant) and kw["write_num"].value is True


def test_decimate_shares_remap_sets_without_filtering_out_of_range():
    src = (PKG / "_decimate.py").read_text(encoding="utf-8")
    assert "def _remap_sets(" not in src
    assert "from ._set_remap import remap_sets" in src
    assert "remap_sets(mesh, out, point_map, cell_maps, unique=True)" in src


def test_remap_sets_raises_on_out_of_range_cell_index():
    np = pytest.importorskip("numpy")
    spec = importlib.util.spec_from_file_location(
        "meshioplusplus._set_remap", PKG / "_set_remap.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    class _Dummy:
        point_sets = None
        cell_sets = {"bad": [np.array([0, 99], dtype=np.int64)]}

    out = types.SimpleNamespace()
    with pytest.raises(IndexError):
        mod.remap_sets(
            _Dummy(),
            out,
            np.array([0, 1], dtype=np.int64),
            [np.array([0, 1], dtype=np.int64)],
            unique=True,
        )
