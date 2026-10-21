"""fit_cad golden digests through Python.

The STEP bytes, the returned mesh and the report of each fixture are pinned in
``fit_cad_goldens.json``, so a change meant to leave the default output alone
(every gated skin-to-CAD stage) provably does. The native backend's rows are
keyed ``"native"``; the OpenCASCADE backend's by the kernel version the core
links (``_core.__occt_version__``), since its bytes depend on it.

The digests rest on the kernel's build and on the C library, so the check is
opt-in:
``MESHIOPLUSPLUS_FIT_CAD_GOLDENS=1``, which the CI OpenCASCADE legs set.
``MESHIOPLUSPLUS_RECORD_GOLDENS=1`` prints the table to paste instead.
"""

import hashlib
import json
import os
import struct
from pathlib import Path

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fit_cad

from .analytic_surfaces import box_skin, jitter

GOLDENS = Path(__file__).with_name("fit_cad_goldens.json")
CHECK = os.environ.get("MESHIOPLUSPLUS_FIT_CAD_GOLDENS", "0") not in ("", "0")
RECORD = os.environ.get("MESHIOPLUSPLUS_RECORD_GOLDENS", "0") not in ("", "0")

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - no compiled core
    _core = None

pytestmark = [
    pytest.mark.skipif(_core is None, reason="needs the C++ core"),
    pytest.mark.skipif(
        not (CHECK or RECORD),
        reason="set MESHIOPLUSPLUS_FIT_CAD_GOLDENS=1 to check fit_cad goldens",
    ),
]


def _fitted(skin, curvature_tolerance=0.0, tolerance=1e-3, min_patch_cells=1):
    seg = mio.segment(
        skin,
        curvature_tolerance=curvature_tolerance,
        min_patch_cells=min_patch_cells,
    )
    return mio.fit_primitives(seg, tolerance=tolerance)


def _saddle_box(n):
    box = box_skin(n)
    p = np.array(box.points, dtype=float)
    p[:, 2] *= 1.0 + 0.8 * (p[:, 0] - 0.5) * (p[:, 1] - 0.5)
    box.points = p
    return box


# Trig-free skins only: NumPy dispatches float64 sin/cos to SIMD kernels by
# CPU, so a fixture built with them can differ in its last bits between
# machines. The C++ table (test_fit_cad_goldens.cpp) covers the curved ones.
FIXTURES = {
    "box": lambda: _fitted(box_skin(3)),
    "saddle_box": lambda: _fitted(_saddle_box(24)),
    "noisy_box": lambda: _fitted(jitter(box_skin(6), 1e-4), tolerance=2e-3),
}


def _mesh_digest(mesh):
    h = hashlib.sha256()

    def array(a):
        a = np.ascontiguousarray(a)
        h.update(str((a.dtype.str, a.shape)).encode())
        h.update(a.tobytes())

    array(mesh.points)
    for block in mesh.cells:
        h.update(block.type.encode())
        array(block.data)
    for name in sorted(mesh.cell_data):
        h.update(name.encode())
        for a in mesh.cell_data[name]:
            array(a)
    for name in sorted(mesh.field_data):
        h.update(name.encode())
        array(mesh.field_data[name])
    return h.hexdigest()


# Report keys added after the digests were pinned, each checked by its own
# tests: the digests guard what the report already said.
_LATER_KEYS = {"closed"}


def _report_digest(report):
    h = hashlib.sha256()
    for key in sorted(report):
        if key in _LATER_KEYS:
            continue
        v = report[key]
        h.update(key.encode())
        if isinstance(v, float):
            h.update(struct.pack("<d", v))
        else:
            h.update(repr(v).encode())
    return h.hexdigest()


def _run(skin, backend, tmp_path):
    path = tmp_path / "fit_cad_golden.step"
    out, report = mio.fit_cad(skin, path, backend=backend, return_report=True)
    return {
        "step": hashlib.sha256(path.read_bytes()).hexdigest(),
        "mesh": _mesh_digest(out),
        "report": _report_digest(report),
    }


def _keys():
    keys = ["native"]
    if _fit_cad.has_occt():
        keys.append(_core.__occt_version__)
    return keys


@pytest.mark.parametrize("key", ["native", "occt"])
def test_fit_cad_goldens(key, tmp_path, capsys):
    if key == "occt":
        if not _fit_cad.has_occt():
            pytest.skip("core built without OpenCASCADE")
        key = _core.__occt_version__
    backend = "native" if key == "native" else "occt"
    table = json.loads(GOLDENS.read_text()) if GOLDENS.exists() else {}
    got = {name: _run(make(), backend, tmp_path) for name, make in FIXTURES.items()}
    if RECORD:
        with capsys.disabled():
            print(json.dumps({key: got}, indent=2, sort_keys=True))
        return
    if key not in table:
        pytest.skip(
            f"no golden rows for {key!r}; record them with MESHIOPLUSPLUS_RECORD_GOLDENS=1"
        )
    assert got == table[key]
