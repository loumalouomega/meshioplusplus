"""Smoke test of the ParaView plugin (tools/paraview-meshioplusplus-plugin.py).

Runs ParaView's own ``pvpython`` in a subprocess: load the plugin, read a file
through "meshio++ reader", check what ParaView got, then write it back through
"meshio++ Writer" and read that with meshio++. Skipped when there is no
``pvpython``, or when its Python cannot import this build (the compiled
``_core`` is tied to one Python version -- which is the plugin's own
requirement, doc/paraview_plugin.md).
"""

import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import pytest

import meshioplusplus

REPO = pathlib.Path(__file__).resolve().parents[2]
PLUGIN = REPO / "tools" / "paraview-meshioplusplus-plugin.py"
PVPYTHON = shutil.which("pvpython")

# The CI job sets this so a missing pvpython fails rather than skips.
_REQUIRED = os.environ.get("MESHIOPLUSPLUS_REQUIRE_PVPYTHON", "") not in ("", "0")

pytestmark = pytest.mark.skipif(
    PVPYTHON is None and not _REQUIRED, reason="pvpython not installed"
)

_SCRIPT = """
import json, sys
import paraview.simple as pvs
from paraview import servermanager

plugin, src, dst = sys.argv[1:4]
pvs.LoadPlugin(plugin, ns=globals())
reader = pvs.meshioreader(FileName=src)
reader.UpdatePipeline()
grid = servermanager.Fetch(reader)
writer = pvs.meshioWriter(Input=reader, FileName=dst)
writer.UpdatePipeline()
print(json.dumps({
    "points": grid.GetNumberOfPoints(),
    "cells": grid.GetNumberOfCells(),
    "point_data": [grid.GetPointData().GetArrayName(i)
                   for i in range(grid.GetPointData().GetNumberOfArrays())],
    "cell_data": [grid.GetCellData().GetArrayName(i)
                  for i in range(grid.GetCellData().GetNumberOfArrays())],
}))
"""


def _env():
    env = dict(os.environ)
    package = pathlib.Path(meshioplusplus.__file__).resolve().parent.parent
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package)] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p]
    )
    return env


def _pvpython_mismatch():
    """Why pvpython cannot run this build, or ``None`` when it can."""
    out = subprocess.run(
        [PVPYTHON, "-c", "import sys; print(sys.version_info[:2])"],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if out.returncode != 0:
        # Not a version question: pvpython itself did not start (a missing
        # shared library, say), and its stderr says why.
        return f"pvpython exited with {out.returncode}: {out.stderr.strip()[-2000:]}"
    ours = str(tuple(sys.version_info[:2]))
    if not out.stdout.strip().endswith(ours):
        return f"pvpython runs Python {out.stdout.strip()!r}, this build {ours}"
    return None


def test_plugin_reads_and_writes_through_paraview(tmp_path):
    assert (
        PVPYTHON is not None
    ), "MESHIOPLUSPLUS_REQUIRE_PVPYTHON is set but there is no pvpython"
    mismatch = _pvpython_mismatch()
    if mismatch:
        if _REQUIRED:
            pytest.fail(mismatch)
        pytest.skip(mismatch)
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3]]), ("triangle", [[1, 2, 4], [0, 1, 4]])],
        point_data={"u": np.arange(5.0)},
        cell_data={"c": [np.array([1.0]), np.array([2.0, 3.0])]},
    )
    src = tmp_path / "in.vtu"
    dst = tmp_path / "out.vtk"
    meshioplusplus.write(src, mesh)
    script = tmp_path / "smoke.py"
    script.write_text(_SCRIPT)
    run = subprocess.run(
        [PVPYTHON, str(script), str(PLUGIN), str(src), str(dst)],
        capture_output=True,
        text=True,
        env=_env(),
        timeout=300,
    )
    assert run.returncode == 0, run.stderr[-3000:]
    seen = json.loads(run.stdout.strip().splitlines()[-1])
    assert seen == {"points": 5, "cells": 3, "point_data": ["u"], "cell_data": ["c"]}

    back = meshioplusplus.read(dst)
    assert len(back.points) == 5
    assert sorted((b.type, len(b.data)) for b in back.cells) == [
        ("tetra", 1),
        ("triangle", 2),
    ]
    np.testing.assert_array_equal(back.point_data["u"], np.arange(5.0))
