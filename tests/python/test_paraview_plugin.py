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

_PRELUDE = """
import json, sys
import paraview.simple as pvs
from paraview import servermanager

pvs.LoadPlugin(sys.argv[1], ns=globals())
args = sys.argv[2:]
"""

_SCRIPT = """
src, dst = args
reader = meshioreader(FileName=src)
reader.UpdatePipeline()
grid = servermanager.Fetch(reader)
writer = meshioWriter(Input=reader, FileName=dst)
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


def _run_pvpython(tmp_path, script, *args):
    """Run ``script`` in pvpython with the plugin loaded; its last line, as JSON.

    The script sees ``pvs`` (paraview.simple), ``servermanager``, the plugin's
    proxies in its globals, and ``args`` (the strings passed here).
    """
    assert (
        PVPYTHON is not None
    ), "MESHIOPLUSPLUS_REQUIRE_PVPYTHON is set but there is no pvpython"
    mismatch = _pvpython_mismatch()
    if mismatch:
        if _REQUIRED:
            pytest.fail(mismatch)
        pytest.skip(mismatch)
    path = tmp_path / "probe.py"
    path.write_text(_PRELUDE + script)
    run = subprocess.run(
        [PVPYTHON, str(path), str(PLUGIN), *map(str, args)],
        capture_output=True,
        text=True,
        env=_env(),
        timeout=300,
    )
    assert run.returncode == 0, run.stderr[-3000:]
    return json.loads(run.stdout.strip().splitlines()[-1])


def test_plugin_reads_and_writes_through_paraview(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3]]), ("triangle", [[1, 2, 4], [0, 1, 4]])],
        point_data={"u": np.arange(5.0)},
        cell_data={"c": [np.array([1.0]), np.array([2.0, 3.0])]},
    )
    src = tmp_path / "in.vtu"
    dst = tmp_path / "out.vtk"
    meshioplusplus.write(src, mesh)
    seen = _run_pvpython(tmp_path, _SCRIPT, src, dst)
    assert seen == {"points": 5, "cells": 3, "point_data": ["u"], "cell_data": ["c"]}

    back = meshioplusplus.read(dst)
    assert len(back.points) == 5
    # The writer keeps the cell order: one tetra block, then the triangles.
    assert [(b.type, len(b.data)) for b in back.cells] == [
        ("tetra", 1),
        ("triangle", 2),
    ]
    np.testing.assert_array_equal(back.point_data["u"], np.arange(5.0))


_VOLUMES = """
reader = meshioreader(FileName=args[0])
sizes = pvs.CellSize(Input=reader)
sizes.UpdatePipeline()
grid = servermanager.Fetch(sizes)
v = grid.GetCellData().GetArray("Volume")
values = [v.GetValue(i) for i in range(v.GetNumberOfTuples())]
writer = meshioWriter(Input=reader, FileName=args[1])
writer.UpdatePipeline()
solids = [x for x in values if x != 0.0]
print(json.dumps({"cells": grid.GetNumberOfCells(), "solids": len(solids),
                  "min": min(solids) if solids else None}))
"""


@pytest.mark.parametrize("case", ["linear", "toroid"])
def test_wedges_are_not_inside_out(tmp_path, case):
    """roadmap 1.1: every wedge ParaView shows has a positive volume, and the
    file the writer writes back has no inverted cell."""
    if case == "linear":
        src = tmp_path / "wedges.msh"
        meshioplusplus.write(
            src,
            meshioplusplus.Mesh(
                np.array(
                    [[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]
                    + [[1, 1, 0], [1, 1, 1]]
                ),
                [("wedge", [[0, 1, 2, 3, 4, 5], [1, 6, 2, 4, 7, 5]])],
            ),
        )
    else:
        src = REPO / "tests" / "python" / "meshes" / "mfem" / "toroid-wedge.mesh"
    dst = tmp_path / "out.vtu"
    seen = _run_pvpython(tmp_path, _VOLUMES, src, dst)
    mesh = meshioplusplus.read(src)
    assert meshioplusplus.compute_quality(mesh)["num_inverted"] == 0
    assert seen["cells"] == sum(len(b) for b in mesh.cells)
    assert seen["solids"] > 0
    assert seen["min"] > 0
    assert meshioplusplus.compute_quality(meshioplusplus.read(dst))["num_inverted"] == 0


_COUNTS = """
reader = meshioreader(FileName=args[0])
reader.UpdatePipeline()
grid = servermanager.Fetch(reader)
print(json.dumps({"points": grid.GetNumberOfPoints(), "cells": grid.GetNumberOfCells()}))
"""


@pytest.mark.parametrize("case", ["polygon5", "line_1d"])
def test_polygon_n_and_1d_meshes_are_read(tmp_path, case):
    """roadmap 1.2: a polygon5 block and a 1-D mesh keep their cells."""
    if case == "polygon5":
        mesh = meshioplusplus.Mesh(
            np.array(
                [[0.0, 0, 0], [1, 0, 0], [1.5, 1, 0], [0.5, 1.6, 0], [-0.5, 1, 0]]
            ),
            [("polygon5", [[0, 1, 2, 3, 4]]), ("triangle", [[0, 1, 2]])],
        )
        src = tmp_path / "in.vtu"
    else:
        mesh = meshioplusplus.Mesh(
            np.array([[0.0], [1.0], [3.0]]), [("line", [[0, 1], [1, 2]])]
        )
        src = tmp_path / "in.msh"
    meshioplusplus.write(src, mesh)
    expected = meshioplusplus.read(src)
    seen = _run_pvpython(tmp_path, _COUNTS, src)
    assert seen == {
        "points": len(expected.points),
        "cells": sum(len(b) for b in expected.cells),
    }
