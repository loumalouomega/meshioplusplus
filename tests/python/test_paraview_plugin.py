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


_NAMES = """
src, state = args
pxm = servermanager.ProxyManager()
names = {}
for group, name in (
    ("sources", "MeshioPlusPlusReader"),
    ("sources", "meshio++ reader"),
    ("writers", "MeshioPlusPlusWriter"),
    ("writers", "meshio++ Writer"),
):
    names[name] = pxm.GetPrototypeProxy(group, name) is not None
# A state file saved with the old reader name loads into the alias, which
# paraview.simple names after its label, "meshio++ reader (old name)".
reader = meshioreaderoldname(FileName=src)
pvs.SaveState(state)
pvs.Delete(reader)
pvs.LoadState(state)
loaded = [p for p in pvs.GetSources().values() if p.GetXMLName() == "meshio++ reader"]
loaded[0].UpdatePipeline()
names["state_cells"] = servermanager.Fetch(loaded[0]).GetNumberOfCells()
# Opening a .vtu in a script picks ParaView's own reader, not this one.
names["vtu_reader"] = pvs.OpenDataFile(src).GetXMLName()
print(json.dumps(names))
"""


def test_stable_names_and_old_aliases(tmp_path):
    """roadmap 6.10.2: stable proxy names; the old ones still load a state."""
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", [[0, 1, 2, 3]])],
    )
    src = tmp_path / "in.vtu"
    meshioplusplus.write(src, mesh)
    seen = _run_pvpython(tmp_path, _NAMES, src, tmp_path / "state.pvsm")
    assert seen.pop("vtu_reader") != "MeshioPlusPlusReader"
    assert seen == {
        "MeshioPlusPlusReader": True,
        "meshio++ reader": True,
        "MeshioPlusPlusWriter": True,
        "meshio++ Writer": True,
        "state_cells": 1,
    }


def _capability(flag):
    from meshioplusplus import _core

    return bool(getattr(_core, flag, False))


@pytest.mark.parametrize(
    "ext, flag",
    [(".xdmf", "__has_hdf5__"), (".exo", "__has_netcdf__"), (".cgns", "__has_hdf5__")],
)
def test_hdf5_and_netcdf_formats_through_the_plugin(tmp_path, ext, flag):
    """roadmap 6.10.1: XDMF, Exodus and CGNS go through the plugin in CI."""
    if not _capability(flag):
        if _REQUIRED:
            pytest.fail(f"this build has no {flag.strip('_')}")
        pytest.skip(f"this build has no {flag.strip('_')}")
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]),
        [("tetra", [[0, 1, 2, 3], [1, 2, 3, 4]])],
        point_data={"u": np.arange(5.0)},
    )
    src = tmp_path / f"in{ext}"
    meshioplusplus.write(src, mesh)
    expected = meshioplusplus.read(src)
    seen = _run_pvpython(tmp_path, _COUNTS, src)
    assert seen == {
        "points": len(expected.points),
        "cells": sum(len(b) for b in expected.cells),
    }


def test_polyhedra_reach_paraview(tmp_path):
    """roadmap 6.2: a VTU of polyhedra has the cell count meshioplusplus.read
    gives, and the writer writes the polyhedra back."""
    src = REPO / "tests" / "python" / "meshes" / "vtu" / "two_pieces_polyhedron.vtu"
    if src.read_bytes().startswith(b"version https://git-lfs"):
        pytest.skip("LFS content not fetched")
    dst = tmp_path / "out.vtu"
    seen = _run_pvpython(tmp_path, _SCRIPT, src, dst)
    mesh = meshioplusplus.read(src)
    assert seen["cells"] == sum(len(b) for b in mesh.cells)
    back = meshioplusplus.read(dst)
    assert sorted((b.type, len(b)) for b in back.cells) == sorted(
        (b.type, len(b)) for b in mesh.cells
    )


_TIME = """
src, dst = args
reader = meshioreader(FileName=src)
reader.UpdatePipelineInformation()
times = list(reader.TimestepValues)
reader.UpdatePipeline(2.0)
grid = servermanager.Fetch(reader)
u = grid.GetPointData().GetArray("u")
writer = meshioWriter(Input=reader, FileName=dst, WriteTimeSteps=1)
writer.UpdatePipeline()
print(json.dumps({"times": times, "u": [u.GetValue(i) for i in range(u.GetNumberOfTuples())]}))
"""


def test_time_steps_in_and_out(tmp_path):
    """roadmap 6.4: a three-step frd shows TimestepValues [0, 0.5, 2], time
    2.0 is the last step, and the writer writes every step back."""
    points = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    steps = [
        (
            t,
            meshioplusplus.Mesh(
                points,
                [("tetra", [[0, 1, 2, 3]])],
                point_data={"u": np.full(4, 10.0 * t)},
            ),
        )
        for t in (0.0, 0.5, 2.0)
    ]
    src = tmp_path / "in.frd"
    dst = tmp_path / "out.frd"
    meshioplusplus.write_sequence(src, steps)
    seen = _run_pvpython(tmp_path, _TIME, src, dst)
    assert seen["times"] == [0.0, 0.5, 2.0]
    assert seen["u"] == [20.0] * 4
    assert meshioplusplus.read_metadata(dst)["time_values"] == [0.0, 0.5, 2.0]


_ARRAYS = """
src = args[0]
reader = meshioreader(FileName=src)
reader.UpdatePipelineInformation()
offered = sorted(reader.PointArrays.Available)
reader.PointArrays = ["v"]
reader.UpdatePipeline()
grid = servermanager.Fetch(reader)
names = [grid.GetPointData().GetArrayName(i)
         for i in range(grid.GetPointData().GetNumberOfArrays())]
print(json.dumps({"offered": offered, "read": names,
                  "addons": list(reader.GetProperty("AddOns"))}))
"""


def test_an_unticked_array_is_absent(tmp_path):
    """roadmap 6.5: the array checkboxes come from read_metadata, and an
    unticked array is absent from the fetched dataset."""
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", [[0, 1, 2, 3]])],
        point_data={"u": np.arange(4.0), "v": np.ones(4)},
    )
    src = tmp_path / "in.vtu"
    meshioplusplus.write(src, mesh)
    seen = _run_pvpython(tmp_path, _ARRAYS, src)
    assert seen["offered"] == ["u", "v"]
    assert seen["read"] == ["v"]
    assert sorted(a.split(":")[0] for a in seen["addons"]) == ["mmg", "netgen", "occt"]
