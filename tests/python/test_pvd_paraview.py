"""Files meshio++ writes, opened by ParaView itself (the closed roadmap §1.1's "done when").

Vanilla VTK has no collection reader (``vtkPVDReader`` ships with ParaView only),
so ``test_pvtu_vtk.py`` can check the parallel indices but not a ``.pvd``. These
drive a real ``pvpython`` -- a separate interpreter that cannot import meshio++ --
over the files this package writes and assert what ParaView reports: the time
steps, the geometry of each, the ghost cells it recognises, and the field data.

Skipped where ``pvpython`` is not on ``PATH``. A second, independent check through
PyVista's pure-Python ``PVDReader`` runs where only PyVista is installed.
"""

import json
import shutil
import subprocess
import textwrap

import numpy as np
import pytest

import meshioplusplus

PVPYTHON = shutil.which("pvpython")
needs_paraview = pytest.mark.skipif(PVPYTHON is None, reason="pvpython not found")

# Runs inside ParaView. Reads each path and prints one JSON line of what it saw.
_PROBE = textwrap.dedent(
    """
    import json, sys
    from paraview.simple import OpenDataFile
    import paraview.servermanager as sm
    from vtkmodules.util import numpy_support as ns

    out = {}
    for path in sys.argv[1:]:
        reader = OpenDataFile(path)
        reader.UpdatePipelineInformation()
        times = list(getattr(reader, "TimestepValues", []) or [])
        steps = []
        for t in times or [None]:
            reader.UpdatePipeline(t) if t is not None else reader.UpdatePipeline()
            d = sm.Fetch(reader)
            ghost = d.GetCellGhostArray()
            fd = d.GetFieldData()
            steps.append({
                "points": d.GetNumberOfPoints(),
                "cells": d.GetNumberOfCells(),
                "cell_arrays": [d.GetCellData().GetArrayName(i)
                                for i in range(d.GetCellData().GetNumberOfArrays())],
                "field": {fd.GetArrayName(i):
                          [float(x) for x in ns.vtk_to_numpy(fd.GetArray(i)).ravel()]
                          for i in range(fd.GetNumberOfArrays())},
                "ghost_flagged": (int((ns.vtk_to_numpy(ghost) != 0).sum())
                                  if ghost is not None else None),
            })
        out[path.rsplit("/", 1)[-1]] = {"times": times, "steps": steps}
    print("JSON:" + json.dumps(out))
    """
)


def _paraview(tmp_path, *paths):
    script = tmp_path / "probe.py"
    script.write_text(_PROBE)
    proc = subprocess.run(
        [PVPYTHON, str(script), *map(str, paths)],
        capture_output=True,
        text=True,
        timeout=240,
        cwd=tmp_path,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.startswith("JSON:")]
    assert lines, f"pvpython produced no result:\n{proc.stdout}\n{proc.stderr}"
    return json.loads(lines[-1][len("JSON:") :])


def _grid(n=6, shift=0.0):
    xs, ys = np.meshgrid(np.arange(n + 1.0), np.arange(n + 1.0))
    points = np.c_[xs.ravel(), ys.ravel(), np.full((n + 1) ** 2, shift)]
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b, c = a + 1, a + n + 1
            tris += [[a, b, c + 1], [a, c + 1, c]]
    mesh = meshioplusplus.Mesh(points, [("triangle", np.array(tris))])
    mesh.point_data["u"] = np.full(len(points), shift)
    mesh.cell_data["c"] = [np.arange(len(tris), dtype=np.float64)]
    return mesh


def _write_pvd(path, entries):
    rows = "".join(f'<DataSet timestep="{t}" part="0" file="{f}"/>' for t, f in entries)
    path.write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
        f'byte_order="LittleEndian"><Collection>{rows}</Collection></VTKFile>'
    )


@needs_paraview
def test_paraview_opens_a_pvd_of_vtu_with_the_written_time_steps(tmp_path):
    times = [0.0, 0.5, 2.0]
    meshioplusplus.write_sequence(
        tmp_path / "run.pvd", ((t, _grid(shift=i)) for i, t in enumerate(times))
    )
    seen = _paraview(tmp_path, tmp_path / "run.pvd")["run.pvd"]
    assert seen["times"] == times
    for step in seen["steps"]:
        assert (step["points"], step["cells"]) == (49, 72)


@needs_paraview
def test_paraview_opens_a_pvd_of_pvtu_with_halos_recognised_as_ghost_cells(tmp_path):
    """The roadmap's "done when": a `.pvd` of `.pvtu` plays in ParaView."""
    if not hasattr(meshioplusplus._core, "partition"):
        pytest.skip("ghost layers need the C++ core")
    entries, owned, halo = [], 0, 0
    for step in range(3):
        mesh = _grid(shift=step)
        mesh.field_data["TimeValue"] = np.array([0.5 * step])
        pieces = meshioplusplus.partition(mesh, 3, ghost_layers=1)
        meshioplusplus.pvtu.write_pieces(tmp_path / f"s{step}.pvtu", pieces)
        entries.append((0.5 * step, f"s{step}.pvtu"))
        layers = [p.cell_data["partition:ghost"][0] for p in pieces]
        owned = sum(int(np.sum(x == 0)) for x in layers)
        halo = sum(int(np.sum(x > 0)) for x in layers)
    _write_pvd(tmp_path / "run.pvd", entries)

    seen = _paraview(tmp_path, tmp_path / "run.pvd")["run.pvd"]
    assert seen["times"] == [0.0, 0.5, 1.0]
    assert owned == 72 and halo > 0
    for step, k in zip(seen["steps"], range(3)):
        assert step["cells"] == owned + halo
        # ParaView recognises our vtkGhostType as its ghost array, flag for flag
        assert step["ghost_flagged"] == halo
        assert {"c", "vtkGhostType"} <= set(step["cell_arrays"])
        # ... and reads the dataset's field data, one value per step
        assert step["field"]["TimeValue"] == [0.5 * k]


@needs_paraview
def test_paraview_opens_a_pvtu_and_reads_its_pieces(tmp_path):
    mesh = _grid()
    mesh.cell_data["partition:part"] = meshioplusplus.partition_labels(mesh, 3)
    meshioplusplus.pvtu.write(tmp_path / "a.pvtu", mesh)
    seen = _paraview(tmp_path, tmp_path / "a.pvtu")["a.pvtu"]
    (step,) = seen["steps"]
    assert step["cells"] == 72
    assert step["points"] > 49  # interface points are duplicated per piece


@needs_paraview
@pytest.mark.parametrize("ext", [".vtu", ".vtp"])
def test_a_time_value_meshio_writes_is_the_datasets_time_in_paraview(tmp_path, ext):
    """VTK's "time in field data" convention, read back by ParaView."""
    mesh = _grid()
    mesh.field_data["TimeValue"] = np.array([0.25])
    meshioplusplus.write(tmp_path / f"timed{ext}", mesh)
    seen = _paraview(tmp_path, tmp_path / f"timed{ext}")[f"timed{ext}"]
    assert seen["times"] == [0.25]
    assert seen["steps"][0]["field"]["TimeValue"] == [0.25]


# --- the second, independent implementation ------------------------------------


def test_pyvistas_independent_pvd_reader_agrees(tmp_path):
    """PyVista parses the same four attributes in pure Python, with the same
    defaults and the same relative-path rule, and dispatches `.pvtu` entries."""
    pv = pytest.importorskip("pyvista")
    times = [0.0, 0.5, 2.0]
    meshioplusplus.write_sequence(
        tmp_path / "run.pvd", ((t, _grid(shift=i)) for i, t in enumerate(times))
    )
    reader = pv.PVDReader(str(tmp_path / "run.pvd"))
    assert reader.time_values == times
    assert [d.path for d in reader.datasets] == [
        f"run/run_{i:04d}.vtu" for i in range(3)
    ]
    reader.set_active_time_value(2.0)
    block = reader.read()[0]
    assert (block.n_points, block.n_cells) == (49, 72)
    np.testing.assert_allclose(block.points[:, 2], 2.0)


def test_pyvista_follows_a_pvd_entry_that_is_a_pvtu(tmp_path):
    pv = pytest.importorskip("pyvista")
    mesh = _grid()
    mesh.cell_data["partition:part"] = meshioplusplus.partition_labels(mesh, 3)
    meshioplusplus.pvtu.write(tmp_path / "s0.pvtu", mesh)
    _write_pvd(tmp_path / "run.pvd", [(0.0, "s0.pvtu")])
    reader = pv.PVDReader(str(tmp_path / "run.pvd"))
    assert reader.time_values == [0.0]
    assert reader.read()[0].n_cells == 72
