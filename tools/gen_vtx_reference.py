"""Record ParaView's reading of the VTX fixtures as the ``vtx`` reader's oracle.

Run with ParaView's ``pvpython`` (6.1 was used)::

    pvpython tools/gen_vtx_reference.py tests/python/meshes/vtx

For every ``<name>.bp`` it writes ``<name>.reference.npz`` holding, per step
``k``, what ``ADIOS2VTXReader`` produced after ``MergeBlocks``: ``k<k>/points``,
``k<k>/connectivity``, ``k<k>/offsets``, ``k<k>/types`` and one
``k<k>/point/<name>`` or ``k<k>/cell/<name>`` per array, plus ``times`` and
``steps`` (the step indices recorded).

Each step is read in a child process: ParaView 6.1 aborts (an ADIOS2 assertion)
on every step after the first of a file written with ``VTXMeshPolicy.reuse``,
where the mesh variables exist in step 0 only. Such steps are left out of the
reference; the tests check them against the raw ADIOS2 arrays instead.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def _one(path, k, out):
    from paraview import servermanager
    from paraview import simple as pv
    from vtkmodules.util.numpy_support import vtk_to_numpy

    reader = pv.ADIOS2VTXReader(FileName=str(path))
    reader.UpdatePipelineInformation()
    times = reader.TimestepValues
    times = [float(t) for t in times] if hasattr(times, "__len__") else [float(times)]
    merged = pv.MergeBlocks(Input=reader)
    merged.UpdatePipeline(float(times[k]))
    grid = servermanager.Fetch(merged)
    res = {"times": np.asarray(times, dtype=np.float64)}
    res["points"] = vtk_to_numpy(grid.GetPoints().GetData()).astype(np.float64)
    cells = grid.GetCells()
    res["connectivity"] = vtk_to_numpy(cells.GetConnectivityArray()).astype(np.int64)
    res["offsets"] = vtk_to_numpy(cells.GetOffsetsArray()).astype(np.int64)
    res["types"] = vtk_to_numpy(grid.GetCellTypes()).astype(np.int64)
    for kind, data in (("point", grid.GetPointData()), ("cell", grid.GetCellData())):
        for i in range(data.GetNumberOfArrays()):
            arr = data.GetArray(i)
            if arr is not None:
                res[f"{kind}/{arr.GetName()}"] = vtk_to_numpy(arr)
    np.savez(out, **res)


def _dump(path):
    merged = {}
    steps = []
    with tempfile.TemporaryDirectory() as tmp:
        k = 0
        while True:
            out = Path(tmp) / f"k{k}.npz"
            proc = subprocess.run(
                [
                    shutil.which("pvpython") or sys.executable,
                    __file__,
                    "--one",
                    str(path),
                    str(k),
                    str(out),
                ],
                capture_output=True,
            )
            if proc.returncode == 0 and out.exists():
                z = np.load(out)
                merged["times"] = z["times"]
                for key in z.files:
                    if key != "times":
                        merged[f"k{k}/{key}"] = z[key]
                steps.append(k)
            if "times" in merged and k + 1 >= len(merged["times"]):
                break
            if "times" not in merged:
                raise SystemExit(f"ParaView cannot read step 0 of {path}")
            k += 1
    merged["steps"] = np.asarray(steps, dtype=np.int64)
    return merged


if __name__ == "__main__":
    if sys.argv[1] == "--one":
        _one(Path(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
        sys.exit(0)
    root = Path(sys.argv[1])
    for bp in sorted(root.glob("*.bp")):
        ref = _dump(bp)
        np.savez_compressed(bp.with_suffix(".reference.npz"), **ref)
        print("wrote", bp.with_suffix(".reference.npz"), "steps", list(ref["steps"]))
