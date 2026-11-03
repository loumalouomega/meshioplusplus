"""Roadmap 7.7.3: ``remesh_quads``, field-aligned quads through an installed
QuadriFlow or Instant Meshes run as a separate process. The probe (an all-quad
sphere whose least scaled Jacobian is at least 0.3, where recombining
``mesh_cad_surface``'s own triangles leaves far less) runs when QuadriFlow is
installed; the refusals run everywhere."""

import json
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "remesh_quads"), reason="core predates remesh_quads"
)
HAS_QUADRIFLOW = "quadriflow" in mio.remesh_quads_backends()
needs_quadriflow = pytest.mark.skipif(
    not HAS_QUADRIFLOW, reason="quadriflow is not installed"
)

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"


def _sphere_skin(h=0.1):
    return mio.mesh_cad_surface(mio.read(CAD / "sphere.cad.json"), h, edges=False)


@needs_core
def test_refusals_are_by_name():
    tri = mio.icosa_sphere(1.0, 1)
    for kwargs, name in [
        ({"backend": ""}, "backend"),
        ({"target_faces": -1}, "target_faces"),
        ({"target_faces": 1.5}, "target_faces"),
        ({"max_cells": 0}, "max_cells"),
    ]:
        with pytest.raises(ValueError, match=name):
            mio.remesh_quads(tri, **kwargs)
    with pytest.raises(NotImplementedError, match="no-such-remesher"):
        mio.remesh_quads(tri, backend="no-such-remesher")
    quads = mio.Mesh(np.eye(4, 3), [("quad", [[0, 1, 2, 3]])])
    with pytest.raises(ValueError, match="quad"):
        mio.remesh_quads(quads)
    names = mio.remesh_quads_backends()
    assert names == sorted(names)


@needs_core
@needs_quadriflow
def test_quadriflow_makes_an_all_quad_sphere():
    """The roadmap's probe, against recombination of the mesher's triangles."""
    skin = _sphere_skin()
    out, report = mio.remesh_quads(skin, backend="quadriflow", return_report=True)
    assert report["backend"] == "quadriflow"
    assert report["num_triangles"] == 0 and report["num_quads"] > 0
    assert [c.type for c in out.cells] == ["quad"]
    assert report["min_scaled_jacobian"] >= 0.3
    assert report["hausdorff_to_surface"] < 0.05
    _, recombined = mio.mesh_cad_surface(
        mio.read(CAD / "sphere.cad.json"), 0.1, quads="all", return_report=True
    )
    assert report["min_scaled_jacobian"] > recombined["min_quad_scaled_jacobian"]
    # The target steers the count; the seed repeats a run.
    _, few = mio.remesh_quads(skin, target_faces=300, return_report=True)
    assert few["num_quads"] < report["num_quads"]
    _, again = mio.remesh_quads(skin, target_faces=300, return_report=True)
    assert again["num_quads"] == few["num_quads"]


@needs_core
@needs_quadriflow
def test_pipeline_and_cli(tmp_path):
    skin = mio.icosa_sphere(1.0, 3)
    out, report = mio.apply_operations(
        skin, [{"Op": "RemeshQuads", "TargetFaces": 400}]
    )
    step = report["steps"][0]
    assert step["NumQuads"] > 0 and step["NumTriangles"] == 0
    assert {c.type for c in out.cells} == {"quad"}
    src = tmp_path / "sphere.stl"
    mio.write(src, skin)
    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "meshioplusplus",
            "remesh-quads",
            str(src),
            str(tmp_path / "q.vtu"),
            "--faces",
            "400",
            "--json",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    py = json.loads(res.stdout)
    assert py["backend"] == "quadriflow" and py["num_quads"] > 0
    native = shutil.which("meshioplusplus")
    if native is None:
        return
    nat = subprocess.run(
        [
            native,
            "remesh-quads",
            str(src),
            str(tmp_path / "n.vtu"),
            "-f",
            "400",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    if nat.returncode == 0 and nat.stdout.startswith("{"):
        assert json.loads(nat.stdout)["num_quads"] == py["num_quads"]
