"""The probe of roadmap 3.2.1: on Netgen's tetrahedra, before and after Netgen's
own optimisation passes, ``optimize_volume`` (the default ``worst`` objective)
leaves the least dihedral angle no lower than it found it. Netgen comes from its
Python package (``pip install netgen-mesher``); the test skips without it.
``tools/optvol_probe.py`` prints the full table."""

import numpy as np
import pytest

import meshioplusplus as mp

netgen_meshing = pytest.importorskip("netgen.meshing")

try:
    from meshioplusplus import _core  # noqa: F401
except ImportError:  # pragma: no cover - a pure-Python build
    pytest.skip("needs the compiled core", allow_module_level=True)


def _netgen_volume(skin, maxh, optsteps):
    nm = netgen_meshing
    mesh = nm.Mesh(dim=3)
    ids = [mesh.Add(nm.MeshPoint(nm.Pnt(*map(float, p)))) for p in skin.points]
    fd = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    for tri in np.asarray(skin.cells[0].data):
        mesh.Add(nm.Element2D(fd, [ids[int(v)] for v in tri]))
    mesh.GenerateVolumeMesh(mp=nm.MeshingParameters(maxh=maxh, optsteps3d=optsteps))
    points = np.array([[p[0], p[1], p[2]] for p in mesh.Points()])
    tets = np.array([[v.nr - 1 for v in el.vertices] for el in mesh.Elements3D()])
    return mp.Mesh(points, [("tetra", tets)])


def _least(mesh):
    q = mp.attach_quality(mesh).cell_data["quality:min_dihedral"]
    return float(np.nanmin(np.concatenate([np.asarray(a).ravel() for a in q])))


@pytest.mark.parametrize("optsteps", [0, 3])
def test_worst_never_lowers_the_least_dihedral_angle_on_netgen_output(optsteps):
    vol = _netgen_volume(mp.icosa_sphere(1.0, 3), 0.2, optsteps)
    out, report = mp.optimize_volume(vol, return_report=True)
    assert report["min_dihedral_before"] == pytest.approx(_least(vol), abs=1e-9)
    assert report["min_dihedral_after"] == pytest.approx(_least(out), abs=1e-9)
    assert report["min_dihedral_after"] >= report["min_dihedral_before"]
    # The boundary is the skin, and nothing inverts.
    assert len(out.points) == len(vol.points)
    assert mp.compute_stats(out)["num_inverted"] == 0
