"""Smoke test of a built wheel's bundled OpenCASCADE (wheels.yml's
CIBW_TEST_COMMAND; standalone, no pytest): the core carries the kernel,
``fit_cad`` builds a valid solid with it on a raw cube skin, and the STEP it
writes reads back through the same kernel."""

import os
import tempfile

import numpy as np

import meshioplusplus as mp
from meshioplusplus import _core

assert _core.__has_occt__, "this wheel was built without OpenCASCADE"

points = np.array(
    [[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)]
)
quads = np.array(
    [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [1, 3, 7, 5],
        [3, 2, 6, 7],
        [2, 0, 4, 6],
    ]
)
skin = mp.Mesh(points, [("quad", quads)])

path = os.path.join(tempfile.mkdtemp(), "cube.step")
solid, report = mp.fit_cad(skin, path, return_report=True)
assert report["backend"] == "occt", report["backend"]
assert report["valid"] and report["num_faces"] == 6, report
assert abs(report["volume"] - 1.0) < 1e-9, report["volume"]

back = mp.read(path)
assert len(back.info["brep"]["faces"]) == 6
print(
    "OpenCASCADE wheel smoke OK:",
    report["num_faces"],
    "faces, volume",
    report["volume"],
)
