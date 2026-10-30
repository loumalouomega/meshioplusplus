"""Smoke test of the base wheel alone (wheels.yml's CIBW_TEST_COMMAND;
standalone, no pytest): it carries no OpenCASCADE, so STEP reading raises an
error naming the add-on, ``fit_cad`` falls back to the native backend, and the
shared core and its CMake package are in the wheel."""

import os
import tempfile

import numpy as np

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus._fit_cad import has_occt

assert not has_occt() and not _core.__has_occt__, "the base wheel must carry no OpenCASCADE"

path = os.path.join(tempfile.mkdtemp(), "cube.step")
points = np.array([[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)])
quads = np.array([[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [1, 3, 7, 5], [3, 2, 6, 7], [2, 0, 4, 6]])
mp.write(path, mp.Mesh(points, [("quad", quads)]))
try:
    mp.read(path)
except NotImplementedError as exc:
    assert "meshioplusplus-occt" in str(exc), exc
else:
    raise AssertionError("STEP reading must need the add-on")

cmake_dir = mp.get_cmake_dir()
assert os.path.isfile(os.path.join(cmake_dir, "meshioplusplusConfig.cmake")), cmake_dir
print("base wheel smoke OK:", mp.__version__)
