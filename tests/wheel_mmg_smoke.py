"""Smoke test of the meshioplusplus-mmg add-on wheel installed beside the base
wheel (wheels.yml's CIBW_TEST_COMMAND; standalone, no pytest): the add-on
registers its remesher with the base wheel's shared core, ``adapt`` remeshes a
tetrahedral sphere, and ``move`` (which needs ELAS) raises naming it."""

import numpy as np

import meshioplusplus as mp
from meshioplusplus import _core, _plugins

assert "mmg" in _plugins.load_plugins(), _plugins.plugin_errors()
assert _core.__has_mmg__ and mp.has_mmg(), "the add-on registered no remesher"
assert not mp.has_mmg_elas(), "the wheel must be built without ELAS (GPL-3.0)"

sphere = mp.tetra_sphere(radius=1.0, subdivisions=2)
adapted, report = mp.adapt(sphere, hsiz=0.4, return_report=True)
assert len(adapted.points) > 0 and len(adapted.cells) > 0, report
try:
    mp.move(sphere, np.zeros((len(sphere.points), 3)), "boundary")
except NotImplementedError as exc:
    assert "ELAS" in str(exc), exc
else:
    raise AssertionError("move must not run without ELAS")
print("MMG wheel smoke OK:", len(adapted.points), "points")
