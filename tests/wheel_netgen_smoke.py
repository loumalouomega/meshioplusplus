"""Smoke test of the meshioplusplus-netgen add-on wheel installed beside the base
wheel (wheels.yml's CIBW_TEST_COMMAND; standalone, no pytest): the add-on
registers its tetrahedral mesher with the base wheel's shared core, and
``tetrahedralize`` fills a cube skin with tetrahedra, keeping its triangles."""

import meshioplusplus as mp
from meshioplusplus import _core, _plugins

assert "netgen" in _plugins.load_plugins(), _plugins.plugin_errors()
assert _core.__has_netgen__ and mp.has_netgen(), "the add-on registered no mesher"

skin = mp.box_surface(cell_type="triangle")
solid, report = mp.tetrahedralize(
    skin, backend="netgen", max_size=0.5, return_report=True
)
assert report["backend"] == "netgen", report["backend"]
assert report["boundary_preserved"] and report["num_tetrahedra"] > 0, report
assert abs(report["volume"] - 1.0) < 1e-9, report["volume"]
print("Netgen wheel smoke OK:", report["num_tetrahedra"], "tetrahedra")
