"""boundary_layer (roadmap 6.9): prismatic layers grown from a surface.

Through the Python wrapper: the sphere and the unit cube grown outward and
inward, the cube's cap at the asked thickness off every face, the face data
carried, an all-quad cube growing hexahedra, the refusals, the CLI verb and
the pipeline step through both engines. The algorithm itself is tested in the
core (``tests/cpp/test_boundary_layer.cpp``).
"""

import json
import pathlib

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

needs_core = pytest.mark.skipif(
    not hasattr(_core, "boundary_layer"), reason="core predates boundary_layer"
)

CAD = pathlib.Path(__file__).parent / "meshes" / "cad_json"


def _skin(name, h, **kw):
    return mio.mesh_cad_surface(mio.read(CAD / name), h, **kw)


def _volume_cells(mesh):
    q = mio.compute_quality(mesh)
    vol = 0.0
    sj = 1.0
    for block, cells in enumerate(mesh.cells):
        if cells.type in ("wedge", "hexahedron", "pyramid", "tetra"):
            vol += float(np.nansum(q["cell_arrays"]["quality:volume"][block]))
            sj = min(
                sj, float(np.nanmin(q["cell_arrays"]["quality:scaled_jacobian"][block]))
            )
    return vol, sj


@needs_core
@pytest.mark.parametrize("reverse", [False, True])
def test_the_sphere_grows_both_ways(reverse):
    skin = _skin("sphere.cad.json", 0.2)
    out, report = mio.boundary_layer(
        skin, 0.02, num_layers=3, reverse=reverse, return_report=True
    )
    assert report["num_wedges"] == report["num_cells"] > 0
    assert report["num_collapsed"] == 0 and report["cap_closed"]
    assert report["min_scaled_jacobian"] > 0.1
    t = 0.02 * (1 + 1.2 + 1.44)
    shell = 4 / 3 * np.pi * ((1 + t) ** 3 - 1 if not reverse else 1 - (1 - t) ** 3)
    vol, sj = _volume_cells(out)
    assert vol == pytest.approx(shell, rel=0.03)
    assert sj > 0.1
    layers = np.concatenate(
        [np.asarray(d) for d in out.cell_data["boundary_layer:layer"]]
    )
    assert set(layers.tolist()) == {0, 1, 2, 3, 4}
    assert "cad:face" in out.cell_data


@needs_core
def test_the_cubes_cap_keeps_its_thickness():
    skin = _skin("unit_cube.cad.json", 0.25)
    out, report = mio.boundary_layer(
        skin, 0.01, num_layers=4, growth=1.3, return_report=True
    )
    total = 0.01 * (1 + 1.3 + 1.3**2 + 1.3**3)
    assert report["min_thickness"] == pytest.approx(total, abs=1e-12)
    layer = out.cell_data["boundary_layer:layer"]
    cap = [i for i, d in enumerate(layer) if len(d) and int(d[0]) == 5][0]
    p = np.asarray(out.points)[np.unique(out.cells[cap].data)]
    far = np.maximum(-p, p - 1.0).max(axis=1)
    assert np.allclose(far, total, atol=1e-12)


@needs_core
def test_quads_grow_hexahedra():
    skin = _skin("unit_cube.cad.json", 0.25, quads="all")
    out, report = mio.boundary_layer(skin, 0.02, num_layers=2, return_report=True)
    assert report["num_hexahedra"] > 0 and report["num_wedges"] == 0
    assert {c.type for c in out.cells} == {"hexahedron", "quad"}


@needs_core
def test_refusals_are_by_name():
    skin = _skin("unit_cube.cad.json", 0.5)
    with pytest.raises(ValueError, match="first_height"):
        mio.boundary_layer(skin, 0.0)
    with pytest.raises(ValueError, match="num_layers"):
        mio.boundary_layer(skin, 0.01, num_layers=0)
    with pytest.raises(ValueError, match="growth"):
        mio.boundary_layer(skin, 0.01, growth=0.5)
    tet = mio.Mesh(np.eye(4, 3), [("tetra", [[0, 1, 2, 3]])])
    with pytest.raises(ValueError, match="tetra"):
        mio.boundary_layer(tet, 0.01)


@needs_core
def test_the_cli_verb(tmp_path, capsys):
    skin = tmp_path / "skin.vtu"
    mio.write(skin, _skin("sphere.cad.json", 0.3))
    out = tmp_path / "layers.vtu"
    code = mio._cli.main(
        ["boundary-layer", str(skin), str(out), "--first-height", "0.01", "--json"]
    )
    assert code in (0, None)
    report = json.loads(capsys.readouterr().out)
    assert report["num_wedges"] > 0 and report["cap_closed"]
    assert "wedge" in {c.type for c in mio.read(out).cells}


@needs_core
@pytest.mark.parametrize("version", [1, 2])
def test_the_pipeline_step(tmp_path, version):
    doc = {
        "Version": version,
        "Input": {"Path": str(CAD / "sphere.cad.json")},
        "Operations": [
            {"Op": "MeshCadSurface", "TargetSize": 0.3},
            {"Op": "BoundaryLayer", "FirstHeight": 0.01, "NumLayers": 2},
        ],
    }
    engines = [("py", mio.run_pipeline)]
    if hasattr(_core, "run_pipeline_json"):
        engines.append(("cpp", lambda d: _core.run_pipeline_json(json.dumps(d))))
    counts = []
    for name, run in engines:
        path = tmp_path / f"{name}.vtu"
        report = run(dict(doc, Output={"Path": str(path)}))
        step = report["steps"][-1]
        assert step["NumCells"] > 0 and step["NumCollapsed"] == 0, name
        counts.append(step["NumCells"])
        assert "wedge" in {c.type for c in mio.read(path).cells}
    assert len(set(counts)) == 1
    bad = dict(doc, Operations=[doc["Operations"][0], {"Op": "BoundaryLayer"}])
    for _, run in engines:
        with pytest.raises(ValueError, match="FirstHeight"):
            run(dict(bad, Output={"Path": str(tmp_path / "bad.vtu")}))
