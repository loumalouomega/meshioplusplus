"""Tests for ``select`` (roadmap 3.1.3): the Python face over the core, its
argument checks, the CLI verb and the pipeline step in both engines. The C++
behaviour is pinned in ``tests/cpp/test_select.cpp``."""

import json

import numpy as np
import pytest

import meshioplusplus as mio

from .analytic_surfaces import box_skin

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def _cube():
    """A 2 x 2 x 2 hexahedral cube of side 2, with side regions xmin ... zmax."""
    return mio.box(dims=(2, 2, 2), size=(2.0, 2.0, 2.0))


def _region(mesh, name, kind):
    found = [r for r in mesh.regions if r.name == name and r.kind == kind]
    assert len(found) == 1
    return found[0]


@needs_core
def test_box_sides_match_the_primitive_region():
    cube = _cube()
    out, report = mio.select(
        cube, "left", "side", box=(-1, -1, -1, 1e-9, 3, 3), return_report=True
    )
    assert report == {
        "num_candidates": 24,
        "num_selected": 4,
        "num_grown": -1,
        "num_components": -1,
    }
    left = _region(out, "left", "side")
    assert left.dim == 2
    np.testing.assert_array_equal(left.entries, _region(cube, "xmin", "side").entries)


@needs_core
def test_is_data_only():
    cube = _cube()
    out = mio.select(cube, "half", plane=((1, 0, 0), (1, 0, 0)), mode="centroid")
    np.testing.assert_array_equal(out.points, cube.points)
    assert len(out.cells) == len(cube.cells)
    np.testing.assert_array_equal(out.cells[0].data, cube.cells[0].data)
    assert len(out.regions) == len(cube.regions) + 1
    assert len(_region(out, "half", "cell").entries) == 4
    # The input is untouched, and a rerun replaces rather than duplicates.
    assert not any(r.name == "half" for r in cube.regions)
    again = mio.select(out, "half", box=((-1, -1, -1), (3, 3, 3)), tag=4)
    assert len(again.regions) == len(out.regions)
    assert _region(again, "half", "cell").tag == 4
    assert len(_region(again, "half", "cell").entries) == 8


@needs_core
def test_where_quality_component_and_points():
    cube = _cube()
    cube.cell_data["mat"] = [np.array([1, 2, 1, 2, 1, 2, 1, 2], dtype=float)]
    out = mio.select(cube, "two", where=("mat", "==", 2))
    np.testing.assert_array_equal(_region(out, "two", "cell").entries, [1, 3, 5, 7])
    _, report = mio.select(
        cube, "good", quality=("scaled_jacobian", ">", 0.5), return_report=True
    )
    assert report["num_selected"] == 8
    _, report = mio.select(cube, "c0", component=0, return_report=True)
    assert report["num_components"] == 1 and report["num_selected"] == 8
    _, report = mio.select(
        cube, "ball", "point", sphere=((0, 0, 0), 1.0), return_report=True
    )
    assert report["num_selected"] == 4


@needs_core
def test_grow_stops_at_feature_edges():
    skin = box_skin(3)
    _, report = mio.select(skin, "top", seed=(0.5, 0.5, 1.2), return_report=True)
    assert report["num_grown"] == 9 and report["num_selected"] == 9
    _, report = mio.select(
        skin, "all", seed=(0.5, 0.5, 1.2), feature_angle=95, return_report=True
    )
    assert report["num_selected"] == 54


@needs_core
def test_feeds_find_interface():
    cube = _cube()
    left = mio.select(cube, "left", plane=((1, 0, 0), (-1, 0, 0)), mode="centroid")
    both = mio.select(left, "right", plane=((1, 0, 0), (1, 0, 0)), mode="centroid")
    _, report = mio.find_interface(both, "left", "right", return_report=True)
    assert report["num_pairs"] == 4


def test_argument_checks():
    cube = mio.box(dims=(1, 1, 1))
    with pytest.raises(ValueError, match="name"):
        mio.select(cube, "", box=(0, 0, 0, 1, 1, 1))
    with pytest.raises(ValueError, match="criterion"):
        mio.select(cube, "x")
    with pytest.raises(ValueError, match="kind"):
        mio.select(cube, "x", "edge", component=0)
    with pytest.raises(ValueError, match="mode"):
        mio.select(cube, "x", component=0, mode="most")
    with pytest.raises(ValueError, match="comparison"):
        mio.select(cube, "x", where=("a", "~", 1))
    with pytest.raises(ValueError, match="box"):
        mio.select(cube, "x", box=(0, 0, 1))
    with pytest.raises(ValueError, match="seed"):
        mio.select(cube, "x", seed=(0, 0))
    with pytest.raises(ValueError, match="component"):
        mio.select(cube, "x", component=-2)


@needs_core
def test_core_errors_are_value_errors():
    cube = _cube()
    with pytest.raises(ValueError, match="no_such"):
        mio.select(cube, "x", where=("no_such", "<", 1))
    with pytest.raises(ValueError, match="metric"):
        mio.select(cube, "x", quality=("no_such", "<", 1))
    with pytest.raises(ValueError, match="points"):
        mio.select(cube, "x", "point", quality=("scaled_jacobian", "<", 1))
    with pytest.raises(ValueError, match="surface"):
        mio.select(cube, "x", seed=(0, 0, 0))


@needs_core
def test_cli_verb(tmp_path, capsys):
    src = str(tmp_path / "cube.vtu")
    mio.write(src, _cube())
    out = str(tmp_path / "sel.vtu")
    args = ["select", src, out, "--name", "top", "--kind", "side", "--seed", "1,1,5"]
    assert mio._cli.main(args + ["--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["num_selected"] == 4 and report["num_grown"] == 4
    back = mio.read(out)
    assert len(_region(back, "top", "side").entries) == 4
    assert (
        mio._cli.main(
            ["select", src, out, "-n", "hot", "--quality", "scaled_jacobian > 0.5"]
        )
        == 0
    )
    assert "8 of 8 cell(s)" in capsys.readouterr().out


@needs_core
def test_pipeline_step_in_both_engines(tmp_path):
    src = tmp_path / "cube.vtu"
    mio.write(src, _cube(), compression=None)
    settings = {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [
            {"Op": "Select", "Name": "top", "Kind": "side", "Seed": [1, 1, 5]},
            {
                "Op": "Select",
                "Name": "near",
                "Kind": "point",
                "Center": [0, 0, 0],
                "Radius": 1.0,
            },
        ],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }
    reports = [mio.run_pipeline(settings)]
    if getattr(_core, "__has_json__", False):
        settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
        reports.append(_core.run_pipeline_json(json.dumps(settings)))
        assert mio.meshes_equal(
            mio.read(str(tmp_path / "out.vtu")), mio.read(settings["Output"]["Path"])
        )
    for report in reports:
        top, near = report["steps"]
        assert top["NumSelected"] == 4 and top["NumGrown"] == 4
        assert near["NumSelected"] == 4
