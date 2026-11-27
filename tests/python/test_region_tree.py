"""Tests for roadmap §8.5 hierarchical regions (/ paths + meshio:tree)."""

import h5py
import numpy as np

import meshioplusplus as mp


def _nested_mesh():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    cells = [
        ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
        ("triangle", np.array([[3, 4, 5]], dtype=np.int64)),
    ]
    regions = [
        mp.Region("Vehicle/Wing/Skin", "cell", np.array([0], dtype=np.int64)),
        mp.Region("Vehicle/Fuselage/Core", "cell", np.array([1], dtype=np.int64)),
    ]
    return mp.Mesh(
        pts,
        cells,
        regions=regions,
        field_data={"meshio:tree": ["Vehicle/EmptySubassembly"]},
    )


def test_vtm_nested_hierarchy_roundtrip(tmp_path):
    m = _nested_mesh()
    out = tmp_path / "nested.vtm"
    mp.write(out, m)

    xml = out.read_text()
    assert 'name="Vehicle"' in xml
    assert 'name="Wing"' in xml
    assert 'name="Skin"' in xml
    assert 'name="EmptySubassembly"' in xml

    back = mp.read(out)
    names = sorted(r.name for r in back.regions if r.kind == "cell")
    assert names == ["Vehicle/Fuselage/Core", "Vehicle/Wing/Skin"]
    assert back.field_data["meshio:tree"] == ["Vehicle/EmptySubassembly"]


def test_vtkhdf_nested_assembly_roundtrip(tmp_path):
    m = _nested_mesh()
    out = tmp_path / "nested.vtkhdf"
    mp.vtkhdf.write(out, m, dataset_type="MultiBlockDataSet")

    with h5py.File(out, "r") as f:
        asm = f["VTKHDF/Assembly"]
        assert "Vehicle/Wing/Skin" in asm
        assert "Vehicle/Fuselage/Core" in asm
        assert "Vehicle/EmptySubassembly" in asm

    back = mp.read(out)
    names = sorted(r.name for r in back.regions if r.kind == "cell")
    assert names == ["Vehicle/Fuselage/Core", "Vehicle/Wing/Skin"]
    assert back.field_data["meshio:tree"] == ["Vehicle/EmptySubassembly"]
