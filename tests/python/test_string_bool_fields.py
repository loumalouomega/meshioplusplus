"""Tests for roadmap §8.2 string, string-list, and boolean field conventions."""

import numpy as np

import meshioplusplus as mp


def _tri_mesh(**kwargs):
    pts = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
    )
    cells = [("triangle", np.array([[0, 1, 2]], dtype=np.int64))]
    return mp.Mesh(pts, cells, **kwargs)


def test_exodus_string_title_roundtrip_strict_core(tmp_path):
    """Roadmap 8.2 Done when: Mesh(..., field_data={'title': 'wing'}) writes
    through the C++ core to Exodus and reads back equal."""
    m = _tri_mesh(field_data={"title": "wing"})
    out = tmp_path / "wing.exo"
    mp.write(out, m)
    back = mp.read(out)
    assert back.field_data["title"] == "wing"


def test_vtu_string_list_and_bool_roundtrip_strict_core(tmp_path):
    mask_pt = np.array([True, False, True], dtype=bool)
    mask_cell = [np.array([True], dtype=bool)]
    mask_fd = np.array([True, False], dtype=bool)
    m = _tri_mesh(
        point_data={"active": mask_pt},
        cell_data={"flag": mask_cell},
        field_data={
            "title": "wing-section",
            "notes": ["alpha", "beta", "gamma-λ"],
            "single_list": ["only-one"],
            "enabled": mask_fd,
        },
    )
    out = tmp_path / "strings_bools.vtu"
    mp.write(out, m)
    back = mp.read(out)

    assert back.field_data["title"] == "wing-section"
    assert back.field_data["notes"] == ["alpha", "beta", "gamma-λ"]
    assert back.field_data["single_list"] == ["only-one"]
    # Hidden metadata keys must not leak into Python field_data
    assert not any(k.endswith(":offsets") for k in back.field_data)
    assert not any(k.startswith("meshio:bool:") for k in back.field_data)

    assert back.point_data["active"].dtype == np.bool_
    np.testing.assert_array_equal(back.point_data["active"], mask_pt)
    assert back.cell_data["flag"][0].dtype == np.bool_
    np.testing.assert_array_equal(back.cell_data["flag"][0], mask_cell[0])
    assert back.field_data["enabled"].dtype == np.bool_
    np.testing.assert_array_equal(back.field_data["enabled"], mask_fd)


def test_cgns_string_and_bool_roundtrip_strict_core(tmp_path):
    m = _tri_mesh(
        point_data={"mask": np.array([False, True, False], dtype=bool)},
        field_data={"case_name": "airfoil", "tags": ["inlet", "outlet"]},
    )
    out = tmp_path / "str.cgns"
    mp.write(out, m)
    back = mp.read(out)
    assert back.field_data["case_name"] == "airfoil"
    assert back.field_data["tags"] == ["inlet", "outlet"]
    assert back.point_data["mask"].dtype == np.bool_
    np.testing.assert_array_equal(back.point_data["mask"], [False, True, False])
