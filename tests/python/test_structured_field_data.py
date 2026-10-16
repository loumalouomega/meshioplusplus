"""Roadmap 3.3.7: field_data and regions in the structured-grid writers
(`.vti`, `.vtr`, `.vts`), both engines."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.vti import _vti
from meshioplusplus.vtr import _vtr
from meshioplusplus.vts import _vts

TWINS = {"vti": _vti, "vtr": _vtr, "vts": _vts}


def _lattice():
    n = 3
    points = np.array(
        [[i, j, k] for k in range(n) for j in range(n) for i in range(n)], float
    )

    def idx(i, j, k):
        return i + n * (j + n * k)

    cells = [
        [
            idx(i, j, k),
            idx(i + 1, j, k),
            idx(i + 1, j + 1, k),
            idx(i, j + 1, k),
            idx(i, j, k + 1),
            idx(i + 1, j, k + 1),
            idx(i + 1, j + 1, k + 1),
            idx(i, j + 1, k + 1),
        ]
        for k in range(2)
        for j in range(2)
        for i in range(2)
    ]
    return meshioplusplus.Mesh(
        points,
        [("hexahedron", cells)],
        field_data={"g": np.array([9.81]), "v": np.array([1, 2, 3])},
        regions=[
            meshioplusplus.Region("walls", "cell", [0, 3, 5], dim=3, tag=7),
            meshioplusplus.Region("inlet", "point", [1, 2]),
        ],
    )


def _check(back):
    assert sorted(back.field_data) == ["g", "v"]
    assert back.field_data["g"][0] == 9.81
    np.testing.assert_array_equal(back.field_data["v"], [1, 2, 3])
    regions = {r.name: r for r in back.regions}
    assert sorted(regions) == ["inlet", "walls"]
    assert list(regions["walls"].entries) == [0, 3, 5]
    assert (regions["walls"].dim, regions["walls"].tag) == (3, 7)
    assert list(regions["inlet"].entries) == [1, 2]


@pytest.mark.parametrize("fmt", ["vti", "vtr", "vts"])
@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("engine", ["core", "python"])
def test_field_data_and_regions_round_trip(fmt, binary, engine, tmp_path):
    path = tmp_path / f"m.{fmt}"
    mesh = _lattice()
    if engine == "core":
        getattr(_core, f"{fmt}_write")(str(path), mesh, binary, True)
    else:
        TWINS[fmt].write(path, mesh, binary=binary)
    _check(meshioplusplus.read(path))
    _check(TWINS[fmt].read(path))  # the Python reader agrees with the core's


@pytest.mark.parametrize("fmt", ["vti", "vtr", "vts"])
def test_engines_write_the_same_binary_bytes(fmt, tmp_path):
    a, b = tmp_path / f"core.{fmt}", tmp_path / f"twin.{fmt}"
    mesh = _lattice()
    getattr(_core, f"{fmt}_write")(str(a), mesh, True, True)
    TWINS[fmt].write(b, mesh, binary=True)
    assert a.read_bytes() == b.read_bytes()
