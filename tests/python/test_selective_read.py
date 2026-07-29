"""Selective/partial reads and ``read_metadata``.

The contract these pin: the options mean the same thing for *every* format.
Readers with a native selective path (vtu/vtp/xdmf/gmsh) skip the unwanted
arrays outright; everything else is read whole and trimmed afterwards. Same
answer either way -- only the cost differs, and ``read_metadata`` says which
happened via ``fell_back_to_full_read``.
"""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus

# Formats with a native selective/metadata path, and some without.
NATIVE = ["vtu", "vtp", "xdmf"]
FALLBACK = ["stl", "obj"]

EXT = {
    "vtu": ".vtu",
    "vtp": ".vtp",
    "xdmf": ".xdmf",
    "gmsh": ".msh",
    "stl": ".stl",
    "obj": ".obj",
}


def _mesh():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
        point_data={"u": np.arange(4.0), "v": np.arange(4.0) * 2},
    )


def _write(tmp_path, fmt):
    path = tmp_path / f"m{EXT[fmt]}"
    # ascii for vtu/vtp: these tests are about selective-read/metadata
    # plumbing, not codecs, and the default `compression="zlib"` would make
    # the native metadata/read path depend on zlib being compiled in (it
    # isn't on Windows CI -- see CMakeLists.txt's MESHIOPLUSPLUS_WITH_ZLIB
    # comment), silently forcing the full-read fallback these tests assert
    # against.
    kwargs = {"binary": False} if fmt in ("vtu", "vtp") else {}
    meshioplusplus.write(path, _mesh(), **kwargs)
    return path


@pytest.mark.parametrize("fmt", NATIVE + FALLBACK)
def test_default_read_is_unchanged(tmp_path, fmt):
    """The whole point: defaults must behave exactly as before."""
    path = _write(tmp_path, fmt)
    mesh = meshioplusplus.read(path)
    assert len(mesh.points) == 4
    assert sum(len(c.data) for c in mesh.cells) == 2


@pytest.mark.parametrize("fmt", NATIVE)
def test_points_only_keeps_geometry_drops_data(tmp_path, fmt):
    path = _write(tmp_path, fmt)
    full = meshioplusplus.read(path)
    bare = meshioplusplus.read(path, points_only=True)

    assert bare.point_data == {}
    assert bare.cell_data == {}
    # Connectivity survives -- points_only narrows data, not topology.
    np.testing.assert_allclose(np.asarray(bare.points), np.asarray(full.points))
    assert [c.type for c in bare.cells] == [c.type for c in full.cells]


@pytest.mark.parametrize("fmt", NATIVE)
def test_arrays_subset_returns_exactly_those(tmp_path, fmt):
    path = _write(tmp_path, fmt)
    got = meshioplusplus.read(path, arrays=["u"])
    assert sorted(got.point_data) == ["u"]
    np.testing.assert_allclose(got.point_data["u"], np.arange(4.0))


@pytest.mark.parametrize("fmt", NATIVE)
def test_empty_arrays_list_means_none_not_all(tmp_path, fmt):
    """``arrays=[]`` is "no arrays"; ``arrays=None`` is "every array"."""
    path = _write(tmp_path, fmt)
    assert meshioplusplus.read(path, arrays=[]).point_data == {}
    assert sorted(meshioplusplus.read(path, arrays=None).point_data) == ["u", "v"]


@pytest.mark.parametrize("fmt", FALLBACK)
def test_options_still_apply_to_formats_without_a_native_path(tmp_path, fmt):
    """Correctness must not depend on whether the format can be fast."""
    path = _write(tmp_path, fmt)
    bare = meshioplusplus.read(path, points_only=True)
    assert bare.point_data == {}
    assert bare.cell_data == {}
    assert len(bare.points) == 4


def test_unknown_requested_names_are_ignored(tmp_path):
    path = _write(tmp_path, "vtu")
    got = meshioplusplus.read(path, arrays=["no-such-array"])
    assert got.point_data == {}
    assert len(got.points) == 4


@pytest.mark.parametrize("fmt", NATIVE + FALLBACK)
def test_read_metadata_matches_a_full_read(tmp_path, fmt):
    path = _write(tmp_path, fmt)
    meta = meshioplusplus.read_metadata(path)
    mesh = meshioplusplus.read(path)

    assert meta["num_points"] == len(mesh.points)
    assert meta["num_cells"] == sum(len(c.data) for c in mesh.cells)
    assert [b["type"] for b in meta["cell_blocks"]] == [c.type for c in mesh.cells]
    assert meta["point_data_names"] == sorted(mesh.point_data)
    assert meta["format"] == fmt


@pytest.mark.parametrize("fmt", NATIVE)
def test_native_metadata_is_not_a_fallback(tmp_path, fmt):
    path = _write(tmp_path, fmt)
    meta = meshioplusplus.read_metadata(path)
    assert meta["fell_back_to_full_read"] is False
    # A native summary skips the coordinates, so it cannot report a bbox.
    assert "bbox_min" not in meta


@pytest.mark.parametrize("fmt", FALLBACK)
def test_fallback_metadata_says_so_and_can_afford_a_bbox(tmp_path, fmt):
    path = _write(tmp_path, fmt)
    meta = meshioplusplus.read_metadata(path)
    assert meta["fell_back_to_full_read"] is True
    assert "bbox_min" in meta and "bbox_max" in meta


def test_metadata_bbox_absent_rather_than_none(tmp_path):
    """ "Not computed" must not be mistakable for a real box at the origin."""
    meta = meshioplusplus.read_metadata(_write(tmp_path, "vtu"))
    assert "bbox_min" not in meta
    assert meta.get("bbox_min") is None  # only via .get, never a stored None


def test_metadata_regions_always_present_key(tmp_path):
    path = _write(tmp_path, "vtu")
    meta = meshioplusplus.read_metadata(path)
    # Always a key, even for a format (VTU) that carries no regions at all.
    assert meta["regions"] == []


def test_metadata_reports_regions_when_the_mesh_carries_them(tmp_path):
    mesh = meshioplusplus.Mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ],
        [("tetra", [[0, 1, 2, 4], [0, 2, 3, 4]])],
        regions=[meshioplusplus.Region("solid", "cell", [0, 1], dim=3, tag=7)],
    )
    path = str(tmp_path / "regions.msh")
    meshioplusplus.write(path, mesh, file_format="gmsh22")

    meta = meshioplusplus.read_metadata(path, file_format="gmsh")
    # 2.2 stores a type per element, so there is no cheap header-only summary to
    # give and it always falls back. (4.1 does have one, and reports its regions
    # natively -- see test_gmsh_41_metadata_reports_regions_without_a_full_read.)
    assert meta["fell_back_to_full_read"] is True
    assert meta["regions"] == [
        {"name": "solid", "kind": "cell", "dim": 3, "tag": 7, "num_entries": 2}
    ]


def test_cli_convert_points_only(tmp_path):
    from meshioplusplus._cli import main

    src = _write(tmp_path, "vtu")
    dst = tmp_path / "out.vtu"
    assert main(["convert", "--points-only", str(src), str(dst)]) is None or True
    assert meshioplusplus.read(dst).point_data == {}


def test_cli_convert_arrays(tmp_path):
    from meshioplusplus._cli import main

    src = _write(tmp_path, "vtu")
    dst = tmp_path / "out.vtu"
    main(["convert", "--arrays", "u", str(src), str(dst)])
    assert sorted(meshioplusplus.read(dst).point_data) == ["u"]


def test_cli_convert_rejects_selective_with_set_conversion(tmp_path):
    from meshioplusplus._cli import main

    src = _write(tmp_path, "vtu")
    dst = tmp_path / "out.vtu"
    with pytest.raises(ValueError, match="cannot be combined"):
        main(["convert", "--points-only", "-s", str(src), str(dst)])


def test_cli_info_fast(tmp_path, capsys):
    from meshioplusplus._cli import main

    main(["info", "--fast", str(_write(tmp_path, "vtu"))])
    out = capsys.readouterr().out
    assert "mesh summary" in out
    assert "Number of points: 4" in out
    assert "triangle: 2" in out
    # vtu has a native path, so no "was read in full" note.
    assert "read in full" not in out


def test_cli_info_fast_reports_fallback(tmp_path, capsys):
    from meshioplusplus._cli import main

    main(["info", "--fast", str(_write(tmp_path, "stl"))])
    assert "read in full" in capsys.readouterr().out


def test_gmsh_41_metadata_reports_regions_without_a_full_read():
    # $Entities and $PhysicalNames are both small and both sit ahead of
    # $Elements, so a 4.1 summary can name the physical groups -- and count
    # their cells -- from block headers alone. It must agree with a real read:
    # a summary that named different groups would be worse than none.
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[2]
    path = str(root / "tests/python/meshes/msh/insulated-4.1.msh")

    meta = meshioplusplus.read_metadata(path)
    assert meta["fell_back_to_full_read"] is False
    assert "gmsh:physical" in meta["cell_data_names"]

    mesh = meshioplusplus.read(path)
    assert sorted(
        (r["name"], r["dim"], r["tag"], r["num_entries"]) for r in meta["regions"]
    ) == sorted((r.name, r.dim, r.tag, len(r.entries)) for r in mesh.regions)
