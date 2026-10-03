"""Mapped and buffered reader views must never escape into a returned mesh."""

import numpy as np
import pytest

import meshioplusplus as pp
from meshioplusplus import _core

from . import helpers
from .text_io_helpers import mesh_arrays, native_snapshot

NATIVE_READERS = {
    "off": _core.off_read,
    "medit": _core.medit_read_ascii,
    "xyz": _core.xyz_read,
    "unv": _core.unv_read,
}
SUFFIXES = {"off": ".off", "medit": ".mesh", "xyz": ".xyz", "unv": ".unv"}
STREAM_READERS = {
    "ensight": _core.ensight_read,
    "vtk": _core.vtk_read,
    "ansys": _core.ansys_read,
    "patran": _core.patran_read,
    "tecplot": _core.tecplot_read,
    "femap": _core.femap_read,
    "pcd": _core.pcd_read,
    "radioss": _core.radioss_read,
    "dex": _core.dex_read,
    "ip": _core.ip_read,
    "flux": _core.flux_read,
    "permas": _core.permas_read,
    "abaqus": _core.abaqus_read,
    "ansysInp": _core.ansysinp_read,
    "marc": _core.marc_read,
    "nastran": _core.nastran_read,
    "netgen": _core.netgen_read,
    "lsdyna": _core.lsdyna_read,
    "gid": _core.gid_read,
    "gmsh": _core.gmsh_read,
    "mfem": _core.mfem_read,
    "mdpa": _core.mdpa_read,
    "flac3d": _core.flac3d_read,
    "obj": _core.obj_read,
    "ply": _core.ply_read,
}
STREAM_SUFFIXES = {
    "ensight": ".case",
    "vtk": ".vtk",
    "ansys": ".msh",
    "patran": ".pat",
    "tecplot": ".dat",
    "femap": ".neu",
    "pcd": ".pcd",
    "radioss": "_0000.rad",
    "dex": ".dex",
    "ip": ".ip",
    "flux": ".pf3",
    "permas": ".dat",
    "abaqus": ".inp",
    "ansysInp": ".cdb",
    "marc": ".dat",
    "nastran": ".bdf",
    "netgen": ".vol",
    "lsdyna": ".k",
    "gid": ".post.msh",
    "gmsh": ".msh",
    "mfem": ".mesh",
    "mdpa": ".mdpa",
    "flac3d": ".f3grid",
    "obj": ".obj",
    "ply": ".ply",
}


@pytest.mark.parametrize("fmt", NATIVE_READERS)
@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_native_view_ownership_and_line_endings(
    tmp_path, monkeypatch, fmt, threshold, ending
):
    monkeypatch.setenv("MESHIOPLUSPLUS_MMAP_THRESHOLD", threshold)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    path = tmp_path / ("mesh" + SUFFIXES[fmt])
    source = helpers.point_cloud_mesh if fmt == "xyz" else helpers.tri_mesh
    pp.write(path, source, file_format=fmt)
    # Compare with the canonical LF reference. NumPy's fromfile behind the
    # Python OFF reader does not reliably mix CRLF with buffered text reads.
    reference = getattr(pp, fmt)._py_read(path)
    text = path.read_text().rstrip("\n")
    if ending == "\r\n":
        text = text.replace("\n", "\r\n")
    path.write_bytes((text + ending).encode())
    native = native_snapshot(path, NATIVE_READERS[fmt].__name__, threshold, tmp_path)
    expected = mesh_arrays(reference)
    assert sorted(native) == sorted(expected)
    for name in native:
        np.testing.assert_array_equal(native[name], expected[name])


def test_xyz_native_embedded_nul_keeps_c_string_prefix_semantics(tmp_path):
    path = tmp_path / "nul.xyz"
    path.write_bytes(b"1\0suffix 2 3\n")
    mesh = _core.xyz_read(str(path))
    np.testing.assert_array_equal(mesh.points, [[1, 2, 3]])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_gmsh_float_spelled_node_ids_and_hex_coordinates(tmp_path, threshold):
    path = tmp_path / "mesh.msh"
    path.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n"
        "2.0 1e0 -0 0x1p1\n1e3 3 4 5\n$EndNodes\n"
        "$Elements\n0\n$EndElements"
    )
    expected = mesh_arrays(_core.gmsh_read(str(path)))
    actual = native_snapshot(path, "gmsh_read", threshold, tmp_path)
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], [[1, 0, 2], [3, 4, 5]])
    assert np.signbit(actual["points"][0, 1])


@pytest.mark.parametrize("delimiter", ["::", "||"])
def test_xyz_native_multichar_delimiter_and_trailing_empty(tmp_path, delimiter):
    path = tmp_path / "delimiter.xyz"
    path.write_text(delimiter.join([" 1 ", " +2 ", " 3e0 ", ""]))
    mesh = _core.xyz_read(str(path), [], delimiter)
    np.testing.assert_array_equal(mesh.points, [[1, 2, 3]])


@pytest.mark.parametrize("fmt", STREAM_READERS)
@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
@pytest.mark.parametrize("ending", ["\n", "\r\n", ""])
def test_remaining_stream_readers_ownership(
    tmp_path, monkeypatch, fmt, threshold, ending
):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    path = tmp_path / ("mesh" + STREAM_SUFFIXES[fmt])
    kwargs = (
        {"binary": False}
        if fmt in ("ensight", "vtk", "ansys", "pcd", "gmsh", "ply")
        else {}
    )
    source = helpers.point_cloud_mesh if fmt == "pcd" else helpers.tri_mesh
    if fmt == "dex":
        source = pp.Mesh(
            source.points,
            [],
            point_data={"field": np.arange(len(source.points), dtype=float)},
        )
    pp.write(path, source, file_format=fmt, **kwargs)
    expected = mesh_arrays(STREAM_READERS[fmt](str(path)))
    for file in tmp_path.iterdir():
        if file.is_file():
            text = file.read_text().rstrip("\n")
            if ending == "\r\n":
                text = text.replace("\n", "\r\n")
            file.write_bytes((text + ending).encode())
    native = native_snapshot(path, STREAM_READERS[fmt].__name__, threshold, tmp_path)
    assert sorted(native) == sorted(expected)
    for name in native:
        np.testing.assert_array_equal(native[name], expected[name])


def _pcd_numeric_file(path, token, kind="I"):
    path.write_bytes(
        (
            "VERSION .7\nFIELDS x y z tag\nSIZE 8 8 8 8\nTYPE F F F "
            + kind
            + "\nCOUNT 1 1 1 1\nWIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA ascii\n"
        ).encode()
        + b"1."
        + b"0" * 80
        + b"e0 2 3 "
        + token
    )
    return path


@pytest.mark.parametrize(
    "token,kind,value",
    [
        (b"+9223372036854775807", "I", (1 << 63) - 1),
        (b"-9223372036854775808", "I", -(1 << 63)),
        (b"-1", "U", (1 << 64) - 1),
        (b"7\0ignored", "I", 7),
    ],
)
def test_pcd_bounded_numbers_keep_native_c_parse_semantics(
    tmp_path, token, kind, value
):
    path = _pcd_numeric_file(tmp_path / "number.pcd", token, kind)
    mesh = _core.pcd_read(str(path))
    np.testing.assert_array_equal(mesh.points, [[1, 2, 3]])
    assert int(mesh.point_data["tag"][0]) == value


@pytest.mark.parametrize(
    "token", [b"9223372036854775808", b"-9223372036854775809", b"7suffix", b"nope"]
)
def test_pcd_bounded_numbers_retain_integer_errors(tmp_path, token):
    path = _pcd_numeric_file(tmp_path / "number.pcd", token)
    with pytest.raises(pp.ReadError, match="non-numeric"):
        _core.pcd_read(str(path))


def test_radioss_engine_bounded_number_tokens(tmp_path):
    path = tmp_path / "numbers_0001.rad"
    path.write_text(
        "#RADIOSS ENGINE\n/ANIM/DT\n1." + "0" * 80 + "e0 +2e0 junk 3suffix\n/END"
    )
    mesh = _core.radioss_read(str(path))
    np.testing.assert_array_equal(mesh.field_data["radioss:engine:ANIM/DT"], [1, 2])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_marc_nested_include_sources_stay_alive(tmp_path, threshold):
    path = tmp_path / "deck.dat"
    path.write_text("title\nmodel\nend\ninclude first.dat\nend option\n")
    first = tmp_path / "first.dat"
    first.write_text("coordinates\n3,1,\n10,1,2,3,\ninclude second.dat\n")
    second = tmp_path / "second.dat"
    second.write_text("20,4,5,6,\n")
    expected = mesh_arrays(_core.marc_read(str(path)))
    actual = native_snapshot(path, "marc_read", threshold, tmp_path, [first, second])
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], [[1, 2, 3], [4, 5, 6]])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_abaqus_nested_include_tokens_and_names_remain_owned(tmp_path, threshold):
    path = tmp_path / "deck.inp"
    first, second = tmp_path / "first.inp", tmp_path / "second.inp"
    path.write_text("*INCLUDE, INPUT=first.inp\n")
    first.write_text(
        "*NODE\n10,1." + "0" * 80 + ",2,3\n"
        "*INCLUDE, INPUT=second.inp\n"
        "*NSET, NSET=a_long_named_node_set_after_the_include\n10,20,\n"
    )
    second.write_text("*NODE\n20,4,5,6\n")
    # The native reader already resolves ids across nested include sources;
    # the Python reference's include-local id map does not handle this deck.
    expected = mesh_arrays(_core.abaqus_read(str(path)))
    actual = native_snapshot(path, "abaqus_read", threshold, tmp_path, [first, second])
    assert actual.keys() == expected.keys()
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], [[1, 2, 3], [4, 5, 6]])
    assert actual["regions/meta"][0, 1] == "a_long_named_node_set_after_the_include"
    np.testing.assert_array_equal(actual["regions/0"], [0, 1])


def test_abaqus_tokens_keep_native_numeric_prefix_behavior(tmp_path):
    path = tmp_path / "prefix.inp"
    path.write_bytes(b"*NODE\n+10suffix,1." + b"0" * 80 + b"suffix,2\0ignored,junk,")
    mesh = _core.abaqus_read(str(path))
    np.testing.assert_array_equal(mesh.points, [[1, 2, 0]])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_nastran_mixed_large_free_continuation_storage(tmp_path, threshold):
    from meshioplusplus.nastran import _nastran as reference

    path = tmp_path / "mixed.bdf"
    path.write_text(
        "BEGIN BULK\n"
        + "".join(
            f"GRID*   {index:<16}{0:<16}{1.25:<16}{2:<16}\n*CONT,3,.\n"
            for index in range(1, 301)
        )
        + "ENDDATA"
    )
    expected = mesh_arrays(reference.read(path))
    actual = native_snapshot(path, "nastran_read", threshold, tmp_path)
    assert actual.keys() == expected.keys()
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], np.tile([1.25, 2, 3], (300, 1)))


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_netgen_gzip_token_views_do_not_escape_the_inflated_buffer(tmp_path, threshold):
    path = tmp_path / "mesh.vol.gz"
    pp.netgen.write(path, helpers.tet_mesh)
    expected = mesh_arrays(_core.netgen_read(str(path)))
    actual = native_snapshot(path, "netgen_read", threshold, tmp_path)
    assert actual.keys() == expected.keys()
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_dex_exponent_normalization_and_long_prefixes(tmp_path, threshold):
    path = tmp_path / "normalization.dex"
    path.write_bytes(
        b"# NB_POINTS = 1\n# NB_COMP = 1\n"
        b"1D0 +2d0 -0D0 3." + b"0" * 100 + b"D0suffix\0ignored"
    )
    expected = mesh_arrays(_core.dex_read(str(path)))
    actual = native_snapshot(path, "dex_read", threshold, tmp_path)
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], [[1, 2, 0]])
    assert np.signbit(actual["points"][0, 2])
    np.testing.assert_array_equal(actual["point_data/dex:field/0"], [3])


@pytest.mark.parametrize("threshold", ["0", str(1 << 40)])
def test_gid_quoted_names_prefix_numbers_and_source_ownership(tmp_path, threshold):
    path = tmp_path / "mesh.post.msh"
    path.write_text(
        'MESH "a long quoted mesh name" dimension 3 ElemType Point Nnode 1\n'
        "Coordinates\n+1suffix +1suffix 2e0suffix 3\nEnd Coordinates\n"
        "Elements\n+1suffix 1suffix +7suffix\nEnd Elements"
    )
    expected = mesh_arrays(_core.gid_read(str(path)))
    actual = native_snapshot(path, "gid_read", threshold, tmp_path)
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name])
    np.testing.assert_array_equal(actual["points"], [[1, 2, 3]])
