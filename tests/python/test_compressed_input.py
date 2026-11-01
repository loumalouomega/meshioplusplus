"""Compressed mesh files read transparently (roadmap 7.1.1): ``cube.msh.gz``,
``.bz2``, ``.xz``, ``.zst`` and a single-member ``.zip``, through the Python
twin (``_compressed.py``) and the native registry (``_core.read_metadata``,
both CLIs)."""

import bz2
import gzip
import importlib.util
import lzma
import os
import subprocess
import tempfile
import zipfile

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _compressed

from . import helpers

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None


def _zstd_compress(data):
    try:
        from compression import zstd

        return zstd.compress(data)
    except ImportError:
        pass
    if importlib.util.find_spec("zstandard"):
        import zstandard

        return zstandard.ZstdCompressor().compress(data)
    pytest.skip("no zstd compressor in this Python")


COMPRESSORS = {
    ".gz": gzip.compress,
    ".bz2": bz2.compress,
    ".xz": lzma.compress,
    ".zst": _zstd_compress,
}


def _zip(data, name, method=zipfile.ZIP_DEFLATED):
    path = tempfile.mktemp(suffix=".zip")
    with zipfile.ZipFile(path, "w", method) as z:
        z.writestr(name, data)
    with open(path, "rb") as fh:
        out = fh.read()
    os.unlink(path)
    return out


@pytest.fixture
def cube(tmp_path):
    mesh = helpers.tet_mesh
    paths = {}
    for name, fmt in (("cube.msh", "gmsh"), ("cube.vtu", "vtu"), ("cube.stl", None)):
        target = mesh if fmt else helpers.tri_mesh
        mp.write(str(tmp_path / name), target, file_format=fmt or "stl")
        paths[name] = (tmp_path / name).read_bytes()
    return tmp_path, paths


def _same(a, b):
    np.testing.assert_allclose(a.points, b.points, atol=1e-12)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)


@pytest.mark.parametrize("suffix", sorted(COMPRESSORS))
@pytest.mark.parametrize("name", ["cube.msh", "cube.vtu", "cube.stl"])
def test_read_compressed(cube, name, suffix):
    tmp, raw = cube
    path = tmp / (name + suffix)
    path.write_bytes(COMPRESSORS[suffix](raw[name]))
    plain = mp.read(
        str(tmp / name), file_format="gmsh" if name.endswith("msh") else None
    )
    _same(
        mp.read(str(path), file_format="gmsh" if name.endswith("msh") else None), plain
    )
    if not name.endswith("msh"):
        _same(mp.read(str(path)), plain)  # format inferred from the inner name


@pytest.mark.parametrize("member", ["cube.vtu", "inner/dir/cube.vtu"])
@pytest.mark.parametrize("archive", ["case.zip", "cube.vtu.zip"])
@pytest.mark.parametrize("method", [zipfile.ZIP_DEFLATED, zipfile.ZIP_STORED])
def test_read_single_member_zip(cube, member, archive, method):
    tmp, raw = cube
    path = tmp / archive
    path.write_bytes(_zip(raw["cube.vtu"], member, method))
    _same(mp.read(str(path)), mp.read(str(tmp / "cube.vtu")))


def test_zip_with_several_files_is_refused(cube):
    tmp, raw = cube
    path = tmp / "two.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("a.vtu", raw["cube.vtu"])
        z.writestr("b.vtu", raw["cube.vtu"])
    with pytest.raises(mp.ReadError, match="holds 2 files"):
        mp.read(str(path))


def test_corrupt_data_raises_read_error(cube):
    tmp, raw = cube
    path = tmp / "cube.vtu.gz"
    path.write_bytes(gzip.compress(raw["cube.vtu"])[:40])
    with pytest.raises(mp.ReadError, match="corrupt or truncated"):
        mp.read(str(path))


def test_multi_file_formats_are_refused_by_name(tmp_path):
    path = tmp_path / "case.pvtu.gz"
    path.write_bytes(gzip.compress(b"<VTKFile/>"))
    with pytest.raises(mp.ReadError, match="pvtu.*several files"):
        mp.read(str(path))
    with pytest.raises(mp.ReadError, match="tetgen.*several files"):
        mp.read(str(path), file_format="tetgen")


def test_compound_extensions_a_format_owns_are_not_stripped():
    assert _compressed.compression_of("mesh.vol.gz") == ""
    assert _compressed.compression_of("mesh.xda.gz") == ""
    assert _compressed.compression_of("deck.dato.gz") == ""
    assert _compressed.compression_of("cube.msh.gz") == ".gz"
    assert _compressed.compression_of("CUBE.VTU.BZ2") == ".bz2"
    assert _compressed.compression_of("cube.vtu") == ""


def test_no_temporary_copies_are_left(cube):
    tmp, raw = cube
    path = tmp / "cube.vtu.gz"
    path.write_bytes(gzip.compress(raw["cube.vtu"]))

    def copies():
        return {
            p
            for p in os.listdir(tempfile.gettempdir())
            if p.startswith("meshioplusplus-")
        }

    before = copies()
    mp.read(str(path))
    assert copies() == before


@pytest.mark.skipif(_core is None, reason="needs the C++ core")
@pytest.mark.parametrize("suffix", [".gz", ".zip"])
def test_native_registry_reads_them_too(cube, suffix):
    """The core's registry (what both CLIs, the C API and the flat bindings
    read through) decompresses the same files."""
    tmp, raw = cube
    path = tmp / ("cube.vtu" + suffix)
    data = raw["cube.vtu"]
    path.write_bytes(gzip.compress(data) if suffix == ".gz" else _zip(data, "cube.vtu"))
    meta = _core.read_metadata(str(path))
    assert meta["format"] == "vtu"
    assert meta["num_points"] == len(helpers.tet_mesh.points)
    assert _core.decompress(".gz", gzip.compress(data)) == data
    with pytest.raises(ValueError, match="unknown suffix"):
        _core.decompress(".rar", data)


@pytest.mark.skipif(_core is None, reason="needs the C++ core")
def test_python_cli_converts_a_compressed_file(cube, capsys):
    tmp, raw = cube
    src = tmp / "cube.vtu.xz"
    src.write_bytes(lzma.compress(raw["cube.vtu"]))
    out = tmp / "out.vtu"
    assert mp._cli.main(["convert", str(src), str(out)]) in (0, None)
    _same(mp.read(str(out)), mp.read(str(tmp / "cube.vtu")))


def _native_cli():
    from .test_cli_json import NATIVE

    return NATIVE


def test_native_cli_converts_a_compressed_file(cube):
    native = _native_cli()
    if native is None:
        pytest.skip("no native CLI binary built")
    tmp, raw = cube
    src = tmp / "cube.vtu.gz"
    src.write_bytes(gzip.compress(raw["cube.vtu"]))
    out = tmp / "native.vtu"
    result = subprocess.run(
        [native, "convert", str(src), str(out)], capture_output=True
    )
    assert result.returncode == 0, result.stderr
    _same(mp.read(str(out)), mp.read(str(tmp / "cube.vtu")))
