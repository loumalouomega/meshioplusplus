"""Reading compressed mesh files transparently (roadmap 7.1.1).

``read("cube.msh.gz")`` strips the compression suffix to infer the format,
decompresses the file into a private temporary directory under the name it had
inside (``cube.msh``), reads that and removes the directory -- the Python twin
of the core's ``detail/compressed_input.hpp``, which the native registry (both
CLIs, the C API and every flat binding) uses for the same files.

``.gz``, ``.bz2`` and ``.xz`` decompress through the standard library, a
single-member ``.zip`` through :mod:`zipfile`, and ``.zst`` through
``compression.zstd`` (Python 3.14), the optional ``zstandard`` package or the C++
core, in that order. A compound extension a format registers itself
(Netgen's ``.vol.gz``, libMesh's ``.xda.gz``, PERMAS's ``.dato.gz``) is that
format's own business and is never stripped, and a multi-file format (one
whose data lives in several files) is refused by name.
"""

from __future__ import annotations

import io
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path

from ._exceptions import ReadError

SUFFIXES = (".gz", ".bz2", ".xz", ".zst", ".zip")


def compression_of(path) -> str:
    """The compression suffix of ``path`` (``".gz"``...), or ``""`` when it has
    none or a format registers the whole compound extension itself."""
    from ._helpers import _ensure_extension, extension_to_filetypes

    name = Path(path).name.lower()
    suffix = Path(name).suffix
    if suffix not in SUFFIXES:
        return ""
    parts = name.split(".")
    for k in range(1, len(parts) - 1):
        compound = "." + ".".join(parts[k:])
        _ensure_extension(compound)
        if dict.__contains__(extension_to_filetypes, compound):
            return ""
    return suffix


def _zip_member(data: bytes):
    import zipfile

    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise ReadError(f"meshio++: zip: {exc}") from exc
    files = [i for i in archive.infolist() if not i.is_dir()]
    if len(files) != 1:
        raise ReadError(
            f"meshio++: zip: the archive holds {len(files)} files; only a single-file "
            "archive is read as a mesh"
        )
    return archive, files[0]


def _unzstd(data: bytes) -> bytes:
    try:
        from compression import zstd  # Python 3.14

        return zstd.decompress(data)
    except ImportError:
        pass
    try:
        import zstandard

        out = io.BytesIO()
        zstandard.ZstdDecompressor().copy_stream(io.BytesIO(data), out)
        return out.getvalue()
    except ImportError:
        pass
    try:
        from . import _core

        if _core.compression_available(".zst"):
            return _core.decompress(".zst", data)
    except ImportError:
        pass
    raise ReadError(
        "meshio++: reading a .zst file needs Python 3.14, the 'zstandard' package or a "
        "core built with MESHIOPLUSPLUS_WITH_ZSTD"
    )


def decompress(suffix: str, data: bytes) -> bytes:
    """The bytes ``data`` decompress to under ``suffix``; a zip's one member."""
    import bz2
    import gzip
    import lzma
    import zlib

    try:
        if suffix == ".gz":
            return gzip.decompress(data)
        if suffix == ".bz2":
            return bz2.decompress(data)
        if suffix == ".xz":
            return lzma.decompress(data)
        if suffix == ".zst":
            return _unzstd(data)
        if suffix == ".zip":
            archive, member = _zip_member(data)
            return archive.read(member)
    except (OSError, EOFError, ValueError, lzma.LZMAError, zlib.error) as exc:
        raise ReadError(
            f"meshio++: {suffix}: corrupt or truncated data ({exc})"
        ) from exc
    raise ValueError(f"meshio++: unknown compression suffix '{suffix}'")


def decompressed_name(path, suffix: str) -> Path:
    """The name ``path``'s content had before compression; for ``case.zip``
    (no inner extension) the member's own name beside it."""
    path = Path(path)
    stripped = path.with_name(path.name[: -len(suffix)])
    if suffix == ".zip" and not stripped.suffix:
        _, member = _zip_member(path.read_bytes())
        stripped = path.with_name(Path(member.filename).name)
    return stripped


@contextmanager
def decompressed(path, suffix: str):
    """A decompressed copy of ``path`` in a temporary directory, under the name
    it had before compression; removed on exit."""
    data = decompress(suffix, Path(path).read_bytes())
    name = decompressed_name(path, suffix).name or "decompressed"
    with tempfile.TemporaryDirectory(prefix="meshioplusplus-") as tmp:
        inner = Path(tmp) / name
        inner.write_bytes(data)
        yield inner


def refuse_multi_file(file_format: str, path) -> None:
    from ._helpers import _MULTIFILE_FORMATS

    if file_format in _MULTIFILE_FORMATS or file_format == "openfoam":
        raise ReadError(
            f"meshio++: '{os.fspath(path)}' is a compressed {file_format} file, whose data "
            "spans several files; decompress it first"
        )
