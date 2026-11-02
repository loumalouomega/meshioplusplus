"""Reading meshes from URLs (roadmap 7.1.2).

``read("https://example.org/cube.vtu")``, ``s3://bucket/run/step_010.vtu``,
``gs://...`` and every other protocol `fsspec <https://filesystem-spec.readthedocs.io>`_
knows: the file is fetched into a private temporary directory under its own
name and read as a local path, so format inference, the compressed-file
handling (``cube.msh.gz``) and every reader work unchanged. ``fsspec`` is an
optional dependency (the ``remote`` extra); the object stores need their own
filesystem package (``s3fs``, ``gcsfs``, ``adlfs``). A chained URL such as
``simplecache::https://...`` keeps a local copy between reads.

Multi-file formats (one whose data lives in files beside the one named) are
refused by name, as they are for buffers. ``file://`` URLs are read too, which
is mostly useful in tests.
"""

from __future__ import annotations

import os
import re
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import PurePosixPath
from urllib.parse import unquote, urlparse

from ._exceptions import ReadError

# A URL scheme: a letter, then letters, digits, `+`, `-` or `.`, then `://`.
# Needs two characters at least, so a Windows drive (`C:\`, `C://x`) never
# matches; chained fsspec URLs (`simplecache::https://`) start with one too.
_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]+://")


def is_url(source) -> bool:
    """Whether ``source`` is a URL rather than a local path."""
    if isinstance(source, os.PathLike):
        return False
    return isinstance(source, str) and bool(_SCHEME.match(source.split("::")[-1]))


def _fsspec():
    try:
        import fsspec
    except ImportError:
        raise ReadError(
            "meshio++: reading a URL needs the optional 'fsspec' package "
            "(pip install fsspec, or meshioplusplus[remote])"
        ) from None
    return fsspec


def url_name(url: str) -> str:
    """The file name a URL's path ends in (query and fragment dropped)."""
    path = urlparse(url.split("::")[-1]).path
    name = PurePosixPath(unquote(path)).name
    return name or "download"


def refuse_multi_file(url: str, file_format) -> None:
    """Refuse a multi-file format by name: a fetched copy has no siblings."""
    from pathlib import Path

    from ._compressed import compression_of, decompressed_name
    from ._helpers import _MULTIFILE_FORMATS, _filetypes_from_path

    def multi(fmt):
        return fmt in _MULTIFILE_FORMATS or fmt == "openfoam"

    if file_format:
        candidates = [file_format]
    else:
        name = Path(url_name(url))
        suffix = compression_of(name)
        if suffix:
            try:
                name = decompressed_name(name, suffix)
            except OSError:
                pass  # a zip named by its member; inferred after the fetch
        try:
            candidates = _filetypes_from_path(name)
        except ReadError:
            candidates = []
    if candidates and all(multi(f) for f in candidates):
        raise ReadError(
            f"meshio++: '{url}' is a {candidates[0]} file, whose data spans several "
            "files; a URL is fetched as one file, so download the whole case first"
        )


@contextmanager
def fetched(url: str, file_format=None):
    """The file at ``url`` downloaded into a temporary directory under its own
    name; removed on exit."""
    fsspec = _fsspec()
    refuse_multi_file(url, file_format)
    with tempfile.TemporaryDirectory(prefix="meshioplusplus-") as tmp:
        local = os.path.join(tmp, url_name(url))
        try:
            with fsspec.open(url, "rb") as src, open(local, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
        except FileNotFoundError as exc:
            raise ReadError(f"meshio++: '{url}' not found") from exc
        except (ImportError, ValueError) as exc:
            # fsspec names the filesystem package a protocol needs (s3fs...),
            # or says it does not know the protocol at all.
            raise ReadError(f"meshio++: cannot read '{url}': {exc}") from exc
        yield local


def exists(url: str) -> bool:
    """Whether ``url`` names an existing file."""
    fs, path = _fsspec().core.url_to_fs(url)
    return fs.exists(path)


def glob(pattern: str) -> list:
    """The URLs a ``*``/``?`` pattern matches, as full URLs, sorted."""
    fs, path = _fsspec().core.url_to_fs(pattern)
    return sorted(fs.unstrip_protocol(p) for p in fs.glob(path))
