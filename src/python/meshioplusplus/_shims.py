"""The read and write shims many format packages share (roadmap §7).

A format package whose shims only try one ``_core`` entry point and fall back
to its Python twin builds them here::

    from .. import _core
    from .._helpers import register_format
    from .._shims import plain_shims
    from ._wkt import read as _py_read
    from ._wkt import write as _py_write

    read, write = plain_shims(__name__, "wkt", "Read ...", "Write ...")
    register_format("wkt", [".wkt"], read, {"wkt": write})

The shims look ``_core``, ``_py_read`` and ``_py_write`` up on the package
when they run, not when they are built, so a test that patches those names on
the package reaches them. ``register_format`` stays in each package:
``tools/gen_format_table.py`` attributes a format to the module that calls it.
Each shim keeps the signature and docstring its package always had.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from typing import Any, TypeVar

from ._fallback import core_declined
from ._files import is_buffer
from ._vtk_codecs import CPP_CODECS, cpp_codec_ok

_F = TypeVar("_F", bound=Callable[..., Any])


def _named(fn: _F, module: str, name: str, doc: str) -> _F:
    fn.__module__ = module
    fn.__name__ = fn.__qualname__ = name
    fn.__doc__ = doc
    return fn


def plain_shims(
    module: str,
    fmt: str,
    read_doc: str,
    write_doc: str,
    read_mode: str = "r",
    write_mode: str = "w",
) -> tuple[Callable[..., Any], Callable[..., Any]]:
    """``read(filename)`` and ``write(filename, mesh)`` for package ``module``.

    Each tries ``_core.<fmt>_read`` / ``_core.<fmt>_write`` on a real path (not
    a buffer opened with ``read_mode`` / ``write_mode``) and falls back to
    ``_py_read`` / ``_py_write`` when :func:`core_declined` says so.
    """

    def read(filename: Any) -> Any:
        pkg = sys.modules[module]
        if not is_buffer(filename, read_mode):
            try:
                return getattr(pkg._core, f"{fmt}_read")(str(filename))
            except Exception as exc:
                if not core_declined(exc, fmt, "read", filename):
                    raise
        return pkg._py_read(filename)

    def write(filename: Any, mesh: Any) -> Any:
        pkg = sys.modules[module]
        if not is_buffer(filename, write_mode):
            try:
                getattr(pkg._core, f"{fmt}_write")(str(filename), mesh)
                return None
            except Exception as exc:
                if not core_declined(exc, fmt, "write", filename):
                    raise
        return pkg._py_write(filename, mesh)

    return (
        _named(read, module, "read", read_doc),
        _named(write, module, "write", write_doc),
    )


def time_values_shim(
    module: str, fmt: str, doc: str, read_mode: str = "r"
) -> Callable[[Any], list[float]]:
    """``time_values(filename)`` for package ``module``.

    Guards the native ``_core.<fmt>_time_values`` call with :func:`is_buffer`
    so file-like objects go directly to ``_py_time_values``, and falls back
    when :func:`core_declined` approves.
    """

    def time_values(filename: Any) -> list[float]:
        pkg = sys.modules[module]
        if not is_buffer(filename, read_mode):
            try:
                return list(getattr(pkg._core, f"{fmt}_time_values")(str(filename)))
            except Exception as exc:
                if not core_declined(exc, fmt, "read", filename):
                    raise
        return list(pkg._py_time_values(filename))

    return _named(time_values, module, "time_values", doc)


def vtk_xml_shims(
    module: str, fmt: str, read_doc: str, write_doc: str
) -> tuple[Callable[..., Any], Callable[..., Any]]:
    """The VTK XML family's ``read(filename, points_only, arrays)`` and
    ``write(filename, mesh, binary, compression, header_type)``.

    The core reads with the narrowing options and writes through
    ``_core.<fmt>_write_codec`` when :func:`cpp_codec_ok` accepts the codec and
    header type; the Python twin takes the rest.
    """

    def read(filename: Any, points_only: bool = False, arrays: Any = None) -> Any:
        pkg = sys.modules[module]
        # points_only/arrays reach the C++ reader, which skips the unwanted
        # <DataArray>/section bodies outright. The Python fallback below has no
        # selective support, so _helpers.read trims its result instead -- same
        # answer, just without the saving.
        if not is_buffer(filename, "r"):
            try:
                return getattr(pkg._core, f"{fmt}_read")(
                    str(filename), points_only=points_only, arrays=arrays
                )
            except Exception as exc:
                if not core_declined(exc, fmt, "read", filename):
                    raise
        return pkg._py_read(filename)

    def write(
        filename: Any,
        mesh: Any,
        binary: bool = True,
        compression: str | None = "zlib",
        header_type: str | None = None,
    ) -> Any:
        pkg = sys.modules[module]
        # zstd/lz4 reach the C++ writer when this build has them; otherwise the
        # Python reference path handles it (needing the `codecs` extra), so a
        # request never silently degrades to a different codec.
        if cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
            try:
                getattr(pkg._core, f"{fmt}_write_codec")(
                    str(filename), mesh, binary, CPP_CODECS[compression]
                )
                return None
            except Exception as exc:
                if not core_declined(exc, fmt, "write", filename):
                    raise
        return pkg._py_write(
            filename,
            mesh,
            binary=binary,
            compression=compression,
            header_type=header_type,
        )

    return (
        _named(read, module, "read", read_doc),
        _named(write, module, "write", write_doc),
    )
