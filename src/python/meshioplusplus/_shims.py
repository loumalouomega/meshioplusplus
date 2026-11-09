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

import sys

from ._fallback import core_declined
from ._files import is_buffer
from ._vtk_codecs import CPP_CODECS, cpp_codec_ok


def _named(fn, module, name, doc):
    fn.__module__ = module
    fn.__name__ = fn.__qualname__ = name
    fn.__doc__ = doc
    return fn


def plain_shims(module, fmt, read_doc, write_doc, read_mode="r", write_mode="w"):
    """``read(filename)`` and ``write(filename, mesh)`` for package ``module``.

    Each tries ``_core.<fmt>_read`` / ``_core.<fmt>_write`` on a real path (not
    a buffer opened with ``read_mode`` / ``write_mode``) and falls back to
    ``_py_read`` / ``_py_write`` when :func:`core_declined` says so.
    """

    def read(filename):
        pkg = sys.modules[module]
        if not is_buffer(filename, read_mode):
            try:
                return getattr(pkg._core, f"{fmt}_read")(str(filename))
            except Exception as exc:
                if not core_declined(exc, fmt, "read", filename):
                    raise
        return pkg._py_read(filename)

    def write(filename, mesh):
        pkg = sys.modules[module]
        if not is_buffer(filename, write_mode):
            try:
                getattr(pkg._core, f"{fmt}_write")(str(filename), mesh)
                return
            except Exception as exc:
                if not core_declined(exc, fmt, "write", filename):
                    raise
        return pkg._py_write(filename, mesh)

    return (
        _named(read, module, "read", read_doc),
        _named(write, module, "write", write_doc),
    )


def vtk_xml_shims(module, fmt, read_doc, write_doc):
    """The VTK XML family's ``read(filename, points_only, arrays)`` and
    ``write(filename, mesh, binary, compression, header_type)``.

    The core reads with the narrowing options and writes through
    ``_core.<fmt>_write_codec`` when :func:`cpp_codec_ok` accepts the codec and
    header type; the Python twin takes the rest.
    """

    def read(filename, points_only=False, arrays=None):
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

    def write(filename, mesh, binary=True, compression="zlib", header_type=None):
        pkg = sys.modules[module]
        # zstd/lz4 reach the C++ writer when this build has them; otherwise the
        # Python reference path handles it (needing the `codecs` extra), so a
        # request never silently degrades to a different codec.
        if cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
            try:
                getattr(pkg._core, f"{fmt}_write_codec")(
                    str(filename), mesh, binary, CPP_CODECS[compression]
                )
                return
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
