"""Content-based mesh-format detection (conservative magic-byte sniffing).

:func:`sniff_format` reads the leading bytes of a file and returns a meshio++
format name on a confident signature match, or ``""`` otherwise. It is used as
a fallback by :func:`meshioplusplus.read` when the format cannot be inferred
from the extension. Mirrors ``src/cpp/src/operations/sniff.cpp``.
"""

from __future__ import annotations

from pathlib import Path


def _sniff_format_py(path) -> str:
    try:
        with open(path, "rb") as f:
            head = f.read(512)
    except OSError:
        return ""
    if not head:
        return ""
    stripped = head.lstrip()

    if b"VTKFile" in head:
        if b"UnstructuredGrid" in head:
            return "vtu"
        if b"PolyData" in head:
            return "vtp"
        # Checked last of the three, matching the C++ twin -- see sniff.cpp.
        if b"ImageData" in head:
            return "vti"
    if b"<Xdmf" in head:
        return "xdmf"
    if stripped.startswith(b"# vtk DataFile"):
        return "vtk"
    if stripped.startswith(b"$MeshFormat"):
        return "gmsh"
    # GiD postprocess -- see the C++ twin for why the geometry match includes
    # the opening quote and why .post.bin/.post.h5 are deliberately absent.
    if stripped.startswith(b"GiD Post Results File"):
        return "gid"
    if stripped.startswith(b'MESH "'):
        return "gid"
    if stripped.startswith((b"ply\n", b"ply\r")) or stripped == b"ply":
        return "ply"
    if stripped.startswith((b"OFF", b"COFF", b"NOFF", b"STOFF")):
        return "off"
    if stripped.startswith(b"solid "):
        return "stl"
    upper = stripped[:8].upper()
    if upper.startswith((b"*HEADING", b"*NODE")):
        return "abaqus"
    return ""


def sniff_format(path) -> str:
    """Guess a mesh file's format from its contents (``""`` if unsure)."""
    try:
        from . import _core

        return _core.sniff_format(str(path))
    except Exception:
        pass
    return _sniff_format_py(Path(path))
