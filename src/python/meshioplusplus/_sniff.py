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
        # The parallel indices and the collection come first, and match the
        # quoted `type=` value: `PUnstructuredGrid` contains `UnstructuredGrid`
        # and `PPolyData` contains `PolyData`, so the loose substring checks
        # below would call an index a piece; and a bare `Collection` would also
        # match `vtkPartitionedDataSetCollection`. A parallel image, structured
        # or rectilinear index is refused outright rather than mistaken for its
        # serial twin.
        for quote in (b'"', b"'"):
            for value, name in (
                (b"Collection", "pvd"),
                (b"PUnstructuredGrid", "pvtu"),
                (b"PPolyData", "pvtp"),
                (b"PImageData", ""),
                (b"PStructuredGrid", ""),
                (b"PRectilinearGrid", ""),
            ):
                if b"type=" + quote + value + quote in head:
                    return name
        if b"UnstructuredGrid" in head:
            return "vtu"
        if b"PolyData" in head:
            return "vtp"
        # Checked last of the four, matching the C++ twin -- see sniff.cpp.
        if b"ImageData" in head:
            return "vti"
        if b"StructuredGrid" in head:  # v11.6.0, roadmap §1 tier B4
            return "vts"
        if b"RectilinearGrid" in head:
            return "vtr"
        if b"MultiBlockDataSet" in head:
            return "vtm"
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
    if stripped.startswith(b"# .PCD") or (
        stripped.startswith(b"VERSION") and b"\nFIELDS" in head
    ):
        return "pcd"
    # CalculiX results: a lone "    1C" record, then the "1U" user header (or the "2C"
    # node block of a file without one); see the C++ twin.
    if head.startswith(b"    1C"):
        rest = head[6:].lstrip(b" \r")
        if rest.startswith(b"\n") and rest[1:7] in (b"    1U", b"    2C"):
            return "frd"
    if stripped.startswith(b"solid "):
        return "stl"
    # LS-DYNA decks open with "*KEYWORD" after any `$` comment lines; checked before
    # the Abaqus rule, as in the C++ twin.
    for line in stripped.split(b"\n"):
        if not line or line[:1] in (b"$", b"\r"):
            continue
        if line[:8].upper().startswith(b"*KEYWORD"):
            return "lsdyna"
        break
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
