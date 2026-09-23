"""Content-based mesh-format detection (conservative magic-byte sniffing).

:func:`sniff_format` reads the leading bytes of a file and returns a meshio++
format name on a confident signature match, or ``""`` otherwise. It is used as
a fallback by :func:`meshioplusplus.read` when the format cannot be inferred
from the extension. Mirrors ``src/cpp/src/operations/sniff.cpp``. A directory
is sniffed by the files it holds (Elmer, OpenFOAM).
"""

from __future__ import annotations

import re
import struct
from pathlib import Path

# Dataset numbers that open an I-DEAS universal file (see sniff.cpp's kUnvIds).
_UNV_IDS = {
    "15",
    "18",
    "55",
    "56",
    "57",
    "58",
    "82",
    "151",
    "164",
    "780",
    "781",
    "2400",
    "2411",
    "2412",
    "2414",
    "2417",
    "2420",
    "2429",
    "2430",
    "2432",
    "2435",
    "2452",
    "2467",
    "2477",
}


def _is_mphbin(head: bytes) -> bool:
    """COMSOL binary: int32 0 1 (version), a tag count, then the first tag's
    length and first character, one int32 each; see sniff.cpp."""
    if len(head) < 20:
        return False
    v = struct.unpack_from("<5i", head)
    return (
        v[0] == 0
        and v[1] == 1
        and 1 <= v[2] <= 4096
        and 1 <= v[3] <= 1024
        and 32 < v[4] < 127
    )


def _is_mphtxt(head: bytes) -> bool:
    """COMSOL text: the same values as tokens, ``#`` comments skipped."""
    text = re.sub(rb"#[^\n]*", b" ", head)
    tokens = text.split()[:5]

    def count(t):
        return t.isdigit() and len(t) < 6 and t != b"0"

    return (
        len(tokens) == 5
        and tokens[0] == b"0"
        and tokens[1] == b"1"
        and count(tokens[2])
        and count(tokens[3])
        and tokens[4][:1].isalpha()
    )


def _has_polymesh(poly: Path) -> bool:
    return (poly / "owner").is_file() and (poly / "faces").is_file()


def _sniff_directory(path: Path) -> str:
    """An Elmer mesh directory or an OpenFOAM case, by the files the readers
    look for; a directory that looks like both, or like neither, is ``""``."""

    def is_partitioning(p):
        return p.name.startswith("partitioning.") and (p / "part.1.header").is_file()

    elmer = (path / "mesh.header").is_file() or is_partitioning(path)
    if not elmer:
        try:
            elmer = any(is_partitioning(child) for child in path.iterdir())
        except OSError:
            pass
    openfoam = (
        (path.name == "polyMesh" and _has_polymesh(path))
        or _has_polymesh(path / "constant" / "polyMesh")
        or _has_polymesh(path / "polyMesh")
        or (path / "processor0" / "constant" / "polyMesh").is_dir()
        or (path / "constant" / "regionProperties").is_file()
    )
    if elmer == openfoam:
        return ""
    return "elmer" if elmer else "openfoam"


def _sniff_format_py(path) -> str:
    path = Path(path)
    if path.is_dir():
        return _sniff_directory(path)
    # The header of an Elmer mesh directory stands for the directory.
    if path.name == "mesh.header" and path.is_file():
        return "elmer"
    try:
        with open(path, "rb") as f:
            head = f.read(512)
    except OSError:
        return ""
    if not head:
        return ""
    stripped = head.lstrip()
    if _is_mphbin(head):
        return "mphbin"
    # FEBio plot file: the magic 0x00464542, in either byte order.
    if head[:4] in (b"BEF\x00", b"\x00FEB"):
        return "xplt"
    # Ansys MAPDL results: a 100-word integer record (length 100, flags
    # 0x80000000) whose first value is the file number, 12.
    if head[:12] == b"d\x00\x00\x00\x00\x00\x00\x80\x0c\x00\x00\x00":
        return "ansys_rst"
    # FEBio input: XML whose root is <febio_spec>.
    if b"<febio_spec" in head:
        return "febio"

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
    if _is_mphtxt(head):
        return "mphtxt"
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
    # I-DEAS universal file: a lone "-1" line, then a known dataset number (a trailing
    # "b" marks a binary dataset such as 58b); see the C++ twin.
    if stripped.startswith(b"-1"):
        rest = stripped[2:].lstrip(b" \r")
        if rest.startswith(b"\n"):
            first = rest[1:].split(None, 1)[:1]
            if first and first[0].rstrip(b"bB").decode("latin-1") in _UNV_IDS:
                return "unv"
    if stripped.startswith(b"solid "):
        return "stl"
    # Code_Aster .mail meshes open, after any `%` comment lines, with a TITRE or
    # COOR_1D/2D/3D block keyword; see the C++ twin.
    for line in stripped.split(b"\n"):
        line = line.lstrip(b" \t\r")
        if not line or line[:1] == b"%":
            continue
        word = re.split(rb"[ \t\r,%]", line.upper(), maxsplit=1)[0]
        if word in (b"TITRE", b"COOR_1D", b"COOR_2D", b"COOR_3D"):
            return "code_aster"
        break
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
