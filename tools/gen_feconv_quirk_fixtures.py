#!/usr/bin/env python3
"""Regenerate the fixtures that reproduce the quirks of FEconv's samples.

FEconv's ``examples/`` (https://github.com/victorsndvg/FEconv) exposed reader
bugs in several formats. Those files are GPL-3 and are not committed; these
small files reproduce each quirk instead. They are written byte by byte here,
not by meshio++, so a reader is never only checked against its own writer.

* ``flux/solids.pf3`` -- one ``tetra``, ``tetra10``, ``pyramid``, ``wedge``,
  ``wedge15``, ``hexahedron`` and ``hexahedron20`` in FLUX's own layout: the
  corners of each solid listed with the base face clockwise seen from inside
  (VTK's order of the mirrored element), then the mid-edge nodes on FLUX's edge
  list (VTK's for the mirrored corners; ``tetra10``'s own (0,1) (0,2) (0,3)
  (1,2) (2,3) (1,3)). Read correctly, every solid is positively oriented with
  its mid-edge nodes on edge midpoints.
* ``flux/truncated.pf3`` -- a header declaring more elements and points than the
  file holds, with sparse node ids (as FLUX excerpts are).
* ``gmsh/indented.msh`` -- a Gmsh 2.2 file with every line indented and one
  element line past the declared count.
* ``medit/no_dimension.mesh`` -- a Medit file without the ``Dimension`` keyword
  and with an empty ``Triangles`` block.
* ``vtu/raw_bigendian.vtu`` -- raw appended data, BigEndian, offsets padded
  with blanks; ``vtu/raw_zlib.vtu`` -- raw appended zlib blocks, UInt64
  headers; ``vtu/base64_appended.vtu`` -- base64 appended data with each
  header encoded apart from its body; ``vtu/binary_bigendian.vtu`` -- inline
  base64, BigEndian; ``vtu/two_pieces_polyhedron.vtu`` -- two pieces, each a
  pyramid and a polyhedron, with cell data in both pieces and point data in
  only one. All five hold the same geometry.
* ``ansys/cells3d.msh`` -- an ASCII Fluent mesh whose cells exist only through
  their faces: two tetrahedra sharing an interior face, a hexahedron, a wedge,
  a pyramid and an 8-node polyhedron (a cube with one side split in two), in a
  mixed cell zone with its type list; wall faces in a mixed face zone; zone
  names from ``(39 ...)`` and ``(45 ...)``.
* ``ansys/tgrid2d.msh`` -- a TGrid-style 2-D mesh (a triangle and a quad):
  node zones out of order, a blank before each body, a mixed face zone.
* ``ansys/gambit2d.msh`` -- a GAMBIT-style 2-D mesh: ``(13(`` with no blank,
  boundary faces listing their cell as ``c1`` with ``c0 = 0``.
* ``ansys/binary3d.msh`` -- binary sections: float64 nodes, an int32 interior
  face zone, an int64 mixed wall zone, and a cell declaration followed by
  ``End of Binary Section`` instead of a body.

    python tools/gen_feconv_quirk_fixtures.py
"""

import base64
import pathlib
import struct
import zlib

MESHES = pathlib.Path(__file__).resolve().parent.parent / "tests" / "python" / "meshes"

# ---------------------------------------------------------------------------
# FLUX .pf3
# ---------------------------------------------------------------------------

# Corners with the base face clockwise seen from inside the element.
_FLUX_CORNERS = {
    "tetra": [(0, 0, 0), (0, 1, 0), (1, 0, 0), (0, 0, 1)],
    "pyramid": [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0), (0.5, 0.5, 1)],
    "wedge": [(0, 0, 0), (0, 1, 0), (1, 0, 0), (0, 0, 1), (0, 1, 1), (1, 0, 1)],
    "hexahedron": [
        (0, 0, 0),
        (0, 1, 0),
        (1, 1, 0),
        (1, 0, 0),
        (0, 0, 1),
        (0, 1, 1),
        (1, 1, 1),
        (1, 0, 1),
    ],
}
_FLUX_CORNERS["tetra10"] = _FLUX_CORNERS["tetra"]
_FLUX_CORNERS["wedge15"] = _FLUX_CORNERS["wedge"]
_FLUX_CORNERS["hexahedron20"] = _FLUX_CORNERS["hexahedron"]

# Mid-edge nodes, as pairs of the corners above, in file order.
_FLUX_EDGES = {
    "tetra10": [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (1, 3)],
    "wedge15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "hexahedron20": [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ],
}

# meshio++ type -> (desc1, desc2, desc3), FLUX's element descriptors.
_FLUX_DESC = {
    "tetra": (5, 4, 10),
    "tetra10": (5, 15, 11),
    "wedge": (6, 207, 12),
    "wedge15": (6, 307, 13),
    "hexahedron": (7, 2202, 15),
    "hexahedron20": (7, 3303, 16),
    "pyramid": (8, 4202, 17),
    "quad8": (4, 303, 8),
}


def _flux_header(nel, nvol, nsur, nnod):
    rows = [
        (3, "NOMBRE DE DIMENSIONS DU DECOUPAGE"),
        (nel, "NOMBRE  D'ELEMENTS"),
        (nvol, "NOMBRE  D'ELEMENTS VOLUMIQUES"),
        (nsur, "NOMBRE  D'ELEMENTS SURFACIQUES"),
        (0, "NOMBRE  D'ELEMENTS LINEIQUES"),
        (0, "NOMBRE  D'ELEMENTS PONCTUELS"),
        (0, "NOMBRE DE MACRO-ELEMENTS"),
        (nnod, "NOMBRE DE POINTS"),
        (1, "NOMBRE DE REGIONS"),
        (1, "NOMBRE DE REGIONS VOLUMIQUES"),
        (0, "NOMBRE DE REGIONS SURFACIQUES"),
        (0, "NOMBRE DE REGIONS LINEIQUES"),
        (0, "NOMBRE DE REGIONS PONCTUELLES"),
        (0, "NOMBRE DE REGIONS MACRO-ELEMENTAIRES"),
        (20, "NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)"),
        (20, "NOMBRE DE POINTS D'INTEGRATION / ELEMENT (MAX)"),
    ]
    out = " Fichier cree par gen_feconv_quirk_fixtures.py\n"
    out += "".join(f"{v:8d}           {label}\n" for v, label in rows)
    # Latin-1 region name, as FLUX writes them.
    out += " NOMS DES REGIONS\n     REGIONS VOLUMIQUES\n SOLIDE_\xc9L\xc9MENTS\n"
    return out


def _flux_element(eid, cell_type, ref, nodes):
    d1, d2, d3 = _FLUX_DESC[cell_type]
    n = len(nodes)
    head = f"{eid:13d}{d1:7d}{d2:7d}{ref:7d}{n:7d}{0:7d}{d3:7d}{n:7d}" + "      0" * 4
    return head + "\n" + "".join(f"{v:7d}" for v in nodes) + "\n"


def write_flux():
    points = []
    elements = []
    for k, cell_type in enumerate(
        [
            "tetra",
            "tetra10",
            "pyramid",
            "wedge",
            "wedge15",
            "hexahedron",
            "hexahedron20",
        ]
    ):
        corners = [(x + 2.0 * k, y, z) for x, y, z in _FLUX_CORNERS[cell_type]]
        mids = [
            tuple((a + b) / 2 for a, b in zip(corners[i], corners[j]))
            for i, j in _FLUX_EDGES.get(cell_type, [])
        ]
        first = len(points) + 1
        points += corners + mids
        elements.append((cell_type, list(range(first, len(points) + 1))))

    text = _flux_header(len(elements), len(elements), 0, len(points))
    text += " DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n"
    for eid, (cell_type, nodes) in enumerate(elements, start=1):
        text += _flux_element(eid, cell_type, 1, nodes)
    text += " COORDONNEES DES NOEUDS\n"
    for i, p in enumerate(points, start=1):
        text += f"{i:8d}" + "".join(f"  {v:.7E}" for v in p) + "\n"
    text += " ==== DECOUPAGE  TERMINE\n"
    (MESHES / "flux").mkdir(parents=True, exist_ok=True)
    (MESHES / "flux" / "solids.pf3").write_bytes(text.encode("latin-1"))

    # An excerpt: 40 elements and 900 points declared, one quad8 on sparse ids.
    ids = [3, 121, 749, 112, 5521, 5525, 5329, 5522]
    coords = [
        (0, 1, 0),
        (1, 1, 0),
        (1, 0, 0),
        (0, 0, 0),
        (0.5, 1, 0),
        (1, 0.5, 0),
        (0.5, 0, 0),
        (0, 0.5, 0),
    ]
    text = _flux_header(40, 0, 40, 900)
    text += " DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n"
    text += _flux_element(3050, "quad8", 4, ids)
    text += " COORDONNEES DES NOEUDS\n"
    for i, p in zip(ids, coords):
        text += f"{i:8d}" + "".join(f"  {v:.7E}" for v in p) + "\n"
    text += " ==== DECOUPAGE  TERMINE\n"
    (MESHES / "flux" / "truncated.pf3").write_bytes(text.encode("latin-1"))


# ---------------------------------------------------------------------------
# Gmsh and Medit
# ---------------------------------------------------------------------------


def write_gmsh():
    lines = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$Nodes",
        "6",
        "1 0.0 0.0 0.0",
        "2 1.0 0.0 0.0",
        "3 1.0 1.0 0.0",
        "4 0.0 1.0 0.0",
        "5 2.0 0.0 0.0",
        "7 2.0 1.0 0.0",
        "$EndNodes",
        "$Elements",
        "2",
        "1 3 2 99 2 1 2 3 4",
        "2 3 2 99 2 2 5 7 3",
        "2 3 2 98 2 2 5 7 3",
        "$EndElements",
    ]
    (MESHES / "gmsh").mkdir(parents=True, exist_ok=True)
    text = "".join(f"     {line}   \n" for line in lines)
    (MESHES / "gmsh" / "indented.msh").write_text(text)


def write_medit():
    text = """MeshVersionFormatted 1

Vertices
5
0.000000000000000E+000 0.000000000000000E+000 0.000000000000000E+000 0
1.00000000000000 0.000000000000000E+000 0.000000000000000E+000 0
0.000000000000000E+000 1.00000000000000 0.000000000000000E+000 0
0.000000000000000E+000 0.000000000000000E+000 1.00000000000000 0
1.00000000000000 1.00000000000000 1.00000000000000 0

Tetrahedra
2
1 2 3 4 1
2 5 3 4 2

Triangles
0

End
"""
    (MESHES / "medit").mkdir(parents=True, exist_ok=True)
    (MESHES / "medit" / "no_dimension.mesh").write_text(text)


# ---------------------------------------------------------------------------
# VTU
# ---------------------------------------------------------------------------

# A unit cube as a polyhedron (points 0-7) and a pyramid on its top face
# (points 4-8), the layout of FEconv's polyhedron2pieces.vtu.
_VTU_POINTS = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
    (0.5, 0.5, 2),
]
_VTU_CONN = [4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7]
_VTU_OFFSETS = [5, 13]
_VTU_TYPES = [14, 42]
_VTU_FACES = [
    6,
    4,
    3,
    0,
    4,
    7,
    4,
    1,
    2,
    6,
    5,
    4,
    0,
    1,
    5,
    4,
    4,
    2,
    3,
    7,
    6,
    4,
    3,
    2,
    1,
    0,
]
_VTU_FACES += [4, 4, 5, 6, 7]
_VTU_FACEOFFSETS = [-1, 31]
_VTU_MATERIAL = [7, 9]
_VTU_HEAT = [float(i) for i in range(9)]

# (name, section, vtk type, struct code, values, components)
_VTU_ARRAYS = [
    ("Points", "Points", "Float64", "d", [v for p in _VTU_POINTS for v in p], 3),
    ("connectivity", "Cells", "Int32", "i", _VTU_CONN, 0),
    ("offsets", "Cells", "Int32", "i", _VTU_OFFSETS, 0),
    ("types", "Cells", "UInt8", "B", _VTU_TYPES, 0),
    ("faces", "Cells", "Int32", "i", _VTU_FACES, 0),
    ("faceoffsets", "Cells", "Int32", "i", _VTU_FACEOFFSETS, 0),
    ("material", "CellData", "Int32", "i", _VTU_MATERIAL, 0),
    ("heat", "PointData", "Float64", "d", _VTU_HEAT, 0),
]


def _vtu_document(byte_order, header_type, compressor, arrays_xml, appended):
    comp = f' compressor="{compressor}"' if compressor else ""
    head = (
        '<?xml version="1.0"?>\n'
        f'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="{byte_order}"'
        f' header_type="{header_type}"{comp}>\n'
        "  <UnstructuredGrid>\n"
        f'    <Piece NumberOfPoints="{len(_VTU_POINTS):10d}"'
        f' NumberOfCells="{len(_VTU_TYPES):10d}">\n'
    )
    body = ""
    for section in ("Points", "Cells", "PointData", "CellData"):
        body += f"      <{section}>\n"
        body += "".join(xml for sec, xml in arrays_xml if sec == section)
        body += f"      </{section}>\n"
    tail = "    </Piece>\n  </UnstructuredGrid>\n"
    if appended is not None:
        enc, payload = appended
        return (
            (head + body + tail + f'  <AppendedData encoding="{enc}">\n   _').encode()
            + payload
            + b"\n  </AppendedData>\n</VTKFile>\n"
        )
    return (head + body + tail + "</VTKFile>\n").encode()


def _da(name, vtk_type, nc, attrs):
    ncomp = f' NumberOfComponents="{nc}"' if nc else ""
    return f'        <DataArray type="{vtk_type}" Name="{name}"{ncomp} {attrs}>\n'


def write_vtu():
    out = MESHES / "vtu"
    out.mkdir(parents=True, exist_ok=True)

    def packed(code, values, endian):
        return struct.pack(f"{endian}{len(values)}{code}", *values)

    # Raw appended, BigEndian, uncompressed, offsets padded with blanks.
    payload = b""
    xml = []
    for name, section, vtk_type, code, values, nc in _VTU_ARRAYS:
        body = packed(code, values, ">")
        xml.append(
            (
                section,
                _da(name, vtk_type, nc, f'format="appended" offset="{len(payload):8d}"')
                + "        </DataArray>\n",
            )
        )
        payload += struct.pack(">I", len(body)) + body
    (out / "raw_bigendian.vtu").write_bytes(
        _vtu_document("BigEndian", "UInt32", None, xml, ("raw", payload))
    )

    # Raw appended, zlib, UInt64 headers, two blocks for the points.
    payload = b""
    xml = []
    for name, section, vtk_type, code, values, nc in _VTU_ARRAYS:
        body = packed(code, values, "<")
        block = 64
        chunks = [body[i : i + block] for i in range(0, len(body), block)] or [b""]
        comp = [zlib.compress(c) for c in chunks]
        last = len(chunks[-1]) if len(body) % block else block
        header = struct.pack(
            f"<{3 + len(comp)}Q", len(comp), block, last, *[len(c) for c in comp]
        )
        xml.append(
            (
                section,
                _da(name, vtk_type, nc, f'format="appended" offset="{len(payload)}"')
                + "        </DataArray>\n",
            )
        )
        payload += header + b"".join(comp)
    (out / "raw_zlib.vtu").write_bytes(
        _vtu_document(
            "LittleEndian", "UInt64", "vtkZLibDataCompressor", xml, ("raw", payload)
        )
    )

    # Base64 appended; each array's header encoded apart from its body.
    payload = ""
    xml = []
    for name, section, vtk_type, code, values, nc in _VTU_ARRAYS:
        body = packed(code, values, "<")
        xml.append(
            (
                section,
                _da(name, vtk_type, nc, f'format="appended" offset="{len(payload)}"')
                + "        </DataArray>\n",
            )
        )
        payload += base64.b64encode(struct.pack("<I", len(body))).decode()
        payload += base64.b64encode(body).decode()
    (out / "base64_appended.vtu").write_bytes(
        _vtu_document("LittleEndian", "UInt32", None, xml, ("base64", payload.encode()))
    )

    # Inline base64, BigEndian.
    xml = []
    for name, section, vtk_type, code, values, nc in _VTU_ARRAYS:
        body = packed(code, values, ">")
        text = base64.b64encode(struct.pack(">I", len(body)) + body).decode()
        xml.append(
            (
                section,
                _da(name, vtk_type, nc, 'format="binary"')
                + f"          {text}\n        </DataArray>\n",
            )
        )
    (out / "binary_bigendian.vtu").write_bytes(
        _vtu_document("BigEndian", "UInt32", None, xml, None)
    )

    # Two ASCII pieces, each the pyramid and the polyhedron.
    def piece(with_heat):
        def arr(name, vtk_type, values, nc=0):
            ncomp = f' NumberOfComponents="{nc}"' if nc else ""
            vals = " ".join(str(v) for v in values)
            return (
                f'        <DataArray type="{vtk_type}" Name="{name}"{ncomp}'
                f' format="ascii">\n          {vals}\n        </DataArray>\n'
            )

        heat = arr("heat", "Float64", _VTU_HEAT) if with_heat else ""
        return (
            f'    <Piece NumberOfPoints="{len(_VTU_POINTS)}"'
            f' NumberOfCells="{len(_VTU_TYPES)}">\n'
            "      <Points>\n"
            + arr("Points", "Float32", [v for p in _VTU_POINTS for v in p], 3)
            + "      </Points>\n      <Cells>\n"
            + arr("connectivity", "Int32", _VTU_CONN)
            + arr("offsets", "Int32", _VTU_OFFSETS)
            + arr("types", "UInt8", _VTU_TYPES)
            + arr("faces", "Int32", _VTU_FACES)
            + arr("faceoffsets", "Int32", _VTU_FACEOFFSETS)
            + "      </Cells>\n      <PointData>\n"
            + heat
            + "      </PointData>\n      <CellData>\n"
            + arr("material", "Int32", _VTU_MATERIAL)
            + "      </CellData>\n    </Piece>\n"
        )

    text = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
        "  <UnstructuredGrid>\n"
        + piece(True)
        + piece(False)
        + "  </UnstructuredGrid>\n</VTKFile>\n"
    )
    (out / "two_pieces_polyhedron.vtu").write_text(text)


# ---------------------------------------------------------------------------
# ANSYS Fluent .msh
# ---------------------------------------------------------------------------


def _newell(pts):
    n = [0.0, 0.0, 0.0]
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        n[0] += (a[1] - b[1]) * (a[2] + b[2])
        n[1] += (a[2] - b[2]) * (a[0] + b[0])
        n[2] += (a[0] - b[0]) * (a[1] + b[1])
    return n


def _fluent_faces(points, cells):
    """Every face of every cell, deduplicated, wound so that its right-hand
    normal points into ``c0`` (Fluent's rule). Cells are ``(faces, centroid)``
    with 1-based ids; returns ``[(nodes, c0, c1)]``, interior faces first."""
    seen = {}
    for cid, (cell_faces, _) in enumerate(cells, start=1):
        for f in cell_faces:
            seen.setdefault(frozenset(f), []).append((cid, f))
    interior, boundary = [], []
    for owners in seen.values():
        (c0, f), rest = owners[0], owners[1:]
        c1 = rest[0][0] if rest else 0
        centre = cells[c0 - 1][1]
        fc = [sum(points[i][k] for i in f) / len(f) for k in range(3)]
        n = _newell([points[i] for i in f])
        if sum(n[k] * (centre[k] - fc[k]) for k in range(3)) < 0:
            f = f[::-1]
        (interior if c1 else boundary).append((list(f), c0, c1))
    return interior, boundary


def _cells3d():
    """(points, cells) with cells as (1-based faces, centroid, fluent type)."""
    pts = []

    def add(p):
        pts.append(p)
        return len(pts)

    cells = []

    def cell(faces, nodes, ftype):
        centre = [sum(pts[i - 1][k] for i in nodes) / len(nodes) for k in range(3)]
        cells.append((faces, centre, ftype))

    # two tetrahedra sharing the face (a, b, c)
    a, b, c = add((0, 0, 0)), add((1, 0, 0)), add((0, 1, 0))
    d, e = add((0, 0, 1)), add((0.3, 0.3, -1))
    cell([(a, b, c), (a, b, d), (b, c, d), (c, a, d)], (a, b, c, d), 2)
    cell([(a, b, c), (a, b, e), (b, c, e), (c, a, e)], (a, b, c, e), 2)
    # a hexahedron and a polyhedron (the same cube, one side split)
    for x0, poly in ((3, False), (6, True)):
        h = [add((x0 + x, y, z)) for z in (0, 1) for y in (0, 1) for x in (0, 1)]
        b0, b1, b3, b2, t0, t1, t3, t2 = h
        faces = [
            (b0, b1, b2, b3),
            (t0, t1, t2, t3),
            (b0, b1, t1, t0),
            (b1, b2, t2, t1),
            (b2, b3, t3, t2),
        ]
        if poly:
            faces += [(b3, b0, t0), (b3, t0, t3)]
        else:
            faces += [(b3, b0, t0, t3)]
        cell(faces, h, 7 if poly else 4)
    # a wedge and a pyramid
    w = [add((9 + x, y, z)) for z in (0, 1) for x, y in ((0, 0), (1, 0), (0, 1))]
    cell(
        [
            (w[0], w[1], w[2]),
            (w[3], w[4], w[5]),
            (w[0], w[1], w[4], w[3]),
            (w[1], w[2], w[5], w[4]),
            (w[2], w[0], w[3], w[5]),
        ],
        w,
        6,
    )
    q = [add((12 + x, y, 0)) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))]
    apex = add((12.5, 0.5, 1))
    cell([tuple(q)] + [(q[i], q[(i + 1) % 4], apex) for i in range(4)], q + [apex], 5)
    points = [tuple(float(v) for v in p) for p in pts]
    return points, cells


def _fluent_ascii_faces(zone, bc, rows, mixed):
    first, last = rows[0][0], rows[-1][0]
    ftype = 0 if mixed else len(rows[0][1])
    out = f"(13 ({zone:x} {first:x} {last:x} {bc:x} {ftype:x})(\n"
    for _, f, c0, c1 in rows:
        lead = f"{len(f):x} " if mixed else ""
        out += lead + " ".join(f"{v:x}" for v in f) + f" {c0:x} {c1:x}\n"
    return out + "))\n"


def write_fluent():
    out = MESHES / "ansys"
    out.mkdir(parents=True, exist_ok=True)

    points, cells = _cells3d()
    # (0-based faces, centroid) per cell
    lookup = [([tuple(v - 1 for v in f) for f in fs], ctr) for fs, ctr, _ in cells]
    interior, boundary = _fluent_faces(points, lookup)
    n_int = len(interior)
    rows_int = [
        (i + 1, [v + 1 for v in f], c0, c1) for i, (f, c0, c1) in enumerate(interior)
    ]
    rows_wall = [
        (n_int + i + 1, [v + 1 for v in f], c0, c1)
        for i, (f, c0, c1) in enumerate(boundary)
    ]
    nf = len(rows_int) + len(rows_wall)
    text = '(0 "gen_feconv_quirk_fixtures.py: cells known only by their faces")\n'
    text += '(1 "fixture")\n(2 3)\n'
    text += f"(10 (0 1 {len(points):x} 0 3))\n(13 (0 1 {nf:x} 0))\n"
    text += f"(12 (0 1 {len(cells):x} 0))\n"
    text += f"(10 (1 1 {len(points):x} 1 3)(\n"
    text += "".join(" ".join(repr(v) for v in p) + "\n" for p in points) + "))\n"
    text += _fluent_ascii_faces(2, 2, rows_int, mixed=False)
    text += _fluent_ascii_faces(3, 3, rows_wall, mixed=True)
    types = " ".join(f"{c[2]:x}" for c in cells)
    text += f"(12 (4 1 {len(cells):x} 1 0)(\n{types}\n))\n"
    text += "(39 (4 fluid block-of-cells)())\n"
    text += "(45 (2 interior shared-face)())\n(45 (3 wall outer-walls)())\n"
    (out / "cells3d.msh").write_text(text)

    # 2-D: a quad (cell 2) left of x = 1 and a triangle (cell 1) right of it,
    # sharing the edge 2-3. Points 1..3 are node zone 7 and 4..5 zone 6, which
    # the file lists first.
    p2 = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (2.0, 0.0)]
    # edges (a, b, c0, c1): walking a -> b leaves c0 on the left
    edges = [
        (2, 3, 2, 1),  # interior: quad on the left going up x=1
        (1, 2, 2, 0),
        (3, 4, 2, 0),
        (4, 1, 2, 0),
        (2, 5, 1, 0),
        (5, 3, 1, 0),
    ]
    text = '(1 "TGrid 2D fixture")\n(2 2)\n(10 (0 1 5 0))\n(13 (0 1 6 0))\n'
    text += "(12 (0 1 2 0))\n"
    text += "(10 (6 4 5 1 2) (\n0.0 1.0\n2.0 0.0\n))\n"
    text += "(10 (7 1 3 1 2) (\n0.0 0.0\n1.0 0.0\n1.0 1.0\n))\n"
    text += "(13 (8 1 1 2 0) (\n2 2 3 2 1\n))\n"
    text += (
        "(13 (9 2 6 3 0) (\n"
        + "".join(f"2 {a:x} {b:x} {c0:x} {c1:x}\n" for a, b, c0, c1 in edges[1:])
        + "))\n"
    )
    text += "(12 (a 1 2 1 0)(\n1 3\n))\n"
    text += "(45 (10 fluid plate)())\n(45 (9 wall rim)())\n"
    (out / "tgrid2d.msh").write_text(text)

    # GAMBIT: the same 2-D mesh, boundary edges with c0 = 0 and reversed.
    text = '(0 "GAMBIT to Fluent File")\n(2 2)\n(10 (0 1 5 0 2))\n'
    text += "(10 (1 1 5 1 2)(\n" + "".join(f"{x} {y}\n" for x, y in p2) + "))\n"
    text += "(13(0 1 6 0))\n(13(8 1 1 2 2)(\n2 3 2 1\n))\n"
    text += (
        "(13(9 2 6 3 2)(\n"
        + "".join(f"{b:x} {a:x} 0 {c0:x}\n" for a, b, c0, _ in edges[1:])
        + "))\n"
    )
    text += "(12 (0 1 2 0))\n(12 (a 1 2 1 0))\n(45 (10 fluid plate)())\n"
    (out / "gambit2d.msh").write_text(text)

    # Binary: the two tetrahedra of cells3d.
    tp = points[:5]
    interior, boundary = _fluent_faces(tp, lookup[:2])
    body = b'(1 "binary fixture")\n(2 3)\n(10 (0 1 5 0 3))\n'
    body += b"(3010 (1 1 5 1 3)\n(" + struct.pack("<15d", *[v for p in tp for v in p])
    body += b")\nEnd of Binary Section 3010)\n"
    f, c0, c1 = interior[0]
    body += b"(2013 (2 1 1 2 3)\n(" + struct.pack("<5i", *[v + 1 for v in f], c0, c1)
    body += b")\nEnd of Binary Section 2013)\n"
    rows = []
    for f, c0, c1 in boundary:
        rows += [len(f)] + [v + 1 for v in f] + [c0, c1]
    body += f"(3013 (3 2 {len(boundary) + 1:x} 3 0)\n(".encode()
    body += struct.pack(f"<{len(rows)}q", *rows)
    body += b")\nEnd of Binary Section 3013)\n"
    body += b"(2012 (4 1 2 1 2)\nEnd of Binary Section 2012)\n"
    body += b"(45 (4 fluid tets)())\n"
    (out / "binary3d.msh").write_bytes(body)


if __name__ == "__main__":
    write_flux()
    write_gmsh()
    write_medit()
    write_vtu()
    write_fluent()
