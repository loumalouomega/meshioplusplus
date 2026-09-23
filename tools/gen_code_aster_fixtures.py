#!/usr/bin/env python3
"""Regenerate the Code_Aster ``.mail`` fixtures under ``tests/python/meshes/code_aster/``.

The files are written here, not by meshio++, so the reader is never only checked
against its own writer. Each element's nodes are listed in Code_Aster's own
numbering, spelled out below edge by edge rather than taken from
``meshioplusplus._node_order`` (the registry was derived from the same Code_Aster
sources, its gmsh reader ``inigms.F90`` and MED reader ``lrmtyp.F90``, but by
composing permutation tables, not by naming edges):

* corners first, in meshio++'s (VTK's) corner order, which Code_Aster shares;
* HEXA20/27 and PENTA15/18: the bottom-ring mid-edges, then the vertical ones,
  then the top ring, then (HEXA27) the face centres ``1234, 1265, 2376, 3487,
  1485, 5678`` and the body centre, or (PENTA18) the quadrilateral face centres;
* TETRA10, PYRAM13, TRIA6/7, QUAD8/9, SEG3: mid-edges in the order of the edges
  listed below, then any face or body centre.

The tests check the reader puts every mid-edge node on its edge's midpoint.

    python tools/gen_code_aster_fixtures.py
"""

import pathlib

import numpy as np

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "code_aster"
)

# Code_Aster keyword -> corner coordinates, mid-edge corner pairs, face-centre
# corner lists and whether there is a body centre, all in Code_Aster numbering.
_HEX = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]
_HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]
_HEX_EDGES += [(0, 4), (1, 5), (2, 6), (3, 7)]
_HEX_EDGES += [(4, 5), (5, 6), (6, 7), (7, 4)]
_HEX_FACES = [(0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (0, 3, 7, 4)]
_HEX_FACES += [(4, 5, 6, 7)]
_WED = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
_WED_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (5, 3)]
_WED_FACES = [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
_TET = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
_TET_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
_PYR = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)]
_PYR_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)]
_TRI = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
_TRI_EDGES = [(0, 1), (1, 2), (2, 0)]
_QUA = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
_QUA_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]

SHAPES = {
    "POI1": ([(0, 0, 0)], [], [], False),
    "SEG2": ([(0, 0, 0), (1, 0, 0)], [], [], False),
    "SEG3": ([(0, 0, 0), (1, 0, 0)], [(0, 1)], [], False),
    "SEG4": ([(0, 0, 0), (1, 0, 0)], [], [], False),
    "TRIA3": (_TRI, [], [], False),
    "TRIA6": (_TRI, _TRI_EDGES, [], False),
    "TRIA7": (_TRI, _TRI_EDGES, [(0, 1, 2)], False),
    "QUAD4": (_QUA, [], [], False),
    "QUAD8": (_QUA, _QUA_EDGES, [], False),
    "QUAD9": (_QUA, _QUA_EDGES, [(0, 1, 2, 3)], False),
    "TETRA4": (_TET, [], [], False),
    "TETRA10": (_TET, _TET_EDGES, [], False),
    "PENTA6": (_WED, [], [], False),
    "PENTA15": (_WED, _WED_EDGES, [], False),
    "PENTA18": (_WED, _WED_EDGES, _WED_FACES, False),
    "PYRAM5": (_PYR, [], [], False),
    "PYRAM13": (_PYR, _PYR_EDGES, [], False),
    "HEXA8": (_HEX, [], [], False),
    "HEXA20": (_HEX, _HEX_EDGES, [], False),
    "HEXA27": (_HEX, _HEX_EDGES, _HEX_FACES, True),
}


def element_nodes(keyword, origin=(0.0, 0.0, 0.0), scale=1.0):
    """The coordinates of one element's nodes, in Code_Aster order."""
    corners, edges, faces, body = SHAPES[keyword]
    c = np.array(corners, dtype=float) * scale + np.array(origin)
    nodes = list(c)
    if keyword == "SEG4":  # the two inner nodes at thirds, from node 1 to node 2
        nodes += [c[0] + (c[1] - c[0]) / 3, c[0] + 2 * (c[1] - c[0]) / 3]
    nodes += [(c[a] + c[b]) / 2 for a, b in edges]
    nodes += [c[list(f)].mean(0) for f in faces]
    if body:
        nodes.append(c.mean(0))
    return nodes


def all_elements():
    """One element of every keyword, spread along x, with commas, comments and
    records wrapped at six nodes a line."""
    lines = ["% One element of every Code_Aster keyword meshio++ reads.", "COOR_3D"]
    elements = []
    count = 0
    for k, keyword in enumerate(SHAPES):
        names = []
        for x in element_nodes(keyword, origin=(2.0 * k, 0.0, 0.0)):
            count += 1
            names.append(f"NO{count}")
            lines.append(f" NO{count}, {x[0]:.4f}, {x[1]:.4f}, {x[2]:.4f}")
        elements.append((keyword, names))
    lines.append("FINSF")
    for keyword, names in elements:
        lines.append(f"{keyword}   % {len(names)} nodes")
        record = f" E_{keyword}"
        for i, name in enumerate(names):
            record += " " + name
            if i % 6 == 5 and i + 1 < len(names):
                lines.append(record)
                record = "        "
        lines.append(record)
        lines.append("FINSF")
    lines.append("GROUP_MA NOM=SOLIDS")
    solids = [f"E_{k}" for k in SHAPES if SHAPES[k][0] in (_HEX, _WED, _TET, _PYR)]
    for j in range(0, len(solids), 5):
        lines.append("  " + " ".join(solids[j : j + 5]))
    lines.append("FINSF")
    lines.append("FIN")
    return "\n".join(lines) + "\n"


def hexa20_block():
    """A 2x1x1 block of HEXA20 with a QUAD8 face and point groups, laid out the
    way Code_Aster's own writer (``PRE_IDEAS``, ``bibfor/stbtrias``) lays a file
    out: a TITRE, a multi-line block header, ``%FORMAT=`` comments, and element
    records of a name and eight nodes a line."""
    coords = {}
    order = []

    def node(p):
        key = tuple(round(v, 6) for v in p)
        if key not in coords:
            coords[key] = f"NO{len(coords) + 1}"
            order.append(key)
        return coords[key]

    cells = []
    for i in range(2):
        cells.append(
            [node(p) for p in element_nodes("HEXA20", origin=(float(i), 0, 0))]
        )
    # The x = 0 face as a QUAD8, corners 1 4 8 5 of the first hexahedron.
    face = element_nodes("QUAD8")
    face = [np.array([0.0, p[0], p[1]]) for p in face]
    face_nodes = [node(p) for p in face]

    lines = [
        "TITRE    NOM=INDEFINI",
        " a 2x1x1 block of HEXA20, written by tools/gen_code_aster_fixtures.py",
        "FINSF",
        "%",
        "COOR_3D    NOM=INDEFINI   NBOBJ=" + str(len(order)),
        "           NUMIN=1                   NUMAX=" + str(len(order)),
        "           AUTEUR=MESHIO++           DATE=23/09/2026",
        "%FORMAT=(1*NOM_DE_NOEUD,3*COORD)",
    ]
    for key in order:
        x, y, z = key
        lines.append(f"  {coords[key]:<8}  {x:21.14E} {y:21.14E} {z:21.14E}")
    lines += ["FINSF", "%", "HEXA20     NOM=INDEFINI   NBOBJ=2"]
    lines.append("%FORMAT=(1*NOM_DE_MAILLE,20*NOM_DE_NOEUD)")
    for m, nodes in enumerate(cells, start=1):
        name = f"MA{m}"
        lines.append(f"{name:<8}" + "".join(f" {n:<8}" for n in nodes[:8]))
        for j in range(8, 20, 8):
            lines.append(" " * 8 + "".join(f" {n:<8}" for n in nodes[j : j + 8]))
    lines += ["FINSF", "%", "QUAD8      NOM=INDEFINI   NBOBJ=1"]
    lines.append("MA3     " + "".join(f" {n:<8}" for n in face_nodes))
    lines += ["FINSF", "%"]
    lines += ["GROUP_MA   NOM=VOLUME", "  MA1 MA2", "FINSF", "%"]
    lines += ["GROUP_MA   NOM=X0", "  MA3", "FINSF", "%"]
    lines += ["GROUP_NO   NOM=X0", "  " + " ".join(face_nodes[:4]), "FINSF", "%"]
    lines += ["FIN"]
    return "\n".join(lines) + "\n"


def plate_2d():
    """A COOR_2D plate: TRIA6, QUAD8 and SEG3 with groups named by their first
    token, and a point listed twice in a group (Code_Aster warns and keeps one)."""
    lines = [
        "COOR_2D",
        " P1 0.0 0.0",
        " P2 1.0 0.0",
        " P3 1.0 1.0",
        " P4 0.0 1.0",
        " P5 2.0 0.0",
        " Q1 0.5 0.0",
        " Q2 1.0 0.5",
        " Q3 0.5 1.0",
        " Q4 0.0 0.5",
        " Q5 1.5 0.0",
        " Q6 1.5 0.5",
        "FINSF",
        "QUAD8",
        " SQUARE P1 P2 P3 P4 Q1 Q2 Q3 Q4",
        "FINSF",
        "TRIA6",
        " WING P2 P5 P3 Q5 Q6 Q2",
        "FINSF",
        "SEG3",
        " BOTTOM1 P1 P2 Q1",
        " BOTTOM2 P2 P5 Q5",
        "FINSF",
        "GROUP_MA",
        " SURFACE SQUARE WING",
        "FINSF",
        "GROUP_MA",
        " EDGE BOTTOM1, BOTTOM2",
        "FINSF",
        "GROUP_NO",
        " CORNERS P1 P5 P1",
        "FINSF",
        "FIN",
    ]
    return "\n".join(lines) + "\n"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, text in (
        ("all_elements.mail", all_elements()),
        ("hexa20_block.mail", hexa20_block()),
        ("plate_2d.mail", plate_2d()),
    ):
        with open(OUT / name, "w", newline="\n") as f:
            f.write(text)
        print(f"wrote {OUT / name}")


if __name__ == "__main__":
    main()
