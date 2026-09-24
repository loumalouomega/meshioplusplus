#!/usr/bin/env python3
"""Regenerate the OpenRadioss starter-deck fixtures under ``tests/python/meshes/radioss/``.

OpenRadioss's own QA decks are CC BY-NC and cannot be copied, so the decks are
written here from the card layouts of OpenRadioss's ``hm_cfg_files`` (``%10d``
integers, ``%20lg`` reals from format 51 on; ``%8d``/``%16lg`` before), each
element's nodes spelled out by position:

* ``/BRIC20``: corners 1-8, then the bottom-ring mid-edges 1-2, 2-3, 3-4, 4-1,
  the vertical ones 1-5 ... 4-8, then the top ring 5-6, 6-7, 7-8, 8-5;
* ``/TETRA10``: corners, then 1-2, 2-3, 3-1, 1-4, 2-4, 3-4 (a separate id line);
* ``/PENTA6``: bottom triangle 1-2-3, top 4-5-6;
* degenerate ``/BRICK`` rows: ``1 2 3 3 5 6 7 7`` and ``1 2 3 1 5 6 7 5``
  (wedges), ``1 2 3 3 4 4 4 4`` (tetra), ``1 2 3 4 5 5 5 5`` (pyramid).

``deck_0000.rad`` (format 2019) holds one element of every read card, parts
with titles, property and material ids, a nested ``/SUBSET``, every group form
the reader resolves (entity ids with a negative removal, ``PART``, ``SUBSET``,
another group, ``/GRNOD/SURF``) plus a ``/GRBRIC/BOX`` it does not, a
``/SURF/SEG`` over solid faces, a shell and a tet face, a comma-separated
``/NODE`` line, a ``/TETRA4`` written with the inverted winding, a single-node
``/SPRING``, an ``#include`` (``deck_include.inc``, which ends at ``#enddata``)
and text after ``/END`` that must not be read.

``old_0000.rad`` is a format-44 deck (8/16 columns): one brick, its part and a
node group.

    python tools/gen_radioss_fixtures.py
"""

import pathlib

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "radioss"
)


class Nodes:
    def __init__(self, start=1):
        self.ids = {}
        self.order = []
        self.next = start

    def __call__(self, x, y, z):
        key = (round(x, 9), round(y, 9), round(z, 9))
        if key not in self.ids:
            self.ids[key] = self.next
            self.order.append((self.next, key))
            self.next += 1
        return self.ids[key]


def box(nodes, x0, y0=0.0, z0=0.0, size=1.0):
    """The 8 corner ids of a cube, Radioss/VTK order (bottom CCW, then top)."""
    x1, y1, z1 = x0 + size, y0 + size, z0 + size
    return [
        nodes(x0, y0, z0), nodes(x1, y0, z0), nodes(x1, y1, z0), nodes(x0, y1, z0),
        nodes(x0, y0, z1), nodes(x1, y0, z1), nodes(x1, y1, z1), nodes(x0, y1, z1),
    ]  # fmt: skip


def mid(nodes, coords, a, b):
    pa, pb = coords[a], coords[b]
    return nodes(*[(p + q) / 2 for p, q in zip(pa, pb)])


def i10(*values, width=10):
    return "".join(f"{v:{width}d}" for v in values) + "\n"


def deck():
    nodes = Nodes()
    lines = []
    # Part 1: bricks, incl. degenerate ones.
    hexa = box(nodes, 0.0)
    w1 = box(nodes, 2.0)
    w1 = [w1[0], w1[1], w1[2], w1[2], w1[4], w1[5], w1[6], w1[6]]  # n3 == n4, n7 == n8
    w2 = box(nodes, 4.0)
    w2 = [w2[0], w2[1], w2[2], w2[0], w2[4], w2[5], w2[6], w2[4]]  # n4 == n1, n8 == n5
    t = box(nodes, 6.0)
    tet_brick = [t[0], t[1], t[2], t[2], t[4], t[4], t[4], t[4]]
    p = box(nodes, 8.0)
    apex = nodes(8.5, 0.5, 1.0)
    pyr = [p[0], p[1], p[2], p[3], apex, apex, apex, apex]
    # Part 2: quadratic solids and a penta.
    coords = {}

    def at(i):
        for nid, key in nodes.order:
            if nid == i:
                return key

    b20 = box(nodes, 10.0)
    for n in b20:
        coords[n] = at(n)

    def ring(a):
        return [mid(nodes, coords, a[k], a[(k + 1) % 4]) for k in range(4)]

    bottom = ring(b20[:4])
    top = ring(b20[4:])
    vertical = [mid(nodes, coords, b20[k], b20[k + 4]) for k in range(4)]
    bric20 = b20 + bottom + vertical + top
    tc = [
        nodes(12.0, 0.0, 0.0),
        nodes(13.0, 0.0, 0.0),
        nodes(12.0, 1.0, 0.0),
        nodes(12.0, 0.0, 1.0),
    ]
    for n in tc:
        coords[n] = at(n)
    tet10 = tc + [
        mid(nodes, coords, tc[a], tc[b])
        for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
    ]
    tet4 = [
        nodes(14.0, 0.0, 0.0),
        nodes(14.0, 1.0, 0.0),
        nodes(15.0, 0.0, 0.0),
        nodes(14.0, 0.0, 1.0),
    ]  # inverted
    penta = [nodes(16.0, 0.0, 0.0), nodes(17.0, 0.0, 0.0), nodes(16.0, 1.0, 0.0),
             nodes(16.0, 0.0, 1.0), nodes(17.0, 0.0, 1.0), nodes(16.0, 1.0, 1.0)]  # fmt: skip
    # Part 3: shells on the top of the hexahedron, and a separate triangle.
    shell_q = hexa[4:]
    tri_a, tri_b, tri_c = (
        nodes(0.0, 0.0, 2.0),
        nodes(1.0, 0.0, 2.0),
        nodes(0.0, 1.0, 2.0),
    )
    tri_d = nodes(1.0, 1.0, 2.0)
    # Part 4: lines.
    la, lb, lc = nodes(0.0, 3.0, 0.0), nodes(1.0, 3.0, 0.0), nodes(2.0, 3.0, 0.0)

    lines.append("#RADIOSS STARTER\n")
    lines.append("# meshio++ synthetic deck (tools/gen_radioss_fixtures.py)\n")
    lines.append("/BEGIN\n")
    lines.append("meshio_deck\n")
    lines.append(i10(2019, 0))
    lines.append("                  kg                  mm                  ms\n")
    lines.append("                  kg                  mm                  ms\n")
    lines.append("/TITLE\nmeshio++ Radioss fixture\n")
    lines.append("#include deck_include.inc\n")
    lines.append(
        "/NODE\n#  node_ID                  Xc                  Yc                  Zc\n"
    )
    all_nodes = list(nodes.order)
    for k, (nid, (x, y, z)) in enumerate(all_nodes):
        if nid == 1:
            lines.append(f"{nid},{x},{y},{z}\n")  # free format
        else:
            lines.append(f"{nid:10d}{x:20.6f}{y:20.6f}{z:20.6f}\n")
    lines.append("/PART/1\nbricks\n" + i10(11, 21, 1))
    lines.append("/BRICK/1\n")
    for eid, conn in ((1, hexa), (2, w1), (3, w2), (4, tet_brick), (5, pyr)):
        lines.append(i10(eid, *conn))
    lines.append("/PART/2\nquadratic\n" + i10(12, 21, 2))
    lines.append(
        "/BRIC20/2\n" + i10(6, *bric20[:8]) + i10(*bric20[8:16]) + i10(*bric20[16:])
    )
    lines.append("/TETRA10/2\n" + i10(7) + i10(*tet10))
    lines.append("/TETRA4/2\n" + i10(8, *tet4))
    lines.append("/PENTA6/2\n" + i10(9, *penta))
    lines.append("/PART/3\nshells\n" + i10(13, 22, 2))
    lines.append("/SHELL/3\n" + i10(1, *shell_q) + i10(2, tri_a, tri_b, tri_c, tri_c))
    lines.append("/SH3N/3\n" + i10(3, tri_b, tri_d, tri_c))
    lines.append("/PART/4\n\n" + i10(14, 23, 0))  # no title
    lines.append("/BEAM/4\n" + i10(1, la, lb, lc))
    lines.append("/TRUSS/4\n" + i10(2, lb, lc))
    lines.append("/SPRING/4\n" + i10(3, la, lc) + i10(4, lc, 0))
    lines.append("/SUBSET/1\nsolid parts\n" + i10(2))
    lines.append("/SUBSET/2\nquadratic subset\n")
    lines.append("/GRNOD/NODE/1\nhex nodes\n" + i10(*hexa))
    lines.append("/GRNOD/PART/2\nquadratic nodes\n" + i10(2))
    lines.append("/GRBRIC/BRIC/3\nsome bricks\n" + i10(1, 2, 3, 4, -3))
    lines.append("/GRBRIC/PART/4\npart 2 solids\n" + i10(2))
    lines.append("/GRSHEL/SHEL/5\nquad shell\n" + i10(1))
    lines.append("/GRBRIC/GRBRIC/6\nunion\n" + i10(3, 4))
    lines.append("/GRBRIC/SUBSET/7\nsubset solids\n" + i10(1))
    lines.append("/GRBRIC/BOX/8\nin a box\n" + i10(1))
    lines.append("/SURF/SEG/1\nskin\n")
    lines.append(i10(1, hexa[0], hexa[3], hexa[2], hexa[1]))  # hex bottom
    lines.append(i10(2, *shell_q))  # the quad shell (and hex top)
    lines.append(i10(3, tc[0], tc[1], tc[2], 0))  # a tet face
    lines.append("/GRNOD/SURF/9\nskin nodes\n" + i10(1))
    lines.append(
        "/MAT/ELAST/21\nsteel\n                7.8e-9\n              210000                 0.3\n"
    )
    lines.append("/END\n")
    lines.append(
        "/NODE\n        99                 9.0                 9.0                 9.0\n"
    )
    (OUT / "deck_0000.rad").write_text("".join(lines))
    include = [
        "# included by deck_0000.rad\n",
        "/PART/5\nincluded truss\n" + i10(15, 23, 0),
        "/TRUSS/5\n" + i10(10, la, lc),
        "#enddata\n",
        "/TRUSS/5\n" + i10(11, lb, lc),
    ]
    (OUT / "deck_include.inc").write_text("".join(include))


def old_deck():
    nodes = Nodes()
    hexa = box(nodes, 0.0)
    lines = ["#RADIOSS STARTER\n", "/BEGIN\n", "old_deck\n", i10(44, 0, width=8)]
    lines.append("                  kg                  mm                  ms\n" * 2)
    lines.append("/NODE\n")
    for nid, (x, y, z) in nodes.order:
        lines.append(f"{nid:8d}{x:16.6f}{y:16.6f}{z:16.6f}\n")
    lines.append("/PART/1\nold part\n" + i10(1, 1, 0, width=8))
    lines.append("/BRICK/1\n" + i10(1, *hexa, width=8))
    lines.append("/GRNOD/NODE/1\nbottom\n" + i10(*hexa[:4], width=8))
    lines.append("/END\n")
    (OUT / "old_0000.rad").write_text("".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    deck()
    old_deck()


if __name__ == "__main__":
    main()
