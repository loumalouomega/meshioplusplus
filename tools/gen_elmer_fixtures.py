#!/usr/bin/env python3
"""Regenerate the Elmer mesh-directory fixtures under ``tests/python/meshes/elmer/``.

The files are written here, not by meshio++, so the reader is never only checked
against its own writer. They follow the layout ElmerGrid writes (``mesh.header``
counts, ``id -1 x y z`` nodes, ``id body type nodes`` elements, ``id boundary
parent1 parent2 type nodes`` boundary elements, ``$ name = id`` names), and every
element lists its nodes in Elmer's own numbering, taken from the reference
coordinates of ElmerSolver's ``elements.def`` rather than from
``meshioplusplus._node_order``:

* 510 / 306 / 408 / 203: corners, then the mid-edge nodes of the edges
  ``(1,2) (2,3) (3,1) (1,4) (2,4) (3,4)`` (tetra), ``(1,2) (2,3) (3,1)``
  (triangle), ``(1,2) (2,3) (3,4) (4,1)`` (quad), ``(1,2)`` (line);
* 827: corners, the bottom-ring mid-edges, the *vertical* mid-edges, the top
  ring, then the mid-height face centres at ``v=-1, u=+1, v=+1, u=-1``, the
  bottom and top centres and the body centre; 409: corners, edges, centre.

Four directories:

* ``tet10_two_bodies``: two quadratic tets in two named bodies, a named outer
  face, the interface between them (two parents) and a parentless point
  condition (101).
* ``hex27_block``: one 827 brick and its 409 bottom face.
* ``quad8_2d``: a 2-D mesh of two 408 quads, 203 edges on both ends and no
  ``mesh.names`` (names default to ``body_<id>``/``boundary_<id>``).
* ``partitioned/partitioning.2``: two triangles split over two parts, the way
  ``ElmerGrid -partition ... -halo`` writes them: shared nodes listed by both
  parts, and part 2 holding a halo copy of element 1 as ``1/1``.

Two more are not written by this script: ``tet10_two_bodies_bin`` and
``quad8_2d_sbin`` are ElmerGrid's own binary output of the text fixtures
(elmerfem ``a8a13b5``, ``ElmerGrid 2 2 <dir> -out <dir>_bin -bin`` and
``-sbin``; ``entities.sif`` left out), the layout ElmerSolver reads.

    python tools/gen_elmer_fixtures.py
"""

import pathlib

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "elmer"
)


class Nodes:
    """Coordinates -> 1-based ids, in order of first use."""

    def __init__(self):
        self.ids = {}

    def __call__(self, xyz):
        key = tuple(float(v) for v in xyz)
        if key not in self.ids:
            self.ids[key] = len(self.ids) + 1
        return self.ids[key]

    def text(self):
        return "".join(
            f"{i} -1 {x:g} {y:g} {z:g}\n" for (x, y, z), i in self.ids.items()
        )


def mid(a, b):
    return tuple((p + q) / 2 for p, q in zip(a, b))


def header(n_nodes, elements, boundary):
    """ElmerGrid's layout: counts, the number of types, then `type count` rows."""
    counts = {}
    for line in elements:
        code = int(line.split()[2])
        counts[code] = counts.get(code, 0) + 1
    for line in boundary:
        code = int(line.split()[4])
        counts[code] = counts.get(code, 0) + 1
    rows = [
        f"{n_nodes:<6d} {len(elements):<6d} {len(boundary):<6d}",
        f"{len(counts):<6d}",
    ]
    rows += [f"{code:<6d} {counts[code]:<6d}" for code in sorted(counts)]
    return "\n".join(rows) + "\n"


def write(directory, files):
    directory.mkdir(parents=True, exist_ok=True)
    for name, text in files.items():
        with open(directory / name, "w", newline="\n") as f:
            f.write(text)
    print(f"wrote {directory}")


def tet10(nodes, c):
    """Elmer 510: corners, then mid-edges (1,2) (2,3) (3,1) (1,4) (2,4) (3,4)."""
    ids = [nodes(p) for p in c]
    for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)):
        ids.append(nodes(mid(c[a], c[b])))
    return ids


def tri6(nodes, c):
    ids = [nodes(p) for p in c]
    for a, b in ((0, 1), (1, 2), (2, 0)):
        ids.append(nodes(mid(c[a], c[b])))
    return ids


def tet10_two_bodies():
    nodes = Nodes()
    A, B, C, D, E = (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)
    t1 = tet10(nodes, [A, B, C, D])
    t2 = tet10(nodes, [B, C, D, E])
    elements = [
        "1 1 510 " + " ".join(map(str, t1)),
        "2 2 510 " + " ".join(map(str, t2)),
    ]
    boundary = [
        "1 2 1 0 306 " + " ".join(map(str, tri6(nodes, [A, B, D]))),
        "2 3 1 2 306 " + " ".join(map(str, tri6(nodes, [B, C, D]))),
        f"3 4 0 0 101 {nodes(A)}",
    ]
    names = (
        "! ----- names for bodies -----\n$ left = 1\n$ right = 2\n"
        "! ----- names for boundaries -----\n$ outer = 2\n$ interface = 3\n"
        "$ corner = 4\n"
    )
    return {
        "mesh.header": header(len(nodes.ids), elements, boundary),
        "mesh.nodes": nodes.text(),
        "mesh.elements": "\n".join(elements) + "\n",
        "mesh.boundary": "\n".join(boundary) + "\n",
        "mesh.names": names,
    }


# Elmer 827 reference coordinates (elements.def), node by node.
# fmt: off
_U827 = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 0, 0]
_V827 = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, 0, 0, 0]
_W827 = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 1, 0]
# fmt: on
# Elmer 409: corners, edges (1,2) (2,3) (3,4) (4,1), centre.
_U409 = [-1, 1, 1, -1, 0, 1, 0, -1, 0]
_V409 = [-1, -1, 1, 1, -1, 0, 1, 0, 0]


def hex27_block():
    """One brick on [0,2]x[0,1]x[0,1], its nodes placed from Elmer's reference
    coordinates, and its bottom face."""
    nodes = Nodes()

    def at(u, v, w):
        return (u + 1.0, (v + 1) / 2, (w + 1) / 2)

    brick = [nodes(at(u, v, w)) for u, v, w in zip(_U827, _V827, _W827)]
    bottom = [nodes(at(u, v, -1)) for u, v in zip(_U409, _V409)]
    elements = ["1 1 827 " + " ".join(map(str, brick))]
    boundary = ["1 1 1 0 409 " + " ".join(map(str, bottom))]
    return {
        "mesh.header": header(len(nodes.ids), elements, boundary),
        "mesh.nodes": nodes.text(),
        "mesh.elements": "\n".join(elements) + "\n",
        "mesh.boundary": "\n".join(boundary) + "\n",
        "mesh.names": "! ----- names for bodies -----\n$ block = 1\n"
        "! ----- names for boundaries -----\n$ bottom = 1\n",
    }


def quad8_2d():
    """Two 408 quads side by side in body 1, the left and right edges as 203s."""
    nodes = Nodes()

    def quad8(c):
        ids = [nodes(p) for p in c]
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            ids.append(nodes(mid(c[a], c[b])))
        return ids

    q1 = quad8([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
    q2 = quad8([(1, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0)])
    elements = [
        "1 1 408 " + " ".join(map(str, q1)),
        "2 1 408 " + " ".join(map(str, q2)),
    ]
    left = [nodes((0, 1, 0)), nodes((0, 0, 0)), nodes((0, 0.5, 0))]
    right = [nodes((2, 0, 0)), nodes((2, 1, 0)), nodes((2, 0.5, 0))]
    boundary = [
        "1 1 1 0 203 " + " ".join(map(str, left)),
        "2 2 2 0 203 " + " ".join(map(str, right)),
    ]
    return {
        "mesh.header": header(len(nodes.ids), elements, boundary),
        "mesh.nodes": nodes.text(),
        "mesh.elements": "\n".join(elements) + "\n",
        "mesh.boundary": "\n".join(boundary) + "\n",
    }


def partitioned():
    """Two triangles over two parts; part 2 also holds a halo copy of element 1."""
    part1 = {
        "part.1.header": "3      1      1     \n2     \n202    1     \n303    1     \n2      0     \n",
        "part.1.nodes": "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 0 1 0\n",
        "part.1.elements": "1 1 303 1 2 3\n",
        "part.1.boundary": "1 1 1 0 202 1 2\n",
        "part.1.shared": "2 2 1 2\n3 2 1 2\n",
    }
    part2 = {
        "part.2.header": "4      2      1     \n2     \n202    1     \n303    2     \n2      0     \n",
        "part.2.nodes": "2 -1 1 0 0\n4 -1 1 1 0\n3 -1 0 1 0\n1 -1 0 0 0\n",
        "part.2.elements": "2 1 303 2 4 3\n1/1 1 303 1 2 3\n",
        "part.2.boundary": "1 2 2 0 202 2 4\n",
        "part.2.shared": "2 2 1 2\n3 2 1 2\n",
    }
    return {**part1, **part2}


def main():
    write(OUT / "tet10_two_bodies", tet10_two_bodies())
    write(OUT / "hex27_block", hex27_block())
    write(OUT / "quad8_2d", quad8_2d())
    write(OUT / "partitioned" / "partitioning.2", partitioned())
    write(
        OUT / "partitioned",
        {
            "mesh.names": "! ----- names for bodies -----\n$ plate = 1\n"
            "! ----- names for boundaries -----\n$ bottom = 1\n$ right = 2\n"
        },
    )


if __name__ == "__main__":
    main()
