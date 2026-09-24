#!/usr/bin/env python3
"""Regenerate the libMesh fixtures under ``tests/python/meshes/libmesh/``.

libMesh (LGPL) is not a dependency and its own sample meshes are not copied: the
streams are written here, value by value, the way ``XdrIO::write`` and
``xdr_cxx.C`` lay them out -- in ASCII (``.xda``: ``value\\t # comment`` scalars,
a vector as its length line then its items, bulk data as bare values) and in XDR
(``.xdr``: big-endian, 4-byte lengths, strings padded to 4 bytes, 8-byte header
integers from 1.3.0 on). Each element lists its nodes in libMesh's order, built
from libMesh's reference-element node positions spelled out below
(``cell_hex27.C``, ``cell_tet14.C``, ``cell_prism18.C``, ``face_quad4.C``), so a
test can check the read cells against meshio++'s own edge and face tables.

Five meshes, each as ``.xda`` and ``.xdr``:

* ``hex27`` (libMesh-1.3.0, 8-byte fields): two HEX27 side by side in
  subdomains 1 (named ``left``) and 2 (unnamed), with unique ids, a side set on
  each end face (id 1 named ``inlet``, id 2 unnamed) and a node set (id 5,
  ``wall``) on the x = 0 face.
* ``one_hex`` (libMesh-0.7.0+): one HEX8 with the six side sets libMesh's own
  reference elements carry, in the old 4-byte header layout.
* ``amr_quad`` (libMesh-1.8.0): a QUAD4 refined once into four children, inline
  p-levels, one extra element integer, an unused node id (NaN coordinates) and a
  side set on the parent's bottom edge that must reach the two bottom children.
* ``tet14_prism18`` (libMesh-1.3.0, 4-byte fields): a TET14 and a PRISM18
  (read as ``tetra10`` and ``wedge18``).
* ``edges_shell`` (libMesh-1.8.0): a HEX20 with a QUAD8 shell on its top face,
  two edge sets (one edge named twice) and both shell faces.

Plus ``hex27.xda.gz`` and ``hex27.xdr.bz2``, compressed as libMesh writes them.

    python tools/gen_libmesh_fixtures.py
"""

import bz2
import gzip
import math
import pathlib
import struct

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "libmesh"
)

# Reference node positions, libMesh order.
_C = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
]
_HEX_EDGES = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
]
_HEX_SIDES = [
    (0, 3, 2, 1),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
    (4, 5, 6, 7),
]


def _mid(points, idx):
    return tuple(sum(points[i][d] for i in idx) / len(idx) for d in range(3))


def hex27_reference():
    p = list(_C)
    p += [_mid(_C, e) for e in _HEX_EDGES]  # 8..19
    p += [_mid(_C, s) for s in _HEX_SIDES]  # 20..25
    p.append((0.5, 0.5, 0.5))  # 26
    return p


def tet14_reference():
    c = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    edges = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
    sides = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]
    return c + [_mid(c, e) for e in edges] + [_mid(c, s) for s in sides]


def prism18_reference():
    c = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (5, 3)]
    quads = [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
    return c + [_mid(c, e) for e in edges] + [_mid(c, q) for q in quads]


class Stream:
    """libMesh's Xdr in WRITE (ASCII) or ENCODE (XDR) mode."""

    def __init__(self, xdr, header_width, type_width):
        self.xdr = xdr
        self.hw = header_width
        self.tw = type_width
        self.out = [] if not xdr else bytearray()

    # -- XDR primitives ---------------------------------------------------------
    def _int(self, v, width):
        self.out += struct.pack(">Q" if width == 8 else ">I", v)

    def _str(self, s):
        b = s.encode()
        self.out += struct.pack(">I", len(b)) + b + b"\0" * ((4 - len(b) % 4) % 4)

    # -- value kinds ------------------------------------------------------------
    def string(self, s, comment):
        if self.xdr:
            self._str(s)
        else:
            self.out.append(f"{s}\t {comment}\n" if comment else f"{s}\n")

    def scalar(self, v, comment, width=None):
        if self.xdr:
            self._int(v, width or self.hw)
        else:
            self.out.append(f"{v}\t {comment}\n")

    def int_vector(self, values, width, comment=""):
        if self.xdr:
            self._int(len(values), 4)
            for v in values:
                self._int(v, width)
            return
        self.out.append(f"{len(values)}\t # vector length\n")
        tail = f"\t {comment}" if comment else ""
        self.out.append("".join(f"{v}\t " for v in values) + tail + "\n")

    def string_vector(self, values, comment=""):
        if self.xdr:
            self._int(len(values), 4)
            for v in values:
                self._str(v)
            return
        self.out.append(f"{len(values)}\t # vector length\n")
        tail = f"\t {comment}" if comment else ""
        self.out.append("".join(f"{v}\t " for v in values) + tail + "\n")

    def ints(self, values):
        if self.xdr:
            for v in values:
                self._int(v, self.tw)
        else:
            self.out.append(" ".join(str(v) for v in values) + "\n")

    def reals(self, values):
        if self.xdr:
            for v in values:
                self.out += struct.pack(">d", v)
        else:
            self.out.append(" ".join(f"{v:.17e}" for v in values) + "\n")

    def data(self):
        return bytes(self.out) if self.xdr else "".join(self.out).encode()


def write_mesh(name, spec):
    for ext, xdr in (("xda", False), ("xdr", True)):
        v = spec["version"]
        v130 = any(x in v for x in ("1.3.0", "1.8.0"))
        v092 = v130 or any(x in v for x in ("0.9.2", "0.9.6", "1.1.0"))
        v096 = v130 or any(x in v for x in ("0.9.6", "1.1.0"))
        v110 = v130 or "1.1.0" in v
        v180 = "1.8.0" in v
        hw = 8 if v130 else 4
        tw = spec.get("type_size", 8) if v130 else 4
        io = Stream(xdr, hw, tw)
        io.string(v, "")
        elems = spec["elements"]
        io.scalar(len(elems), "# number of elements")
        io.scalar(len(spec["coords"]), "# number of nodes")
        bc = "." if spec.get("sides") is not None else "n/a"
        io.string(bc, "# boundary condition specification file")
        io.string(".", "# subdomain id specification file")
        io.string("n/a", "# processor id specification file")
        io.string("." if spec.get("p_level") else "n/a", "# p-level specification file")
        uid = spec.get("uid", False)
        if v092:
            for label, size in (
                ("type", tw),
                ("uid", tw if uid else 0),
                ("pid", 0),
                ("sid", tw),
                ("p-level", tw if spec.get("p_level") else 0),
                ("eid", tw),
                ("side", tw),
                ("bid", tw),
            ):
                io.scalar(size, f"# {label} size")
        if v180:
            extra = spec.get("elem_ints", [])
            io.scalar(tw if extra else 0, "# extra integer size")
            io.string_vector([], "# node integer names")
            io.string_vector(extra, "# elem integer names")
            io.int_vector([], tw, "# elemset codes")
        if v092:
            names = spec.get("subdomain_names", {})
            io.scalar(len(names), "# subdomain id to name map")
            if names:
                io.int_vector(sorted(names), hw)
                io.string_vector([names[k] for k in sorted(names)])
        # Connectivity, one block per level.
        levels = sorted({e["level"] for e in elems})
        for level in levels:
            block = [e for e in elems if e["level"] == level]
            io.scalar(len(block), f"# n_elem at level {level}", width=tw)
            for e in block:
                rec = [e["type"]]
                if uid:
                    rec.append(e["uid"])
                if level > 0:
                    rec.append(e["parent"])
                rec.append(e.get("sid", 0))
                if spec.get("p_level"):
                    rec.append(e.get("p", 0))
                rec += e["nodes"]
                rec += e.get("extra", [])
                io.ints(rec)
        for c in spec["coords"]:
            io.reals(c)
        if v096:
            io.scalar(0, "# presence of unique ids", width=4)
        if bc == "n/a":
            (OUT / f"{name}.{ext}").write_bytes(io.data())
            continue

        def triples(values, names, count_comment):
            if v092:
                io.scalar(len(names), "# sideset id to name map")
                if names:
                    io.int_vector(sorted(names), hw)
                    io.string_vector([names[k] for k in sorted(names)])
            io.scalar(len(values), count_comment)
            for t in values:
                io.ints(t)

        triples(
            spec["sides"],
            spec.get("sideset_names", {}),
            "# number of side boundary conditions",
        )
        if v092:
            names = spec.get("nodeset_names", {})
            io.scalar(len(names), "# nodeset id to name map")
            if names:
                io.int_vector(sorted(names), hw)
                io.string_vector([names[k] for k in sorted(names)])
            io.scalar(len(spec.get("nodesets", [])), "# number of nodesets")
            for t in spec.get("nodesets", []):
                io.ints(t)
        if v110:
            triples(
                spec.get("edges", []),
                spec.get("sideset_names", {}),
                "# number of edge boundary conditions",
            )
            triples(
                spec.get("shellfaces", []),
                spec.get("sideset_names", {}),
                "# number of shellface boundary conditions",
            )
        (OUT / f"{name}.{ext}").write_bytes(io.data())


class Nodes:
    """Global node ids by position."""

    def __init__(self):
        self.ids = {}
        self.coords = []

    def __call__(self, p):
        key = tuple(round(x, 9) for x in p)
        if key not in self.ids:
            self.ids[key] = len(self.coords)
            self.coords.append(tuple(float(x) for x in p))
        return self.ids[key]


def hex27_mesh():
    nodes = Nodes()
    ref = hex27_reference()
    elements = []
    for k in range(2):
        conn = [nodes((x + k, y, z)) for x, y, z in ref]
        elements.append(
            {"type": 12, "uid": 100 + k, "level": 0, "sid": k + 1, "nodes": conn}
        )
    wall = [nodes.ids[key] for key in sorted(nodes.ids) if key[0] == 0.0]
    return {
        "version": "libMesh-1.3.0",
        "elements": elements,
        "coords": nodes.coords,
        "uid": True,
        "subdomain_names": {1: "left"},
        "sides": [(0, 4, 1), (1, 2, 2)],  # element 0's x- face, element 1's x+ face
        "sideset_names": {1: "inlet"},
        "nodesets": [(n, 5) for n in wall],
        "nodeset_names": {5: "wall"},
    }


def one_hex_mesh():
    return {
        "version": "libMesh-0.7.0+",
        "elements": [{"type": 10, "level": 0, "nodes": list(range(8))}],
        "coords": [tuple(2.0 * c - 1.0 for c in p) for p in _C],
        "sides": [(0, s, s) for s in range(6)],
    }


def amr_quad_mesh():
    # Parent corners 0..3 on [0, 2]^2, then the refinement's new nodes, then an
    # id no element uses (its coordinates are NaN, as libMesh writes gaps).
    coords = [
        (0, 0, 0),
        (2, 0, 0),
        (2, 2, 0),
        (0, 2, 0),
        (1, 0, 0),
        (2, 1, 0),
        (1, 2, 0),
        (0, 1, 0),
        (1, 1, 0),
    ]
    coords = [tuple(float(x) for x in c) for c in coords] + [(math.nan,) * 3]
    children = [(0, 4, 8, 7), (4, 1, 5, 8), (8, 5, 2, 6), (7, 8, 6, 3)]
    elements = [{"type": 5, "level": 0, "p": 0, "nodes": [0, 1, 2, 3], "extra": [7]}]
    for k, conn in enumerate(children):
        elements.append(
            {
                "type": 5,
                "level": 1,
                "parent": 0,
                "p": [1, 1, 2, 0][k],
                "nodes": list(conn),
                "extra": [k],
            }
        )
    return {
        "version": "libMesh-1.8.0",
        "elements": elements,
        "coords": coords,
        "p_level": True,
        "elem_ints": ["weight"],
        "sides": [(0, 0, 3)],  # the parent's bottom edge
        "sideset_names": {3: "bottom"},
    }


def tet14_prism18_mesh():
    nodes = Nodes()
    tet = [nodes(p) for p in tet14_reference()]
    prism = [nodes((x + 2.0, y, z)) for x, y, z in prism18_reference()]
    return {
        "version": "libMesh-1.3.0",
        "type_size": 4,
        "elements": [
            {"type": 34, "level": 0, "sid": 1, "nodes": tet},
            {"type": 15, "level": 0, "sid": 2, "nodes": prism},
        ],
        "coords": nodes.coords,
        "sides": None,
    }


def hex20_reference():
    return list(_C) + [_mid(_C, e) for e in _HEX_EDGES]


def edges_shell_mesh():
    # A HEX20 (subdomain 1) with a QUAD8 shell on its top face (subdomain 2).
    # Edge sets: the hex's edge 0 (id 11, named) and its top edge 8, also named
    # by the shell's edge 0 under the same id 12 (one line3 for both). Shell
    # faces: face 0 id 20 (named), face 1 id 21. A side set on the hex bottom.
    nodes = Nodes()
    ref = hex20_reference()
    hexa = [nodes(p) for p in ref]
    top = [4, 5, 6, 7, 16, 17, 18, 19]  # the hex's top corners and edge mids
    shell = [hexa[k] for k in top]
    return {
        "version": "libMesh-1.8.0",
        "elements": [
            {"type": 11, "level": 0, "sid": 1, "nodes": hexa},
            {"type": 6, "level": 0, "sid": 2, "nodes": shell},
        ],
        "coords": nodes.coords,
        "subdomain_names": {2: "skin"},
        "sides": [(0, 0, 1)],
        "sideset_names": {1: "bottom", 11: "axis", 20: "front"},
        "edges": [(0, 0, 11), (0, 8, 12), (1, 0, 12)],
        "shellfaces": [(1, 0, 20), (1, 1, 21)],
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    write_mesh("hex27", hex27_mesh())
    write_mesh("one_hex", one_hex_mesh())
    write_mesh("amr_quad", amr_quad_mesh())
    write_mesh("tet14_prism18", tet14_prism18_mesh())
    write_mesh("edges_shell", edges_shell_mesh())
    # Compressed copies, as libMesh writes them (gzip without a timestamp, so
    # the bytes are reproducible).
    (OUT / "hex27.xda.gz").write_bytes(
        gzip.compress((OUT / "hex27.xda").read_bytes(), mtime=0)
    )
    (OUT / "hex27.xdr.bz2").write_bytes(bz2.compress((OUT / "hex27.xdr").read_bytes()))


if __name__ == "__main__":
    main()
