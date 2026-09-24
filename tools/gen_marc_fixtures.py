"""Write the MSC Marc fixtures of ``tests/python/meshes/marc/``.

No Marc licence was available and no permissively licensed Marc deck or post
file exists, so the fixtures are written here from the layouts of Marc Volume C
(program input) and Volume D (PLDUMP2000, the post file), in the styles Marc
Mentat writes (checked against real Mentat 2020 decks and a real ``.t19``, which
cannot be committed: see ``tests/python/meshes/marc/README.md``).

- ``hex20.dat``: two 20-node bricks (type 21) in fixed 5/10-column fields,
  non-sequential node numbers, and element and node sets using ranges,
  continuation lines, set names, ``EXCEPT`` and ``AND``.
- ``hex20_extended.dat``: the same model in ``EXTENDED`` 10/20-column fields,
  with Mentat's exponent without ``E`` (``-1.5-1``) and fields that touch.
- ``mixed_free.dat``: free (comma) format: a brick, a brick with collapsed
  nodes (a wedge), a 10-node tetrahedron, a shell quad, a beam, an element of a
  type meshio++ does not read, and a face set it skips.
- ``plane_quad8.dat``: 8-node plane-strain quads (type 27) with two coordinates
  per node.
- ``results.t19``: a formatted post file of two 8-node bricks: two increments,
  displacements and reaction forces, and the equivalent stress with the stress
  tensor (codes 17, 311..316) at eight integration points; node and element sets.

Run from the repository root: ``python tools/gen_marc_fixtures.py``.
"""

import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "tests" / "python" / "meshes" / "marc"


# -- the two-brick hex20 model -------------------------------------------------

# Local lattice positions (half-units) of a 20-node brick's nodes, in Marc's
# order: corners, the bottom and top mid-edge rings, the vertical mid-edges.
_HEX20 = [
    (0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
    (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
    (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
    (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
    (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
]  # fmt: skip


def _node_id(i, j, k):
    return 10 * (1 + i + 5 * j + 15 * k)


def hex20_model():
    """(nodes {id: xyz}, elements [(id, type, nodes)]) of two unit bricks side
    by side along x, and the node ids of the faces x = 0 and z = 1."""
    nodes = {}
    for k in range(3):
        for j in range(3):
            for i in range(5):
                if (i % 2) + (j % 2) + (k % 2) <= 1:
                    nodes[_node_id(i, j, k)] = (0.5 * i, 0.5 * j, 0.5 * k)
    elements = []
    for e in range(2):
        ids = [_node_id(2 * e + a, b, c) for a, b, c in _HEX20]
        elements.append((100 * (e + 1), 21, ids))
    x0 = sorted(n for n, p in nodes.items() if p[0] == 0.0)
    top = sorted(n for n, p in nodes.items() if p[2] == 1.0)
    return nodes, elements, x0, top


def _i(v, w):
    return f"{v:>{w}d}"


def _fixed_real(x):
    """A 10-column real with a decimal point (Volume C fixed format)."""
    return f"{x:10.4f}"


def _mentat_real(x):
    """Mentat's 20-column extended real: 16 significant digits and an exponent
    without ``E`` (``6.666666666666666-1``); a negative one fills the field."""
    if x == 0.0:
        text = "0.000000000000000+0"
    else:
        exp = math.floor(math.log10(abs(x)))
        mant = abs(x) / 10.0**exp
        if round(mant, 15) >= 10.0:
            mant, exp = mant / 10.0, exp + 1
        text = f"{'-' if x < 0 else ''}{mant:.15f}{exp:+d}"
    return f"{text:>20}"


def _lines_of(values, per_line, fmt):
    out = []
    for k in range(0, len(values), per_line):
        out.append("".join(fmt(v) for v in values[k : k + per_line]))
    return out


def write_hex20(path, extended):
    nodes, elements, x0, top = hex20_model()
    iw = 10 if extended else 5
    out = [
        "title               two hex20 bricks",
        "$....a fixture written by tools/gen_marc_fixtures.py",
    ]
    if extended:
        out.append("extended")
    out += [
        f"sizing{'':14}{_i(0, iw)}{_i(2, iw)}{_i(len(nodes), iw)}{_i(0, iw)}",
        f"elements{'':12}{_i(21, iw)}",
        "end",
        "$...................",
        "connectivity",
        _i(2, iw) + _i(0, iw) + _i(1, iw),
    ]
    for ident, etype, ids in elements:
        # Id, type and 14 nodes, then the rest on a continuation line.
        out.append(_i(ident, iw) + _i(etype, iw) + "".join(_i(v, iw) for v in ids[:14]))
        out.append("".join(_i(v, iw) for v in ids[14:]))
    out += ["coordinates", _i(3, iw) + _i(len(nodes), iw) + _i(0, iw) + _i(1, iw)]
    # The nodes in an order that is not their numbering; one is defined twice
    # (the later definition wins).
    real = _mentat_real if extended else _fixed_real
    order = sorted(nodes, key=lambda n: (n % 7, n))
    out.append(_i(order[0], iw) + "".join(real(-9.0) for _ in range(3)))
    for ident in order:
        out.append(_i(ident, iw) + "".join(real(v) for v in nodes[ident]))
    # Sets. A continued line ends in C.
    x0_lines = _lines_of(x0, 6, lambda v: _i(v, iw))
    out.append(f"define{'':14}node{'':16}set{'':17}x0_face")
    for line in x0_lines[:-1]:
        out.append(line + "   c")
    out.append(x0_lines[-1])
    out.append(f"define{'':14}node{'':16}set{'':17}top_face")
    for line in _lines_of(top, 6, lambda v: _i(v, iw))[:-1]:
        out.append(line + "   c")
    out.append(_lines_of(top, 6, lambda v: _i(v, iw))[-1])
    out.append(f"define{'':14}node{'':16}set{'':17}top_not_x0")
    out.append("top_face except x0_face")
    out.append(f"define{'':14}element{'':13}set{'':17}left")
    out.append(_i(100, iw))
    out.append(f"define{'':14}element{'':13}set{'':17}right")
    out.append(_i(200, iw))
    out.append(f"define{'':14}element{'':13}set{'':17}both")
    out.append(f"{_i(100, iw)} to {_i(200, iw)} by {_i(100, iw)}")
    out.append(f"define{'':14}element{'':13}set{'':17}both_again")
    out.append("left and right")
    out += [
        "isotropic",
        _i(1, iw),
        "$ a data line of another option that starts with a set name",
        "both",
        "end option",
        "$...................",
        "loadcase            pull",
        "continue",
    ]
    path.write_text("\n".join(out) + "\n")


def write_mixed_free(path):
    out = [
        "title,mixed free-format deck",
        "sizing,0,8,15,0",
        "end",
        "connectivity",
        "7,",
        "1,7,1,2,3,4,5,6,7,8",
        "2,7,2,9,9,3,6,10,10,7",  # nodes 3 = 4 and 7 = 8 collapsed: a wedge
        "3,127,1,2,4,5,11,12,13,14,15,16",
        "4,75,5,6,7,8",
        "5,52,1,2",
        "6,116,1,2,3,4",  # a type meshio++ does not read
        "coordinates",
        "3,16,0,1",
    ]
    coords = {
        1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1),
        9: (2, 0, 0), 10: (2, 0, 1),
    }  # fmt: skip
    # tet10 on corners 1 2 4 5: mid-edges 1-2, 2-4, 4-1, 1-5, 2-5, 4-5.
    tet = [1, 2, 4, 5]
    mids = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    for k, (a, b) in enumerate(mids):
        pa, pb = coords[tet[a]], coords[tet[b]]
        coords[11 + k] = tuple(0.5 * (u + v) for u, v in zip(pa, pb))
    for ident in sorted(coords):
        x, y, z = coords[ident]
        out.append(f"{ident},{float(x)},{float(y)},{float(z)}")
    out += [
        "define,element,set,solids",
        "1,2,c",
        "3",
        "define,node,set,base",
        "1 to 4",
        "define,face,set,skin",
        "1:1",
        "end option",
    ]
    path.write_text("\n".join(out) + "\n")


def write_plane_quad8(path):
    """Two 8-node plane-strain quads (type 27), x-y coordinates only."""
    corner = {1: (0, 0), 2: (1, 0), 3: (2, 0), 4: (0, 1), 5: (1, 1), 6: (2, 1)}
    mid = {
        7: (0.5, 0),
        8: (1.5, 0),
        9: (0.5, 1),
        10: (1.5, 1),
        11: (0, 0.5),
        12: (1, 0.5),
        13: (2, 0.5),
    }
    xy = {**corner, **mid}
    out = ["title               plane strain", "end", "connectivity", "    2    0    1"]
    out.append("    1   27    1    2    5    4    7   12    9   11")
    out.append("    2   27    2    3    6    5    8   13   10   12")
    out += ["coordinates", "    2   13    0    1"]
    for ident in sorted(xy):
        x, y = xy[ident]
        out.append(f"{ident:5d}{float(x):10.3f}{float(y):10.3f}")
    out += [
        "define              element             set                 all",
        "    1 to    2",
        "end option",
    ]
    path.write_text("\n".join(out) + "\n")


# -- the post file -------------------------------------------------------------


def _e13(x):
    """Fortran ``e13.6``: `` 0.123456E+01``."""
    if x == 0.0:
        return " 0.000000E+00"
    exp = math.floor(math.log10(abs(x))) + 1
    mant = round(abs(x) / 10.0**exp, 6)
    if mant >= 1.0:
        mant, exp = mant / 10.0, exp + 1
    digits = f"{mant:.6f}"[2:]
    return f"{'-' if x < 0 else ' '}0.{digits}E{exp:+03d}"


def _ints_block(values):
    return [line for line in _lines_of(values, 6, lambda v: f"{v:13d}")]


def _reals_block(values):
    return [line for line in _lines_of(values, 6, _e13)]


def _block(number, name, body):
    return [f"=beg={number} ({name})".ljust(70)] + body + ["=end="]


# The post file's model: two 8-node bricks (type 7) sharing a face.
T19_NODES = [
    (1, (0, 0, 0)), (2, (1, 0, 0)), (3, (1, 1, 0)), (4, (0, 1, 0)),
    (5, (0, 0, 1)), (6, (1, 0, 1)), (7, (1, 1, 1)), (8, (0, 1, 1)),
    (9, (2, 0, 0)), (10, (2, 1, 0)), (11, (2, 0, 1)), (12, (2, 1, 1)),
]  # fmt: skip
T19_ELEMENTS = [(1, 7, [1, 2, 3, 4, 5, 6, 7, 8]), (2, 7, [2, 9, 10, 3, 6, 11, 12, 7])]
T19_TIMES = [0.5, 1.0]


def t19_stress(inc, element, ip):
    """The stress tensor (xx yy zz xy yz zx) at an integration point."""
    s = (inc + 1) * 10.0
    return [
        s + element + 0.1 * ip,
        -0.5 * s + ip,
        0.25 * s,
        1.5 * element,
        -0.75 * ip,
        0.125 * s,
    ]


def mises(t):
    xx, yy, zz, xy, yz, zx = t
    return math.sqrt(
        0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2)
        + 3.0 * (xy**2 + yz**2 + zx**2)
    )


def t19_displacement(inc, node):
    return [0.001 * (inc + 1) * node, -0.0005 * node, 0.0]


def write_t19(path):
    numnp, numel, nstres = len(T19_NODES), len(T19_ELEMENTS), 8
    codes = [17, 311, 312, 313, 314, 315, 316]
    lm = [len(codes), numnp, numel, 3, nstres, 9, 1, 0, 3, 8, 2, 0, 0, 12, 0, 3, 0, 0]
    lm += [0] * 12
    out = []
    out += _block(50100, "Analysis Title", ["          generated".ljust(70)])
    out += _block(50200, "Analysis Verification Data", _ints_block(lm))
    out += _block(50400, "Dummy", _ints_block([0]))
    out += _block(50500, "Domain Decomposition Information", _ints_block([0, 0]))
    out += _block(
        50600, "Element Variable Postcodes", [f"{c:13d}".ljust(37) for c in codes]
    )
    body = []
    for ident, etype, nodes in T19_ELEMENTS:
        body += _ints_block([ident, etype, len(nodes)] + nodes)
    out += _block(50700, "Element Connectivities", body)
    body = [
        f"{ident:13d}" + "".join(_e13(float(v)) for v in xyz)
        for ident, xyz in T19_NODES
    ]
    out += _block(50800, "Nodal Coordinates", body)
    out += _block(51000, "Nodal Codes and Transformation ID", _ints_block([3] * numnp))
    body = _ints_block([3])
    for name, kind, members in (
        ("fixed", 1, [1, 4, 5, 8]),
        ("loaded", 1, [9, 10, 11, 12]),
        ("second", 0, [2]),
    ):
        body.append(name.ljust(32))
        body += _ints_block([len(members), kind])
        body += _ints_block(members)
    out += _block(51301, "Set Definitions", body)
    out += _block(51501, "Flow Line Data", _ints_block([0] * 6))
    for inc, time in enumerate(T19_TIMES):
        out.append("****")
        out += _block(51600, "Loadcase Title", ["          pull".ljust(70)])
        out += _block(
            51701,
            "Integer Increment Verification Data",
            _ints_block([0, inc + 1, 0, 102, 2, 0, 0, 1, 0, 0, 0, 0]),
        )
        xlm = [time] + [0.0] * 23
        out += _block(
            51801,
            "Real Increment Verification Data",
            _ints_block([24]) + _reals_block(xlm),
        )
        body = []
        for e in range(numel):
            for ip in range(nstres):
                t = t19_stress(inc, e + 1, ip)
                body += _reals_block([mises(t)] + t)
        out += _block(52300, "Element Integration Point Values", body)
        body = _ints_block([2, 6])
        for name, ident, values in (
            (
                "Displacement",
                1,
                [v for n in range(numnp) for v in t19_displacement(inc, n + 1)],
            ),
            (
                "Reaction Force",
                5,
                [
                    (-1.0 if n in (0, 3, 4, 7) else 0.0) * (inc + 1) * (d + 1)
                    for n in range(numnp)
                    for d in range(3)
                ],
            ),
        ):
            body.append(name.ljust(48))
            body += _ints_block([ident, 0, 0, 3, 0, 0, -1, 0, 0, 0, 0, 0])
            body += _reals_block(values)
        out += _block(52401, "Nodal Results", body)
        out.append("----")
    out.append("++++")
    path.write_text("\n".join(out) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    write_hex20(OUT / "hex20.dat", extended=False)
    write_hex20(OUT / "hex20_extended.dat", extended=True)
    write_mixed_free(OUT / "mixed_free.dat")
    write_plane_quad8(OUT / "plane_quad8.dat")
    write_t19(OUT / "results.t19")
    print(f"wrote the Marc fixtures to {OUT}")


if __name__ == "__main__":
    main()
