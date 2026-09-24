#!/usr/bin/env python3
"""Regenerate the Abaqus results-file fixtures under ``tests/python/meshes/abaqus_fil/``.

No free solver writes ``.fil``, so the synthetic files are written here record by
record from the Abaqus Analysis User's Guide, "Results file output format"
(record keys and attribute layouts) and its ASCII item encoding:

* a record is ``[length, key, attributes...]`` of 8-byte words;
* binary: 512-word blocks, each one Fortran record (4-byte length markers), the
  last record of an increment (2001) padded with zeros to the block end;
  integers are 8-byte, in the file's byte order;
* ASCII: ``*`` starts a record; items are ``I`` + a two-digit digit count +
  the digits, ``D`` + a 22-character ``D22.15`` real, ``A`` + 8 characters, run
  together on 80-column lines; after 2001 the line is blank-padded to 80 columns
  and one blank line follows.

``model.fil`` (ASCII), ``model_le.fil`` and ``model_be.fil`` (binary, little- and
big-endian) hold the same model and results:

* nodes 1..27 of a 2 x 1 x 1 block (a C3D20R and a C3D8R) plus four nodes of an
  S4R on top; element 4 is a user element ``U1`` (no meshio++ cell, skipped);
  the C3D20R's node list is split over a 1900 and a 1990 record;
* sets: the node set ``1`` (a label longer than 8 characters, resolved through
  1940 as ``ASSEMBLY_PART-1-1_CLAMPED_NODES``) continued by a 1932 record, and
  the element set ``SOLIDS``;
* two increments (total times 0.5 and 1.0), each with ``S`` and ``E`` at the
  integration points (8 for the C3D20R, 1 for the C3D8R), the S4R's ``S`` at
  section points 1 and 5, ``SINV`` at the C3D8R's centroid, ``NFORC`` at its
  nodes, ``PEEQ`` averaged at nodes 1 and 2, then ``U`` at every node and ``RF``
  at the clamped nodes.

Every value is a closed-form function of its (node, element, point, component,
increment), repeated in ``expected`` below, so tests check them exactly.

``pybaqus/*.fil`` are real Abaqus 2023 ASCII results from pybaqus's test suite
(MIT, Cristóbal Tapia Camú), copied unmodified; see README.md.

    python tools/gen_abaqus_fil_fixtures.py
"""

import pathlib
import struct

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "abaqus_fil"
)


def hex20_nodes(corner_ids, edge_ids):
    return list(corner_ids) + list(edge_ids)


def model():
    """(nodes {label: xyz}, elements [(label, type, [node labels])])."""
    nodes = {}
    label = 0
    grid = {}
    # A 5 x 3 x 3 lattice on [0, 2] x [0, 1] x [0, 1] (quadratic along x in the
    # first cell); only the nodes the elements use are written.
    for k in range(3):
        for j in range(3):
            for i in range(5):
                grid[(i, j, k)] = (0.5 * i, 0.5 * j, 0.5 * k)

    used = {}

    def node(i, j, k):
        nonlocal label
        if (i, j, k) not in used:
            label += 1
            used[(i, j, k)] = label
            nodes[label] = grid[(i, j, k)]
        return used[(i, j, k)]

    c = [
        (0, 0, 0),
        (2, 0, 0),
        (2, 2, 0),
        (0, 2, 0),
        (0, 0, 2),
        (2, 0, 2),
        (2, 2, 2),
        (0, 2, 2),
    ]
    corners = [node(*p) for p in c]
    edges = [
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
    ]
    mids = [node(*[(c[a][d] + c[b][d]) // 2 for d in range(3)]) for a, b in edges]
    hex20 = hex20_nodes(corners, mids)
    c8 = [
        (2, 0, 0),
        (4, 0, 0),
        (4, 2, 0),
        (2, 2, 0),
        (2, 0, 2),
        (4, 0, 2),
        (4, 2, 2),
        (2, 2, 2),
    ]
    hex8 = [node(*p) for p in c8]
    shell = [node(2, 0, 2), node(4, 0, 2), node(4, 2, 2), node(2, 2, 2)]
    elements = [
        (1, "C3D20R", hex20),
        (2, "C3D8R", hex8),
        (3, "S4R", shell),
        (4, "U1", [hex8[0], hex8[1]]),
    ]
    return nodes, elements


def u(node, comp, inc):
    return 0.001 * node + 0.1 * comp + inc


def s(elem, point, sp, comp, inc):
    return 100.0 * elem + 10.0 * point + comp + 0.25 * sp + inc


def e(elem, point, comp, inc):
    return 1e-3 * (elem + point + comp) * inc


def rf(node, comp, inc):
    return -(node + comp) * inc


def sinv(comp, inc):
    return 7.0 * comp + inc


def nforc(node, comp, inc):
    return node * 0.5 + comp + 10.0 * inc


def peeq(node, inc):
    return 0.01 * node * inc


CLAMPED = [1, 2, 3, 4, 5, 6, 7, 8, 9]


def records():
    """Every record as ([tagged item, ...]) with its key first."""
    nodes, elements = model()
    I = lambda v: ("I", v)  # noqa: E731,E741
    D = lambda v: ("D", float(v))  # noqa: E731
    A = lambda v: ("A", v.ljust(8)[:8])  # noqa: E731
    out = []

    def rec(key, *items):
        out.append([I(key)] + list(items))

    rec(
        1921,
        A("6.23-1"),
        A("24-Sep-2"),
        A("026"),
        A("12:00:00"),
        I(len(elements)),
        I(len(nodes)),
        D(1.0),
    )
    heading = "meshio++ synthetic Abaqus results file".ljust(80)
    rec(1922, *[A(heading[8 * k : 8 * k + 8]) for k in range(10)])
    for label, etype, conn in elements:
        if etype == "C3D20R":
            rec(1900, I(label), A(etype), *[I(n) for n in conn[:10]])
            rec(1990, *[I(n) for n in conn[10:]])
        else:
            rec(1900, I(label), A(etype), *[I(n) for n in conn])
    for label in sorted(nodes):
        rec(1901, I(label), *[D(x) for x in nodes[label]])
    rec(1931, A("1"), *[I(n) for n in CLAMPED[:5]])
    rec(1932, *[I(n) for n in CLAMPED[5:]])
    rec(1933, A("SOLIDS"), I(1), I(2))
    long = "ASSEMBLY_PART-1-1_CLAMPED_NODES".ljust(80)
    rec(1940, I(1), *[A(long[8 * k : 8 * k + 8]) for k in range(10)])
    for inc in (1, 2):
        time = 0.5 * inc
        sub = "".ljust(80)
        rec(
            2000,
            D(time),
            D(time),
            D(0.0),
            D(0.0),
            I(1),
            I(1),
            I(inc),
            I(0),
            D(0.0),
            D(0.0),
            D(0.5),
            *[A(sub[8 * k : 8 * k + 8]) for k in range(10)],
        )
        rec(1911, I(0), A(""), A(""))
        for elem, points in ((1, 8), (2, 1)):
            for p in range(1, points + 1):
                rec(1, I(elem), I(p), I(0), I(0), A(""), I(3), I(3), I(0), I(0))
                rec(11, *[D(s(elem, p, 0, c, inc)) for c in range(6)])
                rec(21, *[D(e(elem, p, c, inc)) for c in range(6)])
        for sp in (1, 5):
            rec(1, I(3), I(1), I(sp), I(0), A(""), I(2), I(1), I(0), I(0))
            rec(11, *[D(s(3, 1, sp, c, inc)) for c in range(3)])
        rec(1, I(2), I(0), I(0), I(1), A(""), I(3), I(3), I(0), I(0))
        rec(12, *[D(sinv(c, inc)) for c in range(7)])
        for n in elements[1][2]:
            rec(1, I(2), I(n), I(0), I(2), A(""), I(3), I(3), I(0), I(0))
            rec(15, *[D(nforc(n, c, inc)) for c in range(3)])
        for n in (1, 2):
            rec(1, I(n), I(0), I(0), I(4), A(""), I(0), I(0), I(0), I(0))
            rec(73, D(peeq(n, inc)))
        rec(1911, I(1), A(""))
        for n in sorted(nodes):
            rec(101, I(n), *[D(u(n, c, inc)) for c in range(3)])
        for n in CLAMPED:
            rec(104, I(n), *[D(rf(n, c, inc)) for c in range(3)])
        rec(2001)
    return out


def ascii_bytes(recs):
    text = []
    line = ""

    def put(s):
        nonlocal line
        line += s
        while len(line) >= 80:
            text.append(line[:80])
            line = line[80:]

    for r in recs:
        items = [("I", len(r) + 1)] + r
        put("*")
        for tag, v in items:
            if tag == "I":
                digits = str(v)
                put(f"I{len(digits):2d}{digits}")
            elif tag == "D":
                put("D" + f"{v: .15E}".replace("E", "D"))
            else:
                put("A" + v)
        if r[0] == ("I", 2001):
            text.append(line.ljust(80))
            text.append(" " * 80)
            line = ""
    if line:
        text.append(line)
    return ("\n".join(text) + "\n").encode()


def binary_bytes(recs, order):
    words = []
    for r in recs:
        body = [len(r) + 1] + r
        rec_words = []
        for item in body:
            tag, v = item if isinstance(item, tuple) else ("I", item)
            if tag == "I":
                rec_words.append(struct.pack(order + "q", v))
            elif tag == "D":
                rec_words.append(struct.pack(order + "d", v))
            else:
                rec_words.append(v.encode())
        if r[0] == ("I", 2001):
            # Pad the increment's last record to the end of its 512-word block.
            pad = (512 - (len(words) + len(rec_words)) % 512) % 512
            rec_words[0] = struct.pack(order + "q", len(rec_words) + pad)
            rec_words += [b"\0" * 8] * pad
        words += rec_words
    words += [b"\0" * 8] * ((512 - len(words) % 512) % 512)
    assert len(words) > 1024, "the fixture must span several blocks"
    out = bytearray()
    marker = struct.pack(order + "i", 4096)
    for b in range(0, len(words), 512):
        out += marker + b"".join(words[b : b + 512]) + marker
    return bytes(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    recs = records()
    (OUT / "model.fil").write_bytes(ascii_bytes(recs))
    (OUT / "model_le.fil").write_bytes(binary_bytes(recs, "<"))
    (OUT / "model_be.fil").write_bytes(binary_bytes(recs, ">"))


if __name__ == "__main__":
    main()
