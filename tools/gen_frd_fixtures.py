#!/usr/bin/env python3
"""Regenerate the CalculiX ``.frd`` fixtures under ``tests/python/meshes/frd/``.

The ``.inp`` decks are written from the code below and solved with ``ccx`` (CalculiX
2.23 wrote the committed files); the ``.frd`` results and the ``.dat`` node prints are
copied next to them. ``ccx`` never writes the cgx shell types 7-10 nor the short
(I5) format, so ``cgx_short.frd`` is written by hand here, from the cgx manual's
record layouts. Nothing of CalculiX's own test suite (GPL) is used.

    python tools/gen_frd_fixtures.py [--ccx /usr/bin/ccx]
"""

import argparse
import pathlib
import shutil
import subprocess
import tempfile

import numpy as np

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "frd"
)

_MID = lambda p, q: (p + q) / 2  # noqa: E731


def _with_mids(corners, edges):
    c = np.asarray(corners, float)
    return np.vstack([c] + [_MID(c[a], c[b]) for a, b in edges])


HEX = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
]
WEDGE = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]
TET = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
TRI = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
QUAD = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
LINE = [[0, 0, 0], [1, 0, 0]]

# Abaqus/ccx element -> reference nodes in the element's own order.
ELEMENTS = {
    "C3D8": np.asarray(HEX, float),
    "C3D20": _with_mids(
        HEX,
        [
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
    ),
    "C3D6": np.asarray(WEDGE, float),
    "C3D15": _with_mids(
        WEDGE, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
    ),
    "C3D4": np.asarray(TET, float),
    "C3D10": _with_mids(TET, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]),
    "S3": np.asarray(TRI, float),
    "S6": _with_mids(TRI, [(0, 1), (1, 2), (2, 0)]),
    "S4": np.asarray(QUAD, float),
    "S8": _with_mids(QUAD, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "B31": np.asarray(LINE, float),
    "B32": np.vstack([LINE, [[0.5, 0, 0]]]).astype(float),
    "T3D2": np.asarray(LINE, float),
    "T3D3": np.vstack([LINE, [[0.5, 0, 0]]]).astype(float),
}
SOLIDS = {"C3D8", "C3D20", "C3D6", "C3D15", "C3D4", "C3D10"}


def _wrap_ids(ids):
    """Element card lines: at most 16 entries each, continuation lines end the previous with a comma."""
    out, first = [], True
    ids = [str(i) for i in ids]
    while ids:
        chunk = ids[: 15 if first else 16]
        ids = ids[len(chunk) :]
        out.append(("1, " if first else "") + ", ".join(chunk) + ("," if ids else ""))
        first = False
    return out


def single_element_deck(name):
    pts = ELEMENTS[name]
    n = len(pts)
    truss = name.startswith("T3D")
    lines = ["*HEADING", f"one {name} element, clamped at x=0, loaded at x=1", "*NODE"]
    lines += [f"{i + 1}, {p[0]:.6f}, {p[1]:.6f}, {p[2]:.6f}" for i, p in enumerate(pts)]
    lines.append(f"*ELEMENT, TYPE={name}, ELSET=E")
    lines += _wrap_ids(range(1, n + 1))
    lines += ["*MATERIAL, NAME=M", "*ELASTIC", "210000., 0.3", "*DENSITY", "7.8e-9"]
    if name in SOLIDS:
        lines.append("*SOLID SECTION, ELSET=E, MATERIAL=M")
    elif name.startswith("S"):
        lines += ["*SHELL SECTION, ELSET=E, MATERIAL=M", "0.1"]
    elif name.startswith("B"):
        lines += [
            "*BEAM SECTION, ELSET=E, MATERIAL=M, SECTION=RECT",
            "0.1, 0.1",
            "0.,0.,-1.",
        ]
    else:
        lines += ["*SOLID SECTION, ELSET=E, MATERIAL=M", "0.01"]
    left = [i + 1 for i, p in enumerate(pts) if abs(p[0]) < 1e-9]
    right = [i + 1 for i, p in enumerate(pts) if abs(p[0] - 1) < 1e-9]
    if name in ("S3", "S6"):
        right = [2]
    lines += ["*NSET, NSET=FIX", ", ".join(map(str, left))]
    lines += ["*NSET, NSET=LOAD", ", ".join(map(str, right))]
    lines += ["*BOUNDARY", f"FIX, 1, {3 if name in SOLIDS or truss else 6}"]
    if truss:
        lines += ["*BOUNDARY", f"{right[0]}, 2, 3"]
    lines += ["*STEP", "*STATIC", "*CLOAD", "LOAD, 1, 100."]
    if not truss:
        lines += ["LOAD, 2, 50.", "LOAD, 3, 25."]
    lines += [
        "*NODE FILE",
        "U",
        "*EL FILE",
        "S, E",
        "*NODE PRINT, NSET=LOAD",
        "U",
    ]
    if name in SOLIDS:
        # Integration-point stresses in the .dat print: a different tensor order
        # (sxx syy szz sxy sxz syz) and per-Gauss-point rather than nodal.
        lines += ["*EL PRINT, ELSET=E", "S"]
    lines.append("*END STEP")
    return "\n".join(lines) + "\n"


TWO_HEX_NODES = (
    "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 2, 0, 0\n4, 0, 1, 0\n5, 1, 1, 0\n6, 2, 1, 0\n"
    "7, 0, 0, 1\n8, 1, 0, 1\n9, 2, 0, 1\n10, 0, 1, 1\n11, 1, 1, 1\n12, 2, 1, 1\n"
)
TWO_HEX = "1, 1, 2, 5, 4, 7, 8, 11, 10\n2, 2, 3, 6, 5, 8, 9, 12, 11\n"

DECKS = {
    "mixed": """*HEADING
hex8 + wedge6 sharing a face, two static steps
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
9, 2, 0.5, 0
10, 2, 0.5, 1
*ELEMENT, TYPE=C3D8, ELSET=HEX
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D6, ELSET=WED
2, 2, 9, 3, 6, 10, 7
*MATERIAL, NAME=M
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=HEX, MATERIAL=M
*SOLID SECTION, ELSET=WED, MATERIAL=M
*NSET, NSET=FIX
1, 4, 5, 8
*NSET, NSET=LOAD
9, 10
*BOUNDARY
FIX, 1, 3
*STEP
*STATIC
0.5, 1.
*CLOAD
LOAD, 1, 100.
LOAD, 2, 40.
*NODE FILE
U
*EL FILE
S, E
*NODE PRINT, NSET=LOAD
U
*END STEP
*STEP
*STATIC
1., 1.
*CLOAD
LOAD, 3, 60.
*NODE FILE
U
*EL FILE
S
*NODE PRINT, NSET=LOAD
U
*END STEP
""",
    "freq": (
        "*HEADING\ntwo C3D8, cantilever, three modes\n"
        + TWO_HEX_NODES
        + "*ELEMENT, TYPE=C3D8, ELSET=E\n"
        + TWO_HEX
        + "*MATERIAL, NAME=M\n*ELASTIC\n210000., 0.3\n*DENSITY\n7.8e-9\n"
        "*SOLID SECTION, ELSET=E, MATERIAL=M\n*NSET, NSET=FIX\n1, 4, 7, 10\n*BOUNDARY\nFIX, 1, 3\n"
        "*STEP\n*FREQUENCY\n3\n*NODE FILE\nU\n*END STEP\n"
    ),
    "heat": (
        "*HEADING\ntwo DC3D8, transient heat transfer\n"
        + TWO_HEX_NODES
        + "*ELEMENT, TYPE=DC3D8, ELSET=E\n"
        + TWO_HEX
        + "*MATERIAL, NAME=M\n*CONDUCTIVITY\n50.\n*SPECIFIC HEAT\n5.e8\n*DENSITY\n7.8e-9\n"
        "*SOLID SECTION, ELSET=E, MATERIAL=M\n*NSET, NSET=ALL\n1,2,3,4,5,6,7,8,9,10,11,12\n"
        "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nALL, 20.\n*NSET, NSET=HOT\n1, 4, 7, 10\n"
        "*STEP\n*HEAT TRANSFER\n0.5, 2.\n*BOUNDARY\nHOT, 11, 11, 100.\n*NODE FILE\nNT\n*EL FILE\nHFL\n*END STEP\n"
    ),
}


def cantilever_deck(kind):
    """A 6 x 1 x 1 beam of 12 x 2 x 2 C3D8I elements, clamped at x=0.

    ``static`` is two load steps (the tip load doubles in the second); ``modes`` is the
    first four eigenmodes. The example notebook renders both.
    """
    nx, ny, nz = 12, 2, 2
    node = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)  # noqa: E731
    lines = ["*HEADING", f"cantilever, {nx}x{ny}x{nz} C3D8I, {kind}", "*NODE"]
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                lines.append(
                    f"{node(i, j, k)}, {6.0 * i / nx:.6f}, {j / ny:.6f}, {k / nz:.6f}"
                )
    lines.append("*ELEMENT, TYPE=C3D8I, ELSET=E")
    eid = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                n = [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k),
                ]
                n += [
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i + 1, j + 1, k + 1),
                    node(i, j + 1, k + 1),
                ]
                lines.append(f"{eid}, " + ", ".join(map(str, n)))
    lines += ["*MATERIAL, NAME=STEEL", "*ELASTIC", "210000., 0.3", "*DENSITY", "7.8e-9"]
    lines += ["*SOLID SECTION, ELSET=E, MATERIAL=STEEL"]
    fixed = [node(0, j, k) for k in range(nz + 1) for j in range(ny + 1)]
    tip = [node(nx, j, k) for k in range(nz + 1) for j in range(ny + 1)]
    lines += [
        "*NSET, NSET=FIX",
        ", ".join(map(str, fixed)),
        "*NSET, NSET=TIP",
        ", ".join(map(str, tip)),
    ]
    lines += ["*BOUNDARY", "FIX, 1, 3"]
    if kind == "static":
        for load, op in ((-20.0, "NEW"), (-40.0, "NEW")):
            lines += ["*STEP", "*STATIC", f"*CLOAD, OP={op}", f"TIP, 3, {load}"]
            lines += ["*NODE FILE", "U", "*EL FILE", "S, E", "*END STEP"]
    else:
        lines += ["*STEP", "*FREQUENCY", "4", "*NODE FILE", "U", "*END STEP"]
    return "\n".join(lines) + "\n"


def _e12(x):
    return f"{x:12.5E}"


def cgx_short():
    """A short-format (flag 0) file with the cgx shell types 7-10 and a 8-component result."""
    shapes = [
        (9, ELEMENTS["S4"], 0.0, 2, 3),  # qu4, group 2, material 3
        (10, ELEMENTS["S8"], 2.0, 0, 1),  # qu8
        (7, ELEMENTS["S3"], 4.0, 0, 1),  # tr3
        (8, ELEMENTS["S6"], 6.0, 0, 1),  # tr6
    ]
    nodes, elems = [], []
    for k, (etype, ref, x0, grp, mat) in enumerate(shapes):
        ids = []
        for p in ref:
            nodes.append((len(nodes) + 1, p[0] + x0, p[1], p[2]))
            ids.append(len(nodes))
        elems.append((k + 1, etype, grp, mat, ids))
    out = ["    1C", "    1UCGX SHORT-FORMAT FIXTURE (hand written)".ljust(80).rstrip()]
    out.append(f"    2C{len(nodes):>30}{'':37}0")
    out += [f" -1{i:5d}{_e12(x)}{_e12(y)}{_e12(z)}" for i, x, y, z in nodes]
    out.append(" -3")
    out.append(f"    3C{len(elems):>30}{'':37}0")
    for eid, etype, grp, mat, ids in elems:
        out.append(f" -1{eid:5d}{etype:5d}{grp:5d}{mat:5d}")
        for k in range(0, len(ids), 15):
            out.append(" -2" + "".join(f"{i:5d}" for i in ids[k : k + 15]))
    out.append(" -3")

    def frame(idx, value, analysis, count, blocks):
        """One ``1PSTEP`` + ``100CL`` header, then a result block, per result."""
        rows = []
        for block in blocks:
            rows += [
                f"    1PSTEP{'':25}{idx:1d}{'':11}1{'':11}{idx:1d}{'':10}",
                f"  100CL {100 + idx:4d}{value:12.9f}{count:12d}{'':20}{analysis:2d}{idx:5d}{'':10}{0:2d}",
            ]
            rows += block
        return rows

    def disp(ids):
        rows = [" -4  DISP        4    1"]
        rows += [f" -5  D{c}          1    2    {c}    0" for c in (1, 2, 3)]
        rows.append(" -5  ALL         1    2    0    0    1ALL")
        for i in ids:
            rows.append(
                f" -1{i:5d}{_e12(0.001 * i)}{_e12(-0.002 * i)}{_e12(0.5 - 0.01 * i)}"
            )
        rows.append(" -3")
        return rows

    def wide(ids):
        rows = [" -4  WIDE        8    1"]
        rows += [f" -5  W{c}          1    1    0    0" for c in range(1, 9)]
        for i in ids:
            v = [i * 10.0 + c - 3.5 for c in range(8)]
            rows.append(f" -1{i:5d}" + "".join(_e12(x) for x in v[:6]))
            rows.append(" -2" + " " * 5 + "".join(_e12(x) for x in v[6:]))
        rows.append(" -3")
        return rows

    def temp(ids):
        rows = [" -4  NDTEMP      1    1", " -5  T           1    1    0    0"]
        rows += [f" -1{i:5d}{_e12(20.0 + i)}" for i in ids]
        rows.append(" -3")
        return rows

    every = [n[0] for n in nodes]
    out += frame(1, 0.5, 1, len(every), [disp(every), wide(every)])
    out += frame(2, 1.5, 1, 6, [temp(every[:6])])
    out.append("  9999")
    return "\n".join(out) + "\n"


# Decks that also get a binary (*NODE OUTPUT / *ELEMENT OUTPUT) variant: a linear
# hex (c3d8, cross-checked byte-for-byte against its ASCII sibling), a quadratic hex
# (c3d20, so the element-record node-count skip is exercised beyond the simplest
# case) and the real multi-element cantilever (more than one node/element record, a
# non-trivial numnod). *NODE OUTPUT/*ELEMENT OUTPUT mirror *NODE FILE/*EL FILE's own
# syntax -- verified against ccx 2.23 -- except the keyword is "*ELEMENT OUTPUT", not
# "*EL OUTPUT" (ccx 2.23 does not recognise that spelling and silently drops the
# card with only a warning, not a refusal).
BINARY_DECKS = ("c3d8", "c3d20", "cantilever_static")


def _binary_variant(text):
    return text.replace("*NODE FILE\n", "*NODE OUTPUT\n").replace(
        "*EL FILE\n", "*ELEMENT OUTPUT\n"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ccx", default=shutil.which("ccx"), help="path to the ccx solver")
    args = ap.parse_args()
    decks = OUT / "decks"
    decks.mkdir(parents=True, exist_ok=True)
    named = {name.lower(): single_element_deck(name) for name in ELEMENTS}
    named.update(DECKS)
    named["cantilever_static"] = cantilever_deck("static")
    named["cantilever_modes"] = cantilever_deck("modes")
    for name, text in named.items():
        (decks / f"{name}.inp").write_text(text)
    for name in BINARY_DECKS:
        (decks / f"{name}_bin.inp").write_text(_binary_variant(named[name]))
    (OUT / "cgx_short.frd").write_text(cgx_short())
    if not args.ccx:
        print("ccx not found: decks and cgx_short.frd written, results not regenerated")
        return
    for name in named:
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(decks / f"{name}.inp", tmp)
            run = subprocess.run(
                [args.ccx, name], cwd=tmp, capture_output=True, text=True
            )
            frd = pathlib.Path(tmp) / f"{name}.frd"
            if run.returncode != 0 or not frd.exists():
                raise SystemExit(
                    f"ccx failed on {name}:\n{run.stdout[-800:]}{run.stderr[-800:]}"
                )
            shutil.copy(frd, OUT / f"{name}.frd")
            dat = pathlib.Path(tmp) / f"{name}.dat"
            if "NODE PRINT" in named[name] and dat.exists():
                shutil.copy(dat, OUT / f"{name}.dat")
            elif (OUT / f"{name}.dat").exists():
                (OUT / f"{name}.dat").unlink()
    for name in BINARY_DECKS:
        bname = f"{name}_bin"
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(decks / f"{bname}.inp", tmp)
            run = subprocess.run(
                [args.ccx, bname], cwd=tmp, capture_output=True, text=True
            )
            frd = pathlib.Path(tmp) / f"{bname}.frd"
            if run.returncode != 0 or not frd.exists():
                raise SystemExit(
                    f"ccx failed on {bname}:\n{run.stdout[-800:]}{run.stderr[-800:]}"
                )
            shutil.copy(frd, OUT / f"{bname}.frd")


if __name__ == "__main__":
    main()
