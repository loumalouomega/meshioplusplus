#!/usr/bin/env python3
"""Regenerate the Tecplot binary fixtures under ``tests/python/meshes/tecplot/plt/``.

Each fixture is a pair: an ASCII ``<name>.dat`` written here and a binary
``<name>.plt`` with the same content, written by Tecplot's own TecIO library
(the classic ``TECINI142``/``TECZNE142``/``TECDAT142``/``TECNOD142`` API, which
is what ``preplot`` is built on). The reader's acceptance test is that each
``.plt`` converts identically to its ``.dat``, so the binary reader is never
checked against itself.

TecIO is not a dependency of meshio++. Build it once from a Tecplot
distribution or from SU2's ``externals/tecio/teciosrc`` (Tecplot's royalty-free
licence, see its ``tecio_license_agreement.txt``)::

    cmake -S teciosrc -B build -DCMAKE_POLICY_VERSION_MINIMUM=3.5
    cmake --build build --target tecio -j2

then point this script at the source directory (for ``TECIO.h``) and the
static library::

    TECIO_SRC=.../teciosrc TECIO_LIB=.../build/libtecio.a \\
        python tools/gen_tecplot_plt_fixtures.py

The ASCII files alone need no TecIO (``--ascii-only``).

Fixtures:

* ``fe_mixed``: a 3-D tetrahedron zone and a brick zone, nodal ``T`` and a
  cell-centred ``P``, double precision, plus every header record the reader
  must skip (dataset, zone and variable auxiliary data, text, a polyline
  geometry, custom labels and a user record).
* ``fe_2d``: single precision triangle, quadrilateral and line-segment zones; a
  passive variable, variables shared from an earlier zone, shared connectivity
  and a user-defined face-neighbour connection; written in the foreign
  (big-endian) byte order.
* ``ordered``: I-, IJ- and IJK-ordered zones with nodal and cell-centred
  variables.
* ``transient``: two strands over three solution times, later zones sharing
  the coordinates and connectivity of the first.
"""

import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "tests" / "python" / "meshes" / "tecplot" / "plt"

ZONETYPE = {
    "ORDERED": 0,
    "FELINESEG": 1,
    "FETRIANGLE": 2,
    "FEQUADRILATERAL": 3,
    "FETETRAHEDRON": 4,
    "FEBRICK": 5,
}
NODES_PER_CELL = {
    "FELINESEG": 2,
    "FETRIANGLE": 3,
    "FEQUADRILATERAL": 4,
    "FETETRAHEDRON": 4,
    "FEBRICK": 8,
}


def _values(n, seed, single):
    rng = np.random.default_rng(seed)
    v = np.round(rng.uniform(-10.0, 10.0, n), 3)
    return v.astype(np.float32).astype(np.float64) if single else v


def _zone(title, ztype, nvar, **kw):
    z = {
        "title": title,
        "type": ztype,
        "loc": [1] * nvar,  # 1 nodal, 0 cell-centred (TecIO's convention)
        "passive": [0] * nvar,
        "share": [0] * nvar,  # 1-based source zone, 0 = own data
        "connshare": 0,
        "time": None,
        "strand": 0,
        "data": {},
        "conn": None,
        "faces": [],  # local one-to-one face neighbours (cz, fz, nz), 1-based
        "aux": [],
    }
    z.update(kw)
    return z


def _n_nodes(z):
    if z["type"] == "ORDERED":
        return z["I"] * z["J"] * z["K"]
    return z["nodes"]


def _n_cells(z):
    if z["type"] == "ORDERED":
        n = 1
        for d in (z["I"], z["J"], z["K"]):
            if d > 1:
                n *= d - 1
        return n
    return len(z["conn"])


def _fill(ds, seed=0):
    """Gives every owned, active variable of every zone its values."""
    single = not ds["double"]
    for zi, z in enumerate(ds["zones"]):
        for v in range(len(ds["variables"])):
            if z["share"][v] or z["passive"][v] or v in z["data"]:
                continue
            n = _n_nodes(z) if z["loc"][v] else _n_cells(z)
            z["data"][v] = _values(n, seed + 100 * zi + v, single)


# --- datasets -----------------------------------------------------------------


def fe_mixed():
    variables = ["X", "Y", "Z", "T", "P"]
    tet_xyz = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    hex_xyz = np.array(
        [[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (2.0, 3.0, 4.0)]
    )
    # nodes ordered x fastest: 6 per layer -> two hexes side by side
    hconn = []
    for c in range(2):
        b = [c, c + 1, c + 4, c + 3]
        hconn.append([n + 1 for n in b] + [n + 7 for n in b])
    tet = _zone(
        "tets",
        "FETETRAHEDRON",
        5,
        nodes=5,
        conn=[[1, 2, 3, 4], [2, 3, 4, 5]],
        loc=[1, 1, 1, 1, 0],
        aux=[("Material", "steel")],
    )
    tet["data"] = {0: tet_xyz[:, 0], 1: tet_xyz[:, 1], 2: tet_xyz[:, 2]}
    brick = _zone(
        "bricks",
        "FEBRICK",
        5,
        nodes=12,
        conn=hconn,
        loc=[1, 1, 1, 1, 0],
    )
    brick["data"] = {0: hex_xyz[:, 0], 1: hex_xyz[:, 1], 2: hex_xyz[:, 2]}
    ds = {
        "title": "FE mixed",
        "variables": variables,
        "zones": [tet, brick],
        "double": True,
        "foreign": False,
        "extras": True,
    }
    _fill(ds, 1)
    return ds


def fe_2d():
    variables = ["X", "Y", "U", "Q"]
    xy = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]], dtype=float)
    tri = _zone(
        "triangles",
        "FETRIANGLE",
        4,
        nodes=6,
        conn=[[1, 2, 3], [1, 3, 4], [2, 5, 6], [2, 6, 3]],
        loc=[1, 1, 1, 0],
        faces=[(1, 1, 2), (2, 3, 1)],
    )
    tri["data"] = {0: xy[:, 0], 1: xy[:, 1]}
    quad = _zone(
        "quads",
        "FEQUADRILATERAL",
        4,
        nodes=6,
        conn=[[1, 2, 3, 4], [2, 5, 6, 3]],
        loc=[1, 1, 1, 0],
        share=[1, 1, 0, 0],
        passive=[0, 0, 0, 1],
    )
    quad2 = _zone(
        "quads_again",
        "FEQUADRILATERAL",
        4,
        nodes=6,
        conn=None,
        loc=[1, 1, 1, 0],
        share=[1, 1, 0, 0],
        connshare=2,
    )
    quad2["conn"] = quad["conn"]  # only for the cell count
    line = _zone(
        "edges",
        "FELINESEG",
        4,
        nodes=3,
        conn=[[1, 2], [2, 3]],
        loc=[1, 1, 1, 1],
    )
    line["data"] = {0: np.array([0.0, 1.0, 2.0]), 1: np.array([-1.0, -1.0, -1.0])}
    ds = {
        "title": "FE 2D",
        "variables": variables,
        "zones": [tri, quad, quad2, line],
        "double": False,
        "foreign": True,
        "extras": False,
    }
    _fill(ds, 2)
    return ds


def ordered():
    variables = ["X", "Y", "Z", "V", "C"]
    zones = []
    for title, (ni, nj, nk) in (
        ("i_line", (5, 1, 1)),
        ("ij_plane", (4, 3, 1)),
        ("ijk_block", (3, 3, 2)),
    ):
        i, j, k = np.meshgrid(
            np.arange(ni), np.arange(nj), np.arange(nk), indexing="ij"
        )
        # I fastest, then J, then K
        x = (1.5 * i).transpose(2, 1, 0).ravel() + 0.25 * j.transpose(2, 1, 0).ravel()
        y = (1.0 * j).transpose(2, 1, 0).ravel()
        z = (0.5 * k).transpose(2, 1, 0).ravel()
        zone = _zone(title, "ORDERED", 5, I=ni, J=nj, K=nk, loc=[1, 1, 1, 1, 0])
        zone["data"] = {0: x, 1: y, 2: z}
        zones.append(zone)
    ds = {
        "title": "Ordered",
        "variables": variables,
        "zones": zones,
        "double": True,
        "foreign": False,
        "extras": False,
    }
    _fill(ds, 3)
    return ds


def transient():
    variables = ["X", "Y", "S"]
    xy = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
    conn = [[1, 2, 5, 4], [2, 3, 6, 5]]
    zones = []
    for step, t in enumerate((0.0, 0.5, 1.25)):
        for strand in (1, 2):
            first = len(zones) == 0
            z = _zone(
                f"strand{strand}_t{step}",
                "FEQUADRILATERAL",
                3,
                nodes=6,
                conn=conn,
                time=t,
                strand=strand,
                share=[0, 0, 0] if first else [1, 1, 0],
                connshare=0 if first else 1,
            )
            if first:
                z["data"] = {0: xy[:, 0], 1: xy[:, 1]}
            zones.append(z)
    ds = {
        "title": "Transient",
        "variables": variables,
        "zones": zones,
        "double": True,
        "foreign": False,
        "extras": False,
    }
    _fill(ds, 4)
    return ds


DATASETS = {
    "fe_mixed": fe_mixed,
    "fe_2d": fe_2d,
    "ordered": ordered,
    "transient": transient,
}


# --- ASCII --------------------------------------------------------------------


def _ranges(idx):
    """1-based variable numbers -> ``[1-3,5]``."""
    idx = sorted(idx)
    parts, start = [], None
    for k, v in enumerate(idx):
        if start is None:
            start = v
        if k + 1 == len(idx) or idx[k + 1] != v + 1:
            parts.append(str(start) if start == v else f"{start}-{v}")
            start = None
    return "[" + ",".join(parts) + "]"


def _fmt(v):
    return repr(float(v))


def write_ascii(ds, path):
    nv = len(ds["variables"])
    out = [f'TITLE = "{ds["title"]}"']
    out.append("VARIABLES = " + " ".join(f'"{v}"' for v in ds["variables"]))
    for z in ds["zones"]:
        head = [f'ZONE T="{z["title"]}"']
        if z["type"] == "ORDERED":
            head.append(f'I={z["I"]}, J={z["J"]}, K={z["K"]}, ZONETYPE=ORDERED')
        else:
            head.append(
                f'NODES={z["nodes"]}, ELEMENTS={_n_cells(z)}, ZONETYPE={z["type"]}'
            )
        head.append("DATAPACKING=BLOCK")
        cc = [v + 1 for v in range(nv) if z["loc"][v] == 0]
        if cc:
            head.append(f"VARLOCATION=({_ranges(cc)}=CELLCENTERED)")
        shared = {}
        for v in range(nv):
            if z["share"][v]:
                shared.setdefault(z["share"][v], []).append(v + 1)
        if shared:
            head.append(
                "VARSHARELIST=("
                + ", ".join(f"{_ranges(vs)}={src}" for src, vs in shared.items())
                + ")"
            )
        passive = [v + 1 for v in range(nv) if z["passive"][v]]
        if passive:
            head.append(f"PASSIVEVARLIST=({_ranges(passive)})")
        if z["connshare"]:
            head.append(f'CONNECTIVITYSHAREZONE={z["connshare"]}')
        if z["time"] is not None:
            head.append(f'SOLUTIONTIME={_fmt(z["time"])}, STRANDID={z["strand"]}')
        if z["faces"]:
            head.append(
                f'FACENEIGHBORMODE=LOCALONETOONE, FACENEIGHBORCONNECTIONS={len(z["faces"])}'
            )
        out.append(", ".join(head))
        for v in range(nv):
            if v in z["data"]:
                vals = z["data"][v]
                for s in range(0, len(vals), 5):
                    out.append(" ".join(_fmt(x) for x in vals[s : s + 5]))
        if z["type"] != "ORDERED" and not z["connshare"]:
            for c in z["conn"]:
                out.append(" ".join(str(n) for n in c))
            for f in z["faces"]:
                out.append(" ".join(str(n) for n in f))
    path.write_text("\n".join(out) + "\n")


# --- TecIO driver ---------------------------------------------------------------


def _c_array(ctype, name, vals):
    body = ", ".join(
        repr(float(v)) if ctype in ("double", "float") else str(int(v)) for v in vals
    )
    return f"    static const {ctype} {name}[] = {{{body}}};\n"


def driver_source(ds, plt_name):
    nv = len(ds["variables"])
    vis_double = 1 if ds["double"] else 0
    ctype = "double" if ds["double"] else "float"
    src = [
        '#include <cstdio>\n#include <cstdlib>\n#include "TECIO.h"\n',
        '#define CHECK(x) do { if ((x) != 0) { std::fprintf(stderr, "fail: %s\\n", #x); std::exit(1); } } while (0)\n',
        "int main() {\n",
        f"    INTEGER4 fileFormat = 0, fileType = 0, debug = 0, visDouble = {vis_double};\n",
    ]
    if ds["foreign"]:
        src.append("    INTEGER4 foreign = 1; TECFOREIGN142(&foreign);\n")
    src.append(
        f'    CHECK(TECINI142("{ds["title"]}", "{" ".join(ds["variables"])}", "{plt_name}", ".", '
        "&fileFormat, &fileType, &debug, &visDouble));\n"
    )
    if ds["extras"]:
        src.append('    CHECK(TECAUXSTR142("Common.Author", "meshio++ fixture"));\n')
    for zi, z in enumerate(ds["zones"]):
        zt = ZONETYPE[z["type"]]
        if z["type"] == "ORDERED":
            imx, jmx, kmx = z["I"], z["J"], z["K"]
        else:
            imx, jmx, kmx = z["nodes"], _n_cells(z), 0
        t = z["time"] if z["time"] is not None else 0.0
        src.append("    {\n")
        src.append(
            f"        INTEGER4 zt = {zt}, imx = {imx}, jmx = {jmx}, kmx = {kmx}, zero = 0, one = 1;\n"
        )
        src.append(
            f"        INTEGER4 strand = {z['strand']}, parent = 0, connshare = {z['connshare']};\n"
        )
        src.append(f"        INTEGER4 nfc = {len(z['faces'])}, fnmode = 0;\n")
        src.append(f"        double t = {t!r};\n")
        src.append(
            "        INTEGER4 passive[] = {"
            + ", ".join(map(str, z["passive"]))
            + "};\n"
        )
        src.append(
            "        INTEGER4 loc[] = {" + ", ".join(map(str, z["loc"])) + "};\n"
        )
        src.append(
            "        INTEGER4 share[] = {" + ", ".join(map(str, z["share"])) + "};\n"
        )
        src.append(
            f'        CHECK(TECZNE142("{z["title"]}", &zt, &imx, &jmx, &kmx, &zero, &zero, &zero, &t, '
            "&strand, &parent, &one, &nfc, &fnmode, &zero, &zero, &zero, passive, loc, share, "
            "&connshare));\n"
        )
        for name, value in z["aux"]:
            src.append(f'        CHECK(TECZAUXSTR142("{name}", "{value}"));\n')
        for v in range(nv):
            if v not in z["data"]:
                continue
            vals = z["data"][v]
            src.append("    " + _c_array(ctype, f"d{zi}_{v}", vals))
            src.append(
                f"        {{ INTEGER4 n = {len(vals)}; CHECK(TECDAT142(&n, d{zi}_{v}, &visDouble)); }}\n"
            )
        if z["type"] != "ORDERED" and not z["connshare"]:
            flat = [n for c in z["conn"] for n in c]
            src.append("    " + _c_array("INTEGER4", f"c{zi}", flat))
            src.append(f"        CHECK(TECNOD142(c{zi}));\n")
            if z["faces"]:
                flat = [n for f in z["faces"] for n in f]
                src.append("    " + _c_array("INTEGER4", f"f{zi}", flat))
                src.append(f"        CHECK(TECFACE142(f{zi}));\n")
        src.append("    }\n")
    if ds["extras"]:
        src.append(
            """    {
        INTEGER4 var = 4;
        CHECK(TECVAUXSTR142(&var, "Units", "Pa"));
        double x = 0.5, y = 0.5, z = 0.0, h = 12.0, margin = 20.0, lt = 0.1, angle = 0.0, ls = 1.0;
        INTEGER4 pos = 0, attach = 0, zone = 1, font = 1, units = 2, box = 0, bc = 0, bf = 7, anchor = 0, tc = 0, scope = 1, clip = 0;
        CHECK(TECTXT142(&x, &y, &z, &pos, &attach, &zone, &font, &units, &h, &box, &margin, &lt,
                        &bc, &bf, &angle, &anchor, &ls, &tc, &scope, &clip, "fixture text", ""));
        INTEGER4 color = 0, fill = 7, filled = 0, gtype = 0, pattern = 0, nell = 72, astyle = 0, aattach = 0;
        double plen = 2.0, thick = 0.1, asize = 5.0, aangle = 12.0;
        INTEGER4 nseg = 1, segpts[] = {3};
        static const float gx[] = {0.0f, 1.0f, 2.0f}, gy[] = {0.0f, 1.0f, 0.0f}, gz[] = {0.0f, 0.0f, 0.0f};
        INTEGER4 gpos = 0;
        CHECK(TECGEO142(&x, &y, &z, &gpos, &attach, &zone, &color, &fill, &filled, &gtype, &pattern,
                        &plen, &thick, &nell, &astyle, &aattach, &asize, &aangle, &scope, &clip,
                        &nseg, segpts, gx, gy, gz, ""));
        CHECK(TECLAB142("\\"low\\" \\"high\\""));
        CHECK(TECUSR142("user record"));
    }
"""
        )
    src.append("    CHECK(TECEND142());\n    return 0;\n}\n")
    return "".join(src)


def write_plt(ds, path, tecio_src, tecio_lib):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        (tmp / "drv.cpp").write_text(driver_source(ds, path.name))
        cxx = os.environ.get("CXX", "c++")
        subprocess.run(
            [
                cxx,
                "-O0",
                "-I",
                str(tecio_src),
                "drv.cpp",
                str(tecio_lib),
                "-lpthread",
                "-o",
                "drv",
            ],
            cwd=tmp,
            check=True,
        )
        subprocess.run(["./drv"], cwd=tmp, check=True, stdout=subprocess.DEVNULL)
        shutil.copyfile(tmp / path.name, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ascii-only", action="store_true")
    ap.add_argument("--out", type=pathlib.Path, default=OUT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    tecio_src = os.environ.get("TECIO_SRC")
    tecio_lib = os.environ.get("TECIO_LIB")
    if not args.ascii_only and not (tecio_src and tecio_lib):
        sys.exit("set TECIO_SRC and TECIO_LIB, or pass --ascii-only")
    for name, make in DATASETS.items():
        ds = make()
        write_ascii(ds, args.out / f"{name}.dat")
        if not args.ascii_only:
            write_plt(ds, args.out / f"{name}.plt", tecio_src, tecio_lib)
        print("wrote", name)


if __name__ == "__main__":
    main()
