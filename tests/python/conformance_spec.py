"""The format conformance matrix: what each writable format keeps.

One canonical mesh -- every linear cell type, float64/int32/vector point and
cell data, field data and one region of each kind -- is written to and read
back from every format in ``meshioplusplus.formats()["writable"]``. What
survives is *observed* (:func:`observe`) and compared, exactly, with the
declaration in :data:`SPEC`; ``tests/python/test_conformance.py`` fails in both
directions, so a format that starts losing something and a format that starts
keeping something both show up as a diff to review. The same table renders
``doc/conformance.md`` and the "Round trip" column of ``doc/formats.md``
(``tools/gen_conformance_table.py``), and ``benchmark/bench.py`` uses
:func:`canonical_mesh` restricted to a format's cells as its input.

Observed values, per format:

``cells``
    ``{cell type: outcome}`` for every cell type the writer accepts on its own
    (a type the writer refuses is absent). ``"exact"``: the same cells, node
    for node (compared by coordinates, so point renumbering is fine);
    ``"reordered"``: the same node sets in another order; ``"as <types>"``:
    written as other cell types; ``"lost"``: silently dropped.
``points``
    ``"exact"``, ``"approx"`` (within 1e-6: single precision), ``"2d"`` (the
    z coordinate dropped), ``"subset"`` (the points of lost cells are gone,
    the rest exact), ``"none"`` (a field-only format) or ``"changed"``.
``point_data`` / ``cell_data``
    ``{array: outcome}``: ``"missing"``, or the dtype read back, suffixed
    ``"~"`` when the values are only approximately equal.
``field_data``
    ``True`` when ``fd`` survives.
``regions``
    The region kinds (``point``/``cell``/``side``) whose region came back
    under its name and kind; ``doc/regions.md`` and
    ``test_region_roundtrip.py`` check membership in depth.
``error``
    Set instead of the rest when the full mesh cannot be written or read back
    at all (with the exception type), or ``"write-only"``.

To refresh after an intended change: ``python tools/gen_conformance_table.py
--observe`` prints the observed table in :data:`SPEC`'s syntax.
"""

from __future__ import annotations

import pathlib
import tempfile

import numpy as np

import meshioplusplus

# ----------------------------------------------------------------------------
# The canonical mesh
# ----------------------------------------------------------------------------

# Disjoint cells, each positively oriented in VTK node order, on coordinates
# that are exact in binary -- except the vertex, whose 0.1 shows single
# precision up as "approx".
_BLOCKS = {
    "vertex": [[0.1, 0.2, 0.3]],
    "line": [[2, 0, 0], [3, 0, 0]],
    "triangle": [[0, 2, 0], [1, 2, 0], [0, 3, 0]],
    "quad": [[2, 2, 0], [3, 2, 0], [3, 3, 0], [2, 3, 0]],
    "tetra": [[0, 0, 4], [1, 0, 4], [0, 1, 4], [0, 0, 5]],
    "hexahedron": [
        [2, 0, 4],
        [3, 0, 4],
        [3, 1, 4],
        [2, 1, 4],
        [2, 0, 5],
        [3, 0, 5],
        [3, 1, 5],
        [2, 1, 5],
    ],
    "wedge": [[0, 2, 4], [1, 2, 4], [0, 3, 4], [0, 2, 5], [1, 2, 5], [0, 3, 5]],
    "pyramid": [[2, 2, 4], [3, 2, 4], [3, 3, 4], [2, 3, 4], [2.5, 2.5, 5]],
}
CELL_TYPES = tuple(_BLOCKS)
DATA_NAMES = {
    "point_data": ("p_f64", "p_i32", "p_vec"),
    "cell_data": ("c_f64", "c_i32"),
}


_PLANAR = ("vertex", "line", "triangle", "quad")


def canonical_mesh(cell_types=CELL_TYPES, data=True, regions=True, planar=False):
    """The canonical mesh restricted to ``cell_types`` (in canonical order).

    ``planar`` drops the z column, for the formats that only take 2-D points
    (every planar cell lies in z = 0 except the vertex, which moves onto it).
    """
    points, cells = [], []
    for ctype in CELL_TYPES:
        if ctype not in cell_types:
            continue
        corners = _BLOCKS[ctype]
        start = len(points)
        points.extend(corners)
        cells.append((ctype, np.array([range(start, start + len(corners))])))
    points = np.array(points, dtype=float)
    if planar:
        points = points[:, :2].copy()
    n = len(points)
    kw = {}
    if data:
        kw["point_data"] = {
            "p_f64": np.arange(n, dtype=np.float64) * 0.5,
            "p_i32": np.arange(n, dtype=np.int32) + 7,
            "p_vec": np.column_stack([np.arange(n) * 0.25, np.ones(n), -np.arange(n)]),
        }
        kw["cell_data"] = {
            "c_f64": [np.array([1.5 + i]) for i in range(len(cells))],
            "c_i32": [np.array([10 + i], dtype=np.int32) for i in range(len(cells))],
        }
        kw["field_data"] = {"fd": np.array([1.5, 2.5])}
    if regions:
        cell_ids = {ctype: i for i, (ctype, _) in enumerate(cells)}
        regs = [meshioplusplus.Region("rp", "point", [0, n - 1])]
        top = [c for c in ("tetra", "hexahedron", "wedge", "pyramid") if c in cell_ids]
        solid = top[0] if top else (cells[-1][0] if cells else None)
        if solid is not None:
            regs.append(meshioplusplus.Region("rc", "cell", [cell_ids[solid]]))
        if "tetra" in cell_ids:
            regs.append(
                meshioplusplus.Region("rs", "side", [[cell_ids["tetra"], 0]], dim=2)
            )
        kw["regions"] = regs
    return meshioplusplus.Mesh(points, cells, **kw)


# ----------------------------------------------------------------------------
# Observation
# ----------------------------------------------------------------------------

_EXTENSION = {
    "abaqus": ".inp",
    "ansysInp": ".inp",
    "dolfin-xml": ".xml",
    "elmer": "",
    "gmsh22": ".msh",
    "lsdyna": ".k",
    "openfoam": "",
    "vtk42": ".vtk",
    "vtk51": ".vtk",
    "z88": "",
}
# Writers whose files the registry reads under another name.
_READ_AS = {"vtk42": "vtk", "vtk51": "vtk", "gmsh22": "gmsh"}


def _target(tmp: pathlib.Path, fmt: str) -> pathlib.Path:
    from meshioplusplus._helpers import extension_to_filetypes

    ext = _EXTENSION.get(fmt)
    if ext is None:
        ext = next((e for e, t in extension_to_filetypes.items() if fmt in t), ".dat")
    if fmt == "z88":
        return tmp / "z88i1.txt"
    return tmp / f"mesh{ext}"


def round_trip(fmt: str, mesh, tmp: pathlib.Path):
    """Write ``mesh`` as ``fmt`` into ``tmp`` and read it back."""
    path = _target(tmp, fmt)
    if path.suffix == "" and fmt in ("elmer", "openfoam"):
        path.mkdir(parents=True, exist_ok=True)
    meshioplusplus.write(path, mesh, file_format=fmt)
    return meshioplusplus.read(path, file_format=_READ_AS.get(fmt, fmt))


def _key(p, nd=6):
    return tuple(round(float(x), nd) for x in p)


def _point_index(points):
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] == 0:
        # A field-only format (MFF, DEX without geometry) reads no coordinates.
        pts = np.zeros((0, 3))
    if pts.ndim == 2 and pts.shape[1] == 2:
        pts = np.column_stack([pts, np.zeros(len(pts))])
    return {_key(p, 5): i for i, p in enumerate(pts)}, pts


def _values_outcome(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    tag = str(b.dtype)
    if b.size != a.size:
        return tag + "?"
    b = b.reshape(a.shape)
    if np.array_equal(a.astype(np.float64), b.astype(np.float64)):
        return tag
    if np.allclose(a.astype(np.float64), b.astype(np.float64), atol=1e-6):
        return tag + "~"
    return tag + "?"


def observe(fmt: str, cell_types=None, planar=None) -> dict:
    """Round-trip the canonical mesh through ``fmt`` and describe what survives.

    ``cell_types`` is the set to write; ``None`` probes each type on its own
    first (what :data:`SPEC` records as the format's accepted cells). A format
    that takes no 3-D mesh at all is retried with the planar cells in 2-D
    (``planar``; the result then carries ``"input": "2d"``).
    """
    formats = meshioplusplus.formats()
    if (
        fmt not in formats["readable"]
        and _READ_AS.get(fmt, fmt) not in formats["readable"]
    ):
        return {"error": "write-only"}
    if planar is None:
        out = observe(fmt, cell_types, planar=False)
        if out.get("error") == "no cell type round-trips":
            flat = observe(fmt, cell_types, planar=True)
            if "points" in flat:
                return flat
        return out
    candidates = [t for t in CELL_TYPES if not planar or t in _PLANAR]
    errors = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        wrote_something = False
        if cell_types is None:
            cell_types = []
            for i, ctype in enumerate(candidates):
                # Bare first; then with data, which some writers require (DEX
                # writes a nodal field and nothing else).
                for data in (False, True):
                    sub = tmp / f"probe{i}{int(data)}"
                    sub.mkdir()
                    try:
                        back = round_trip(
                            fmt, canonical_mesh([ctype], data, False, planar), sub
                        )
                    except Exception as exc:
                        errors.append(type(exc).__name__)
                        continue
                    wrote_something = True
                    if sum(len(b.data) for b in back.cells):
                        cell_types.append(ctype)
                    break
            if not cell_types and wrote_something:
                # A points-only format: it takes the cells and keeps none.
                cell_types = list(candidates)
        if not cell_types:
            if errors and all(e == "ImportError" for e in errors):
                return {"error": "ImportError"}
            return {"error": "no cell type round-trips"}
        mesh = canonical_mesh(cell_types, planar=planar)
        full = tmp / "full"
        full.mkdir()
        try:
            back = round_trip(fmt, mesh, full)
        except Exception as exc:
            return {"cells": {c: "?" for c in cell_types}, "error": type(exc).__name__}
    out = _describe(mesh, back)
    if planar:
        out["input"] = "2d"
    return out


def _describe(mesh, back) -> dict:
    index, pts = _point_index(back.points)
    _, src = _point_index(mesh.points)
    # Points: nearest read point for every written point.
    if pts.shape[0] == 0:
        points = "none"
    else:
        devs = np.array([np.min(np.linalg.norm(pts - p, axis=1)) for p in src])
        flat = np.asarray(back.points).shape[1] == 2 or np.allclose(pts[:, 2], 0)
        if devs.max() == 0:
            points = "exact"
        elif devs.max() < 1e-6:
            points = "approx"
        elif flat and not np.allclose(src[:, 2], 0):
            points = "2d"
        elif (devs < 1e-6).any():
            points = "subset"
        else:
            points = "changed"

    def keys_of(m, pts_):
        out = {}
        for bi, block in enumerate(m.cells):
            for ci, row in enumerate(np.asarray(block.data)):
                out.setdefault(block.type, []).append(
                    (tuple(_key(pts_[j], 5) for j in row), bi, ci)
                )
        return out

    written = keys_of(mesh, src)
    read_p = pts if pts.size else np.zeros((0, 3))
    read = keys_of(back, read_p)
    where = {}  # written cell key (sorted) -> (read block, read cell)
    for t, rows in read.items():
        for k, bi, ci in rows:
            where[tuple(sorted(k))] = (t, k, bi, ci)
    cells = {}
    for t, rows in written.items():
        k = rows[0][0]
        hit = where.get(tuple(sorted(k)))
        if hit is None:
            # Split into other cells (e.g. quads as triangles)?
            corners = set(k)
            parts = sorted(
                {
                    rt
                    for rt, rk, _, _ in where.values()
                    if set(rk) <= corners and rt != t
                }
            )
            cells[t] = "as " + "+".join(parts) if parts else "lost"
        elif hit[0] != t:
            cells[t] = "as " + hit[0]
        else:
            cells[t] = "exact" if hit[1] == k else "reordered"

    point_data = {}
    for name in DATA_NAMES["point_data"]:
        if name not in back.point_data:
            point_data[name] = "missing"
            continue
        if pts.size == 0:  # no coordinates to match on: compare in order
            point_data[name] = _values_outcome(
                mesh.point_data[name], back.point_data[name]
            )
            continue
        order = [index.get(_key(p, 5)) for p in src]
        if any(o is None for o in order):
            point_data[name] = "missing"
            continue
        point_data[name] = _values_outcome(
            mesh.point_data[name], np.asarray(back.point_data[name])[order]
        )

    cell_data = {}
    for name in DATA_NAMES["cell_data"]:
        if name not in back.cell_data:
            cell_data[name] = "missing"
            continue
        a, b = [], []
        for t, rows in written.items():
            k, bi, ci = rows[0]
            hit = where.get(tuple(sorted(k)))
            if hit is None or hit[2] >= len(back.cell_data[name]):
                continue
            a.append(np.asarray(mesh.cell_data[name][bi])[ci])
            b.append(np.asarray(back.cell_data[name][hit[2]])[hit[3]])
        cell_data[name] = _values_outcome(a, np.array(b)) if a else "missing"

    regions = sorted(
        r.kind
        for r in mesh.regions
        if any(q.name == r.name and q.kind == r.kind for q in back.regions)
    )
    return {
        "cells": cells,
        "points": points,
        "point_data": point_data,
        "cell_data": cell_data,
        "field_data": "fd" in back.field_data,
        "regions": regions,
    }


# ----------------------------------------------------------------------------
# The declarations
# ----------------------------------------------------------------------------

SPEC: dict[str, dict] = {
    "abaqus": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point", "side"],
    },
    "ansys": {
        "cells": {
            "triangle": "lost",
            "quad": "lost",
            "tetra": "reordered",
            "hexahedron": "reordered",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
        "note": "A Fluent mesh stores volume cells and the faces that bound them; "
        "lower-dimensional cells that are not boundary faces of a volume cell are "
        "dropped with a warning.",
    },
    "ansysInp": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "avsucd": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "cae": {
        "cells": {
            "triangle": "lost",
            "quad": "lost",
            "tetra": "as triangle",
            "hexahedron": "as triangle",
            "wedge": "as triangle",
            "pyramid": "as triangle",
        },
        "points": "subset",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "The physics-ML `.npz` stores the surface of the volume cells.",
    },
    "cgns": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "code_aster": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "dex": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "lost",
            "hexahedron": "lost",
            "wedge": "lost",
            "pyramid": "lost",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "DEX stores one nodal field over a node set: cells are not kept.",
    },
    "dolfin-xml": {
        "cells": {"triangle": "lost", "tetra": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int64", "p_vec": "missing"},
        "cell_data": {"c_f64": "float64", "c_i32": "int64"},
        "field_data": False,
        "regions": [],
    },
    "elmer": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
    },
    "ensight": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "approx",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": False,
        "regions": [],
    },
    "exodus": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "febio": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64?", "c_i32": "float64?"},
        "field_data": False,
        "regions": ["cell", "point", "side"],
    },
    "femap": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "missing"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "flac3d": {
        "cells": {
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "reordered",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "flux": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "freefem": {
        "cells": {"triangle": "exact", "tetra": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "gid": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": False,
        "regions": [],
    },
    "gltf": {
        "error": "write-only",
        "note": "Written for viewers; meshio++ reads glTF only as the scene it wrote, not as "
        "a mesh round trip.",
    },
    "gmsh": {
        "cells": {
            "vertex": "?",
            "line": "?",
            "triangle": "?",
            "quad": "?",
            "tetra": "?",
            "hexahedron": "?",
            "wedge": "?",
            "pyramid": "?",
        },
        "error": "WriteError",
        "note": "The Gmsh 4.1 writer needs `gmsh:dim_tags` point data to place more than one "
        "cell type into entities, and refuses a mixed mesh without it; `gmsh22` "
        "writes the same mesh.",
    },
    "gmsh22": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": True,
        "regions": ["cell"],
    },
    "h5m": {
        "cells": {"line": "exact", "triangle": "exact", "tetra": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "hmf": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "ip": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "lost",
            "hexahedron": "lost",
            "wedge": "lost",
            "pyramid": "lost",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "An integration-point cloud: points and nodal values, no cells.",
    },
    "libmesh": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point", "side"],
    },
    "lsdyna": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point", "side"],
    },
    "mdpa": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "med": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "medit": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "mfem": {
        "cells": {
            "line": "lost",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
    },
    "mff": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "lost",
            "hexahedron": "lost",
            "wedge": "lost",
            "pyramid": "lost",
        },
        "points": "none",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "A field without geometry: one array of values, no points or cells.",
    },
    "mfm": {
        "cells": {
            "line": "?",
            "triangle": "?",
            "quad": "?",
            "tetra": "?",
            "hexahedron": "?",
            "wedge": "?",
        },
        "error": "WriteError",
        "note": "MFM holds one element type per file, so the mixed canonical mesh is refused; "
        "each type alone round-trips.",
    },
    "mphbin": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
    },
    "mphtxt": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
    },
    "nastran": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell"],
    },
    "netgen": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "neuroglancer": {
        "cells": {"triangle": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "obj": {
        "cells": {"triangle": "exact", "quad": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "off": {
        "cells": {"triangle": "exact", "quad": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "openfoam": {
        "cells": {
            "tetra": "reordered",
            "hexahedron": "reordered",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point", "side"],
    },
    "patran": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "pcd": {
        "cells": {
            "vertex": "exact",
            "line": "as vertex",
            "triangle": "as vertex",
            "quad": "as vertex",
            "tetra": "as vertex",
            "hexahedron": "as vertex",
            "wedge": "as vertex",
            "pyramid": "as vertex",
        },
        "points": "approx",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "A point cloud: every node is kept as a vertex, cells are not.",
    },
    "permas": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "ply": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "as triangle",
            "hexahedron": "as quad",
            "wedge": "as quad+triangle",
            "pyramid": "as quad+triangle",
        },
        "points": "subset",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "PLY stores faces: volume cells are written as their skin.",
    },
    "pmsh": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "exact",
            "hexahedron": "as tetra",
            "wedge": "as tetra",
            "pyramid": "as tetra",
        },
        "points": "approx",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
        "note": "The physics-ML mesh stores tetrahedra: other volume cells are simplexified, "
        "and non-volume blocks dropped.",
    },
    "pvd": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "pvtp": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "pvtu": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "stl": {
        "cells": {
            "triangle": "lost",
            "tetra": "as triangle",
            "hexahedron": "as triangle",
            "wedge": "as triangle",
            "pyramid": "as triangle",
        },
        "points": "subset",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "STL stores triangles: volume cells are written as their skin, and other "
        "blocks are dropped when volume cells are present.",
    },
    "su2": {
        "cells": {
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "svg": {"error": "write-only", "note": "A 2-D drawing, write-only."},
    "tecplot": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "as hexahedron",
            "pyramid": "as hexahedron",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "missing"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": False,
        "regions": ["cell"],
        "note": "Tecplot has no wedge or pyramid zone type: both are written as degenerate "
        "bricks and come back as hexahedra.",
    },
    "tetgen": {
        "cells": {"tetra": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "tikz": {"error": "write-only", "note": "A 2-D drawing, write-only."},
    "triangle": {
        "cells": {"triangle": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "input": "2d",
        "note": "Triangle reads and writes 2-D triangulations only, so this row is observed "
        "on the planar part of the canonical mesh.",
    },
    "ugrid": {
        "cells": {
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "unv": {
        "cells": {
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "float64"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "usd": {
        "cells": {
            "triangle": "lost",
            "quad": "lost",
            "tetra": "as triangle",
            "hexahedron": "as quad",
            "wedge": "as quad+triangle",
            "pyramid": "as quad+triangle",
        },
        "points": "subset",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "A UsdGeom mesh stores faces: volume cells are written as their skin.",
    },
    "vti": {
        "error": "no cell type round-trips",
        "note": "ImageData holds one regular hexahedral lattice; the canonical mesh is not "
        "one ([VTI](./formats/vti.md)).",
    },
    "vtk": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "vtk42": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "vtk51": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "vtkhdf": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "vtm": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "vtp": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "vtr": {
        "error": "no cell type round-trips",
        "note": "RectilinearGrid holds one axis-aligned lattice; the canonical mesh is not "
        "one.",
    },
    "vts": {
        "error": "no cell type round-trips",
        "note": "StructuredGrid holds one curvilinear lattice; the canonical mesh is not one.",
    },
    "vtu": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
    },
    "wkt": {
        "cells": {"triangle": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
    },
    "xdmf": {
        "cells": {
            "vertex": "exact",
            "line": "exact",
            "triangle": "exact",
            "quad": "exact",
            "tetra": "exact",
            "hexahedron": "exact",
            "wedge": "exact",
            "pyramid": "exact",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": False,
        "regions": [],
    },
    "xyz": {
        "cells": {
            "vertex": "exact",
            "line": "as vertex",
            "triangle": "as vertex",
            "quad": "as vertex",
            "tetra": "as vertex",
            "hexahedron": "as vertex",
            "wedge": "as vertex",
            "pyramid": "as vertex",
        },
        "points": "exact",
        "point_data": {"p_f64": "float64", "p_i32": "float64", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": [],
        "note": "A point cloud: every node is kept as a vertex, cells are not.",
    },
    "z88": {
        "cells": {"line": "exact", "tetra": "exact", "hexahedron": "exact"},
        "points": "exact",
        "point_data": {"p_f64": "missing", "p_i32": "missing", "p_vec": "missing"},
        "cell_data": {"c_f64": "missing", "c_i32": "missing"},
        "field_data": False,
        "regions": ["cell", "point"],
    },
    "zarr": {
        "cells": {
            "vertex": "lost",
            "line": "lost",
            "triangle": "lost",
            "quad": "lost",
            "tetra": "exact",
            "hexahedron": "as tetra",
            "wedge": "as tetra",
            "pyramid": "as tetra",
        },
        "points": "approx",
        "point_data": {"p_f64": "float64", "p_i32": "int32", "p_vec": "float64"},
        "cell_data": {"c_f64": "float64", "c_i32": "int32"},
        "field_data": True,
        "regions": [],
        "note": "As `pmsh`: tetrahedra only, other volume cells simplexified.",
    },
}
