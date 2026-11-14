"""Conduit Mesh Blueprint in memory (roadmap 10.1).

:func:`to_conduit` fills a ``conduit.Node`` with a mesh as the Mesh Blueprint
describes it (one domain: ``coordsets/coords``, ``topologies/mesh``,
``fields/<name>`` and ``state``), and :func:`from_conduit` reads one back. That
node is what Catalyst 2 and Ascent receive, so a simulation that already holds
one reaches every meshio++ writer through it, and the other way round.

Conduit's Python module is not on PyPI (the ``conduit`` package there is an
unrelated project); it ships with Conduit, Ascent and ParaView builds, and with
conda-forge's ``conduit``. Everything here also works on a plain nested ``dict``
of the same paths, which is what the tests use, so nothing but the import of a
real ``conduit.Node`` needs the library. See ``doc/interop.md``.
"""

import numpy as np

from ._common import warn

__all__ = ["to_conduit", "from_conduit", "has_conduit"]

#: Blueprint ``elements/shape`` names of the fixed-size meshio++ cell types; the
#: node order of each is VTK's, which is meshio++'s for these linear cells.
_SHAPES = {
    "vertex": "point",
    "line": "line",
    "triangle": "tri",
    "quad": "quad",
    "tetra": "tet",
    "hexahedron": "hex",
    "wedge": "wedge",
    "pyramid": "pyramid",
}
_TYPES = {v: k for k, v in _SHAPES.items()}
#: ``shape_map`` ids of a ``mixed`` topology: the VTK cell type numbers.
_SHAPE_IDS = {
    "point": 1,
    "line": 3,
    "tri": 5,
    "quad": 9,
    "tet": 10,
    "hex": 12,
    "wedge": 13,
    "pyramid": 14,
}
_AXES = ("x", "y", "z")
TIME_KEY = "meshio:time"


def has_conduit():
    """Whether Conduit's Python module (with ``conduit.blueprint``) is importable."""
    try:
        import conduit.blueprint  # noqa: F401
    except Exception:
        return False
    return True


def _conduit_node(op):
    try:
        import conduit
        import conduit.blueprint  # noqa: F401
    except ImportError as e:
        raise ImportError(
            f"meshio++: {op}: Conduit's Python module is not installed; it ships "
            "with Conduit, Ascent and ParaView (or conda-forge's `conduit`), or "
            "pass a dict to fill instead"
        ) from e
    return conduit.Node()


# --------------------------------------------------------------------------- #
# path access over a conduit.Node or a nested dict                            #
# --------------------------------------------------------------------------- #
def _set(node, path, value):
    if isinstance(node, dict):
        *parents, leaf = path.split("/")
        for part in parents:
            node = node.setdefault(part, {})
        node[leaf] = value
    else:
        node[path] = value


def _has(node, path):
    if isinstance(node, dict):
        for part in path.split("/"):
            if not isinstance(node, dict) or part not in node:
                return False
            node = node[part]
        return True
    return bool(node.has_path(path))


def _get(node, path):
    if isinstance(node, dict):
        for part in path.split("/"):
            node = node[part]
        return node
    return node[path]


def _names(node, path=None):
    obj = node if path is None else _get(node, path)
    if isinstance(obj, dict):
        return list(obj.keys())
    return list(obj.child_names())


def _is_object(value):
    return isinstance(value, dict) or hasattr(value, "child_names")


def _text(node, path, default=None):
    if not _has(node, path):
        return default
    value = _get(node, path)
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray) and value.dtype.kind in "SU":
        return str(value.tolist())
    return str(value)


def _array(value):
    return np.asarray(value).reshape(-1)


# --------------------------------------------------------------------------- #
# to_conduit                                                                  #
# --------------------------------------------------------------------------- #
def _plan(mesh, notes):
    """Which blocks are written and as what topology: ``(kind, blocks)``."""
    kept, polygons, polyhedra = [], [], []
    for k, block in enumerate(mesh.cells):
        if block.type in _SHAPES:
            kept.append(k)
        elif block.type.startswith("polyhedron"):
            polyhedra.append(k)
        elif block.type.startswith("polygon"):
            polygons.append(k)
        else:
            notes.append(
                f"cell block {k} ('{block.type}') has no Blueprint shape; dropped"
            )
    if polyhedra and not kept and not polygons:
        return "polyhedral", polyhedra
    if polygons and not kept and not polyhedra:
        return "polygonal", polygons
    for k in polygons + polyhedra:
        notes.append(
            f"cell block {k} ('{mesh.cells[k].type}') cannot share a topology with "
            "the other blocks; dropped"
        )
    if len({mesh.cells[k].type for k in kept}) > 1:
        return "mixed", kept
    return "fixed", kept


def _rows(block):
    """A polygon block's rows as int64 arrays."""
    return [np.asarray(r, dtype=np.int64).reshape(-1) for r in block.data]


def _offsets(sizes):
    out = np.zeros(len(sizes), dtype=np.int64)
    if len(sizes) > 1:
        np.cumsum(sizes[:-1], out=out[1:])
    return out


def _write_topology(node, mesh, kind, blocks):
    base = "topologies/mesh"
    _set(node, f"{base}/type", "unstructured")
    _set(node, f"{base}/coordset", "coords")
    if kind == "fixed":
        block = mesh.cells[blocks[0]] if blocks else None
        shape = _SHAPES[block.type] if block is not None else "point"
        conn = (
            np.concatenate(
                [np.asarray(mesh.cells[k].data, dtype=np.int64).ravel() for k in blocks]
            )
            if blocks
            else np.zeros(0, dtype=np.int64)
        )
        _set(node, f"{base}/elements/shape", shape)
        _set(node, f"{base}/elements/connectivity", conn)
        return
    if kind == "mixed":
        shapes, sizes, conn = [], [], []
        for k in blocks:
            data = np.asarray(mesh.cells[k].data, dtype=np.int64)
            shape = _SHAPES[mesh.cells[k].type]
            shapes.append(np.full(len(data), _SHAPE_IDS[shape], dtype=np.int64))
            sizes.append(np.full(len(data), data.shape[1], dtype=np.int64))
            conn.append(data.ravel())
        sizes = np.concatenate(sizes)
        _set(node, f"{base}/elements/shape", "mixed")
        for name in sorted({_SHAPES[mesh.cells[k].type] for k in blocks}):
            _set(node, f"{base}/elements/shape_map/{name}", np.int64(_SHAPE_IDS[name]))
        _set(node, f"{base}/elements/shapes", np.concatenate(shapes))
        _set(node, f"{base}/elements/sizes", sizes)
        _set(node, f"{base}/elements/offsets", _offsets(sizes))
        _set(node, f"{base}/elements/connectivity", np.concatenate(conn))
        return
    if kind == "polygonal":
        rows = [r for k in blocks for r in _rows(mesh.cells[k])]
        sizes = np.array([len(r) for r in rows], dtype=np.int64)
        _set(node, f"{base}/elements/shape", "polygonal")
        _set(node, f"{base}/elements/sizes", sizes)
        _set(node, f"{base}/elements/offsets", _offsets(sizes))
        _set(
            node,
            f"{base}/elements/connectivity",
            np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64),
        )
        return
    # polyhedral: each cell lists its faces, each face its points
    from ._ragged import polyhedron_csr

    face_conn, face_sizes, cell_faces, cell_sizes = [], [], [], []
    face_base = 0
    for k in blocks:
        flat, offsets, face_offsets = polyhedron_csr(mesh.cells[k].data)
        face_conn.append(flat)
        face_sizes.append(np.diff(offsets))
        per_cell = np.diff(face_offsets)
        cell_sizes.append(per_cell)
        cell_faces.append(face_base + np.arange(len(offsets) - 1, dtype=np.int64))
        face_base += len(offsets) - 1
    face_sizes = np.concatenate(face_sizes)
    cell_sizes = np.concatenate(cell_sizes)
    _set(node, f"{base}/elements/shape", "polyhedral")
    _set(node, f"{base}/elements/connectivity", np.concatenate(cell_faces))
    _set(node, f"{base}/elements/sizes", cell_sizes)
    _set(node, f"{base}/elements/offsets", _offsets(cell_sizes))
    _set(node, f"{base}/subelements/shape", "polygonal")
    _set(node, f"{base}/subelements/connectivity", np.concatenate(face_conn))
    _set(node, f"{base}/subelements/sizes", face_sizes)
    _set(node, f"{base}/subelements/offsets", _offsets(face_sizes))


def _write_field(node, name, association, values):
    base = f"fields/{name}"
    _set(node, f"{base}/association", association)
    _set(node, f"{base}/topology", "mesh")
    values = np.ascontiguousarray(values)
    if values.ndim == 1:
        _set(node, f"{base}/values", values)
        return
    flat = values.reshape(len(values), -1)
    labels = _AXES if flat.shape[1] <= 3 else [f"c{i}" for i in range(flat.shape[1])]
    for i in range(flat.shape[1]):
        _set(node, f"{base}/values/{labels[i]}", np.ascontiguousarray(flat[:, i]))


def to_conduit(mesh, node=None):
    """``mesh`` as a Conduit Mesh Blueprint node (one domain).

    The points are the ``explicit`` coordset ``coords``; the cells are the
    ``unstructured`` topology ``mesh``: one shape when the mesh has one fixed-size
    cell type, ``mixed`` (``shape_map`` ids are VTK's) for several, ``polygonal``
    or ``polyhedral`` (faces as ``subelements``) for a mesh of only polygons or
    only polyhedra. ``point_data`` and ``cell_data`` become ``vertex`` and
    ``element`` fields, a vector as one component child per axis;
    ``field_data["meshio:time"]`` is ``state/time``. Blocks with no Blueprint
    shape (higher-order cells, polygons beside other cells) are dropped, with a
    warning, and so are regions and other field data, which the Blueprint has no
    place for.

    ``node`` is the ``conduit.Node`` to fill (a new one by default) or a ``dict``,
    filled with the same paths as nested dicts. Returns it.
    """
    if node is None:
        node = _conduit_node("to_conduit")
    notes = []
    points = np.asarray(mesh.points, dtype=np.float64)
    points = points.reshape(len(points), -1)
    for d in range(min(3, points.shape[1])):
        _set(
            node,
            f"coordsets/coords/values/{_AXES[d]}",
            np.ascontiguousarray(points[:, d]),
        )
    _set(node, "coordsets/coords/type", "explicit")

    kind, blocks = _plan(mesh, notes)
    _write_topology(node, mesh, kind, blocks)

    from ._interop import _flatten_cell_data

    for name in sorted(mesh.point_data):
        values = np.asarray(mesh.point_data[name])
        if values.dtype.kind not in "biuf" or values.shape[:1] != (len(points),):
            notes.append(f"point_data '{name}' is not numeric per point; dropped")
            continue
        _write_field(node, name, "vertex", values)
    for name in sorted(mesh.cell_data):
        values = _flatten_cell_data(mesh, name, blocks, notes)
        if values is None:
            continue
        if values.dtype.kind not in "biuf":
            notes.append(f"cell_data '{name}' is not numeric; dropped")
            continue
        _write_field(node, name, "element", values)

    time = mesh.field_data.get(TIME_KEY) if mesh.field_data else None
    if time is not None and np.asarray(time).size:
        _set(node, "state/time", float(np.asarray(time).ravel()[0]))
    others = [k for k in (mesh.field_data or {}) if k != TIME_KEY]
    if others:
        notes.append(
            f"{len(others)} field_data entr(ies) have no Blueprint place; dropped"
        )
    if getattr(mesh, "regions", None):
        notes.append(f"{len(mesh.regions)} region(s) have no Blueprint place; dropped")
    for note in notes:
        warn(f"to_conduit: {note}")
    return node


# --------------------------------------------------------------------------- #
# from_conduit                                                                #
# --------------------------------------------------------------------------- #
def _domain(node, domain):
    """The one domain to read: the node itself, or its ``domain``-th child."""
    if _has(node, "coordsets"):
        return node
    children = [c for c in _names(node) if _is_object(_get(node, c))]
    domains = [c for c in children if _has(_get(node, c), "coordsets")]
    if not domains:
        raise ValueError("from_conduit: the node holds no Mesh Blueprint domain")
    if domain is None:
        if len(domains) > 1:
            warn(
                f"from_conduit: the node holds {len(domains)} domains; reading the "
                "first (pass domain= to choose)"
            )
        domain = 0
    return _get(node, domains[int(domain)])


def _coords(node, name):
    base = f"coordsets/{name}"
    kind = _text(node, f"{base}/type", "explicit")
    if kind == "explicit":
        axes = [a for a in _AXES if _has(node, f"{base}/values/{a}")]
        if not axes:
            axes = [a for a in ("r", "z") if _has(node, f"{base}/values/{a}")]
        if not axes:
            raise ValueError(f"from_conduit: coordset '{name}' has no values")
        columns = [
            _array(_get(node, f"{base}/values/{a}")).astype(np.float64) for a in axes
        ]
        return np.column_stack(columns), None
    if kind == "uniform":
        dims = [
            int(_array(_get(node, f"{base}/dims/{d}"))[0])
            for d in ("i", "j", "k")
            if _has(node, f"{base}/dims/{d}")
        ]
        origin = [
            (
                float(_array(_get(node, f"{base}/origin/{a}"))[0])
                if _has(node, f"{base}/origin/{a}")
                else 0.0
            )
            for a in _AXES[: len(dims)]
        ]
        spacing = [
            (
                float(_array(_get(node, f"{base}/spacing/d{a}"))[0])
                if _has(node, f"{base}/spacing/d{a}")
                else 1.0
            )
            for a in _AXES[: len(dims)]
        ]
        axes = [o + s * np.arange(n) for o, s, n in zip(origin, spacing, dims)]
    elif kind == "rectilinear":
        axes = [
            _array(_get(node, f"{base}/values/{a}")).astype(np.float64)
            for a in _AXES
            if _has(node, f"{base}/values/{a}")
        ]
    else:
        raise ValueError(f"from_conduit: coordset type '{kind}' is not read")
    grids = np.meshgrid(*axes, indexing="ij")
    # i fastest, as the Blueprint numbers the points of a grid
    points = np.column_stack([g.transpose().ravel() for g in grids])
    return points, [len(a) for a in axes]


def _grid_cells(dims):
    """Quads or hexahedra of a logically structured grid of ``dims`` points (i fastest)."""
    if len(dims) == 1:
        i = np.arange(dims[0] - 1)
        return "line", np.column_stack([i, i + 1])
    if len(dims) == 2:
        ni, nj = dims
        j, i = np.meshgrid(np.arange(nj - 1), np.arange(ni - 1), indexing="ij")
        a = (j * ni + i).ravel()
        return "quad", np.column_stack([a, a + 1, a + ni + 1, a + ni])
    ni, nj, nk = dims
    k, j, i = np.meshgrid(
        np.arange(nk - 1), np.arange(nj - 1), np.arange(ni - 1), indexing="ij"
    )
    a = (k * ni * nj + j * ni + i).ravel()
    bottom = [a, a + 1, a + ni + 1, a + ni]
    return "hexahedron", np.column_stack(bottom + [b + ni * nj for b in bottom])


def _split(array, sizes):
    return [array[o : o + s] for o, s in zip(_offsets(sizes), sizes)]


def _unstructured(node, base):
    """``(blocks, order)``: the cell blocks, and per block the element ids it holds."""
    shape = _text(node, f"{base}/elements/shape")
    conn = _array(_get(node, f"{base}/elements/connectivity")).astype(np.int64)
    if shape in _TYPES:
        width = {
            "point": 1,
            "line": 2,
            "tri": 3,
            "quad": 4,
            "tet": 4,
            "hex": 8,
            "wedge": 6,
            "pyramid": 5,
        }[shape]
        n = len(conn) // width
        return [(_TYPES[shape], conn.reshape(n, width))], [np.arange(n)]
    sizes = _array(_get(node, f"{base}/elements/sizes")).astype(np.int64)
    offsets = (
        _array(_get(node, f"{base}/elements/offsets")).astype(np.int64)
        if _has(node, f"{base}/elements/offsets")
        else _offsets(sizes)
    )
    rows = [conn[o : o + s] for o, s in zip(offsets, sizes)]
    if shape == "polygonal":
        return [("polygon", rows)], [np.arange(len(rows))]
    if shape == "polyhedral":
        sub = f"{base}/subelements"
        fconn = _array(_get(node, f"{sub}/connectivity")).astype(np.int64)
        fsizes = _array(_get(node, f"{sub}/sizes")).astype(np.int64)
        foffsets = (
            _array(_get(node, f"{sub}/offsets")).astype(np.int64)
            if _has(node, f"{sub}/offsets")
            else _offsets(fsizes)
        )
        faces = [fconn[o : o + s] for o, s in zip(foffsets, fsizes)]
        cells = [[faces[f] for f in row] for row in rows]
        return [("polyhedron", cells)], [np.arange(len(cells))]
    if shape != "mixed":
        raise ValueError(f"from_conduit: element shape '{shape}' is not read")
    names = {
        int(_array(_get(node, f"{base}/elements/shape_map/{n}"))[0]): n
        for n in _names(node, f"{base}/elements/shape_map")
    }
    shapes = _array(_get(node, f"{base}/elements/shapes")).astype(np.int64)
    blocks, order = [], []
    for sid in dict.fromkeys(shapes.tolist()):
        name = names.get(sid)
        if name not in _TYPES:
            raise ValueError(
                f"from_conduit: mixed shape id {sid} ('{name}') is not read"
            )
        ids = np.flatnonzero(shapes == sid)
        blocks.append((_TYPES[name], np.array([rows[i] for i in ids], dtype=np.int64)))
        order.append(ids)
    return blocks, order


def _field_values(node, base):
    values = _get(node, f"{base}/values")
    if not _is_object(values):
        return np.asarray(values).reshape(-1)
    columns = [
        np.asarray(_get(node, f"{base}/values/{c}")).reshape(-1)
        for c in _names(node, f"{base}/values")
    ]
    return np.column_stack(columns)


def from_conduit(node, domain=None, topology=None):
    """A ``Mesh`` from a Conduit Mesh Blueprint node (a ``conduit.Node`` or a dict).

    Reads one domain (``domain=``, the first by default for a multi-domain node)
    and one topology (``topology=``, the first by default): ``explicit``,
    ``uniform`` and ``rectilinear`` coordsets; ``unstructured`` topologies of one
    shape, ``mixed``, ``polygonal`` and ``polyhedral``, and ``uniform``,
    ``rectilinear`` and ``structured`` ones expanded to lines, quads or
    hexahedra. ``vertex`` and ``element`` fields on that topology become point
    and cell data; ``state/time`` becomes ``field_data["meshio:time"]``.
    ``matsets`` and ``adjsets`` are not read.
    """
    from ._mesh import Mesh

    node = _domain(node, domain)
    topologies = _names(node, "topologies")
    if not topologies:
        raise ValueError("from_conduit: the node has no topology")
    topo = topology if topology is not None else topologies[0]
    base = f"topologies/{topo}"
    points, dims = _coords(node, _text(node, f"{base}/coordset", "coords"))
    kind = _text(node, f"{base}/type", "unstructured")
    if kind == "unstructured":
        blocks, order = _unstructured(node, base)
    elif kind in ("uniform", "rectilinear", "structured"):
        if kind == "structured":
            dims = [
                int(_array(_get(node, f"{base}/elements/dims/{d}"))[0]) + 1
                for d in ("i", "j", "k")
                if _has(node, f"{base}/elements/dims/{d}")
            ]
        if dims is None:
            raise ValueError(f"from_conduit: topology '{topo}' has no grid dimensions")
        cell_type, cells = _grid_cells(dims)
        blocks, order = [(cell_type, cells)], [np.arange(len(cells))]
    else:
        raise ValueError(f"from_conduit: topology type '{kind}' is not read")

    point_data, cell_data = {}, {}
    if _has(node, "fields"):
        for name in _names(node, "fields"):
            fbase = f"fields/{name}"
            on = _text(node, f"{fbase}/topology", topo)
            if on != topo:
                continue
            association = _text(node, f"{fbase}/association")
            values = _field_values(node, fbase)
            if association == "vertex":
                point_data[name] = values
            elif association == "element":
                cell_data[name] = [values[ids] for ids in order]
            else:
                warn(
                    f"from_conduit: field '{name}' is not per vertex or element; skipped"
                )
    field_data = {}
    if _has(node, "state/time"):
        field_data[TIME_KEY] = np.array([float(_array(_get(node, "state/time"))[0])])
    return Mesh(
        points,
        blocks,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
