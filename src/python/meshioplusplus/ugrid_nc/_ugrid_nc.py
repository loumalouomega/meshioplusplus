"""UGRID-conventions NetCDF meshes: the Python reference of ``formats/ugrid_nc.cpp``.

See that header (and ``doc/formats/ugrid_nc.md``) for the layout: one 1-D or 2-D
``mesh_topology``, faces grouped by node count (triangle, quad, then each larger
count), edges as a trailing ``line`` block, data by ``location``, a leading time
dimension as steps. Through ``netCDF4``.
"""

from __future__ import annotations

import re

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import CellBlock, Mesh

NAME_ATT = "meshioplusplus_name"
_FILL_INT = -2147483647
_FILL_INT64 = -9223372036854775806


def _netcdf4():
    try:
        import netCDF4
    except ImportError as exc:
        raise ImportError(
            "meshio++: ugrid_nc needs netCDF4 (pip install netCDF4)"
        ) from exc
    return netCDF4


def _att(var, name):
    value = getattr(var, name, None) if name in var.ncattrs() else None
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    if isinstance(value, str):
        return value.rstrip("\x00 \t\r\n")
    return value


def _int_att(var, name, default=None):
    value = _att(var, name)
    if value is None:
        return default
    try:
        return int(np.asarray(value).reshape(-1)[0])
    except (TypeError, ValueError):
        return default


def _var(nc, name):
    if name not in nc.variables:
        raise ReadError(
            f"meshio++: ugrid_nc: the mesh names a variable '{name}' the file does not have"
        )
    return nc.variables[name]


def _connectivity(nc, name, along, num_nodes):
    var = _var(nc, name)
    if var.ndim != 2:
        raise ReadError(f"meshio++: ugrid_nc: '{name}' must be two-dimensional")
    var.set_auto_mask(False)
    raw = np.asarray(var[:], dtype=np.int64)
    if along is not None and var.dimensions[1] == along and var.dimensions[0] != along:
        raw = raw.T
    start = _int_att(var, "start_index", 0)
    fill = _int_att(var, "_FillValue")
    rows = []
    for r, row in enumerate(raw):
        keep = (row != _FILL_INT) & (row != _FILL_INT64) & (row >= start)
        if fill is not None:
            keep &= row != fill
        nodes = row[keep] - start
        if len(nodes) and nodes.max() >= num_nodes:
            raise ReadError(
                f"meshio++: ugrid_nc: '{name}' names node {nodes.max() + start} of "
                f"{num_nodes}"
            )
        rows.append(nodes)
    return rows


def _scan(nc):
    topo = None
    for name, var in nc.variables.items():
        if _att(var, "cf_role") == "mesh_topology":
            topo = name
            break
    if topo is None:
        if "num_nodes" in nc.dimensions and "num_dim" in nc.dimensions:
            raise ReadError(
                "meshio++: ugrid_nc: this netCDF file is an Exodus mesh; read it with "
                "file_format='exodus'"
            )
        raise ReadError(
            'meshio++: ugrid_nc: no variable has cf_role = "mesh_topology"; this netCDF '
            "file is not a UGRID mesh"
        )
    tvar = nc.variables[topo]
    if _int_att(tvar, "topology_dimension", 2) == 3:
        raise ReadError(
            f"meshio++: ugrid_nc: '{topo}' is a volume (3-D) topology; meshio++ reads "
            "1-D and 2-D UGRID meshes"
        )
    coords = (_att(tvar, "node_coordinates") or "").split()
    if len(coords) < 2:
        raise ReadError(
            f"meshio++: ugrid_nc: '{topo}' needs at least two node_coordinates"
        )
    node_dim = _var(nc, coords[0]).dimensions[0]
    faces = _att(tvar, "face_node_connectivity")
    face_dim = None
    if faces:
        face_dim = _att(tvar, "face_dimension") or _var(nc, faces).dimensions[0]
    edges = _att(tvar, "edge_node_connectivity")
    edge_dim = None
    if edges:
        edge_dim = _att(tvar, "edge_dimension") or _var(nc, edges).dimensions[0]
    time_dim = None
    for var in nc.variables.values():
        if _att(var, "mesh") != topo or not _att(var, "location"):
            continue
        dims = var.dimensions
        if len(dims) >= 2 and dims[0] not in (node_dim, face_dim, edge_dim):
            time_dim = dims[0]
            break
    times = []
    if time_dim is not None:
        n = len(nc.dimensions[time_dim])
        if time_dim in nc.variables and n:
            times = [
                float(t) for t in np.asarray(nc.variables[time_dim][:]).reshape(-1)
            ]
        else:
            times = [float(k) for k in range(n)]
    return topo, coords, node_dim, face_dim, edge_dim, time_dim, times


def _resolve_step(times, time_step):
    n = len(times)
    if n == 0:
        if time_step in (0, -1):
            return 0
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    k = time_step + n if time_step < 0 else time_step
    if not 0 <= k < n:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {n} "
            f"step{'' if n == 1 else 's'}"
        )
    return k


def _pad(like, rows):
    like = np.asarray(like)
    out = np.empty((rows,) + like.shape[1:], dtype=like.dtype)
    out[...] = np.nan if like.dtype.kind == "f" else 0
    return out


def time_values(filename):
    nc4 = _netcdf4()
    with nc4.Dataset(str(filename), "r") as nc:
        return _scan(nc)[6]


def read(filename, time_step=0):
    nc4 = _netcdf4()
    try:
        nc = nc4.Dataset(str(filename), "r")
    except OSError as exc:
        raise ReadError(
            f"meshio++: ugrid_nc: cannot open '{filename}': {exc}"
        ) from None
    with nc:
        topo, coords, node_dim, face_dim, edge_dim, time_dim, times = _scan(nc)
        tvar = nc.variables[topo]
        nn = len(nc.dimensions[node_dim])
        points = np.zeros((nn, 3))
        for c, name in enumerate(coords[:3]):
            var = _var(nc, name)
            var.set_auto_mask(False)
            points[:, c] = np.asarray(var[:], dtype=np.float64).reshape(-1)
        cells, face_blocks = [], []
        faces = _att(tvar, "face_node_connectivity")
        rows = _connectivity(nc, faces, face_dim, nn) if faces else []
        by_count = {}
        for r, row in enumerate(rows):
            if len(row) < 3:
                raise ReadError(f"meshio++: ugrid_nc: face {r} has {len(row)} nodes")
            by_count.setdefault(len(row), []).append(r)
        for count in sorted(by_count):
            idx = by_count[count]
            data = np.asarray([rows[r] for r in idx], dtype=np.int64)
            kind = {3: "triangle", 4: "quad"}.get(count, "polygon")
            cells.append(CellBlock(kind, data))
            face_blocks.append(np.asarray(idx, dtype=np.int64))
        num_edges = 0
        edges = _att(tvar, "edge_node_connectivity")
        if edges:
            erows = _connectivity(nc, edges, edge_dim, nn)
            num_edges = len(erows)
            if num_edges:
                if any(len(r) != 2 for r in erows):
                    raise ReadError(
                        "meshio++: ugrid_nc: an edge does not have two nodes"
                    )
                cells.append(CellBlock("line", np.asarray(erows, dtype=np.int64)))
        step = _resolve_step(times, time_step)
        field_data = {}
        if times:
            field_data["meshio:time"] = np.array([times[step]])
        point_data, face_data, edge_data = {}, {}, {}
        locations = {"node": node_dim, "face": face_dim, "edge": edge_dim}
        for vname, var in nc.variables.items():
            if _att(var, "mesh") != topo:
                continue
            location = _att(var, "location")
            loc_dim = locations.get(location)
            if loc_dim is None:
                continue
            name = _att(var, NAME_ATT) or vname
            dims = var.dimensions
            var.set_auto_mask(False)
            k = 1 if (dims and time_dim is not None and dims[0] == time_dim) else 0
            if k >= len(dims) or dims[k] != loc_dim:
                warn(
                    f"ugrid_nc: '{vname}' is not laid out as ([time,] {location}"
                    "[, components]); skipped"
                )
                continue
            data = np.asarray(var[step] if k else var[:])
            width = int(np.prod(data.shape[1:])) if data.ndim > 1 else 1
            data = (
                data.reshape(data.shape[0]) if width == 1 else data.reshape(-1, width)
            )
            {"node": point_data, "face": face_data, "edge": edge_data}[location][
                name
            ] = data
        cell_data = {}
        names = list(face_data) + [n for n in edge_data if n not in face_data]
        for name in sorted(names):
            like = face_data[name] if name in face_data else edge_data[name]
            blocks = [
                face_data[name][idx] if name in face_data else _pad(like, len(idx))
                for idx in face_blocks
            ]
            if num_edges:
                blocks.append(
                    edge_data[name] if name in edge_data else _pad(like, num_edges)
                )
            if blocks:
                cell_data[name] = blocks
    return Mesh(
        points, cells, point_data=point_data, cell_data=cell_data, field_data=field_data
    )


def _sanitize(name):
    out = re.sub(r"[^A-Za-z0-9_.\-]", "_", name)
    if not out or not (out[0].isalpha() or out[0] == "_"):
        out = "_" + out
    return out


def write(filename, mesh):
    nc4 = _netcdf4()
    face_blocks, line_blocks = [], []
    for b, block in enumerate(mesh.cells):
        if block.type in ("triangle", "quad", "polygon"):
            face_blocks.append(b)
        elif block.type == "line":
            line_blocks.append(b)
        else:
            raise WriteError(
                f"meshio++: ugrid_nc: a '{block.type}' block has no UGRID 1-D or 2-D "
                "topology; write triangles, quads, polygons and lines"
            )
    face_rows = []
    for b in face_blocks:
        data = mesh.cells[b].data
        face_rows += [np.asarray(r, dtype=np.int64) for r in data]
    edge_rows = [np.asarray(mesh.cells[b].data, dtype=np.int64) for b in line_blocks]
    points = np.asarray(mesh.points, dtype=np.float64)
    with_z = points.shape[1] > 2 and bool(np.any(points[:, 2] != 0.0))
    with nc4.Dataset(str(filename), "w", format="NETCDF4") as nc:
        nc.Conventions = "CF-1.8 UGRID-1.0"
        nc.createDimension("nmesh_node", len(points))
        topo = nc.createVariable("mesh", "i4")
        topo.cf_role = "mesh_topology"
        topo.long_name = "Topology data of the mesh"
        topo.topology_dimension = np.int32(2 if face_blocks else 1)
        topo.node_coordinates = (
            "mesh_node_x mesh_node_y mesh_node_z"
            if with_z
            else "mesh_node_x mesh_node_y"
        )
        for c, axis in enumerate("xyz"):
            if axis == "z" and not with_z:
                break
            v = nc.createVariable(f"mesh_node_{axis}", "f8", ("nmesh_node",))
            v.standard_name = f"projection_{axis}_coordinate"
            v[:] = points[:, c] if c < points.shape[1] else 0.0
        if face_blocks:
            max_nodes = max(len(r) for r in face_rows)
            nc.createDimension("nmesh_face", len(face_rows))
            nc.createDimension("nmesh_face_max_nodes", max_nodes)
            fv = nc.createVariable(
                "mesh_face_nodes",
                "i4",
                ("nmesh_face", "nmesh_face_max_nodes"),
                fill_value=np.int32(-1),
            )
            fv.cf_role = "face_node_connectivity"
            fv.start_index = np.int32(0)
            conn = np.full((len(face_rows), max_nodes), -1, dtype=np.int32)
            for r, row in enumerate(face_rows):
                conn[r, : len(row)] = row
            fv[:] = conn
            topo.face_node_connectivity = "mesh_face_nodes"
            topo.face_dimension = "nmesh_face"
        if line_blocks:
            edges = np.concatenate(edge_rows)
            nc.createDimension("nmesh_edge", len(edges))
            nc.createDimension("Two", 2)
            ev = nc.createVariable("mesh_edge_nodes", "i4", ("nmesh_edge", "Two"))
            ev.cf_role = "edge_node_connectivity"
            ev.start_index = np.int32(0)
            ev[:] = edges.astype(np.int32)
            topo.edge_node_connectivity = "mesh_edge_nodes"
            topo.edge_dimension = "nmesh_edge"

        def define(name, var_name, location, dim, array):
            array = np.asarray(array)
            dims = (dim,)
            if array.ndim > 1:
                width = int(np.prod(array.shape[1:]))
                array = array.reshape(len(array), width)
                if width != 1:
                    comp = f"n{var_name}_comp"
                    nc.createDimension(comp, width)
                    dims = (dim, comp)
                else:
                    array = array.reshape(-1)
            v = nc.createVariable(var_name, array.dtype, dims)
            v.mesh = "mesh"
            v.location = location
            if var_name != name:
                setattr(v, NAME_ATT, name)
            v[:] = array

        for name in sorted(mesh.point_data):
            define(name, _sanitize(name), "node", "nmesh_node", mesh.point_data[name])
        for name in sorted(mesh.cell_data):
            arrays = mesh.cell_data[name]
            if len(arrays) != len(mesh.cells):
                continue
            fparts = [np.asarray(arrays[b]) for b in face_blocks]
            eparts = [np.asarray(arrays[b]) for b in line_blocks]
            if fparts:
                define(
                    name, _sanitize(name), "face", "nmesh_face", np.concatenate(fparts)
                )
            if eparts:
                define(
                    name,
                    _sanitize(name) + ("_edge" if fparts else ""),
                    "edge",
                    "nmesh_edge",
                    np.concatenate(eparts),
                )
