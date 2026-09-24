"""
The Tecplot zone model both readers decode into, and the one mesh builder
they share (the Python twin of ``tecplot.cpp``'s zone/timeline/step code).

A source (ASCII tokens or binary offsets) fills each zone's header fields and
provides ``source.own_data(zone_index)`` -> ``(columns, connectivity)``, the
variables the zone owns (neither shared nor passive) and its own FE
connectivity (0-based, ``None`` when shared or ordered). Sharing, passive
variables, ordered connectivity, time steps and regions are resolved here,
once for both.
"""

import numpy as np

from .._exceptions import ReadError
from .._mesh import Mesh
from .._regions import Region, block_bases

FE_TYPES = {
    "FELINESEG": ("line", 2),
    "FETRIANGLE": ("triangle", 3),
    "FEQUADRILATERAL": ("quad", 4),
    "FETETRAHEDRON": ("tetra", 4),
    "FEBRICK": ("hexahedron", 8),
}


class Zone:
    def __init__(self):
        self.title = ""
        self.type_name = "ORDERED"
        self.ordered = False
        self.ijk = (1, 1, 1)
        self.num_nodes = 0
        self.num_cells = 0
        self.cell_centered = []
        self.var_share = {}  # variable -> 0-based source zone
        self.passive = set()
        self.conn_share = -1  # 0-based source zone, -1 = own
        self.has_solution_time = False
        self.solution_time = 0.0
        self.has_strand = False
        self.strand = 0

    def finish(self):
        """Derives node/cell counts of an ordered zone from I/J/K."""
        if self.ordered:
            dims = [d for d in self.ijk if d > 1]
            self.num_nodes = int(np.prod(self.ijk))
            # a single-point zone is one vertex
            self.num_cells = int(np.prod([d - 1 for d in dims])) if dims else 1

    def meshio_type(self):
        if self.ordered:
            n = sum(d > 1 for d in self.ijk)
            return {0: "vertex", 1: "line", 2: "quad", 3: "hexahedron"}[n]
        if self.type_name in ("FEPOLYGON", "FEPOLYHEDRON"):
            raise ReadError(
                f"Tecplot: {self.type_name} zones (polygonal/polyhedral) are not supported"
            )
        if self.type_name not in FE_TYPES:
            raise ReadError(f"Tecplot: unsupported zone type {self.type_name}")
        return FE_TYPES[self.type_name][0]

    def owns(self, var):
        return var not in self.var_share and var not in self.passive

    def data_length(self, var):
        return self.num_cells if self.cell_centered[var] else self.num_nodes


def ordered_connectivity(ijk):
    """Cells of an ordered zone over its node index ``i + I*(j + J*k)``: lines,
    quads or hexahedra over the dimensions longer than one, in VTK order."""
    ni, nj, nk = ijk
    if ni * nj * nk == 1:
        return np.zeros((1, 1), dtype=np.int64)
    idx = np.arange(ni * nj * nk, dtype=np.int64).reshape(nk, nj, ni)
    # Collapse unit dimensions so a J- or JK-ordered zone is meshed like an
    # I- or IJ-ordered one.
    grid = idx.transpose(2, 1, 0)  # (i, j, k)
    grid = grid.reshape([d for d in (ni, nj, nk) if d > 1] or [1])
    if grid.ndim == 1:
        return np.column_stack([grid[:-1], grid[1:]])
    if grid.ndim == 2:
        a = grid[:-1, :-1]
        b = grid[1:, :-1]
        c = grid[1:, 1:]
        d = grid[:-1, 1:]
        # i fastest, then j -- the cell-centred value order
        return np.stack([a, b, c, d], axis=-1).transpose(1, 0, 2).reshape(-1, 4)
    corners = [
        grid[:-1, :-1, :-1],
        grid[1:, :-1, :-1],
        grid[1:, 1:, :-1],
        grid[:-1, 1:, :-1],
        grid[:-1, :-1, 1:],
        grid[1:, :-1, 1:],
        grid[1:, 1:, 1:],
        grid[:-1, 1:, 1:],
    ]
    return np.stack(corners, axis=-1).transpose(2, 1, 0, 3).reshape(-1, 8)


def timeline(zones):
    """One entry per step: the zone indices of that step (``tecplot_timeline``)."""
    if not zones[0].has_solution_time:
        return [list(range(len(zones)))]
    by_time = {}
    for k, z in enumerate(zones):
        if not z.has_solution_time:
            continue
        if zones[0].has_strand and z.has_strand and z.strand != zones[0].strand:
            continue
        by_time.setdefault(z.solution_time, []).append(k)
    return [by_time[t] for t in sorted(by_time)]


def xyz_indices(variables):
    xi = yi = zi = -1
    for k, v in enumerate(variables):
        u = v.upper()
        if u == "X":
            xi = k
        elif u == "Y":
            yi = k
        elif u == "Z":
            zi = k
    return xi, yi, zi


def _points_origin(idx, zones, xi, yi, zi):
    """The zone whose points zone ``idx`` literally reuses: itself, unless every
    one of its nodal variables (coordinates and fields alike) is shared from one
    earlier zone, followed down the chain. A zone that shares only its
    coordinates but carries its own nodal values keeps its own copy of the
    points, so no value is lost."""
    seen = set()
    cur = idx
    while cur not in seen:
        seen.add(cur)
        z = zones[cur]
        src = z.var_share.get(xi)
        if src is None:
            return cur
        nodal = [v for v in range(len(z.cell_centered)) if not z.cell_centered[v]]
        if any(z.var_share.get(v) != src for v in nodal):
            return cur
        cur = src
    return idx


class _Decoder:
    def __init__(self, source, zones, variables):
        self.source = source
        self.zones = zones
        self.variables = variables
        self.cache = {}

    def zone(self, idx, depth=0):
        if idx in self.cache:
            return self.cache[idx]
        if depth > len(self.zones):
            raise ReadError("Tecplot: circular VARSHARELIST/CONNECTIVITYSHAREZONE")
        z = self.zones[idx]
        mtype = z.meshio_type()
        own, conn = self.source.own_data(idx)
        cols = []
        for v in range(len(self.variables)):
            if v in z.var_share:
                src = z.var_share[v]
                if not 0 <= src < len(self.zones) or src == idx:
                    raise ReadError(
                        f"Tecplot: zone {idx + 1} shares from bad zone {src + 1}"
                    )
                cols.append(self.zone(src, depth + 1)["cols"][v])
            elif v in z.passive:
                cols.append(np.full(z.data_length(v), np.nan))
            else:
                cols.append(own[v])
        if z.ordered:
            conn = ordered_connectivity(z.ijk)
        elif z.conn_share >= 0:
            if not 0 <= z.conn_share < len(self.zones) or z.conn_share == idx:
                raise ReadError(
                    f"Tecplot: zone {idx + 1} shares connectivity from bad zone {z.conn_share + 1}"
                )
            conn = self.zone(z.conn_share, depth + 1)["conn"]
        out = {
            "cols": cols,
            "cc": z.cell_centered,
            "type": mtype,
            "conn": conn,
            "nodes": z.num_nodes,
            "cells": z.num_cells,
        }
        self.cache[idx] = out
        return out


def build_step(step_zones, zones, variables, source):
    """``tecplot_build_step_mesh``: one block and one Cell region per zone."""
    xi, yi, zi = xyz_indices(variables)
    if xi < 0:
        raise ReadError("Variable 'X' not found")
    if yi < 0:
        raise ReadError("Variable 'Y' not found")
    dec = _Decoder(source, zones, variables)
    decoded = [dec.zone(i) for i in step_zones]

    pos = {z: k for k, z in enumerate(step_zones)}
    offset = [0] * len(decoded)
    owns = [True] * len(decoded)
    total = 0
    for k, idx in enumerate(step_zones):
        origin = _points_origin(idx, zones, xi, yi, zi)
        if origin != idx and pos.get(origin, len(decoded)) < k:
            offset[k] = offset[pos[origin]]
            owns[k] = False
        else:
            offset[k] = total
            total += decoded[k]["nodes"]

    ndim = 3 if zi >= 0 else 2
    points = np.empty((total, ndim))
    for k, d in enumerate(decoded):
        if not owns[k]:
            continue
        sl = slice(offset[k], offset[k] + d["nodes"])
        points[sl, 0] = d["cols"][xi]
        points[sl, 1] = d["cols"][yi]
        if zi >= 0:
            points[sl, 2] = d["cols"][zi]

    cells = [
        (d["type"], np.asarray(d["conn"], dtype=np.int64) + offset[k])
        for k, d in enumerate(decoded)
    ]
    # A variable is cell data where a zone stores it cell-centred and point
    # data where it stores it at the nodes; a variable that is nodal in some
    # zones and cell-centred in others becomes both, NaN where absent.
    point_data, cell_data = {}, {}
    for v, name in enumerate(variables):
        if v in (xi, yi, zi):
            continue
        cc = [bool(d["cc"][v]) for d in decoded]
        if any(cc):
            cell_data[name] = [
                (
                    np.asarray(d["cols"][v], dtype=np.float64)
                    if c
                    else np.full(d["cells"], np.nan)
                )
                for d, c in zip(decoded, cc)
            ]
        if not all(cc):
            col = np.full(total, np.nan)
            for k, d in enumerate(decoded):
                if owns[k] and not cc[k]:
                    col[offset[k] : offset[k] + d["nodes"]] = d["cols"][v]
            point_data[name] = col
    cell_data["tecplot:zone"] = [
        np.full(d["cells"], idx, dtype=np.int64) for idx, d in zip(step_zones, decoded)
    ]
    mesh = Mesh(points, cells, point_data=point_data, cell_data=cell_data)

    bases = block_bases(mesh.cells)
    used = []
    regions = []
    for k, idx in enumerate(step_zones):
        name = zones[idx].title or f"zone_{k}"
        unique = name
        suffix = 2
        while unique in used:
            unique = f"{name}_{suffix}"
            suffix += 1
        used.append(unique)
        entries = np.arange(bases[k], bases[k] + decoded[k]["cells"], dtype=np.int64)
        regions.append(Region(unique, "cell", entries, dim=-1, tag=idx))
    mesh.regions = regions
    return mesh


def metadata(zones, variables):
    """``read_tecplot_metadata``: cell blocks and point count of the first
    step, and the time values of a transient file."""
    steps = timeline(zones)
    xi, yi, zi = xyz_indices(variables)
    first = steps[0]
    pos = {z: k for k, z in enumerate(first)}
    blocks = []
    total = 0
    for k, idx in enumerate(first):
        origin = _points_origin(idx, zones, xi, yi, zi)
        if not (origin != idx and pos.get(origin, len(first)) < k):
            total += zones[idx].num_nodes
        blocks.append((zones[idx].meshio_type(), zones[idx].num_cells))
    times = (
        [zones[s[0]].solution_time for s in steps] if zones[0].has_solution_time else []
    )
    return blocks, total, times
