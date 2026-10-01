"""Conforming adjacency between named cell regions.

The numpy implementation is the reference for ``operations/interfaces.hpp``.
It intentionally matches facets by exact corner-node ids: geometrically close
but separately numbered faces are a proximity-interface problem, not
conforming adjacency.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from ._fallback import core_op_declined
from ._mesh import Mesh, topological_dimension
from ._regions import Region, block_bases
from ._skin import _CELL_FACES
from ._surface import _CELL_EDGES

_PREFIX = "meshio++: region_adjacency: "


def _selector(spec, default_kind="cell"):
    if isinstance(spec, str):
        return {"name": spec, "kind": default_kind}
    sel = {"name": str(spec["name"]), "kind": spec.get("kind") or default_kind}
    for key in ("dim", "tag"):
        if spec.get(key) is not None and int(spec[key]) != -2:
            sel[key] = int(spec[key])
    return sel


def _selected_regions(mesh, regions):
    all_regions = list(getattr(mesh, "regions", None) or [])
    specs = [] if regions is None else list(regions)
    if not specs:
        selected = sorted(
            (r for r in all_regions if r.kind == "cell"), key=lambda region: region.key
        )
        block_groups = len(selected) < 2 and len(mesh.cells) >= 2
        if block_groups:
            selected = [
                Region(f"block:{i}", "cell", []) for i in range(len(mesh.cells))
            ]
    else:
        selected = []
        block_groups = False
        for raw in specs:
            spec = _selector(raw, "cell")
            matches = [
                r
                for r in all_regions
                if r.name == spec["name"]
                and r.kind == spec["kind"]
                and ("dim" not in spec or r.dim == spec["dim"])
                and ("tag" not in spec or r.tag == spec["tag"])
            ]
            if len(matches) != 1:
                state = "does not exist" if not matches else "is ambiguous"
                raise ValueError(
                    f"{_PREFIX}Cell region {spec['name']!r} {state}; "
                    "specify kind, dim and/or tag"
                )
            if matches[0] in selected:
                raise ValueError(f"{_PREFIX}the same region was selected twice")
            selected.append(matches[0])
    if len(selected) < 2:
        raise ValueError(
            f"{_PREFIX}select at least two Cell regions or use at least two cell blocks"
        )
    return selected, block_groups


def _cell_facets(block, row):
    """Yield ``(corner-node-list, local-facet-index)`` for a cell."""
    data = block.data
    if block.type.startswith("polyhedron"):
        for facet, face in enumerate(data[row]):
            yield [int(node) for node in face], facet
        return
    if block.type.startswith("polygon"):
        nodes = [int(node) for node in data[row]]
        for facet, node in enumerate(nodes):
            yield [node, nodes[(facet + 1) % len(nodes)]], facet
        return
    if block.type in _CELL_FACES:
        conn = np.asarray(data[row], dtype=np.int64)
        for facet, (_name, ncorner, local) in enumerate(_CELL_FACES[block.type]):
            yield [int(conn[i]) for i in local[:ncorner]], facet
        return
    if block.type in _CELL_EDGES:
        conn = np.asarray(data[row], dtype=np.int64)
        for facet, (_name, _ncorner, local) in enumerate(_CELL_EDGES[block.type]):
            yield [int(conn[local[0]]), int(conn[local[1]])], facet


def _measure(points, nodes):
    xyz = np.zeros((len(nodes), 3), dtype=np.float64)
    pts = np.asarray(points, dtype=np.float64)
    xyz[:, : min(3, pts.shape[1])] = pts[np.asarray(nodes, dtype=np.int64), :3]
    if len(nodes) == 2:
        return float(np.linalg.norm(xyz[1] - xyz[0]))
    area = 0.0
    for i in range(1, len(nodes) - 1):
        area += 0.5 * float(
            np.linalg.norm(np.cross(xyz[i] - xyz[0], xyz[i + 1] - xyz[0]))
        )
    return area


def _region_adjacency_py(mesh, regions=None):
    selected, block_groups = _selected_regions(mesh, regions)
    bases = block_bases(mesh.cells)
    total = int(bases[-1])
    membership = [[] for _ in range(total)]
    if block_groups:
        for block_id, base in enumerate(bases[:-1]):
            for cell in range(int(base), int(bases[block_id + 1])):
                membership[cell].append(block_id)
    else:
        for rid, region in enumerate(selected):
            entries = np.asarray(region.entries, dtype=np.int64).reshape(-1)
            for cell in entries[(entries >= 0) & (entries < total)]:
                membership[int(cell)].append(rid)

    # Dict insertion order never reaches output: keys, region pairs and owner
    # pairs are explicitly sorted below for backend/thread-independent results.
    facets = {}
    for base, block in zip(bases[:-1], mesh.cells):
        dim = (
            3
            if block.type.startswith("polyhedron")
            else topological_dimension.get(block.type, 0)
        )
        if dim not in (2, 3):
            continue
        for row in range(len(block.data)):
            cell = int(base + row)
            if not membership[cell]:
                continue
            for nodes, facet in _cell_facets(block, row):
                if len(nodes) < 2 or any(n < 0 or n >= len(mesh.points) for n in nodes):
                    continue
                key = tuple(sorted(nodes))
                if len(set(key)) != len(key):
                    continue
                facets.setdefault(key, []).append((nodes, cell, facet))

    pair_facets = {}
    for key in sorted(facets):
        owners = sorted(facets[key], key=lambda item: (item[1], item[2]))
        for i, a in enumerate(owners):
            for b in owners[i + 1 :]:
                if a[1] == b[1]:
                    continue
                for member_a in membership[a[1]]:
                    for member_b in membership[b[1]]:
                        if member_a == member_b:
                            continue
                        if member_a > member_b:
                            region_a, region_b = member_b, member_a
                            owner_a, owner_b = b, a
                        else:
                            region_a, region_b = member_a, member_b
                            owner_a, owner_b = a, b
                        pair_key = (region_a, region_b, key)
                        record = pair_facets.setdefault(
                            pair_key,
                            {
                                "a": owner_a,
                                "b": owner_b,
                                "cells_a": set(),
                                "cells_b": set(),
                            },
                        )
                        if (owner_a[1], owner_a[2]) < (record["a"][1], record["a"][2]):
                            record["a"] = owner_a
                        if (owner_b[1], owner_b[2]) < (record["b"][1], record["b"][2]):
                            record["b"] = owner_b
                        record["cells_a"].add(owner_a[1])
                        record["cells_b"].add(owner_b[1])

    by_type = {"line": [], "triangle": [], "quad": [], "polygon": []}
    for (region_a, region_b, _key), record in sorted(pair_facets.items()):
        owner_a, owner_b = record["a"], record["b"]
        nodes = owner_a[0]
        cell_type = {2: "line", 3: "triangle", 4: "quad"}.get(len(nodes), "polygon")
        by_type[cell_type].append(
            {
                "nodes": nodes,
                "region_a": region_a,
                "region_b": region_b,
                "cell_a": owner_a[1],
                "facet_a": owner_a[2],
                "cell_b": owner_b[1],
                "facet_b": owner_b[2],
                "shared_count": len(record["cells_a"]) + len(record["cells_b"]),
                "measure": _measure(mesh.points, nodes),
            }
        )

    blocks = []
    cell_data = {
        key: []
        for key in (
            "interface:region_a",
            "interface:region_b",
            "interface:parent_cell_a",
            "interface:parent_facet_a",
            "interface:parent_cell_b",
            "interface:parent_facet_b",
            "interface:shared_count",
            "interface:measure",
        )
    }
    pair_entries = {}
    output_cell = 0
    for cell_type, records in by_type.items():
        if not records:
            continue
        cells = [record["nodes"] for record in records]
        if cell_type == "polygon":
            conn = cells
        else:
            conn = np.asarray(cells, dtype=np.int64).reshape(-1, len(cells[0]))
        blocks.append((cell_type, conn))
        for field, key in (
            ("region_a", "interface:region_a"),
            ("region_b", "interface:region_b"),
            ("cell_a", "interface:parent_cell_a"),
            ("facet_a", "interface:parent_facet_a"),
            ("cell_b", "interface:parent_cell_b"),
            ("facet_b", "interface:parent_facet_b"),
            ("shared_count", "interface:shared_count"),
            ("measure", "interface:measure"),
        ):
            dtype = np.float64 if field == "measure" else np.int64
            cell_data[key].append(
                np.asarray([record[field] for record in records], dtype=dtype)
            )
        for record in records:
            pair_name = f"adjacency:{record['region_a']}:{record['region_b']}"
            pair_entries.setdefault(pair_name, []).append(output_cell)
            output_cell += 1

    kept_regions = [
        r.copy() for r in getattr(mesh, "regions", ()) or () if r.kind == "point"
    ]
    kept_regions.extend(
        Region(name, "cell", entries) for name, entries in pair_entries.items()
    )
    kept_regions.sort(key=lambda region: region.key)
    return Mesh(
        np.array(mesh.points, copy=True),
        blocks,
        point_data={
            key: np.array(value, copy=True) for key, value in mesh.point_data.items()
        },
        cell_data=cell_data,
        field_data={
            key: np.array(value, copy=True) for key, value in mesh.field_data.items()
        },
        regions=kept_regions,
    )


def _resolve_region(mesh, raw, kind, operation):
    all_regions = list(getattr(mesh, "regions", None) or [])
    if isinstance(raw, Region):
        if raw.kind != kind:
            raise ValueError(f"meshio++: {operation}: expected a {kind} region")
        return raw
    spec = _selector(raw, kind)
    if spec["kind"] != kind:
        raise ValueError(
            f"meshio++: {operation}: selector {spec['name']!r} must have kind {kind!r}"
        )
    if spec["name"].startswith("block:") and kind == "cell":
        try:
            block = int(spec["name"].split(":", 1)[1])
        except ValueError as exc:
            raise ValueError(
                f"meshio++: {operation}: invalid block selector {spec['name']!r}"
            ) from exc
        if block < 0 or block >= len(mesh.cells):
            raise ValueError(f"meshio++: {operation}: block selector is out of range")
        return Region(
            spec["name"],
            "cell",
            np.arange(
                block_bases(mesh.cells)[block], block_bases(mesh.cells)[block + 1]
            ),
        )
    matches = [
        region
        for region in all_regions
        if region.name == spec["name"]
        and region.kind == kind
        and ("dim" not in spec or spec["dim"] == -2 or region.dim == spec["dim"])
        and ("tag" not in spec or spec["tag"] == -2 or region.tag == spec["tag"])
    ]
    if len(matches) != 1:
        state = "does not exist" if not matches else "is ambiguous"
        raise ValueError(
            f"meshio++: {operation}: {kind.capitalize()} region {spec['name']!r} {state}; "
            "pin dim and/or tag"
        )
    return matches[0]


def _region_cell_mask(mesh, raw, operation):
    region = _resolve_region(mesh, raw, "cell", operation)
    total = int(block_bases(mesh.cells)[-1])
    mask = np.zeros(total, dtype=bool)
    ids = np.asarray(region.entries, dtype=np.int64).reshape(-1)
    if np.any((ids < 0) | (ids >= total)):
        raise ValueError(
            f"meshio++: {operation}: selected Cell region contains an invalid cell id"
        )
    mask[ids] = True
    return region, mask


def _boundary_facets(mesh, mask):
    bases = block_bases(mesh.cells)
    indexed = defaultdict(list)
    for base, block in zip(bases[:-1], mesh.cells):
        for row in range(len(block.data)):
            cell = int(base + row)
            if not mask[cell]:
                continue
            for nodes, facet in _cell_facets(block, row):
                key = tuple(sorted(nodes))
                if len(key) > 1 and len(set(key)) == len(key):
                    indexed[key].append((list(nodes), cell, facet))
    out = []
    for key in sorted(indexed):
        owners = sorted(indexed[key], key=lambda item: (item[1], item[2]))
        owners = list({owner[1]: owner for owner in owners}.values())
        if len(owners) == 1:
            out.append(owners[0])
    return sorted(out, key=lambda item: (item[1], item[2], item[0]))


def _xyz(mesh, point):
    out = np.zeros(3, dtype=np.float64)
    p = np.asarray(mesh.points, dtype=np.float64)[int(point)]
    out[: min(3, p.size)] = p[:3]
    return out


def _normal(mesh, nodes):
    xyz = np.asarray([_xyz(mesh, node) for node in nodes])
    if len(xyz) == 2:
        vector = xyz[1] - xyz[0]
        normal = np.array([vector[1], -vector[0], 0.0])
    else:
        normal = np.zeros(3)
        for i in range(1, len(xyz) - 1):
            normal += np.cross(xyz[i] - xyz[0], xyz[i + 1] - xyz[0])
    length = np.linalg.norm(normal)
    return normal / length if length > 0 else np.zeros(3)


def _closest_triangle(point, a, b, c):
    # Ericson's closest-point regions, matching the core point/triangle kernel.
    ab, ac, ap = b - a, c - a, point - a
    d1, d2 = np.dot(ab, ap), np.dot(ac, ap)
    if d1 <= 0 and d2 <= 0:
        q, weights = a, np.array([1.0, 0.0, 0.0])
    else:
        bp = point - b
        d3, d4 = np.dot(ab, bp), np.dot(ac, bp)
        if d3 >= 0 and d4 <= d3:
            q, weights = b, np.array([0.0, 1.0, 0.0])
        else:
            vc = d1 * d4 - d3 * d2
            if vc <= 0 and d1 >= 0 and d3 <= 0:
                v = d1 / (d1 - d3)
                q, weights = a + v * ab, np.array([1 - v, v, 0.0])
            else:
                cp = point - c
                d5, d6 = np.dot(ab, cp), np.dot(ac, cp)
                if d6 >= 0 and d5 <= d6:
                    q, weights = c, np.array([0.0, 0.0, 1.0])
                else:
                    vb = d5 * d2 - d1 * d6
                    if vb <= 0 and d2 >= 0 and d6 <= 0:
                        w = d2 / (d2 - d6)
                        q, weights = a + w * ac, np.array([1 - w, 0.0, w])
                    else:
                        va = d3 * d6 - d5 * d4
                        if va <= 0 and d4 - d3 >= 0 and d5 - d6 >= 0:
                            w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
                            q, weights = b + w * (c - b), np.array([0.0, 1 - w, w])
                        else:
                            den = va + vb + vc
                            if den <= 0:
                                candidates = []
                                for edge_id, (u, v) in enumerate(
                                    ((a, b), (b, c), (c, a))
                                ):
                                    edge = v - u
                                    t = (
                                        np.clip(
                                            np.dot(point - u, edge)
                                            / np.dot(edge, edge),
                                            0,
                                            1,
                                        )
                                        if np.dot(edge, edge)
                                        else 0
                                    )
                                    q_edge = u + t * edge
                                    edge_weights = (
                                        np.array([1.0 - t, t, 0.0])
                                        if edge_id == 0
                                        else (
                                            np.array([0.0, 1.0 - t, t])
                                            if edge_id == 1
                                            else np.array([t, 0.0, 1.0 - t])
                                        )
                                    )
                                    candidates.append(
                                        (
                                            np.dot(point - q_edge, point - q_edge),
                                            q_edge,
                                            edge_weights,
                                        )
                                    )
                                _d2, q, weights = min(
                                    candidates, key=lambda item: item[0]
                                )
                            else:
                                v, w = vb / den, vc / den
                                q, weights = a + v * ab + w * ac, np.array(
                                    [1 - v - w, v, w]
                                )
    return q, weights, float(np.dot(point - q, point - q))


def _project(mesh, point, facet):
    nodes, cell, local_facet = facet
    xyz = [_xyz(mesh, node) for node in nodes]
    best = None
    subfacet = 0
    if len(nodes) == 2:
        edge = xyz[1] - xyz[0]
        t = (
            np.clip(np.dot(point - xyz[0], edge) / np.dot(edge, edge), 0, 1)
            if np.dot(edge, edge)
            else 0
        )
        q = xyz[0] + t * edge
        best = (float(np.dot(point - q, point - q)), q, np.array([1 - t, t, 0.0]))
    else:
        for i in range(1, len(nodes) - 1):
            q, weights, distance_sq = _closest_triangle(
                point, xyz[0], xyz[i], xyz[i + 1]
            )
            if best is None or distance_sq < best[0]:
                best = (distance_sq, q, weights)
                subfacet = i - 1
    if best is None:
        return None
    distance_sq, q, weights = best
    return {
        "distance_sq": distance_sq,
        "point": q,
        "weights": weights,
        "normal": _normal(mesh, nodes),
        "cell": cell,
        "facet": local_facet,
        "subfacet": subfacet,
        "nodes": nodes,
    }


def _facet_center(mesh, facet):
    return np.mean([_xyz(mesh, node) for node in facet[0]], axis=0)


def _mean_edge(mesh, facets):
    lengths = []
    for nodes, _cell, _facet in facets:
        xyz = [_xyz(mesh, node) for node in nodes]
        lengths.extend(
            np.linalg.norm(xyz[(i + 1) % len(xyz)] - xyz[i]) for i in range(len(xyz))
        )
    return float(np.mean(lengths)) if lengths else 0.0


def _facet_measure(mesh, nodes):
    return _measure(mesh.points, nodes)


def _make_interface_mesh(master, output_facets):
    from ._mesh import Mesh

    order = ("line", "triangle", "quad", "polygon")
    cells, cell_data = [], {
        name: []
        for name in (
            "interface:parent_cell",
            "interface:parent_facet",
            "interface:partner_cell",
            "interface:partner_facet",
            "interface:gap",
            "interface:measure",
        )
    }
    type_names = {2: "line", 3: "triangle", 4: "quad"}
    for cell_type in order:
        records = [
            r
            for r in output_facets
            if type_names.get(len(r["nodes"]), "polygon") == cell_type
        ]
        if not records:
            continue
        conn = [r["nodes"] for r in records]
        cells.append(
            (
                cell_type,
                conn if cell_type == "polygon" else np.asarray(conn, dtype=np.int64),
            )
        )
        for name, key, dtype in (
            ("interface:parent_cell", "cell", np.int64),
            ("interface:parent_facet", "facet", np.int64),
            ("interface:partner_cell", "partner_cell", np.int64),
            ("interface:partner_facet", "partner_facet", np.int64),
            ("interface:gap", "gap", np.float64),
            ("interface:measure", "measure", np.float64),
        ):
            cell_data[name].append(np.asarray([r[key] for r in records], dtype=dtype))
    kept_regions = [
        r.copy() for r in getattr(master, "regions", ()) or () if r.kind == "point"
    ]
    return Mesh(
        np.array(master.points, copy=True),
        cells,
        point_data={k: np.array(v, copy=True) for k, v in master.point_data.items()},
        cell_data=cell_data if cells else {},
        field_data={k: np.array(v, copy=True) for k, v in master.field_data.items()},
        regions=kept_regions,
    )


def _find_interface_py(
    mesh_a,
    region_a,
    region_b,
    mesh_b=None,
    mode="conforming",
    master="a",
    gap_tolerance=0.0,
    angle_tolerance=30.0,
    overlap_tolerance=0.0,
):
    mesh_b = mesh_a if mesh_b is None else mesh_b
    _reg_a, mask_a = _region_cell_mask(mesh_a, region_a, "find_interface")
    _reg_b, mask_b = _region_cell_mask(mesh_b, region_b, "find_interface")
    facets_a, facets_b = _boundary_facets(mesh_a, mask_a), _boundary_facets(
        mesh_b, mask_b
    )
    if not facets_a or not facets_b:
        raise ValueError(
            "meshio++: find_interface: both selected parts must have boundary facets"
        )
    if mode not in ("conforming", "proximity") or master not in ("a", "b"):
        raise ValueError("meshio++: find_interface: invalid mode or master side")
    if (
        not np.isfinite(gap_tolerance)
        or gap_tolerance < 0
        or not np.isfinite(overlap_tolerance)
        or overlap_tolerance < 0
    ):
        raise ValueError(
            "meshio++: find_interface: tolerances must be finite and non-negative"
        )
    if not np.isfinite(angle_tolerance) or not 0 <= angle_tolerance <= 180:
        raise ValueError(
            "meshio++: find_interface: angle_tolerance must be in [0, 180]"
        )
    matches, matched_a, matched_b = [], set(), set()
    if mode == "conforming":
        by_a, by_b = defaultdict(list), defaultdict(list)
        for facet in facets_a:
            by_a[tuple(sorted(facet[0]))].append(facet)
        for facet in facets_b:
            by_b[tuple(sorted(facet[0]))].append(facet)
        for key in sorted(by_a.keys() & by_b.keys()):
            for a in by_a[key]:
                for b in by_b[key]:
                    if mesh_a is mesh_b and a[1] == b[1]:
                        continue
                    if any(
                        not np.array_equal(_xyz(mesh_a, node), _xyz(mesh_b, node))
                        for node in a[0]
                    ):
                        continue
                    matches.append((a, b, 0.0, 0.0))
                    matched_a.add((a[1], a[2]))
                    matched_b.add((b[1], b[2]))
    else:
        tolerance = (
            gap_tolerance
            or 0.01
            * (_mean_edge(mesh_a, facets_a) + _mean_edge(mesh_b, facets_b))
            * 0.5
        )
        tolerance += overlap_tolerance
        limit_sq = tolerance * tolerance
        normal_limit = -np.cos(np.deg2rad(angle_tolerance))
        for a in facets_a:
            center = _facet_center(mesh_a, a)
            projections = [_project(mesh_b, center, b) for b in facets_b]
            projections = [
                p for p in projections if p is not None and p["distance_sq"] <= limit_sq
            ]
            normal_a = _normal(mesh_a, a[0])
            if not np.linalg.norm(normal_a):
                continue
            projections = [
                p
                for p in projections
                if np.linalg.norm(p["normal"])
                and np.dot(normal_a, p["normal"]) <= normal_limit
            ]
            if not projections:
                continue
            best = min(
                projections, key=lambda p: (p["distance_sq"], p["cell"], p["facet"])
            )
            if any(
                min(
                    (
                        hit["distance_sq"]
                        for b in facets_b
                        if (hit := _project(mesh_b, _xyz(mesh_a, node), b)) is not None
                    ),
                    default=np.inf,
                )
                > limit_sq
                for node in a[0]
            ):
                continue
            delta = center - best["point"]
            gap_a = float(np.dot(delta, normal_a))
            gap_b = float(np.dot(delta, best["normal"]))
            b = (best["nodes"], best["cell"], best["facet"])
            matches.append((a, b, gap_a, gap_b))
            matched_a.add((a[1], a[2]))
            matched_b.add((b[1], b[2]))
    matches.sort(key=lambda item: (item[0][1], item[0][2], item[1][1], item[1][2]))
    use_a = master == "a"
    master_mesh = mesh_a if use_a else mesh_b
    output_facets, side_a, side_b = [], [], []
    area = max_gap = 0.0
    for a, b, gap_a, gap_b in matches:
        side_a.append([a[1], a[2]])
        side_b.append([b[1], b[2]])
        mf, partner = (a, b) if use_a else (b, a)
        gap = 0.0 if mode == "conforming" else gap_a if use_a else gap_b
        measure = _facet_measure(master_mesh, mf[0])
        output_facets.append(
            {
                "nodes": mf[0],
                "cell": mf[1],
                "facet": mf[2],
                "partner_cell": partner[1],
                "partner_facet": partner[2],
                "gap": gap,
                "measure": measure,
            }
        )
        area += measure
        max_gap = max(max_gap, abs(gap))
    from ._regions import Region

    result = {
        "mesh": _make_interface_mesh(master_mesh, output_facets),
        "side_a": Region(
            "interface:side_a",
            "side",
            np.asarray(side_a, dtype=np.int64).reshape(-1, 2),
        ),
        "side_b": Region(
            "interface:side_b",
            "side",
            np.asarray(side_b, dtype=np.int64).reshape(-1, 2),
        ),
        "report": {
            "num_pairs": len(matches),
            "area": area,
            "max_gap": max_gap,
            "unmatched_a": len(facets_a) - len(matched_a),
            "unmatched_b": len(facets_b) - len(matched_b),
        },
    }
    return result


def _contact_facets(mesh, mask):
    bases = block_bases(mesh.cells)
    boundary_mask = np.array(mask, copy=True)
    facets = []
    for base, block in zip(bases[:-1], mesh.cells):
        if block.type.startswith("polyhedron"):
            continue
        dim = topological_dimension.get(block.type, 0)
        if dim != 2 or np.asarray(mesh.points).shape[1] < 3:
            continue
        data = block.data
        for row in range(len(data)):
            cell = int(base + row)
            if not mask[cell]:
                continue
            nodes = np.asarray(data[row], dtype=np.int64).reshape(-1)
            if block.type.startswith("VTK_LAGRANGE_"):
                boundary_mask[cell] = False
                continue
            if block.type.startswith("triangle"):
                nodes = nodes[:3]
            elif block.type.startswith("quad"):
                nodes = nodes[:4]
            facets.append((nodes.tolist(), cell, 0))
            boundary_mask[cell] = False
    facets.extend(_boundary_facets(mesh, boundary_mask))
    return sorted(facets, key=lambda item: (item[1], item[2], item[0]))


def _contact_pairs_py(
    slave_mesh,
    slave_points,
    master_mesh,
    master_cells,
    tolerance=0.0,
    require_complete=False,
):
    point_region = _resolve_region(slave_mesh, slave_points, "point", "contact_pairs")
    _master_region, mask = _region_cell_mask(master_mesh, master_cells, "contact_pairs")
    facets = _contact_facets(master_mesh, mask)
    if not facets:
        raise ValueError(
            "meshio++: contact_pairs: the master region has no queryable facets"
        )
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError(
            "meshio++: contact_pairs: tolerance must be finite and non-negative"
        )
    if tolerance == 0:
        tolerance = 0.01 * _mean_edge(master_mesh, facets)
    ids, cells, local_facets, subfacets, coords, closest, gaps, normals, unmatched = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for raw_node in np.asarray(point_region.entries, dtype=np.int64).reshape(-1):
        node = int(raw_node)
        if node < 0 or node >= len(slave_mesh.points):
            raise ValueError(
                "meshio++: contact_pairs: slave Point region contains an invalid point id"
            )
        point = _xyz(slave_mesh, node)
        projections = [_project(master_mesh, point, facet) for facet in facets]
        projections = [p for p in projections if p is not None]
        best = (
            min(projections, key=lambda p: (p["distance_sq"], p["cell"], p["facet"]))
            if projections
            else None
        )
        ids.append(node)
        if best is None or best["distance_sq"] > tolerance * tolerance:
            cells.append(-1)
            local_facets.append(-1)
            subfacets.append(-1)
            coords.append([0.0, 0.0, 0.0])
            closest.append([0.0, 0.0, 0.0])
            gaps.append(0.0)
            normals.append([0.0, 0.0, 0.0])
            unmatched.append(node)
        else:
            cells.append(best["cell"])
            local_facets.append(best["facet"])
            subfacets.append(best["subfacet"])
            coords.append(best["weights"].tolist())
            closest.append(best["point"].tolist())
            gaps.append(float(np.dot(point - best["point"], best["normal"])))
            normals.append(best["normal"].tolist())
    if require_complete and unmatched:
        raise ValueError(
            f"meshio++: contact_pairs: {len(unmatched)} slave point(s) have no master facet within tolerance"
        )
    return {
        "slave_point": np.asarray(ids, dtype=np.int64),
        "master_cell": np.asarray(cells, dtype=np.int64),
        "master_facet": np.asarray(local_facets, dtype=np.int64),
        "master_subfacet": np.asarray(subfacets, dtype=np.int64),
        "local_coordinates": np.asarray(coords, dtype=np.float64).reshape(-1, 3),
        "closest_point": np.asarray(closest, dtype=np.float64).reshape(-1, 3),
        "gap": np.asarray(gaps, dtype=np.float64),
        "normal": np.asarray(normals, dtype=np.float64).reshape(-1, 3),
        "unmatched": np.asarray(unmatched, dtype=np.int64),
    }


def _split_interface_py(mesh, side, add_cohesive=False):
    from ._mesh import Mesh

    side_region = _resolve_region(mesh, side, "side", "split_interface")
    if any(block.type.startswith("polyhedron") for block in mesh.cells):
        raise ValueError(
            "meshio++: split_interface: polyhedron cells are not supported"
        )
    bases = block_bases(mesh.cells)
    total = int(bases[-1])
    cell_nodes = []
    for block in mesh.cells:
        cell_nodes.extend(
            [np.asarray(row, dtype=np.int64).reshape(-1).tolist() for row in block.data]
        )
    if np.asarray(side_region.entries).size == 0:
        raise ValueError("meshio++: split_interface: Side region is empty")
    side_entries = np.asarray(side_region.entries, dtype=np.int64).reshape(-1, 2)
    facet_map = defaultdict(list)
    for base, block in zip(bases[:-1], mesh.cells):
        for row in range(len(block.data)):
            cell = int(base + row)
            for nodes, local_facet in _cell_facets(block, row):
                key = tuple(sorted(nodes))
                facet_map[key].append((cell, local_facet, list(nodes), block.type))
    for key in facet_map:
        owners = sorted(facet_map[key], key=lambda item: (item[0], item[1]))
        facet_map[key] = list({owner[0]: owner for owner in owners}.values())
    cut_keys, seed_cells = set(), set()
    for cell, facet in side_entries:
        cell, facet = int(cell), int(facet)
        if cell < 0 or cell >= total:
            raise ValueError(
                "meshio++: split_interface: Side region has an invalid cell id"
            )
        key = None
        block_id = int(np.searchsorted(bases[1:], cell, side="right"))
        row = cell - int(bases[block_id])
        for nodes, local_facet in _cell_facets(mesh.cells[block_id], row):
            if local_facet == facet:
                key = tuple(sorted(nodes))
                break
        if key is None or key not in facet_map:
            raise ValueError(
                "meshio++: split_interface: Side region has an invalid local facet id"
            )
        cut_keys.add(key)
        seed_cells.add(cell)

    incident = [[] for _ in range(len(mesh.points))]
    for cell, nodes in enumerate(cell_nodes):
        for node in nodes:
            if node < 0 or node >= len(mesh.points):
                raise ValueError(
                    "meshio++: split_interface: input connectivity has an invalid point id"
                )
            incident[node].append(cell)
    graphs = [defaultdict(set) for _ in range(len(mesh.points))]
    for node, cells in enumerate(incident):
        for cell in cells:
            graphs[node][cell]
    for key, owners in facet_map.items():
        if key in cut_keys or len(owners) < 2:
            continue
        for i, a in enumerate(owners):
            for b in owners[i + 1 :]:
                common = set(cell_nodes[a[0]]) & set(cell_nodes[b[0]])
                for node in common:
                    graphs[node][a[0]].add(b[0])
                    graphs[node][b[0]].add(a[0])

    component_for, keeper_for, duplicate_map, duplicate_sources = {}, {}, {}, []
    for node, graph in enumerate(graphs):
        unseen = set(graph)
        components = []
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            stack, comp = [seed], []
            while stack:
                cell = stack.pop()
                comp.append(cell)
                for other in graph[cell]:
                    if other in unseen:
                        unseen.remove(other)
                        stack.append(other)
            components.append(sorted(comp))
        components.sort(key=lambda comp: comp[0])
        if not components:
            continue
        scores = [sum(cell in seed_cells for cell in comp) for comp in components]
        keeper = max(
            range(len(components)), key=lambda i: (scores[i], -components[i][0])
        )
        keeper_for[node] = keeper
        for ci, comp in enumerate(components):
            for cell in comp:
                component_for[(node, cell)] = ci
        for ci in range(len(components)):
            if ci != keeper:
                duplicate_map[(node, ci)] = len(mesh.points) + len(duplicate_sources)
                duplicate_sources.append(node)

    def remap(node, cell):
        component = component_for.get((int(node), int(cell)))
        return duplicate_map.get((int(node), component), int(node))

    cohesive = []
    if add_cohesive:
        seen = set()
        for cell, facet in side_entries:
            cell, facet = int(cell), int(facet)
            block_id = int(np.searchsorted(bases[1:], cell, side="right"))
            row = cell - int(bases[block_id])
            chosen = next(
                (f for f in _cell_facets(mesh.cells[block_id], row) if f[1] == facet),
                None,
            )
            if chosen is None:
                continue
            key = tuple(sorted(chosen[0]))
            for other in facet_map[key]:
                if other[0] == cell:
                    continue
                pair_key = tuple(sorted(((cell, facet), (other[0], other[1]))))
                if pair_key in seen:
                    continue
                seen.add(pair_key)
                if len(chosen[0]) not in (2, 3, 4):
                    raise ValueError(
                        "meshio++: split_interface: cohesive cells support line, triangle and quad facets"
                    )
                trace_a = [remap(n, cell) for n in chosen[0]]
                by_node = {n: n for n in other[2]}
                trace_b = [remap(by_node[n], other[0]) for n in chosen[0]]
                etype = {2: "line", 3: "wedge", 4: "hexahedron"}[len(chosen[0])]
                cohesive.append((etype, trace_a, trace_b))
        cohesive.sort(key=lambda item: (item[0], item[1], item[2]))

    point_sources = list(range(len(mesh.points))) + duplicate_sources
    points = np.asarray(mesh.points)[np.asarray(point_sources, dtype=np.int64)].copy()
    cells = []
    for base, block in zip(bases[:-1], mesh.cells):
        if block.type.startswith("polygon"):
            rows = []
            for row, connectivity in enumerate(block.data):
                rows.append(
                    [remap(int(node), int(base + row)) for node in connectivity]
                )
            cells.append((block.type, rows))
        else:
            conn = np.asarray(block.data, dtype=np.int64).copy()
            for row in range(len(conn)):
                conn[row] = [remap(int(node), int(base + row)) for node in conn[row]]
            cells.append((block.type, conn))
    by_type = defaultdict(list)
    for item in cohesive:
        by_type[item[0]].append(item)
    for cell_type in ("line", "wedge", "hexahedron"):
        if by_type[cell_type]:
            width = 2 if cell_type == "line" else 6 if cell_type == "wedge" else 8
            rows = []
            for _type, trace_a, trace_b in by_type[cell_type]:
                rows.append(trace_a if cell_type == "line" else trace_a + trace_b)
            cells.append(
                (cell_type, np.asarray(rows, dtype=np.int64).reshape(-1, width))
            )

    point_data = {
        key: np.asarray(values)[np.asarray(point_sources)].copy()
        for key, values in mesh.point_data.items()
    }
    cell_data = {
        key: [np.array(values, copy=True) for values in blocks]
        for key, blocks in mesh.cell_data.items()
    }
    if by_type["line"]:
        if "cohesive:trace_b" in cell_data:
            raise ValueError(
                "meshio++: split_interface: input already has cell_data 'cohesive:trace_b'"
            )
        for key, arrays in cell_data.items():
            component_shape = arrays[0].shape[1:] if arrays else ()
            dtype = arrays[0].dtype if arrays else np.float64
            for cell_type in ("line", "wedge", "hexahedron"):
                if by_type[cell_type]:
                    arrays.append(
                        np.zeros(
                            (len(by_type[cell_type]),) + component_shape, dtype=dtype
                        )
                    )
        trace_blocks = [
            np.full((len(block.data), 2), -1, dtype=np.int64) for block in mesh.cells
        ]
        for cell_type in ("line", "wedge", "hexahedron"):
            if by_type[cell_type]:
                values = np.full((len(by_type[cell_type]), 2), -1, dtype=np.int64)
                if cell_type == "line":
                    values[:] = [item[2] for item in by_type[cell_type]]
                trace_blocks.append(values)
        cell_data["cohesive:trace_b"] = trace_blocks
    else:
        for arrays in cell_data.values():
            for cell_type in ("line", "wedge", "hexahedron"):
                if by_type[cell_type]:
                    shape = arrays[0].shape[1:] if arrays else ()
                    dtype = arrays[0].dtype if arrays else np.float64
                    arrays.append(
                        np.zeros((len(by_type[cell_type]),) + shape, dtype=dtype)
                    )

    regions = []
    for region in getattr(mesh, "regions", ()) or ():
        copied = region.copy()
        if copied.kind == "point":
            entries = set(
                np.asarray(copied.entries, dtype=np.int64).reshape(-1).tolist()
            )
            entries.update(
                len(mesh.points) + i
                for i, source in enumerate(duplicate_sources)
                if source in entries
            )
            copied.entries = np.asarray(sorted(entries), dtype=np.int64)
        regions.append(copied)
    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data={
            key: np.array(value, copy=True) for key, value in mesh.field_data.items()
        },
        regions=regions,
    ), {
        "num_duplicated_points": len(duplicate_sources),
        "num_cohesive_cells": len(cohesive),
    }


def region_adjacency(mesh, regions=None):
    """Return a mesh of facets shared by pairs of named Cell regions.

    ``regions`` is an ordered sequence of region names or selector dictionaries
    (``{"name", "kind"?, "dim"?, "tag"?}``). When omitted, every Cell region
    is selected in canonical order; if fewer than two Cell regions exist and
    the mesh has multiple cell blocks, each block is treated as a group. Two
    selected groups are adjacent where different source cells share exactly
    the same corner-node set. This is
    conforming adjacency; separately meshed or merely nearby faces are not
    matched.

    The result retains source points and point/field data. Its cells are shared
    facets, with cell data ``interface:region_a``/``interface:region_b`` (indices
    into the selected-region order), ``interface:parent_cell_a``/``_b``,
    ``interface:parent_facet_a``/``_b``, ``interface:shared_count`` and
    ``interface:measure``. Region indices are also represented by Cell regions
    named ``adjacency:<a>:<b>``. In 3-D the measure is fan-triangulated area; in
    2-D it is edge length.
    """
    try:
        from . import _core

        return _core.region_adjacency(mesh, regions)
    except Exception as exc:
        if not core_op_declined(exc, "region_adjacency"):
            raise
    return _region_adjacency_py(mesh, regions)


def _selector_payload(value, default_kind):
    if isinstance(value, Region):
        return {
            "name": value.name,
            "kind": value.kind,
            "dim": value.dim,
            "tag": value.tag,
        }
    if isinstance(value, str):
        return value
    payload = dict(value)
    payload.setdefault("kind", default_kind)
    return payload


def find_interface(
    mesh,
    region_a,
    region_b,
    *,
    mesh_b=None,
    mode="conforming",
    master="a",
    gap_tolerance=0.0,
    angle_tolerance=30.0,
    overlap_tolerance=0.0,
    return_report=False,
):
    """Find conforming or proximity-matched boundary facets between two parts.

    `region_a` and `region_b` are Cell-region names or selectors. `mesh_b=None`
    selects both parts from `mesh`; otherwise B is taken from a second mesh.
    Proximity mode checks facet vertices and centroid, with opposing normals;
    a zero gap tolerance derives one percent of the mean boundary edge length.
    The returned facet mesh stores parent/partner cell and facet ids, gap, and
    measure. With `return_report=True`, returns `(mesh, report)`; the report's
    `side_a`/`side_b` are Side regions on the original input mesh(es), ready for
    `split_interface`.
    """
    selector_a = _selector_payload(region_a, "cell")
    selector_b = _selector_payload(region_b, "cell")
    result = None
    try:
        from . import _core

        result = _core.find_interface(
            mesh,
            selector_a,
            selector_b,
            mesh_b,
            mode,
            master,
            float(gap_tolerance),
            float(angle_tolerance),
            float(overlap_tolerance),
        )
    except Exception as exc:
        if not core_op_declined(exc, "find_interface"):
            raise
    if result is None:
        result = _find_interface_py(
            mesh,
            region_a,
            region_b,
            mesh_b,
            mode,
            master,
            float(gap_tolerance),
            float(angle_tolerance),
            float(overlap_tolerance),
        )
    report = dict(result["report"])
    for key in ("side_a", "side_b"):
        if not isinstance(result[key], Region):
            result[key] = Region(
                "interface:" + key,
                "side",
                np.asarray(result[key], dtype=np.int64).reshape(-1, 2),
            )
        report[key] = result[key]
    return (result["mesh"], report) if return_report else result["mesh"]


def contact_pairs(
    slave_mesh,
    slave_points,
    master_cells,
    *,
    master_mesh=None,
    tolerance=0.0,
    require_complete=False,
):
    """Project a Point region to the closest facets of a Cell region.

    Returns arrays `slave_point`, `master_cell`, `master_facet`, `master_subfacet`,
    `local_coordinates`, `closest_point`, signed `gap`, `normal` and
    `unmatched`. Triangle local coordinates are barycentric weights in the
    selected fan triangle; `master_subfacet` is its zero-based fan ordinal
    (zero for triangles and edges). Edge coordinates are endpoint weights in
    slots 0/1.
    """
    selector_points = _selector_payload(slave_points, "point")
    selector_cells = _selector_payload(master_cells, "cell")
    actual_master = slave_mesh if master_mesh is None else master_mesh
    try:
        from . import _core

        return _core.contact_pairs(
            slave_mesh,
            selector_points,
            actual_master,
            selector_cells,
            float(tolerance),
            bool(require_complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "contact_pairs"):
            raise
    return _contact_pairs_py(
        slave_mesh,
        slave_points,
        actual_master,
        master_cells,
        float(tolerance),
        bool(require_complete),
    )


def split_interface(mesh, side, *, add_cohesive=False, return_report=False):
    """Duplicate point fans along a Side region, optionally adding cohesive cells.

    Polyhedra are refused. Triangle/quad facets add wedge/hexahedron cohesive
    cells. A 2-D cohesive line uses its connectivity for trace A and
    `cohesive:trace_b` (two Int64 cell-data components) for the opposite trace.
    Pass a Side `Region` returned by `find_interface(..., return_report=True)` or
    a selector for a Side region already on `mesh`.
    """
    temporary = mesh
    selector = side
    if isinstance(side, Region):
        if side.kind != "side":
            raise ValueError("meshio++: split_interface: expected a Side region")
        temporary = mesh.copy()
        name = side.name
        if any(r.name == name and r.kind == "side" for r in temporary.regions):
            name = "__split_interface_side"
            side = Region(name, "side", side.entries, side.dim, side.tag)
        temporary.regions.append(side.copy())
        selector = {"name": name, "kind": "side", "dim": side.dim, "tag": side.tag}
    elif isinstance(side, dict):
        selector = dict(side)
        selector.setdefault("kind", "side")
    result = None
    try:
        from . import _core

        result = _core.split_interface(
            temporary, _selector_payload(selector, "side"), bool(add_cohesive)
        )
    except Exception as exc:
        if not core_op_declined(exc, "split_interface"):
            raise
    if result is None:
        out, report = _split_interface_py(temporary, selector, bool(add_cohesive))
    else:
        out = result["mesh"]
        report = {
            "num_duplicated_points": int(result["num_duplicated_points"]),
            "num_cohesive_cells": int(result["num_cohesive_cells"]),
        }
    return (out, report) if return_report else out


__all__ = ["region_adjacency", "find_interface", "contact_pairs", "split_interface"]
