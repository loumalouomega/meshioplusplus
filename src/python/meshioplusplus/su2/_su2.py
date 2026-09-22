"""
I/O SU2 mesh format
<https://su2code.github.io/docs_v7/Mesh-File/>

Single-file multizone (`NZONE=` then one `IZONE= i` section per zone, each a
complete NDIME/NPOIN/NELEM/NMARK mesh) is read and written alongside the
single-zone layout: zones are never welded (each keeps its own point
numbering, offset when concatenated), `cell_data["su2:zone"]` names each
cell's zone, and a marker whose `MARKER_TAG` is a string (not just a plain
integer) gets its own `Region` -- `"zone_<i>/<name>"` when multizone,
plain `"<name>"` for a single-zone file -- so the name survives instead of
collapsing into an anonymous auto-incremented `su2:tag`. See
`doc/formats/su2.md` and the C++ twin, `formats/su2.cpp`.
"""

import numpy as np

from .._common import _pick_first_int_data, warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._regions import Region

# follows VTK conventions
su2_type_to_numnodes = {
    3: 2,  # line
    5: 3,  # triangle
    9: 4,  # quad
    10: 4,  # tetra
    12: 8,  # hexahedron
    13: 6,  # wedge
    14: 5,  # pyramid
}
su2_to_meshio_type = {
    3: "line",
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
}
meshio_to_su2_type = {
    "line": 3,
    "triangle": 5,
    "quad": 9,
    "tetra": 10,
    "hexahedron": 12,
    "wedge": 13,
    "pyramid": 14,
}


def read(filename):
    with open_file(filename, "r") as f:
        lines = f.read().split("\n")
    return _read_lines(lines)


def _next_key(lines, li):
    """The next non-blank, non-comment line's ``KEY``, or ``None`` at EOF /
    a line with no ``=`` -- without consuming it."""
    n = len(lines)
    while li < n:
        line = lines[li].strip()
        if not line or line[0] == "%":
            li += 1
            continue
        if "=" not in line:
            return None, li
        return line.split("=", 1)[0].strip(), li
    return None, li


def _read_lines(lines):
    li = 0

    # NZONE= (if present) is always the first key: a single-file multizone
    # mesh's own header, before the first (implicit) IZONE= 1.
    nzone = 1
    multizone = False
    key, li0 = _next_key(lines, li)
    if key == "NZONE":
        rest = lines[li0].split("=", 1)[1]
        nzone = int(rest.split()[0])
        multizone = True
        li = li0 + 1

    zones = []
    for z in range(nzone):
        # Skip blank/comment lines and, for a multizone file, the NZONE=/
        # IZONE= markers themselves; the zone body starts at NDIME.
        while True:
            key, li0 = _next_key(lines, li)
            if key in ("NZONE", "IZONE"):
                li = li0 + 1
                continue
            break
        zone, li = _read_zone_body(lines, li, z)
        zones.append(zone)
    if multizone and len(zones) != nzone:
        raise ReadError(f"SU2: NZONE={nzone} but found {len(zones)} IZONE section(s)")

    dim = zones[0]["dim"] if zones else 0
    total_points = sum(len(z["points"]) for z in zones)
    points = (
        np.concatenate([z["points"] for z in zones], axis=0)
        if zones
        else np.zeros((0, dim))
    )
    offset = 0
    blocks = []  # dicts: type, conn (n, k), tag (n,), zone (n,)
    for z in zones:
        for blk in z["blocks"]:
            blk["conn"] = blk["conn"] + offset
        blocks.extend(z["blocks"])
        offset += len(z["points"])
    assert offset == total_points

    blocks = _merge_by_type(blocks)

    cells = []
    tags = []
    zone_ids = []
    for blk in blocks:
        cells.append(CellBlock(blk["type"], blk["conn"]))
        tags.append(blk["tag"].astype(np.int32))
        if multizone:
            zone_ids.append(blk["zone"].astype(np.int32))
    cell_data = {"su2:tag": tags}
    if multizone:
        cell_data["su2:zone"] = zone_ids

    mesh = Mesh(points, cells, cell_data=cell_data)

    # Regions: one "zone_<i>" per zone (multizone only) and one per named
    # marker. A purely numeric marker gets no region: su2:tag already
    # carries it losslessly.
    block_bases = np.cumsum([0] + [len(c.data) for c in cells])
    zone_cells = {}
    marker_cells = {}
    for bi, blk in enumerate(blocks):
        base = int(block_bases[bi])
        for r in range(len(blk["tag"])):
            g = base + r
            zc = int(blk["zone"][r])
            tg = int(blk["tag"][r])
            zone_cells.setdefault(zc, []).append(g)
            marker_cells.setdefault((zc, tg), []).append(g)
    if multizone:
        for zc in sorted(zone_cells):
            mesh.regions.append(
                Region(f"zone_{zc}", "cell", np.asarray(zone_cells[zc], dtype=np.int64))
            )
    for zc, tg in sorted(marker_cells):
        name = zones[zc]["marker_names"].get(tg)
        if name is None:
            continue
        rname = f"zone_{zc}/{name}" if multizone else name
        mesh.regions.append(
            Region(rname, "cell", np.asarray(marker_cells[(zc, tg)], dtype=np.int64))
        )
    return mesh


def _read_zone_body(lines, li, zone):
    """One zone's own NDIME/NPOIN/NELEM/NMARK/MARKER_* body -- everything a
    standalone (single-zone) file also has. Stops at the next ``IZONE=`` (the
    next zone) or end of file, without consuming it."""
    n = len(lines)
    dim = 0
    points = np.zeros((0, 0))
    blocks = []
    marker_names = {}  # tag -> string name, only for a non-numeric MARKER_TAG
    next_tag_id = 0
    current_tag = 0
    expected_nmarkers = 0
    markers_found = 0

    while li < n:
        line = lines[li].strip()
        if not line or line[0] == "%":
            li += 1
            continue
        if "=" not in line:
            li += 1
            continue
        name, rest_of_line = line.split("=", 1)
        name = name.strip()
        if name == "IZONE":
            break  # the next zone: let the caller consume it
        li += 1

        if name == "NDIME":
            dim = int(rest_of_line)
            if dim not in (2, 3):
                raise ReadError(f"Invalid dimension value {line}")
        elif name == "NPOIN":
            first_line = np.array(lines[li].split(), dtype="f8")
            li += 1
            extra_columns = first_line.shape[0] - dim
            num_verts = int(rest_of_line.split()[0]) - 1
            rest_tokens = " ".join(lines[li : li + num_verts]).split()
            li += num_verts
            rest_pts = np.array(rest_tokens, dtype="f8").reshape(
                num_verts, dim + extra_columns
            )
            if extra_columns > 0:
                first_line = first_line[:-extra_columns]
                rest_pts = rest_pts[:, :-extra_columns]
            points = np.vstack([first_line, rest_pts])
        elif name in ("NELEM", "MARKER_ELEMS"):
            num_elems = int(rest_of_line)
            elem_lines = lines[li : li + num_elems]
            li += num_elems
            first_tokens = elem_lines[0].split()
            nnodes = su2_type_to_numnodes[int(first_tokens[0])]
            has_extra_column = len(first_tokens) == nnodes + 2
            cell_array = np.fromiter(" ".join(elem_lines).split(), dtype="i8")
            cells_, _ = _translate_cells(cell_array, has_extra_column)
            tag = 0 if name == "NELEM" else current_tag
            for eltype, conn in cells_.items():
                blocks.append(
                    {
                        "type": eltype,
                        "conn": conn,
                        "tag": np.full(len(conn), tag, dtype=np.int64),
                        "zone": np.full(len(conn), zone, dtype=np.int64),
                    }
                )
        elif name == "NMARK":
            expected_nmarkers = int(rest_of_line)
        elif name == "MARKER_TAG":
            text = rest_of_line.strip()
            try:
                current_tag = int(text)
            except ValueError:
                next_tag_id += 1
                current_tag = next_tag_id
                marker_names[current_tag] = text
            markers_found += 1

    if markers_found != expected_nmarkers:
        warn(
            f"expected {expected_nmarkers} markers according to NMARK value "
            f"but found only {markers_found}"
        )
    return (
        {"dim": dim, "points": points, "blocks": blocks, "marker_names": marker_names},
        li,
    )


def _merge_by_type(blocks):
    """Merges same-type blocks (across every zone, volume or boundary alike)
    into one block per type -- a straight generalization of the single-zone
    reader's own "merge same-type boundary blocks" pass, which becomes a
    no-op there (a single NELEM/MARKER_ELEMS call already groups by type)."""
    by_type = {}
    order = []
    for blk in blocks:
        if blk["type"] not in by_type:
            by_type[blk["type"]] = []
            order.append(blk["type"])
        by_type[blk["type"]].append(blk)
    merged = []
    for t in order:
        parts = by_type[t]
        merged.append(
            {
                "type": t,
                "conn": np.concatenate([p["conn"] for p in parts], axis=0),
                "tag": np.concatenate([p["tag"] for p in parts]),
                "zone": np.concatenate([p["zone"] for p in parts]),
            }
        )
    return merged


def _translate_cells(data, has_extra_column=False):
    # adapted from _vtk.py
    # Translate input array  into the cells dictionary.
    # `data` is a one-dimensional vector with
    # (vtk cell type, p0, p1, ... ,pk, vtk cell type, p10, p11, ..., p1k, ...

    entry_offset = 1
    if has_extra_column:
        entry_offset += 1

    # Collect types into bins.
    # See <https://stackoverflow.com/q/47310359/353337> for better
    # alternatives.
    types = []
    i = 0
    while i < len(data):
        types.append(data[i])
        i += su2_type_to_numnodes[data[i]] + entry_offset

    types = np.array(types)
    bins = {u: np.where(types == u)[0] for u in np.unique(types)}

    # Deduct offsets from the cell types. This is much faster than manually
    # going through the data array. Slight disadvantage: This doesn't work for
    # cells with a custom number of points.
    numnodes = np.empty(len(types), dtype=int)
    for tpe, idx in bins.items():
        numnodes[idx] = su2_type_to_numnodes[tpe]
    offsets = np.cumsum(numnodes + entry_offset) - (numnodes + entry_offset)

    cells = {}
    cell_data = {}
    for tpe, b in bins.items():
        meshio_type = su2_to_meshio_type[tpe]
        nnodes = su2_type_to_numnodes[tpe]
        indices = np.add.outer(offsets[b], np.arange(1, nnodes + 1))
        cells[meshio_type] = data[indices]

    return cells, cell_data


def _write_zone_body(f, mesh, dim, cell_indices, tag_key, marker_names, subset):
    """Writes one zone's own NDIME/NPOIN/NELEM/NMARK/MARKER_* body --
    everything a standalone (single-zone) file also has -- restricted to
    ``cell_indices`` (list of ``(block_index, row)``). ``subset=False`` (a
    whole-mesh, single-zone write): every point is written in its original
    order, connectivity untouched. ``subset=True`` (one zone of a multizone
    write): points are remapped to a dense, zone-local 0-based numbering in
    first-appearance order, so a multizone write never shares a point index
    between zones."""
    if dim == 2:
        vtypes, btypes = ("triangle", "quad"), ("line",)
    else:
        vtypes, btypes = ("tetra", "hexahedron", "wedge", "pyramid"), (
            "triangle",
            "quad",
        )

    if subset:
        remap = {}
        old_ids = []

        def remap_of(old_id):
            old_id = int(old_id)
            idx = remap.get(old_id)
            if idx is None:
                idx = len(old_ids)
                remap[old_id] = idx
                old_ids.append(old_id)
            return idx

        for bi, row in cell_indices:
            for v in mesh.cells[bi].data[row]:
                remap_of(v)
    else:
        old_ids = list(range(len(mesh.points)))

        def remap_of(old_id):
            return int(old_id)

    f.write(f"NDIME= {dim}\n".encode())
    f.write(f"NPOIN= {len(old_ids)}\n".encode())
    np.savetxt(f, mesh.points[old_ids])

    def tag_of(bi, row):
        if tag_key is None:
            return 1
        return int(mesh.cell_data[tag_key][bi][row])

    def write_cells(indices):
        for bi, row in indices:
            conn = mesh.cells[bi].data[row]
            remapped = [remap_of(v) for v in conn]
            f.write(
                (
                    f"{meshio_to_su2_type[mesh.cells[bi].type]} "
                    + " ".join(str(v) for v in remapped)
                    + "\n"
                ).encode()
            )

    volume = [(bi, r) for bi, r in cell_indices if mesh.cells[bi].type in vtypes]
    f.write(f"NELEM= {len(volume)}\n".encode())
    write_cells(volume)

    boundary = [(bi, r) for bi, r in cell_indices if mesh.cells[bi].type in btypes]
    tag_counts = {}
    for bi, r in boundary:
        tag_counts[tag_of(bi, r)] = tag_counts.get(tag_of(bi, r), 0) + 1

    f.write(f"NMARK= {len(tag_counts)}\n".encode())
    for tag in tag_counts:
        name = marker_names.get(tag, str(tag))
        f.write(f"MARKER_TAG= {name}\n".encode())
        f.write(f"MARKER_ELEMS= {tag_counts[tag]}\n".encode())
        write_cells([(bi, r) for bi, r in boundary if tag_of(bi, r) == tag])


def write(filename, mesh):
    dim = mesh.points.shape[1]
    for cell_block in mesh.cells:
        if cell_block.type not in meshio_to_su2_type:
            warn(
                f".su2 does not support elements of type {cell_block.type}.\nSkipping ..."
            )

    # su2:zone drives zone splitting below, not marker tags: exclude it here
    # so it is never reported as an unwritable "other" int array.
    tag_candidates = {k: v for k, v in mesh.cell_data.items() if k != "su2:zone"}
    tag_key, other = _pick_first_int_data(tag_candidates)
    if tag_key and other:
        warn(
            "su2 file format can only write one cell data array. "
            "Picking {}, skipping {}.".format(tag_key, ", ".join(other))
        )

    zone_key = "su2:zone" if "su2:zone" in mesh.cell_data else None
    by_zone = {}
    for bi, cell_block in enumerate(mesh.cells):
        for r in range(len(cell_block.data)):
            zone = int(mesh.cell_data[zone_key][bi][r]) if zone_key else 0
            by_zone.setdefault(zone, []).append((bi, r))
    multizone = len(by_zone) > 1

    # Marker names: a region "zone_<i>/<name>" (multizone) or "<name>"
    # (single-zone) supplies a marker's text where the mesh has one.
    marker_names = {}  # zone -> {tag: name}
    for region in mesh.regions:
        if region.kind != "cell" or len(region.entries) == 0:
            continue
        bi, row = _global_to_block_row(mesh, int(region.entries[0]))
        if bi is None:
            continue
        zone = int(mesh.cell_data[zone_key][bi][row]) if zone_key else 0
        tag = int(mesh.cell_data[tag_key][bi][row]) if tag_key else 1
        name = region.name
        prefix = f"zone_{zone}/"
        if multizone and name.startswith(prefix):
            name = name[len(prefix) :]
        elif multizone:
            continue  # a "zone_<i>" region itself, or one from another zone
        marker_names.setdefault(zone, {})[tag] = name

    with open_file(filename, "wb") as f:
        if not multizone:
            all_cells = next(iter(by_zone.values())) if by_zone else []
            _write_zone_body(
                f, mesh, dim, all_cells, tag_key, marker_names.get(0, {}), subset=False
            )
            return
        f.write(f"NZONE= {len(by_zone)}\n".encode())
        for ordinal, zone in enumerate(sorted(by_zone), start=1):
            f.write(f"\nIZONE= {ordinal}\n".encode())
            _write_zone_body(
                f,
                mesh,
                dim,
                by_zone[zone],
                tag_key,
                marker_names.get(zone, {}),
                subset=True,
            )


def _global_to_block_row(mesh, global_index):
    """``(block_index, row)`` of the global (block-major) cell index
    ``global_index``, or ``(None, 0)`` if out of range."""
    base = 0
    for bi, cell_block in enumerate(mesh.cells):
        count = len(cell_block.data)
        if global_index < base + count:
            return bi, global_index - base
        base += count
    return None, 0


# NOTE: format registration now lives in meshioplusplus/su2/__init__.py, which wraps the
# reader/writer above with the C++-backed fast paths.
