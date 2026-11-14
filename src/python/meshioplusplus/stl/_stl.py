"""
I/O for the STL format, cf.
<https://en.wikipedia.org/wiki/STL_(file_format)>.
"""

from __future__ import annotations

import os
from typing import Union

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._regions import Region


def read(filename):
    with open_file(filename, "rb") as f:
        # Checking if the file is ASCII format is normally done by checking if the
        # first 5 characters of the header is "solid".
        # ```
        # header = f.read(80).decode()
        # ```
        # Unfortunately, there are mesh files out there which are binary and still put
        # "solid" there.
        # A suggested alternative is to pretend the file is binary, read the
        # num_triangles and see if it matches the file size
        # (https://stackoverflow.com/a/7394842/353337).
        filesize_bytes = os.path.getsize(filename)
        if filesize_bytes < 80:
            return _read_ascii(f)

        f.read(80)
        num_triangles = np.fromfile(f, count=1, dtype="<u4")[0]
        # for each triangle, one has 3 float32 (facet normal), 9 float32 (facet),
        # and 1 int16 (attribute count), 50 bytes in total
        if 84 + num_triangles * 50 == filesize_bytes:
            return _read_binary(f, num_triangles)

        # rewind and skip header
        f.seek(0)
        header = f.readline().decode(errors="replace").strip()
        return _read_ascii(f, header)


# np.loadtxt is super slow
# Code adapted from <https://stackoverflow.com/a/8964779/353337>.
def iter_loadtxt(
    infile,
    skiprows: int = 0,
    comments: Union[str, tuple[str, ...]] = "#",
    dtype=float,
    usecols: Union[tuple[int], None] = None,
    solids: Union[list, None] = None,
):
    nrows = 0

    def iter_func():
        nonlocal nrows
        items = None
        for _ in range(skiprows):
            try:
                next(infile)
            except StopIteration:
                raise ReadError("EOF Skipped too many rows")

        for line in infile:
            line = line.decode().strip()
            if line.startswith(comments):
                if solids is not None and line.startswith("solid"):
                    solids.append((line[5:].strip(), nrows // 4))
                continue
            nrows += 1
            # remove all text
            items = line.split()[-3:]
            usecols_ = range(len(items)) if usecols is None else usecols
            for idx in usecols_:
                yield dtype(items[idx])

        if items is None:
            iter_loadtxt.rowlength = 3
            return

        iter_loadtxt.rowlength = len(items) if usecols is None else len(usecols)

    data = np.fromiter(iter_func(), dtype=dtype)
    return data.reshape((-1, iter_loadtxt.rowlength))


def _read_ascii(f, header=None):
    # The file has the form
    # ```
    # solid foo
    #   facet normal 0.455194 -0.187301 -0.870469
    #    outer loop
    #     vertex 266.36 234.594 14.6145
    #     vertex 268.582 234.968 15.6956
    #     vertex 267.689 232.646 15.7283
    #    endloop
    #   endfacet
    #   # [...] more facets [...]
    # endsolid
    # ```
    # In the interest of speed, don't verify the format and instead just skip the text.

    # TODO Pandas is MUCH faster than numpy for i/o, see
    # <https://stackoverflow.com/a/18260092/353337>.
    # import pandas
    # data = pandas.read_csv(
    #     f,
    #     skiprows=lambda row: row == 0 or (row - 1) % 7 in [0, 1, 5, 6],
    #     skipfooter=1,
    #     usecols=(1, 2, 3),
    # )

    # np.loadtxt is super slow
    # data = np.loadtxt(
    #     f,
    #     comments=["solid", "facet", "outer loop", "endloop", "endfacet", "endsolid"],
    #     usecols=(1, 2, 3),
    # )
    solids = []
    if header is not None and header.startswith("solid"):
        solids.append((header[5:].strip(), 0))
    data = iter_loadtxt(
        f,
        comments=("solid", "outer loop", "endloop", "endfacet", "endsolid"),
        # usecols=(1, 2, 3),
        solids=solids,
    )

    if data.shape[0] % 4 != 0:
        raise ReadError()

    # split off the facet normals
    facet_rows = np.zeros(len(data), dtype=bool)
    facet_rows[0::4] = True
    facet_normals = data[facet_rows]
    data = data[~facet_rows]

    if data.shape[0] == 0:
        points = []
        cells = {}
        cell_data = {}
    else:
        facets = np.split(data, data.shape[0] // 3)
        points, cells = data_from_facets(facets)
        cell_data = {"facet_normals": [facet_normals]}

    regions = _solid_regions(solids, len(data) // 3) if len(cells) else []
    return Mesh(points, cells, cell_data=cell_data, regions=regions)


def _solid_regions(solids, num_facets):
    """Several ``solid`` records of an ASCII file as cell regions named after
    them (``solid_<k>`` when unnamed, ``_<k>`` appended to repeated names)."""
    if len(solids) < 2:
        return []
    regions = []
    used = set()
    for k, (name, first) in enumerate(solids):
        last = solids[k + 1][1] if k + 1 < len(solids) else num_facets
        name = name or f"solid_{k}"
        while name in used:
            name = f"{name}_{k}"
        used.add(name)
        regions.append(Region(name, "cell", np.arange(first, last, dtype=np.int64)))
    return regions


def data_from_facets(facets):
    # Now, all facets contain the point coordinate. Try to identify individual points
    # and build the data arrays.
    if len(facets) == 0:
        points = np.empty((0, 3), dtype=float)
        cells = []
    else:
        pts = np.concatenate(facets)
        # TODO equip `unique()` with a tolerance
        # Use return_index so we can use sort on `idx` such that the order is
        # preserved; see <https://stackoverflow.com/a/15637512/353337>.
        _, idx, inv = np.unique(pts, axis=0, return_index=True, return_inverse=True)
        k = np.argsort(idx)
        points = pts[idx[k]]
        inv_k = np.argsort(k)
        cells = [CellBlock("triangle", inv_k[inv].reshape(-1, 3))]
    return points, cells


def _read_binary(f, num_triangles: int):
    # for each triangle, one has 3 float32 (facet normal), 9 float32 (facet), and 1
    # int16 (attribute count)
    out = np.fromfile(
        f,
        count=num_triangles,
        dtype=np.dtype(
            [("normal", "<f4", (3,)), ("facet", "<f4", (3, 3)), ("attr count", "<i2")]
        ),
    )
    # discard normals, attribute count
    facets = out["facet"]
    # if not np.all(out["attr count"] == 0):
    #     print(out["attr count"])
    #     raise ReadError("Nonzero attr count")

    points, cells = data_from_facets(facets)
    cell_data = {}
    attr = out["attr count"]
    if len(attr) and np.any(attr):
        cell_data["stl:attribute"] = [attr.astype("<u2")]
    return Mesh(points, cells, cell_data=cell_data)


def _skin_triangles(mesh):
    """The linearized boundary skin of a volume mesh, quads triangulated."""
    from .._skin import _extract_skin_py

    skin = _extract_skin_py(mesh, linearize=True)
    tri = []
    for block in skin.cells:
        if block.type == "triangle":
            tri.append(np.asarray(block.data, dtype=np.int64))
        else:  # quad -> (0,1,2), (0,2,3)
            quad = np.asarray(block.data, dtype=np.int64)
            tri.append(quad[:, [0, 1, 2]])
            tri.append(quad[:, [0, 2, 3]])
    return Mesh(skin.points, [("triangle", np.concatenate(tri))])


def write(filename, mesh, binary=False, skin=True):
    from .._skin import _has_skinnable_cells, _is_volume_type

    if skin and _has_skinnable_cells(mesh):
        dropped = sum(1 for block in mesh.cells if not _is_volume_type(block.type))
        if dropped > 0:
            warn(
                "STL: writing the extracted skin of the volume cells; "
                f"{dropped} pre-existing non-volume cell block(s) dropped "
                "(pass skin=False for the legacy behavior)."
            )
            _provenance.note(
                "cells-dropped",
                f"{dropped} non-volume cell block(s) dropped in favour of the "
                "extracted skin of the volume cells",
            )
        mesh = _skin_triangles(mesh)

    if "triangle" not in {block.type for block in mesh.cells}:
        warn("STL can only write triangle cells. No triangle cells found.")
    if len(mesh.cells) > 1:
        invalid = {block.type for block in mesh.cells if block.type != "triangle"}
        invalid = ", ".join(invalid)
        warn(f"STL can only write triangle cells. Discarding {invalid}.")

    if mesh.points.shape[1] == 2:
        warn(
            "STL requires 3D points, but 2D points given. Appending 0 third component."
        )
        points = np.column_stack([mesh.points, np.zeros_like(mesh.points[:, 0])])
    else:
        points = mesh.points

    pts = points[mesh.get_cells_type("triangle")]
    if "facet_normals" in mesh.cell_data:
        normals = mesh.get_cell_data("facet_normals", "triangle")
    else:
        normals = np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0])
        nrm = np.sqrt(np.einsum("ij,ij->i", normals, normals))
        normals = (normals.T / nrm).T

    solids = _solids_to_write(mesh) if not binary else None
    if binary:
        attr = None
        if "stl:attribute" in mesh.cell_data:
            attr = mesh.get_cell_data("stl:attribute", "triangle")
        _write_binary(filename, pts, normals, attr)
    else:
        _write_ascii(filename, pts, normals, solids)


def _solids_to_write(mesh):
    """``[(name, triangle positions)]`` when the mesh has two or more cell
    regions over its triangles, with the triangles no region holds in a final
    unnamed solid; else ``None`` (one solid)."""
    regions = [r for r in getattr(mesh, "regions", []) if r.kind == "cell"]
    if len(regions) < 2:
        return None
    position = {}
    n = 0
    base = 0
    for block in mesh.cells:
        if block.type == "triangle":
            for i in range(len(block.data)):
                position[base + i] = n
                n += 1
        base += len(block.data)
    out = []
    held = set()
    for r in regions:
        tris = [position[int(c)] for c in r.entries if int(c) in position]
        held.update(tris)
        out.append((r.name, np.array(sorted(tris), dtype=np.int64)))
    rest = np.array([i for i in range(n) if i not in held], dtype=np.int64)
    if len(rest):
        out.append(("", rest))
    return out


def _write_ascii(filename, pts, normals, solids=None):
    with open_file(filename, "w") as fh:
        if solids is None:
            _write_solid(fh, "", pts, normals)
        else:
            for name, idx in solids:
                _write_solid(fh, name, pts[idx], normals[idx])


def _write_solid(fh, name, pts, normals):
    fh.write(f"solid {name}\n" if name else "solid\n")
    for local_pts, normal in zip(pts, normals):
        out = (
            "\n".join(
                [
                    "facet normal {} {} {}".format(*normal),
                    " outer loop",
                    "  vertex {} {} {}".format(*local_pts[0]),
                    "  vertex {} {} {}".format(*local_pts[1]),
                    "  vertex {} {} {}".format(*local_pts[2]),
                    " endloop",
                    "endfacet",
                ]
            )
            + "\n"
        )
        fh.write(out)
    fh.write(f"endsolid {name}\n" if name else "endsolid\n")


def _write_binary(filename, pts, normals, attr=None):
    with open_file(filename, "wb") as fh:
        # 80-byte header data, NUL-padded to match the C++ writer exactly.
        header = (
            _provenance.lines(_provenance.SlotTier.BOUNDED)[0].encode().ljust(80, b"\0")
        )
        fh.write(header)

        fh.write(np.array(len(pts)).astype("<u4"))

        dtype = np.dtype(
            [
                ("normal", ("<f4", 3)),
                ("points", ("<f4", (3, 3))),
                ("attr", "<u2"),
            ]
        )
        a = np.empty(len(pts), dtype=dtype)
        a["normal"] = normals
        a["points"] = pts
        a["attr"] = 0 if attr is None else attr
        a.tofile(fh)


# NOTE: format registration now lives in meshioplusplus/stl/__init__.py, which wraps the
# reader/writer below with the C++-backed fast paths.
