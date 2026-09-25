# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""A small, standalone ParaView series writer: one ``.vtu`` per step and part,
indexed by a ``.pvd``.

The vendor-route scripts under ``contrib/`` run inside a vendor's own Python
(``abaqus python``, Marc's PyPost), which has numpy but neither h5py nor
meshio++. This module needs numpy only, runs on Python 2.7 and 3, and writes
files that meshio++ reads as a time sequence (``meshioplusplus.read_sequence``)
and ParaView opens directly. Copy it next to the script that imports it.

Layout: ``run.pvd`` plus a sibling ``run/`` directory holding
``run_<step>_<part>.vtu``; each ``.vtu`` is an UnstructuredGrid with base64
binary arrays (UInt64 headers) and its time as ``TimeValue`` field data.
"""

from __future__ import division, print_function

import base64
import os
import struct
from xml.sax.saxutils import quoteattr

import numpy as np

# numpy dtype kind/size -> VTK type name.
_VTK_TYPES = {
    ("f", 4): "Float32",
    ("f", 8): "Float64",
    ("i", 1): "Int8",
    ("i", 2): "Int16",
    ("i", 4): "Int32",
    ("i", 8): "Int64",
    ("u", 1): "UInt8",
    ("u", 2): "UInt16",
    ("u", 4): "UInt32",
    ("u", 8): "UInt64",
}


def _vtk_type(array):
    key = (array.dtype.kind, array.dtype.itemsize)
    if key not in _VTK_TYPES:
        raise TypeError("no VTK type for dtype %s" % array.dtype)
    return _VTK_TYPES[key]


def _encode(array):
    """Base64 of a UInt64 byte count followed by the little-endian bytes."""
    data = np.ascontiguousarray(array)
    data = data.astype(data.dtype.newbyteorder("<"), copy=False).tobytes()
    return base64.b64encode(struct.pack("<Q", len(data)) + data).decode("ascii")


def _data_array(name, array, extra=""):
    array = np.asarray(array)
    if array.dtype == np.bool_:
        array = array.astype(np.uint8)
    ncomp = 1 if array.ndim < 2 else int(np.prod(array.shape[1:]))
    return (
        '<DataArray type="%s" Name=%s NumberOfComponents="%d" format="binary"%s>'
        "%s</DataArray>\n"
        % (_vtk_type(array), quoteattr(name), ncomp, extra, _encode(array.ravel()))
    )


def write_vtu(path, points, blocks, point_data=None, cell_data=None, time=None):
    """Write one UnstructuredGrid.

    ``points`` is (n, 3); ``blocks`` is a list of ``(vtk_type, connectivity)``
    with a 0-based (cells, nodes) connectivity per block; ``cell_data`` values
    are one array over every block's cells, in block order.
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2:
        raise ValueError("points must be (n, dim)")
    if points.shape[1] < 3:
        points = np.hstack([points, np.zeros((len(points), 3 - points.shape[1]))])
    conn, offsets, types = [], [], []
    end = 0
    for vtk_type, block in blocks:
        block = np.asarray(block, dtype=np.int64)
        if block.size == 0:
            continue
        rows, npc = block.shape
        conn.append(block.ravel())
        offsets.append(end + npc * np.arange(1, rows + 1, dtype=np.int64))
        end += rows * npc
        types.append(np.full(rows, vtk_type, dtype=np.uint8))
    conn = np.concatenate(conn) if conn else np.zeros(0, np.int64)
    offsets = np.concatenate(offsets) if offsets else np.zeros(0, np.int64)
    types = np.concatenate(types) if types else np.zeros(0, np.uint8)

    out = [
        '<?xml version="1.0"?>\n',
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian"'
        ' header_type="UInt64">\n',
        "<UnstructuredGrid>\n",
    ]
    if time is not None:
        out.append("<FieldData>\n")
        out.append(
            _data_array(
                "TimeValue", np.array([time], dtype=np.float64), ' NumberOfTuples="1"'
            )
        )
        out.append("</FieldData>\n")
    out.append(
        '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n' % (len(points), len(types))
    )
    out.append("<Points>\n" + _data_array("Points", points) + "</Points>\n")
    out.append("<Cells>\n")
    out.append(_data_array("connectivity", conn))
    out.append(_data_array("offsets", offsets))
    out.append(_data_array("types", types))
    out.append("</Cells>\n")
    for tag, data in (("PointData", point_data), ("CellData", cell_data)):
        if data:
            out.append("<%s>\n" % tag)
            for name in sorted(data):
                out.append(_data_array(name, data[name]))
            out.append("</%s>\n" % tag)
    out.append("</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
    with open(path, "w") as f:
        f.write("".join(out))


def _time_text(t):
    return repr(float(t))


class PvdSeries(object):
    """Collects steps and parts, writes one ``.vtu`` each, then the ``.pvd``.

    >>> series = PvdSeries("run.pvd")
    >>> series.add(0.0, points, [(10, tets)], point_data={"U": u}, part=0, name="PART-1")
    >>> series.close()
    """

    def __init__(self, path):
        self.path = path
        stem = os.path.splitext(os.path.basename(path))[0]
        self.stem = stem
        self.dir = os.path.join(os.path.dirname(os.path.abspath(path)), stem)
        if not os.path.isdir(self.dir):
            os.makedirs(self.dir)
        self.entries = []
        self.steps = {}

    def add(
        self, time, points, blocks, point_data=None, cell_data=None, part=0, name=None
    ):
        step = self.steps.setdefault(float(time), len(self.steps))
        fname = "%s_%04d_%d.vtu" % (self.stem, step, part)
        write_vtu(
            os.path.join(self.dir, fname),
            points,
            blocks,
            point_data=point_data,
            cell_data=cell_data,
            time=time,
        )
        self.entries.append((float(time), part, name, self.stem + "/" + fname))

    def close(self):
        lines = [
            '<?xml version="1.0"?>\n',
            '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">\n',
            "<Collection>\n",
        ]
        for time, part, name, fname in self.entries:
            lines.append(
                '<DataSet timestep="%s" part="%d"%s file=%s/>\n'
                % (
                    _time_text(time),
                    part,
                    "" if name is None else " name=%s" % quoteattr(name),
                    quoteattr(fname),
                )
            )
        lines.append("</Collection>\n</VTKFile>\n")
        with open(self.path, "w") as f:
            f.write("".join(lines))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        return False
