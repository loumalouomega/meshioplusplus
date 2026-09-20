"""ParaView collection ``.pvd``, the pure-Python reference.

A ``.pvd`` is ``<VTKFile type="Collection"><Collection><DataSet timestep= part=
group= file=/>...``: a time-indexed list of serial or parallel XML files, never
legacy ``.vtk``. It has two axes. ``timestep`` selects the *step* (``time_step=``
on read) and, within a step, ``part`` selects the *piece* (``piece=``). ``read``
gives step 0 with every part merged; ``read_sequence`` walks the time axis.

The index holds no geometry: writing means one ``.vtu`` per step in a sibling
directory named after the index's stem, with zero-padded names and relative
paths; reading resolves each ``file=`` against the index's own directory.
"""

from __future__ import annotations

import math
import os
import xml.etree.ElementTree as ET

import numpy as np

from .. import _provenance
from .. import _pvtk_index as _ix
from .._exceptions import ReadError, WriteError

TIME_KEY = "meshio:time"


def _parse(filename):
    """``(entries, times)``: the ``<DataSet>`` entries and the sorted distinct times."""
    filename = str(filename)
    try:
        tree = ET.parse(filename)
    except ET.ParseError as e:
        raise ReadError(f"meshio++: pvd: could not parse {filename}: {e}") from None
    root = tree.getroot()
    if root.tag != "VTKFile":
        raise ReadError(f"meshio++: pvd: expected tag 'VTKFile': {filename}")
    if root.get("type") != "Collection":
        raise ReadError(
            f"meshio++: pvd: expected type Collection, got {root.get('type')!r}: "
            f"{filename}"
        )
    coll = root.find("Collection")
    if coll is None:
        raise ReadError(f"meshio++: pvd: expected tag 'Collection': {filename}")

    entries = []
    for order, el in enumerate(coll.findall("DataSet")):
        try:
            # A missing timestep is 0, as ParaView reads it; a missing part is 0.
            time = float(el.get("timestep", "0"))
            part = int(el.get("part", "0"))
        except ValueError as e:
            raise ReadError(
                f"meshio++: pvd: bad timestep/part in {filename}: {e}"
            ) from None
        entries.append(
            {
                "time": time,
                "part": part,
                "order": order,
                "group": el.get("group", ""),
                "name": el.get("name", ""),
                "file": el.get("file"),
            }
        )
    return entries, sorted({e["time"] for e in entries})


def _region_name(entry):
    if entry["name"]:
        return entry["name"]
    base = f"part_{entry['part']}"
    return f"{entry['group']}/{base}" if entry["group"] else base


def _resolve_step(time_step, count):
    k = int(time_step)
    resolved = count + k if k < 0 else k
    if resolved < 0 or resolved >= count:
        raise ReadError(
            f"meshio++: time step {k} is out of range: this file has {count} "
            f"{'step' if count == 1 else 'steps'}"
        )
    return resolved


def read(filename, time_step=0, piece=None, ghosts="keep"):
    _ix.check_ghosts(ghosts)
    entries, times = _parse(filename)
    if not entries:
        return _ix.empty_mesh()

    step = _resolve_step(time_step, len(times))
    chosen = sorted(
        (e for e in entries if e["time"] == times[step]),
        key=lambda e: (e["part"], e["order"]),
    )

    def one(entry):
        path = _ix.resolve_path(filename, entry["file"], "pvd")
        mesh = _ix.read_child(path, ghosts, _ix.COLLECTION_EXTENSIONS)
        return _ix.drop_ghosts(mesh) if ghosts == "drop" else mesh

    if piece is not None:
        out = one(chosen[_ix.resolve_piece(piece, len(chosen))])
    else:
        out = _ix.merge_pieces(
            [one(e) for e in chosen], [_region_name(e) for e in chosen]
        )
    out.field_data[TIME_KEY] = np.array([times[step]])
    # Every step's time, so `read_metadata` can report them off a full read.
    out.time_values = list(times)
    return out


def _index_text(entries):
    lines = ['<?xml version="1.0"?>']
    lines.append('<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">')
    lines.append(_provenance.render_xml_comment(_provenance.SlotTier.BLOCK))
    lines.append("<Collection>")
    for time, rel in entries:
        lines.append(
            f'<DataSet timestep="{float(time)!r}" part="0" file={_ix.quote(rel)}/>'
        )
    lines.append("</Collection>")
    lines.append("</VTKFile>")
    return "\n".join(lines) + "\n"


class SeriesWriter:
    """Stream ``(time, mesh)`` steps into one ``.pvd`` plus a ``.vtu`` per step.

    One mesh is alive at a time, and the index is rewritten after every step, so
    a run that is killed leaves a collection ParaView opens covering every
    finished step.
    """

    def __init__(self, path, binary=True, compression="zlib", header_type=None):
        self.path = str(path)
        self._binary = binary
        self._compression = compression
        self._header_type = header_type
        stem = os.path.splitext(os.path.basename(self.path))[0]
        parent = os.path.dirname(self.path)
        self._dir = os.path.join(parent, stem) if parent else stem
        self._stem = stem
        self._entries = []
        os.makedirs(self._dir, exist_ok=True)

    def write(self, time, mesh):
        from .. import vtu

        if not math.isfinite(float(time)):
            raise WriteError("meshio++: pvd: a step's time must be finite")
        name = f"{self._stem}_{len(self._entries):04d}.vtu"
        vtu.write(
            os.path.join(self._dir, name),
            mesh,
            binary=self._binary,
            compression=self._compression,
            header_type=self._header_type,
        )
        self._entries.append((float(time), f"{self._stem}/{name}"))
        self._flush()

    def _flush(self):
        with open(self.path, "w") as f:
            f.write(_index_text(self._entries))

    def close(self):
        if not self._entries:
            raise WriteError("meshio++: pvd: a collection needs at least one step")
        self._flush()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.close()
        return False


def _mesh_time(mesh):
    value = mesh.field_data.get(TIME_KEY)
    if value is not None and np.asarray(value).size == 1:
        return float(np.asarray(value).ravel()[0])
    return 0.0


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """A one-step collection; the step's time is ``field_data['meshio:time']`` if set."""
    with SeriesWriter(filename, binary, compression, header_type) as writer:
        writer.write(_mesh_time(mesh), mesh)
