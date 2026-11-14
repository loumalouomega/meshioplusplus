# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Write the meshes Catalyst 2 or Ascent hands an in situ pipeline as a meshio++ series.

Both pass a Conduit node holding a Mesh Blueprint (roadmap 10.2); this script reads
it with ``meshioplusplus.from_conduit`` and writes one step per call, so a run
lands as a file meshio++, ParaView and every registry writer read afterwards.

**Catalyst 2.** The Catalyst API calls ``catalyst_initialize``, ``catalyst_execute``
and ``catalyst_finalize`` with a node in the Catalyst protocol: the time in
``catalyst/state/time``, each channel's Blueprint in
``catalyst/channels/<name>/data`` (its ``type`` ``mesh`` or ``multimesh``). Use the
module functions below from a Catalyst implementation or a simulation's Python
adaptor; the output path comes from ``catalyst/meshioplusplus/output`` in the
initialize node, else ``$MIO_CATALYST_OUTPUT``, else ``catalyst.vtkhdf``, and the
channel from ``catalyst/meshioplusplus/channel`` (the first mesh channel by default).

**Ascent.** As an Ascent Python extract (``type: python``, ``params/file`` this
script), ``ascent_data()`` is the published Blueprint and ``ascent_extract()``
below writes it. Ascent runs the extract afresh at every step, so it writes a file
per step, numbered by the Blueprint's ``state/cycle``: the path is
``$MIO_CATALYST_OUTPUT`` when it holds ``{step}``, else ``ascent_{step}.vtu``.

A ``.vtkhdf`` path is one transient file (the grid must stay the same from step
to step); a path with ``{step}`` writes one file per step in any format meshio++
writes (``run_{step}.vtu``, ``run_{step}.e``). Read the result back with
``meshioplusplus.read_sequence``. Needs ``meshioplusplus`` (and h5py for
``.vtkhdf``) in the Python that runs the pipeline; Conduit comes with the host.
"""

import os

import numpy as np

import meshioplusplus

_ENV = "MIO_CATALYST_OUTPUT"
_DEFAULT = "catalyst.vtkhdf"
TIME_KEY = "meshio:time"


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


def _names(node):
    return list(node.keys()) if isinstance(node, dict) else list(node.child_names())


def _scalar(value):
    return np.asarray(value).reshape(-1)[0]


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


class BlueprintSeries:
    """One output series, a step per :meth:`write`.

    ``path`` is a ``.vtkhdf`` file (fixed grid, one file) or a path with
    ``{step}`` (a file per step, any format). ``file_format`` forces the writer.
    """

    def __init__(self, path, file_format=None):
        self.path = str(path)
        self.file_format = file_format
        self.num_steps = 0
        self._series = None
        self._grid = None

    def write(self, time, mesh, step=None):
        """Write ``mesh`` at ``time``; ``step`` numbers a ``{step}`` file (the count
        of steps written by default)."""
        mesh.field_data[TIME_KEY] = np.array([float(time)])
        if "{step}" in self.path:
            number = self.num_steps if step is None else int(step)
            out = self.path.replace("{step}", "%06d" % number)
            meshioplusplus.write(out, mesh, file_format=self.file_format)
        else:
            self._write_series(time, mesh)
        self.num_steps += 1

    def _write_series(self, time, mesh):
        from meshioplusplus.vtkhdf import TimeSeriesWriter

        grid = (
            np.asarray(mesh.points).shape,
            [(c.type, len(c.data)) for c in mesh.cells],
        )
        if self._series is None:
            self._series = TimeSeriesWriter(self.path)
            self._series.__enter__()
            self._series.write_points_cells(mesh.points, mesh.cells)
            self._grid = grid
        elif grid != self._grid:
            raise ValueError(
                "blueprint_series: the mesh changed at step %d; a .vtkhdf series holds "
                "one grid -- give a path with {step} to write a file per step"
                % self.num_steps
            )
        self._series.write_data(
            time, point_data=mesh.point_data, cell_data=mesh.cell_data
        )

    def close(self):
        if self._series is not None:
            self._series.__exit__(None, None, None)
            self._series = None


# --------------------------------------------------------------------------- #
# Catalyst 2                                                                  #
# --------------------------------------------------------------------------- #
_STATE = {"series": None, "channel": None}


def catalyst_initialize(params):
    """Open the series; ``params`` is the node ``catalyst_initialize`` received."""
    path = os.environ.get(_ENV, _DEFAULT)
    channel = None
    if params is not None and _has(params, "catalyst/meshioplusplus/output"):
        path = _text(_get(params, "catalyst/meshioplusplus/output"))
    if params is not None and _has(params, "catalyst/meshioplusplus/channel"):
        channel = _text(_get(params, "catalyst/meshioplusplus/channel"))
    _STATE["series"] = BlueprintSeries(path)
    _STATE["channel"] = channel


def mesh_from_catalyst(params, channel=None):
    """``(time, mesh)`` from a ``catalyst_execute`` node."""
    time = (
        float(_scalar(_get(params, "catalyst/state/time")))
        if _has(params, "catalyst/state/time")
        else 0.0
    )
    channels = _get(params, "catalyst/channels")
    names = _names(channels)
    if channel is None:
        names = [
            n
            for n in names
            if _text(_get(channels, n + "/type")) in ("mesh", "multimesh")
        ]
        if not names:
            raise ValueError("blueprint_series: the node has no mesh channel")
        channel = names[0]
    data = _get(channels, channel + "/data")
    return time, meshioplusplus.from_conduit(data)


def catalyst_execute(params):
    """Write the step ``params`` describes."""
    if _STATE["series"] is None:
        catalyst_initialize(None)
    time, mesh = mesh_from_catalyst(params, _STATE["channel"])
    step = (
        int(_scalar(_get(params, "catalyst/state/timestep")))
        if _has(params, "catalyst/state/timestep")
        else None
    )
    _STATE["series"].write(time, mesh, step)


def catalyst_finalize(params=None):
    """Close the series."""
    if _STATE["series"] is not None:
        _STATE["series"].close()
    _STATE["series"] = None


# --------------------------------------------------------------------------- #
# Ascent                                                                      #
# --------------------------------------------------------------------------- #
def _state(data, name):
    """``state/<name>`` of a Blueprint domain, or of the first domain of a multi-domain node."""
    if _has(data, "state/" + name):
        return _scalar(_get(data, "state/" + name))
    if not _has(data, "coordsets") and _names(data):
        first = _get(data, _names(data)[0])
        if _has(first, "state/" + name):
            return _scalar(_get(first, "state/" + name))
    return None


def ascent_extract(data=None):
    """Write the Blueprint an Ascent Python extract publishes (``ascent_data()``):
    one file per call, numbered by ``state/cycle``."""
    if data is None:
        data = ascent_data()  # noqa: F821 -- an Ascent extract's builtin
    path = os.environ.get(_ENV, "")
    if "{step}" not in path:
        path = "ascent_{step}.vtu"
    time, cycle = _state(data, "time"), _state(data, "cycle")
    BlueprintSeries(path).write(
        0.0 if time is None else float(time),
        meshioplusplus.from_conduit(data),
        None if cycle is None else int(cycle),
    )


if __name__ == "__main__" and "ascent_data" in globals():
    ascent_extract()
