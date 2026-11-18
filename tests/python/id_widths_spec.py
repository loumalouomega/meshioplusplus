"""What each writable format does with an integer above 2**31 (roadmap 3.19).

One canonical mesh per format (the cell types its row in ``conformance_spec.SPEC``
declares) carries an int64 point array and an int64 cell array whose values are
``2**31 + 5`` and up. :func:`observe` writes it and reads it back and names the
outcome per array:

``exact``
    the same integers came back, as the dtype named after it (``exact int64``);
``double`` / ``rounded``
    the format stores the values as floating point; ``double`` came back equal,
    ``rounded`` came back to within 1e-4 relative but not equal: the format keeps
    single precision, or a few significant digits of text, so an integer above
    2**24 is no longer exact (the format's own limit, not a narrowing);
``dropped``
    the format does not carry that array;
``WRAPPED``
    the integers came back as *other* integers: a silent narrowing, always a
    defect.

A writer that refuses is reported as ``error: <exception>`` for the whole mesh;
a refusal must name the format (the value of :func:`names_the_format`), and an
optional library that is not installed is ``not observed``.
``tools/gen_id_widths_table.py`` renders the same table into ``doc/id_widths.md``,
and ``test_id_widths.py`` checks the declarations in :data:`EXPECTED` against it.
"""

from __future__ import annotations

import pathlib
import re
import signal
import tempfile
import warnings

try:  # imported as part of the test package
    from . import conformance_spec as cs
except ImportError:  # imported by tools/gen_id_widths_table.py
    import conformance_spec as cs
import numpy as np

BIG = 2**31 + 5
_TIMEOUT_S = 30

# Formats whose writer needs an optional library (or build option): not pinned when
# the library is present, because CI legs differ in what they have installed.
LIBRARY_FORMATS = frozenset({"gid", "szplt", "step", "usd", "vtx", "zarr"})

EXPECTED_IDS = {
    "abaqus": {"cell": "dropped", "point": "dropped"},
    "abaqus_fil": {"cell": "exact int64", "point": "exact int64"},
    "exodus": {"cell": "exact int64", "point": "exact int64"},
    "libmesh": {"error": "WriteError"},
    "mdpa": {"cell": "dropped", "point": "exact int64"},
    "med": {"cell": "exact int64", "point": "exact int64"},
    "radioss_anim": {"error": "WriteError"},
}

# Declared outcomes, from ``python tools/gen_id_widths_table.py --observe``.
EXPECTED = {
    "3mf": {"cell": "dropped", "point": "dropped"},
    "abaqus": {"cell": "dropped", "point": "dropped"},
    "abaqus_fil": {"cell": "dropped", "point": "dropped"},
    "adcirc": {"cell": "dropped", "point": "dropped"},
    "ansys": {"cell": "dropped", "point": "dropped"},
    "ansysInp": {"cell": "dropped", "point": "dropped"},
    "ansys_rst": {"cell": "dropped", "point": "dropped"},
    "avsucd": {"cell": "dropped", "point": "double"},
    "cad_json": {"cell": "dropped", "point": "dropped"},
    "cae": {"cell": "rounded", "point": "dropped"},
    "cgns": {"cell": "dropped", "point": "double"},
    "code_aster": {"cell": "dropped", "point": "dropped"},
    "dex": {"cell": "dropped", "point": "double"},
    "dolfin-xml": {"cell": "exact int64", "point": "exact int64"},
    "dxf": {"cell": "dropped", "point": "dropped"},
    "elmer": {"cell": "dropped", "point": "dropped"},
    "ensight": {"cell": "rounded", "point": "rounded"},
    "exodus": {"cell": "exact int64", "point": "exact int64"},
    "febio": {"cell": "double", "point": "double"},
    "femap": {"cell": "double", "point": "double"},
    "flac3d": {"cell": "dropped", "point": "dropped"},
    "flux": {"cell": "dropped", "point": "dropped"},
    "frd": {"cell": "dropped", "point": "rounded"},
    "freefem": {"cell": "dropped", "point": "dropped"},
    "freesurfer": {"cell": "dropped", "point": "dropped"},
    "gid": {"cell": "rounded", "point": "rounded"},
    "gifti": {"error": "WriteError"},
    "gltf": {"cell": "dropped", "point": "rounded"},
    "gmsh": {"cell": "double", "point": "double"},
    "gmsh22": {"cell": "double", "point": "double"},
    "h5m": {"cell": "dropped", "point": "exact int64"},
    "hmf": {"cell": "exact int64", "point": "exact int64"},
    "ip": {"cell": "dropped", "point": "double"},
    "las": {"cell": "dropped", "point": "dropped"},
    "libmesh": {"cell": "dropped", "point": "dropped"},
    "lsdyna": {"cell": "dropped", "point": "dropped"},
    "lsdyna_binout": {"error": "WriteError"},
    "lsdyna_d3plot": {"cell": "dropped", "point": "dropped"},
    "marc": {"cell": "dropped", "point": "dropped"},
    "marc_t19": {"cell": "rounded", "point": "rounded"},
    "mdpa": {"cell": "dropped", "point": "double"},
    "med": {"cell": "exact int64", "point": "exact int64"},
    "medit": {"cell": "dropped", "point": "dropped"},
    "medit_sol": {"error": "WriteError"},
    "mfem": {"cell": "dropped", "point": "dropped"},
    "mff": {"cell": "dropped", "point": "dropped"},
    "mfm": {"error": "WriteError"},
    "mphbin": {"cell": "dropped", "point": "dropped"},
    "mphtxt": {"cell": "dropped", "point": "dropped"},
    "nastran": {"cell": "dropped", "point": "dropped"},
    "nastran_h5": {"cell": "dropped", "point": "double"},
    "nastran_op2": {"cell": "dropped", "point": "dropped"},
    "netgen": {"error": "WriteError"},
    "neuroglancer": {"cell": "dropped", "point": "dropped"},
    "obj": {"cell": "dropped", "point": "dropped"},
    "off": {"cell": "dropped", "point": "dropped"},
    "openfoam": {"cell": "double", "point": "double"},
    "patran": {"cell": "dropped", "point": "dropped"},
    "pcd": {"cell": "dropped", "point": "exact int64"},
    "permas": {"cell": "dropped", "point": "dropped"},
    "ply": {"cell": "dropped", "point": "exact int64"},
    "pmsh": {"cell": "exact int64", "point": "exact int64"},
    "pvd": {"cell": "exact int64", "point": "exact int64"},
    "pvtp": {"cell": "exact int64", "point": "exact int64"},
    "pvtu": {"cell": "exact int64", "point": "exact int64"},
    "radioss": {"cell": "dropped", "point": "dropped"},
    "radioss_anim": {"cell": "rounded", "point": "rounded"},
    "radioss_th": {"cell": "dropped", "point": "dropped"},
    "selafin": {"error": "WriteError"},
    "step": {"error": "not observed"},
    "stl": {"cell": "dropped", "point": "dropped"},
    "su2": {"cell": "dropped", "point": "dropped"},
    "svg": {"cell": "dropped", "point": "dropped"},
    "szplt": {"error": "not observed"},
    "tecplot": {"cell": "double", "point": "double"},
    "tetgen": {"cell": "dropped", "point": "dropped"},
    "tikz": {"cell": "dropped", "point": "dropped"},
    "triangle": {"cell": "dropped", "point": "dropped"},
    "ugrid": {"cell": "dropped", "point": "dropped"},
    "ugrid_nc": {"cell": "exact int64", "point": "exact int64"},
    "unv": {"cell": "double", "point": "double"},
    "usd": {"error": "not observed"},
    "vtk": {"cell": "exact int64", "point": "exact int64"},
    "vtk42": {"cell": "exact int64", "point": "exact int64"},
    "vtk51": {"cell": "exact int64", "point": "exact int64"},
    "vtkhdf": {"cell": "exact int64", "point": "exact int64"},
    "vtm": {"cell": "exact int64", "point": "exact int64"},
    "vtp": {"cell": "exact int64", "point": "exact int64"},
    "vtu": {"cell": "exact int64", "point": "exact int64"},
    "vtx": {"error": "not observed"},
    "wkt": {"cell": "dropped", "point": "dropped"},
    "xdmf": {"cell": "exact int64", "point": "exact int64"},
    "xplt": {"cell": "rounded", "point": "rounded"},
    "xyz": {"cell": "dropped", "point": "double"},
    "z88": {"cell": "dropped", "point": "dropped"},
    "zarr": {"error": "not observed"},
}


class _Timeout(Exception):
    pass


def _alarm(*_):
    raise _Timeout()


def _needs_a_build(message):
    """True for a refusal that only says this build lacks an optional library."""
    return "needs a build with -D" in message or "-DMESHIOPLUSPLUS_WITH_" in message


def names_the_format(fmt, message):
    """True when an error message names the format.

    Compared without punctuation and case, so ``lsdyna_binout`` matches "LS-DYNA
    binout", and a variant (``medit_sol``) matches its family name ("Medit").
    """

    def norm(text):
        return re.sub(r"[^a-z0-9]", "", text.lower())

    text = norm(message)
    return norm(fmt) in text or norm(fmt.split("_")[0]) in text


def _array_outcome(back, want):
    """Classify one array that came back; compares the set of distinct values.

    Writers may renumber points, split blocks or pad a fixed-width array with
    zeros, so values rather than positions are compared: every nonzero value
    that came back must be one that went in (exactly, or to single precision).
    """
    if back is None:
        return "dropped"
    got = np.asarray(back).reshape(-1)
    if got.size == 0:
        return "dropped"
    sent = np.unique(want.astype(np.float64))
    came = np.unique(got.astype(np.float64))
    # Zero padding and NaN fill (a block the format has no value for) are not ids.
    came = came[(came != 0) & ~np.isnan(came)]
    if came.size == 0:
        return "dropped"
    if got.dtype.kind in "iu":
        ok = np.isin(came, sent).all()
        return f"exact {got.dtype}" if ok else "WRAPPED"
    if got.dtype.kind == "f":
        if np.isin(came, sent).all():
            return "double" if got.dtype.itemsize >= 8 else "single (exact)"
        near = (
            np.abs(came[:, None] - sent[None, :]) <= 1e-4 * np.abs(sent[None, :])
        ).any(1)
        return "rounded" if near.all() else "WRAPPED"
    return "WRAPPED"


# The formats that expose the ids they store for nodes and elements through a
# data key of the mesh model: ``{format: (point key, cell key)}``, ``None`` for a
# kind the format has no id for. Other formats renumber, or take an id only
# through their own numbering, so their ids are documented, not probed.
ID_KEYS = {
    "abaqus": ("abaqus:id", "abaqus:id"),
    "abaqus_fil": ("abaqus:id", "abaqus:id"),
    "exodus": ("global:node_id", "global:cell_id"),
    "libmesh": ("libmesh:id", None),
    "mdpa": ("mdpa:id", "mdpa:id"),
    "med": ("med:num", "med:num"),
    "radioss_anim": ("radioss:node_id", "radioss:element_id"),
}


def observe(fmt, keys=None):
    """``{"point": outcome, "cell": outcome}``, ``{"error": ...}`` or ``None``.

    With ``keys`` (from :data:`ID_KEYS`) the arrays are the format's id arrays,
    one distinct id per node and per cell, instead of the generic data arrays.
    """
    spec = cs.SPEC.get(fmt)
    if not spec or "cells" not in spec:
        return None
    pkey, ckey = keys if keys else ("p_big", "c_big")
    ctypes = [c for c in cs.CELL_TYPES if c in spec["cells"]]
    mesh = cs.canonical_mesh(
        ctypes, data=False, regions=False, planar=spec.get("input") == "2d"
    )
    n = len(mesh.points)
    point = np.arange(n, dtype=np.int64) + BIG
    ncells = sum(len(b.data) for b in mesh.cells)
    want_cells = (
        np.arange(ncells, dtype=np.int64) + BIG
        if keys
        else np.full(ncells, BIG, dtype=np.int64)
    )
    if pkey:
        mesh.point_data[pkey] = point.copy()
    if ckey:
        offsets = np.cumsum([0] + [len(b.data) for b in mesh.cells])
        mesh.cell_data[ckey] = [
            want_cells[offsets[i] : offsets[i + 1]].copy()
            for i in range(len(mesh.cells))
        ]
    # SIGALRM does not exist on Windows: there the probe runs without a timeout.
    use_alarm = hasattr(signal, "SIGALRM")
    if use_alarm:
        old = signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(_TIMEOUT_S)
    try:
        with tempfile.TemporaryDirectory() as tmp, warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                back = cs.round_trip(fmt, mesh, pathlib.Path(tmp))
            except _Timeout:
                return {"error": "TIMEOUT"}
            except (ImportError, NotImplementedError):
                # an optional library (or the OpenCASCADE backend) is missing
                return {"error": "not observed"}
            except Exception as exc:  # noqa: BLE001 - every refusal is classified
                if _needs_a_build(str(exc)):
                    return {"error": "not observed"}
                named = names_the_format(fmt, str(exc))
                return {"error": f"{type(exc).__name__}", "named": named}
    finally:
        if use_alarm:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
    cells = None
    if ckey and ckey in back.cell_data:
        cells = np.concatenate(
            [np.asarray(a).reshape(-1) for a in back.cell_data[ckey]]
        )
    pts = back.point_data.get(pkey) if pkey else None
    if pts is not None:
        pts = np.asarray(pts).reshape(-1)
    return {
        "point": _array_outcome(pts, point) if pkey else "no id",
        "cell": _array_outcome(cells, want_cells) if ckey else "no id",
    }


def observe_ids(fmt):
    """:func:`observe` on the format's id arrays (``None`` when it has none)."""
    keys = ID_KEYS.get(fmt)
    return observe(fmt, keys) if keys else None
