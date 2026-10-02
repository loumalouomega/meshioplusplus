from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._med import read as _py_read
from ._med import write as _py_write
from ._medmulti import _resolve_mesh_names
from ._medmulti import read_med_multi as _py_read_multi
from ._medmulti import write_med_multi as _py_write_multi

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)

# The C++ core handles the mesh-representation part of MED exactly (points,
# point/cell tags, families with GRO group names, named regions derived from
# those families (one per group name -- see doc/regions.md), the optional
# NUM global-numbering datasets, an INFOS_GENERALES version check, mesh-level
# metadata, node orientation, and POG/POG2 ragged polygons). It deliberately
# DEFERS to the Python implementation (by raising, caught below) for anything
# the Python reference does that the C++ path does not replicate
# byte-for-byte: enhanced field units, step metadata and ELNO/ELGA data.
# Named meshes and ordinary nodal/element profiles are handled natively.


def read(filename, time_step: int = 0, *, mesh_name=None):
    """Read a MED file (C++ core when built with HDF5, Python/h5py fallback).

    ``time_step`` selects one step of a multi-step ``CHA`` field (0 = first,
    negative counts from the end), resolved the same way the C API/Fortran/
    Julia/R/WASM surfaces already do -- see :func:`meshioplusplus.med.read`'s
    C++ counterpart, ``read_med``. A non-default value forces the C++ path
    (it has no Python-fallback equivalent) and an out-of-range step raises
    :class:`~meshioplusplus.ReadError` naming the field and its step count.
    ``mesh_name`` selects a named mesh from a multi-mesh file.
    """
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            # The C++ path already returns a mesh whose `.regions` (and so
            # `.point_sets`/`.cell_sets`, the compat views over them) were
            # derived from the same family tables the Python reader's
            # `_families_to_point_sets`/`_families_to_cell_sets` build --
            # see `med_attach_point_regions`/`med_attach_cell_regions` in
            # med.cpp. Re-deriving them here would be redundant at best and,
            # since the property setters *replace* all regions of their kind,
            # would silently discard the dim/tag `mesh_to_py` already
            # attached to each one.
            if mesh_name is not None:
                return _core.med_read_named(str(filename), mesh_name, time_step)
            return _core.med_read(str(filename), time_step)
        except Exception as exc:
            if time_step:
                raise
            if not core_declined(exc, "med", "read", filename):
                raise
    if mesh_name is not None:
        meshes, names = _py_read_multi(filename)
        if mesh_name not in names:
            from .._exceptions import ReadError

            raise ReadError(f"MED: no mesh named '{mesh_name}'")
        return meshes[names.index(mesh_name)]
    return _py_read(filename)


def read_med_multi(filename, **kwargs):
    """Read all MED meshes, returning the historical (meshes, mesh_names) pair."""
    if (
        _HAS_HDF5
        and hasattr(_core, "med_read_named")
        and not kwargs
        and not is_buffer(filename, "r")
    ):
        try:
            names = _core.med_mesh_names(str(filename))
            return [_core.med_read_named(str(filename), name) for name in names], names
        except Exception as exc:
            if not core_declined(exc, "med", "read", filename):
                raise
    return _py_read_multi(filename, **kwargs)


def write_med_multi(filename, meshes, mesh_names=None, med_version="4.1.0", **kwargs):
    """Write several named meshes to one MED file, using native HDF5 when supported."""
    enhanced = any(
        "med:field_units" in mesh.field_data
        or "med:step_meta" in mesh.field_data
        or mesh.field_data.get("med:nom")
        or _names_encode_a_timestep(mesh)
        for mesh in meshes
    )
    if (
        _HAS_HDF5
        and hasattr(_core, "med_write_multi")
        and not kwargs
        and not enhanced
        and not is_buffer(filename, "w")
    ):
        try:
            _core.med_write_multi(
                str(filename),
                meshes,
                _resolve_mesh_names(meshes, mesh_names),
                str(med_version),
            )
            return
        except Exception as exc:
            if not core_declined(exc, "med", "write", filename):
                raise
    return _py_write_multi(filename, meshes, mesh_names, med_version, **kwargs)


def _names_encode_a_timestep(mesh):
    """Whether any array name uses the ``"Name[idx] - pdt"`` multi-timestep
    convention ``_med._parse_med_field_name`` recognizes.

    The C++ writer has no notion of this encoding at all -- it would write
    "Temperature[0] - 0.0" and "Temperature[1] - 1.0" as two unrelated
    fields instead of two timesteps of one field named "Temperature". A mesh
    using this convention must defer to the Python writer, which groups them
    correctly.
    """
    from ._med import _parse_med_field_name

    for name in list(mesh.point_data):
        if _parse_med_field_name(name)[1] is not None:
            return True
    for name in list(mesh.cell_data):
        if _parse_med_field_name(name)[1] is not None:
            return True
    return False


def write(filename, mesh, med_version="4.1.0", **kwargs):
    """Write a MED file (C++ core when built with HDF5, Python/h5py fallback).

    The C++ path handles a single-timestep field (CHA) write directly --
    ordinary ``point_data``/``cell_data`` arrays, no units, no component
    names, no multiple timesteps. Three signals mean the caller wants more
    than that, and must go to Python instead:

    - ``med:field_units``/``med:step_meta`` in ``field_data`` -- Python-only
      conventions (dicts, not arrays), so this check has to happen *here*
      rather than in the C++ core: they cannot survive the Python->C++ mesh
      conversion at all (it silently drops any non-numeric ``field_data``
      entry), so a check made there would never see them and could not defer
      correctly. ``med:nom`` is unaffected -- it is passed through explicitly
      below and always was.
    - an array name using the ``"Name[idx] - pdt"`` multi-timestep encoding
      (see ``_names_encode_a_timestep``).

    Named regions need no signal here: when the mesh carries no native
    ``point_tags``/``cell_tags`` of its own, the C++ writer synthesizes them
    from ``mesh.regions`` directly (``med_point_regions_to_tags``/
    ``med_cell_regions_to_tags`` in med.cpp), the same combo-per-name-set
    algorithm ``_ensure_med_families`` below uses for the Python path.
    """
    wants_enhanced_fields = (
        "med:field_units" in mesh.field_data
        or "med:step_meta" in mesh.field_data
        or _names_encode_a_timestep(mesh)
    )
    if (
        _HAS_HDF5
        and not kwargs
        and not wants_enhanced_fields
        and not is_buffer(filename, "w")
    ):
        point_tags = getattr(mesh, "point_tags", None) or {}
        cell_tags = getattr(mesh, "cell_tags", None) or {}
        med_nom = mesh.field_data.get("med:nom", [])
        point_tag_groups = getattr(mesh, "point_tag_groups", None) or {}
        cell_tag_groups = getattr(mesh, "cell_tag_groups", None) or {}
        try:
            _core.med_write(
                str(filename),
                mesh,
                dict(point_tags),
                dict(cell_tags),
                list(med_nom),
                getattr(mesh, "mesh_name", "mesh") or "mesh",
                getattr(mesh, "description", "") or "",
                getattr(mesh, "unit_time", "") or "",
                getattr(mesh, "unit_coords", "") or "",
                dict(point_tag_groups),
                dict(cell_tag_groups),
                str(med_version),
            )
            return
        except Exception as exc:
            if not core_declined(exc, "med", "write", filename):
                raise
    return _py_write(filename, mesh, med_version=med_version, **kwargs)


register_format("med", [".med"], read, {"med": write})

__all__ = ["read", "write", "read_med_multi", "write_med_multi"]
