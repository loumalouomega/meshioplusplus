from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._med import read as _py_read
from ._med import write as _py_write
from ._medmulti import read_med_multi, write_med_multi

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)

# The C++ core handles the mesh-representation part of MED exactly (points,
# point/cell tags, families with GRO group names, named regions derived from
# those families (one per group name -- see doc/regions.md), the optional
# NUM global-numbering datasets, an INFOS_GENERALES version check, mesh-level
# metadata, node orientation, and POG/POG2 ragged polygons). It deliberately
# DEFERS to the Python implementation (by raising, caught below) for anything
# the Python reference does that the C++ path does not replicate
# byte-for-byte: fields (CHA) with the MED-4.1 bitmask / units / step
# metadata, gmsh:physical family bridging, non-default profiles, and
# multi-mesh files. `read_med_multi` / `write_med_multi` stay on Python.


def read(filename):
    """Read a MED file (C++ core when built with HDF5, Python/h5py fallback)."""
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
            return _core.med_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


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
        except Exception:
            pass
    return _py_write(filename, mesh, med_version=med_version, **kwargs)


register_format("med", [".med"], read, {"med": write})

__all__ = ["read", "write", "read_med_multi", "write_med_multi"]
