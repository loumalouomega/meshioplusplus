from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._med import read as _py_read
from ._med import write as _py_write
from ._medmulti import read_med_multi, write_med_multi

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)

# The C++ core handles the mesh-representation part of MED exactly (points,
# point/cell tags, families with GRO group names, mesh-level metadata, node
# orientation, and POG/POG2 ragged polygons). It deliberately DEFERS to the
# Python implementation (by raising, caught below) for anything the Python
# reference does that the C++ path does not replicate byte-for-byte: fields
# (CHA) with the MED-4.1 bitmask / units / step metadata, gmsh:physical family
# bridging, non-default profiles, and multi-mesh files. `read_med_multi` /
# `write_med_multi` stay on Python.


def read(filename):
    """Read a MED file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            mesh = _core.med_read(str(filename))
            # The C++ path returns point/cell tags + families; reconstruct the
            # point_sets/cell_sets exactly as the Python reader does (reusing
            # its helpers) so named families round-trip identically.
            from ._med import _families_to_cell_sets, _families_to_point_sets

            mesh.point_sets = _families_to_point_sets(
                getattr(mesh, "point_tags", {}) or {},
                mesh.point_data.get("point_tags"),
            )
            mesh.cell_sets = _families_to_cell_sets(
                getattr(mesh, "cell_tags", {}) or {},
                mesh.cell_data.get("cell_tags"),
                len(mesh.cells),
            )
            return mesh
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, med_version="4.1.0", **kwargs):
    """Write a MED file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not kwargs and not is_buffer(filename, "w"):
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
