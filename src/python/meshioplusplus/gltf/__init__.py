from typing import Union

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._gltf import write as _py_write
from ._gltf import write_buffer as _py_write_buffer


def write(
    filename,
    mesh,
    split_angle: float = 30.0,
    normal_weight: str = "angle",
    normals: bool = True,
    fields: bool = True,
    color_by: Union[str, None] = None,
    component: Union[int, None] = None,
    cmap: str = "viridis",
    vmin: Union[float, None] = None,
    vmax: Union[float, None] = None,
    nan_color: str = "#808080",
    unlit: bool = True,
    up_axis: str = "auto",
    recenter: bool = True,
    scale: float = 1.0,
    by_region: bool = True,
    container: str = "auto",
):
    """Write a glTF 2.0 file (C++ core for real file paths, Python fallback).

    ``.glb`` writes the binary container; ``.gltf`` writes JSON with a
    ``<name>.bin`` beside it (``container`` overrides the suffix). The surface
    of the mesh is exported: the skin of volume cells, 2-D cells, ``line``
    cells as ``LINES`` and ``vertex`` cells or a cell-less mesh as ``POINTS``,
    one named node per cell region.

    glTF normals are per vertex, so points are duplicated where the surface
    creases by more than ``split_angle`` degrees (``[0, 180]``); ``normal_weight``
    is ``"angle"`` or ``"area"``, and ``normals=False`` writes none.

    Every one-to-four component ``point_data`` array is exported raw as
    ``_NAME`` (``fields=False`` skips them). ``color_by`` names a ``point_data``
    or ``cell_data`` array to bake into ``COLOR_0`` through ``cmap``
    (``viridis``, ``coolwarm`` or ``turbo``) over ``vmin``..``vmax`` (default:
    the finite range of what is exported), the material then being
    ``KHR_materials_unlit`` unless ``unlit=False``; non-finite values take
    ``nan_color``. Multi-component arrays reduce to ``component`` or to their
    magnitude.

    The output is Y-up, right-handed, ``float32``. ``up_axis`` names the source
    axis that points up (``"auto"`` is ``y`` for a flat mesh and ``z``
    otherwise); the rotation, ``scale`` (source unit to metres) and the
    bounding-box centre that ``recenter`` subtracts all live on the root node,
    not in the coordinates. ``by_region=False`` writes one node.
    """
    kwargs = dict(
        up_axis=up_axis,
        normal_weight=normal_weight,
        normals=normals,
        fields=fields,
        recenter=recenter,
        by_region=by_region,
        unlit=unlit,
        split_angle=split_angle,
        scale=scale,
        color_by=color_by,
        component=component,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        nan_color=nan_color,
    )
    if is_buffer(filename, "w"):
        if container not in ("auto", "glb", "binary"):
            raise ValueError(
                "meshio++: gltf: only the binary (glb) container can be written to a buffer"
            )
        return _py_write_buffer(filename, mesh, **kwargs)
    try:
        # By keyword: the names match the py::arg list in bindings/python/_core.cpp.
        _core.gltf_write(
            str(filename),
            mesh,
            container=container,
            up_axis=up_axis,
            normal_weight=normal_weight,
            normals=normals,
            fields=fields,
            recenter=recenter,
            by_region=by_region,
            unlit=unlit,
            split_angle=float(split_angle),
            scale=float(scale),
            color_by=color_by or "",
            component=component,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            nan_color=nan_color,
        )
        return None
    except Exception as exc:
        if not core_declined(exc, "gltf", "write", filename):
            raise
    return _py_write(filename, mesh, container=container, **kwargs)


register_format("gltf", [".gltf", ".glb"], None, {"gltf": write})

__all__ = ["write"]
