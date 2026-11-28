"""Integration-point data operations: ip_to_points, ip_to_cells, points_to_ip."""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from ._mesh import Mesh


def _extract_ip_layout(mesh: Mesh, name: str) -> Optional[tuple[int, int, int]]:
    """Return ``(points, width, rule_id)`` for ``cell_data[name]`` if annotated."""
    fd = mesh.field_data or {}
    key = f"meshio:ip:{name}"
    if key in fd:
        arr = np.asarray(fd[key]).ravel()
        if arr.size >= 2:
            pts = int(arr[0])
            width = int(arr[1])
            rule_id = int(arr[2]) if arr.size >= 3 else -1
            if pts >= 1 and width >= 1:
                return (pts, width, rule_id)
    for prefix in ("abaqus:layout:", "lsdyna_d3plot:layout:"):
        lk = f"{prefix}{name}"
        if lk in fd:
            arr = np.asarray(fd[lk]).ravel()
            if arr.size >= 2:
                pts = int(arr[0])
                width = int(arr[1])
                if pts >= 1 and width >= 1:
                    return (pts, width, -1)
    gid_k = f"gid:gauss_points:{name}"
    if gid_k in fd and name in (mesh.cell_data or {}):
        arr = np.asarray(fd[gid_k]).ravel()
        if arr.size >= 1:
            g = int(arr[0])
            for blk in mesh.cell_data[name]:
                b = np.asarray(blk)
                if b.shape[0] > 0:
                    cols = b.shape[1] if b.ndim > 1 else 1
                    if g >= 1 and cols % g == 0:
                        return (g, cols // g, -1)
                    break
    gmsh_k = f"gmsh:element_node_components:{name}"
    if gmsh_k in fd and name in (mesh.cell_data or {}):
        arr = np.asarray(fd[gmsh_k]).ravel()
        if arr.size >= 1:
            w = int(arr[0])
            for blk in mesh.cell_data[name]:
                b = np.asarray(blk)
                if b.shape[0] > 0:
                    cols = b.shape[1] if b.ndim > 1 else 1
                    if w >= 1 and cols % w == 0:
                        return (cols // w, w, 0)
                    break
    return None


def ip_to_points(
    mesh: Mesh,
    *,
    arrays: Optional[Sequence[str]] = None,
    drop_source: bool = False,
    return_report: bool = False,
):
    """Extrapolate integration-point ``cell_data`` arrays to ``point_data``.

    Uses element shape-function extrapolation (least-squares when the
    integration-point count differs from the node count) averaged across
    incident cells. Element-nodal arrays (``rule_id == 0``) average directly
    at shared nodes.
    """
    from . import _core

    arr_list = list(arrays) if arrays is not None else []
    out, report = _core.ip_to_points(mesh, arr_list, bool(drop_source))
    return (out, report) if return_report else out


def ip_to_cells(
    mesh: Mesh,
    *,
    arrays: Optional[Sequence[str]] = None,
    drop_metadata: bool = True,
    return_report: bool = False,
):
    """Average integration-point ``cell_data`` arrays to one value per cell.

    Uses quadrature weights when the rule is known, or uniform weights when
    ``rule_id == -1``.
    """
    from . import _core

    arr_list = list(arrays) if arrays is not None else []
    out, report = _core.ip_to_cells(mesh, arr_list, bool(drop_metadata))
    return (out, report) if return_report else out


def points_to_ip(
    mesh: Mesh,
    *,
    arrays: Optional[Sequence[str]] = None,
    rule_id: int = 2,
    drop_source: bool = False,
    return_report: bool = False,
):
    """Evaluate ``point_data`` arrays at element integration points.

    Attaches the sampled ``cell_data`` arrays and records
    ``field_data["meshio:ip:<name>"] = [points, width, rule_id]``.
    ``rule_id == 0`` produces element-nodal data; ``rule_id >= 1`` selects
    the Gauss quadrature rule of that exactness order.
    """
    from . import _core

    arr_list = list(arrays) if arrays is not None else []
    out, report = _core.points_to_ip(mesh, arr_list, int(rule_id), bool(drop_source))
    return (out, report) if return_report else out
