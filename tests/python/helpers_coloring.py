"""Fixtures and the shared option matrix for the SVG/TikZ colouring tests.

Both writers implement the same feature over the same resolution layer
(``detail/face_color.hpp`` and its twin ``_facecolor.py``), so they are worth
exercising against one matrix rather than two drifting ones.
"""

import copy

import numpy as np

import meshioplusplus

from . import helpers


def volume_mesh():
    """A tetra mesh carrying point, cell, vector and NaN-bearing arrays."""
    mesh = copy.deepcopy(helpers.tet_mesh)
    n_cells = [len(b.data) for b in mesh.cells]
    n_points = len(mesh.points)
    mesh.cell_data["tag"] = [
        np.arange(n, dtype=np.float64) for n in n_cells
    ]  # distinct per cell
    mesh.cell_data["velocity"] = [
        np.arange(n * 3, dtype=np.float64).reshape(n, 3) for n in n_cells
    ]
    mesh.cell_data["holey"] = [
        np.where(np.arange(n) % 2 == 0, np.arange(n, dtype=np.float64), np.nan)
        for n in n_cells
    ]
    mesh.point_data["T"] = np.arange(n_points, dtype=np.float64)
    mesh.point_data["flat"] = np.full(n_points, 2.5)  # degenerate range
    return mesh


def surface_mesh_3d():
    """A non-flat 3D surface mesh -- the projected path with no skin step."""
    return meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.5], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2]])],
        cell_data={"tag": [np.array([0.0, 1.0])]},
        point_data={"T": np.array([0.0, 1.0, 2.0, 3.0])},
    )


def flat_mesh():
    """A flat 2D mesh -- the legacy path, which must stay untouched when unset."""
    return meshioplusplus.Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        [("triangle", [[0, 1, 2], [1, 3, 2]])],
        cell_data={"tag": [np.array([0.5, np.nan])]},
        point_data={"T": np.array([0.0, 1.0, 2.0, 3.0])},
    )


# (id, mesh factory, writer kwargs) -- the matrix both writers are pinned on.
COLOR_CASES = [
    ("volume-cell", volume_mesh, {"color_by": "tag"}),
    ("volume-point", volume_mesh, {"color_by": "T"}),
    ("volume-vector-magnitude", volume_mesh, {"color_by": "velocity"}),
    ("volume-vector-component", volume_mesh, {"color_by": "velocity", "component": 1}),
    ("volume-explicit-range", volume_mesh, {"color_by": "T", "vmin": 1.0, "vmax": 3.0}),
    ("volume-nan", volume_mesh, {"color_by": "holey"}),
    ("volume-degenerate-range", volume_mesh, {"color_by": "flat"}),
    ("volume-coolwarm", volume_mesh, {"color_by": "tag", "cmap": "coolwarm"}),
    ("volume-turbo", volume_mesh, {"color_by": "T", "cmap": "turbo"}),
    ("volume-colorbar", volume_mesh, {"color_by": "tag", "colorbar": True}),
    ("volume-point-colorbar", volume_mesh, {"color_by": "T", "colorbar": True}),
    ("surface3d-cell", surface_mesh_3d, {"color_by": "tag"}),
    ("surface3d-point", surface_mesh_3d, {"color_by": "T", "colorbar": True}),
    ("flat-cell", flat_mesh, {"color_by": "tag"}),
    ("flat-point", flat_mesh, {"color_by": "T"}),
    ("flat-colorbar", flat_mesh, {"color_by": "tag", "colorbar": True}),
]

CASE_IDS = [c[0] for c in COLOR_CASES]
CASE_ARGS = [(c[1], c[2]) for c in COLOR_CASES]
