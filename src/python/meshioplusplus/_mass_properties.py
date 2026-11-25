"""Mass, centre of mass and inertia of a body (roadmap 9.4).

``mass_properties(mesh, density=...)`` integrates a density over a volume mesh,
a closed surface or a shell and returns the mass, the centre of mass, the
inertia tensor, its principal moments and axes and the radii of gyration, for
the whole mesh and for every named cell region. See ``doc/mass_properties.md``.

**C++-core only, with no numpy fallback**: the integrals reuse the core's
polyhedron decomposition, whose exact arithmetic a second implementation would
have to repeat. This module raises when ``_core`` cannot be used.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

__all__ = ["mass_properties"]


def mass_properties(
    mesh,
    *,
    density: float = 1.0,
    density_array: Optional[str] = None,
    region: Optional[str] = None,
    thickness: float = 0.0,
    about: str = "centroid",
) -> Dict[str, Any]:
    """Mass properties of a mesh.

    :param mesh: a mesh with volume cells, or a closed triangle / polygon
        surface, or any surface when ``thickness`` is given.
    :param density: the density of every cell.
    :param density_array: name of a scalar ``cell_data`` array of per-cell
        densities, replacing ``density`` (volume cells and shells only).
    :param region: restrict the whole-mesh entry to this cell region.
    :param thickness: above zero, surface cells are a shell of this thickness
        (open surfaces allowed); ignored for volume cells.
    :param about: ``"centroid"`` (default) or ``"origin"``, the point the
        inertia tensor is taken about.
    :returns: ``{"kind": "solid" | "shell", "total": entry, "regions":
        [entry, ...], "warnings": [...]}`` where each entry holds ``name``,
        ``num_cells``, ``measure`` (volume, or area for a shell), ``mass``,
        ``centre``, ``inertia`` (``xx yy zz xy yz zx``), ``principal_moments``
        (ascending), ``principal_axes`` (a row-major 3x3, one axis per
        column) and ``gyration_radii``.
    :raises ValueError: on an unknown region or array, a non-finite density,
        an open or inconsistently wound surface without ``thickness``, or no
        cell that can be integrated.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _as_mesh, _core_op

    return _core_op("mass_properties")(
        _as_mesh(mesh),
        float(density),
        "" if density_array is None else str(density_array),
        "" if region is None else str(region),
        float(thickness),
        str(about),
    )
