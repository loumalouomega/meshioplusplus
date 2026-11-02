"""The lowest eigenpairs of a triangle surface's Laplace-Beltrami operator
(roadmap 7.3.2): spectral positional encodings and shape descriptors.

``laplacian_eigenmodes(surface, num_modes=10)`` returns the surface with
``point_data["laplacian:modes"]`` (``(n, num_modes)``, M-orthonormal) and
``field_data["laplacian:eigenvalues"]``, from the cotangent stiffness against
the lumped mass (the operators ``geodesic_distance`` uses). See
``doc/eigenmodes.md``.

**C++-core only, with no numpy fallback** (the ``geodesic_distance``
precedent): the block iteration stops on a rounded test, and within a
repeated eigenvalue the basis it returns depends on every rounding on the
way, so a twin would return another, equally valid, basis.
"""

from __future__ import annotations

__all__ = ["laplacian_eigenmodes", "MODES_NAME", "EIGENVALUES_NAME"]

#: Point data: the eigenvectors. Twins ``kLaplacianModesName``.
MODES_NAME = "laplacian:modes"
#: Field data: the eigenvalues. Twins ``kLaplacianEigenvaluesName``.
EIGENVALUES_NAME = "laplacian:eigenvalues"


def laplacian_eigenmodes(
    mesh,
    num_modes: int = 10,
    *,
    tolerance: float = 1e-8,
    max_iterations: int = 0,
    return_report: bool = False,
):
    """The ``num_modes`` lowest eigenpairs of the surface Laplacian.

    Solves ``K phi = lambda M phi`` over the surface's linear triangles (the
    cotangent stiffness ``K``, the lumped mass ``M``; Neumann on a
    boundary). Returns the mesh with ``point_data["laplacian:modes"]``
    (Float64 ``(n, num_modes)``, ``M``-orthonormal, each signed so its
    largest entry is positive) and ``field_data["laplacian:eigenvalues"]``
    (ascending; a zero per connected piece). On the unit sphere the
    eigenvalues approach ``l (l + 1)``, ``2 l + 1`` times each.

    :param mesh: a surface: its ``triangle`` blocks (never modified). A volume
        cell, a quad, a polygon or a higher-order cell is refused naming the
        operation that makes linear triangles.
    :param num_modes: eigenpairs wanted, the smallest.
    :param tolerance: the relative residual every wanted pair must reach.
    :param max_iterations: block iterations at most; 0 takes 300.
    :param return_report: also return ``{"eigenvalues", "iterations",
        "max_residual", "converged", "num_degenerate"}``.
    :returns: the mesh, or ``(mesh, report)``.
    :raises ValueError: on an out-of-scope cell, no triangle, or
        ``num_modes`` below 1 or above the surface's points.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("laplacian_eigenmodes")(
        mesh, int(num_modes), float(tolerance), int(max_iterations)
    )
    out = res.pop("mesh")
    if not res["converged"]:
        from ._common import warn

        warn(
            f"laplacian_eigenmodes: no convergence in {res['iterations']} iterations "
            f"(residual {res['max_residual']:.3g}); raise max_iterations"
        )
    return (out, res) if return_report else out
