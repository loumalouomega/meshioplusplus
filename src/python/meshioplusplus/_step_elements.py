"""Fitted primitives as STEP elements (roadmap 7.6.5).

Each plane, sphere, cylinder, cone or torus ``fit_primitives`` fitted, written
as exact CAD geometry a user can snap to: a bounded plane face, a ball, a
lidded cylinder, a frustum, a torus. C++-core only (see ``doc/formats/step.md``).
"""

from __future__ import annotations

import os

__all__ = ["write_step_elements"]


def write_step_elements(path, mesh, *, patches=None, segments: int = 64):
    """Write the fitted primitives of ``mesh`` as STEP geometry.

    :param path: the ``.step`` file to create.
    :param mesh: a mesh carrying ``segment:patch`` cell data and the
        ``fit:patch_*`` table (``segment``, then ``fit_primitives``).
    :param patches: the patch labels to write; ``None`` writes every fitted one.
    :param segments: facets round a cylinder, cone or torus in the tessellation
        the B-rep writer is given (the faces are the exact surfaces).
    :returns: the report: ``num_planes``, ``num_spheres``, ``num_cylinders``,
        ``num_cones``, ``num_tori``, ``num_skipped`` and ``warnings``.
    :raises meshioplusplus.WriteError: on a mesh with no fit table, an unknown
        patch, or nothing to write.
    :raises NotImplementedError: when the compiled C++ core is unavailable.
    """
    from ._hull import _core_op

    res = _core_op("step_write_elements")(
        os.fspath(path),
        mesh,
        [int(p) for p in (patches or [])],
        int(segments),
    )
    report = {
        k: int(res[k])
        for k in (
            "num_planes",
            "num_spheres",
            "num_cylinders",
            "num_cones",
            "num_tori",
            "num_skipped",
        )
    }
    report["warnings"] = [str(w) for w in res["warnings"]]
    return report
