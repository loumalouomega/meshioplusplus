"""Turn the fitted patches of a closed skin into a CAD solid.

The last step of the skin-to-CAD chain :func:`segment` -> :func:`fit_primitives`
-> :func:`fit_cad`, through the optional OpenCASCADE backend of the core
(``operations/fit_cad.hpp``). Trimming, sewing and checking analytic surfaces
is a CAD kernel's work, so there is no numpy twin: without a core built with
``-DMESHIOPLUSPLUS_WITH_OCCT=ON`` the function raises ``NotImplementedError``
naming that flag.
"""

from __future__ import annotations

import math
import os

from ._fit_primitives import TYPES

_PREFIX = "meshio++: fit_cad: "
_FLAG = "-DMESHIOPLUSPLUS_WITH_OCCT=ON"

_REPORT_FLOATS = (
    "tolerance",
    "volume",
    "skin_volume",
    "volume_ratio",
    "max_deviation",
)
_REPORT_INTS = (
    "num_faces",
    "num_plane",
    "num_sphere",
    "num_cylinder",
    "num_cone",
    "num_torus",
    "num_faceted_faces",
    "num_patches",
    "num_demoted",
    "num_exact_edges",
    "num_approx_edges",
    "num_solids",
    "num_bspline",
    "num_split_patches",
    "num_fillet",
    "num_helical",
    "num_variable_fillet",
    "num_blend",
    "num_rounded_thread",
)


def has_occt() -> bool:
    """Whether the compiled core carries the OpenCASCADE backend."""
    try:
        from . import _core
    except ImportError:
        return False
    return bool(getattr(_core, "__has_occt__", False))


def _types_mask(types) -> int:
    if types is None:
        return sum(1 << t for t in range(1, 6))
    if isinstance(types, str):
        types = [types]
    mask = 0
    for t in types:
        if t not in TYPES[1:]:
            raise ValueError(
                f"{_PREFIX}unknown primitive type {t!r} "
                f"(expected one of {', '.join(TYPES[1:])})"
            )
        mask |= 1 << TYPES.index(t)
    return mask


def fit_cad(
    mesh,
    path=None,
    *,
    tolerance: float | None = None,
    exact_edges: bool = True,
    types=None,
    mesh_deflection: float | None = None,
    max_demotion_rounds: int = 8,
    freeform: bool = True,
    freeform_degree: int = 3,
    freeform_max_control: int = 64,
    freeform_split: bool = True,
    fillets: bool = True,
    threads: bool = True,
    return_report: bool = False,
):
    """Build a CAD solid from the fitted patches of a closed skin.

    Every patch :func:`fit_primitives` accepted becomes one plane, sphere,
    cylinder, cone or torus face trimmed by the patch boundary. A patch no
    primitive explains is tried, with ``fillets``, as a rolling-ball fillet
    when it runs between two analytic faces (of one radius, else of a radius
    varying along it), and, with ``threads``, as a helical thread face: about
    the axis of a cylinder or cone beside it, or of the screw the patch itself
    fits, with a straight flank or any other profile (a rounded root or
    crest). Failing those it becomes, with ``freeform``, a B-spline face
    fitted over the patch flattened onto its best-fit plane or unrolled around
    a cylinder, and one both fold is, with ``freeform_split``, split into
    pieces that flatten, each its own face under a new patch label. A patch
    between two fitted faces, one of them freeform, is then tried, with
    ``fillets``, as a rolling-ball blend. The rest of the skin becomes planar
    facets. The faces are sewn into a solid that
    OpenCASCADE checks, and a face that fails the check or strays further than
    ``tolerance`` from the skin is demoted to facets, so a feature the fitter
    cannot explain comes back faceted, never silently wrong.

    :param mesh: a closed, consistently wound surface mesh, normally the
        output of :func:`fit_primitives`. Without ``segment:patch`` or the
        ``fit:patch_*`` table, :func:`segment` and :func:`fit_primitives` run
        first with their default options.
    :param path: write the solid here as AP214 STEP (millimetres, a pinned
        header timestamp); ``None`` writes nothing.
    :param tolerance: the largest distance a face may lie from the skin; by
        default three times the largest ``fit:patch_tolerance`` of an accepted
        patch (``1e-3`` of the bounding-box diagonal when none was accepted).
    :param exact_edges: try the exact intersection of two analytic surfaces
        before approximating the edge between them with a B-spline.
    :param types: which accepted types become analytic faces, from ``plane``,
        ``sphere``, ``cylinder``, ``cone``, ``torus`` (default all); an empty
        sequence facets everything.
    :param mesh_deflection: the chordal deflection of the re-meshed solid; by
        default a quarter of the tolerance.
    :param max_demotion_rounds: how many times failed faces are demoted and
        the solid rebuilt before every patch is faceted.
    :param freeform: try a patch no primitive explains as a B-spline face.
    :param freeform_degree: that B-spline's highest degree, 3 to 8.
    :param freeform_max_control: the most grid points per direction its fit
        samples, 4 to 1024.
    :param freeform_split: split a freeform patch that folds over both its
        flattenings into pieces that do not.
    :param fillets: try a patch between two fitted faces as a rolling-ball
        fillet or blend.
    :param threads: try a patch as a helical thread face.
    :param return_report: also return the report.
    :returns: the solid re-meshed as triangles with cell data ``fit_cad:face``,
        ``fit_cad:type`` (0 a facet, 6 a freeform B-spline, 7 a fillet or
        blend, 8 a helical thread face, else as ``fit:type``) and
        ``fit_cad:patch``, and the solid as NURBS faces,
        loops and edges in ``mesh.info["brep"]`` (face ``i`` of
        ``fit_cad:face`` has id ``i + 1``), so writing it as ``.cad.json`` gives
        Kratos CAD JSON. With ``return_report``, ``(mesh, report)``, the report
        holding ``valid``, the tolerance, the solid's and the skin's volume and
        their ratio, the Hausdorff ``max_deviation``, the faces by type
        (``num_bspline`` the freeform ones, ``num_fillet``, ``num_helical``),
        of those the fillets of varying radius (``num_variable_fillet``), the
        blends beside a freeform face (``num_blend``) and the thread faces of
        a profile that is not straight (``num_rounded_thread``),
        ``num_split_patches``, the demoted patches and the exact and
        approximated edges.
    :raises NotImplementedError: when the core was built without OpenCASCADE.
    """
    opts = {
        "path": "" if path is None else os.fspath(path),
        "tolerance": 0.0 if tolerance is None else float(tolerance),
        "exact_edges": bool(exact_edges),
        "types": _types_mask(types),
        "mesh_deflection": 0.0 if mesh_deflection is None else float(mesh_deflection),
        "max_demotion_rounds": int(max_demotion_rounds),
        "freeform": bool(freeform),
        "freeform_degree": int(freeform_degree),
        "freeform_max_control": int(freeform_max_control),
        "freeform_split": bool(freeform_split),
        "fillets": bool(fillets),
        "threads": bool(threads),
    }
    if tolerance is not None and not (
        opts["tolerance"] > 0.0 and math.isfinite(opts["tolerance"])
    ):
        raise ValueError(f"{_PREFIX}the tolerance must be finite and > 0")
    if mesh_deflection is not None and not (
        opts["mesh_deflection"] > 0.0 and math.isfinite(opts["mesh_deflection"])
    ):
        raise ValueError(f"{_PREFIX}mesh_deflection must be finite and > 0")
    if opts["max_demotion_rounds"] < 0:
        raise ValueError(f"{_PREFIX}max_demotion_rounds must be >= 0")
    if not 3 <= opts["freeform_degree"] <= 8:
        raise ValueError(f"{_PREFIX}freeform_degree must be 3 to 8")
    if not 4 <= opts["freeform_max_control"] <= 1024:
        raise ValueError(f"{_PREFIX}freeform_max_control must be 4 to 1024")

    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}needs the compiled core built with {_FLAG}; "
            "there is no pure-Python CAD kernel"
        ) from exc
    if not getattr(_core, "__has_occt__", False):
        raise NotImplementedError(
            f"{_PREFIX}this build has no OpenCASCADE backend; rebuild with {_FLAG} "
            "(OpenCASCADE 7.6 or later)"
        )
    res = _core.fit_cad(
        mesh,
        opts["path"],
        opts["tolerance"],
        opts["exact_edges"],
        opts["types"],
        opts["mesh_deflection"],
        opts["max_demotion_rounds"],
        opts["freeform"],
        opts["freeform_degree"],
        opts["freeform_max_control"],
        opts["freeform_split"],
        opts["fillets"],
        opts["threads"],
    )
    out = res["mesh"]
    info = out.info if isinstance(getattr(out, "info", None), dict) else {}
    info["brep"] = res["brep"]
    out.info = info
    if not return_report:
        return out
    report = {"valid": bool(res["valid"])}
    report.update({k: float(res[k]) for k in _REPORT_FLOATS})
    report.update({k: int(res[k]) for k in _REPORT_INTS})
    return out, report


__all__ = ["fit_cad", "has_occt"]
