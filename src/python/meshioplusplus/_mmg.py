"""The optional MMG backend (roadmap 3.2): ``adapt`` and its metric helper.

`MMG <https://www.mmgtools.org/>`_ -- ``mmg2d``, ``mmgs`` and ``mmg3d`` --
remeshes simplicial meshes. It is LGPL-3.0 and linked into a build only when
configured with ``-DMESHIOPLUSPLUS_WITH_MMG=ON`` (``doc/mmg.md``), never into
the published wheels. MMG's result cannot be reproduced by a second
implementation, so there is no numpy twin (the :func:`remesh_volume`
precedent): without the backend these raise ``NotImplementedError`` naming the
flag. :func:`metric_from_hessian` needs no MMG.

Public API:

* :func:`has_mmg`, :func:`has_mmg_elas` -- what the compiled core carries.
* :func:`adapt` -- remesh to a constant size, a size field or a metric.
* :func:`metric_from_hessian` -- an anisotropic metric from a Hessian.
* :func:`discretize_levelset` -- embed a level set and remesh both sides.
* :func:`mesh_domain` -- triangulate a planar domain from its boundary lines.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "has_mmg",
    "has_mmg_elas",
    "adapt",
    "metric_from_hessian",
    "discretize_levelset",
    "mesh_domain",
]

_PREFIX = "meshio++: "


def _core_or_raise(what):
    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}{what}: needs the compiled C++ core (meshioplusplus._core) "
            "built with -DMESHIOPLUSPLUS_WITH_MMG=ON; it has no numpy twin"
        ) from exc
    return _core


def has_mmg() -> bool:
    """Whether the compiled core carries the MMG backend."""
    try:
        from . import _core
    except ImportError:
        return False
    return bool(getattr(_core, "__has_mmg__", False))


def has_mmg_elas() -> bool:
    """Whether the core's MMG was built with ELAS, which ``move`` needs."""
    try:
        from . import _core
    except ImportError:
        return False
    return bool(getattr(_core, "__has_mmg_elas__", False))


def _frozen_ids(mesh, frozen, what):
    """Point ids to pin: an index array, a boolean mask or a point region's name."""
    if frozen is None:
        return None
    if isinstance(frozen, str):
        for region in getattr(mesh, "regions", None) or []:
            if region.name == frozen and region.kind == "point":
                return np.asarray(region.entries, dtype=np.int64).ravel()
        raise ValueError(f"{_PREFIX}{what}: no point region named {frozen!r}")
    frozen = np.asarray(frozen)
    if frozen.dtype == bool:
        if frozen.shape != (len(mesh.points),):
            raise ValueError(
                f"{_PREFIX}{what}: a frozen mask needs one entry per point "
                f"({len(mesh.points)}), not {frozen.shape}"
            )
        return np.flatnonzero(frozen).astype(np.int64)
    return frozen.astype(np.int64).ravel()


def _parameters(
    mesh,
    what,
    hmin,
    hmax,
    hausd,
    hgrad,
    hgradreq,
    angle_detection,
    optim,
    noinsert,
    noswap,
    nomove,
    nosurf,
    nreg,
    opnbdy,
    mem,
    verbose,
    frozen,
    required,
    local,
):
    """The dict of shared MMG parameters the core reads."""
    if isinstance(required, str):
        required = [required]
    return {
        "hmin": hmin,
        "hmax": hmax,
        "hausd": hausd,
        "hgrad": hgrad,
        "hgradreq": hgradreq,
        "angle_detection": angle_detection,
        "optim": bool(optim),
        "noinsert": bool(noinsert),
        "noswap": bool(noswap),
        "nomove": bool(nomove),
        "nosurf": bool(nosurf),
        "nreg": bool(nreg),
        "opnbdy": bool(opnbdy),
        "mem": mem,
        "verbose": int(verbose),
        "frozen": _frozen_ids(mesh, frozen, what),
        "required": list(required) if required else None,
        "local": dict(local) if local else None,
    }


def _field_name(mesh, field, default, what, kind):
    """A point-data name for ``field``: as given, or an array put on a shallow copy."""
    if field is None:
        return mesh, ""
    if isinstance(field, str):
        return mesh, field
    array = np.asarray(field, dtype=float)
    if len(array) != len(mesh.points):
        raise ValueError(
            f"{_PREFIX}{what}: the {kind} array needs one row per point "
            f"({len(mesh.points)}), not {len(array)}"
        )
    from ._mesh import Mesh

    copy = Mesh(
        mesh.points,
        mesh.cells,
        point_data={**mesh.point_data, default: array},
        cell_data=dict(mesh.cell_data),
        field_data=dict(mesh.field_data),
        regions=list(getattr(mesh, "regions", None) or []),
    )
    return copy, default


def adapt(
    mesh,
    hsiz=None,
    size=None,
    metric=None,
    *,
    hmin=None,
    hmax=None,
    hausd=None,
    hgrad=None,
    hgradreq=None,
    angle_detection=45.0,
    optim=False,
    noinsert=False,
    noswap=False,
    nomove=False,
    nosurf=False,
    nreg=False,
    opnbdy=False,
    mem=None,
    verbose=-1,
    frozen=None,
    required=None,
    local=None,
    return_report=False,
):
    """Remesh ``mesh`` to a size or a metric with MMG (roadmap 3.2.1).

    A mesh with tetrahedra goes to ``mmg3d`` (wedges pass as MMG prisms,
    triangles and quads as boundary faces, lines as ridges), a planar triangle
    mesh (2-D points, or 3-D ones on one z) to ``mmg2d`` and any other
    triangle surface to ``mmgs``.

    :param hsiz: a constant target edge size.
    :param size: the target edge size per point: a ``point_data`` name or an
        array of one value per point.
    :param metric: an anisotropic metric per point: a ``point_data`` name or a
        ``(n, 6)`` (``(n, 3)`` in 2-D) array in the lower-triangle order of a
        Medit ``.sol``, ``m11 m21 m22 m31 m32 m33`` -- the order
        :func:`metric_from_hessian` writes. Give at most one of the three;
        none lets MMG size from the mesh itself.
    :param hmin, hmax: the smallest and largest edge size.
    :param hausd: how far the boundary may move from the input's.
    :param hgrad, hgradreq: the gradation, and the gradation from required
        entities.
    :param angle_detection: the ridge angle in degrees; ``None`` or ``False``
        turns the detection off.
    :param optim, noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem: MMG's
        switches of those names (``nosurf`` and ``opnbdy`` not for mmgs).
    :param verbose: MMG's verbosity (-1, the default, is silent).
    :param frozen: points MMG must keep: ids, a boolean mask, or a point
        region's name.
    :param required: names of regions whose entities MMG must keep.
    :param local: per-region sizes, ``{region: {"hmin": .., "hmax": ..,
        "hausd": ..}}``.
    :param return_report: also return the report dict (``status``
        ``"success"`` or ``"low_failure"``, ``library``, point and element
        counts in and out, the ``quality`` summary of
        ``quality:aspect_ratio``, and ``warnings``).
    :returns: the adapted mesh, with cell and side regions (and a
        ``medit:ref`` array) carried under their names, point regions when
        their points are required, ``field_data`` copied; other point and cell
        data are dropped with a warning (compose with :func:`interpolate`).
    :raises NotImplementedError: without the MMG backend.
    :raises ValueError: for cells MMG does not take (convert_cells first), a
        missing field or region, or a bad parameter.
    :raises RuntimeError: when MMG fails (``MMG5_STRONGFAILURE``).
    """
    _core = _core_or_raise("adapt")
    given = sum(x is not None for x in (hsiz, size, metric))
    if given > 1:
        raise ValueError(f"{_PREFIX}adapt: give one of hsiz, size and metric")
    if hsiz is not None and not hsiz > 0:
        raise ValueError(f"{_PREFIX}adapt: hsiz must be positive")
    work, size_name = _field_name(mesh, size, "adapt:size", "adapt", "size")
    work, metric_name = _field_name(work, metric, "adapt:metric", "adapt", "metric")
    params = _parameters(
        mesh,
        "adapt",
        hmin,
        hmax,
        hausd,
        hgrad,
        hgradreq,
        angle_detection,
        optim,
        noinsert,
        noswap,
        nomove,
        nosurf,
        nreg,
        opnbdy,
        mem,
        verbose,
        frozen,
        required,
        local,
    )
    out = _core.adapt(
        work,
        hsiz=-1.0 if hsiz is None else float(hsiz),
        size=size_name,
        metric=metric_name,
        params=params,
    )
    result = out.pop("mesh")
    return (result, out) if return_report else result


_MMG_KEYWORDS = (
    "hmin",
    "hmax",
    "hausd",
    "hgrad",
    "hgradreq",
    "angle_detection",
    "optim",
    "noinsert",
    "noswap",
    "nomove",
    "nosurf",
    "nreg",
    "opnbdy",
    "mem",
    "verbose",
    "frozen",
    "required",
    "local",
)


def _mmg_keywords(mesh, what, options):
    """The shared MMG parameters from keyword arguments, as :func:`adapt` takes them."""
    unknown = sorted(set(options) - set(_MMG_KEYWORDS))
    if unknown:
        raise TypeError(f"{_PREFIX}{what}: unexpected keyword(s) {', '.join(unknown)}")
    get = options.get
    return _parameters(
        mesh,
        what,
        get("hmin"),
        get("hmax"),
        get("hausd"),
        get("hgrad"),
        get("hgradreq"),
        get("angle_detection", 45.0),
        get("optim", False),
        get("noinsert", False),
        get("noswap", False),
        get("nomove", False),
        get("nosurf", False),
        get("nreg", False),
        get("opnbdy", False),
        get("mem"),
        get("verbose", -1),
        get("frozen"),
        get("required"),
        get("local"),
    )


def discretize_levelset(
    mesh,
    level_set,
    isovalue=0.0,
    *,
    materials=None,
    remove_small_components=None,
    return_report=False,
    **options,
):
    """Embed a level set in ``mesh`` and remesh both sides with MMG (roadmap 3.2.2).

    ``mmg3d`` for tetrahedra, ``mmg2d`` for a planar triangle mesh, ``mmgs``
    for a surface. The cut follows the ``isovalue`` of ``level_set`` exactly
    and conformally, which :func:`isosurface` (an extracted surface) and
    :func:`crop` (whole cells) do not.

    :param level_set: the level set per point: a ``point_data`` name or an
        array (a signed distance, such as :func:`compute_sdf` gives).
    :param isovalue: the value to cut at.
    :param materials: several materials, ``[(ref, split, ref_inside,
        ref_outside), ...]`` over the input's ``medit:ref``: the cells of
        ``ref`` are cut into ``ref_inside`` and ``ref_outside``, or kept whole
        when ``split`` is false (MMG's ``Set_multiMat``).
    :param remove_small_components: remove the parasitic components of the
        cut smaller than this fraction of the volume (MMG's ``rmc``).
    :param options: :func:`adapt`'s MMG parameters (``hmin``, ``hmax``,
        ``hausd``, ..., ``frozen``, ``required``, ``local``).
    :returns: the cut mesh. Its elements carry MMG's references in
        ``cell_data["medit:ref"]`` -- 3 inside (below the isovalue), 2
        outside, or the materials' -- and the cell regions ``inside`` and
        ``outside``; the interface is the side region ``interface``, on the
        inside elements. The input's side regions come back, its cell regions
        do not (name them with ``materials``); other data is dropped. With
        ``return_report``, also :func:`adapt`'s report.
    :raises NotImplementedError: without the MMG backend.
    """
    _core = _core_or_raise("discretize_levelset")
    work, name = _field_name(
        mesh, level_set, "discretize_levelset:ls", "discretize_levelset", "level set"
    )
    params = _mmg_keywords(mesh, "discretize_levelset", options)
    out = _core.discretize_levelset(
        work,
        name,
        isovalue=float(isovalue),
        remove_small_components=(
            -1.0 if remove_small_components is None else float(remove_small_components)
        ),
        materials=[
            (int(r), bool(sp), int(rin), int(rout))
            for r, sp, rin, rout in materials or []
        ],
        params=params,
    )
    result = out.pop("mesh")
    return (result, out) if return_report else result


def mesh_domain(
    boundary, hsiz=None, size=None, *, holes=None, return_report=False, **options
):
    """Triangulate the planar domain bounded by ``boundary``'s lines with MMG
    (``MMG2D_mmg2dmesh``, roadmap 3.2.3).

    MMG meshes every region the lines enclose and gives each its own
    reference; it has no notion of a hole, so ``holes`` names the regions to
    leave out by a point inside each. :func:`triangulate_polygon_2d` meshes a
    polygon with holes natively at one size; this adds MMG's size fields.

    :param boundary: a mesh of ``line`` cells (and vertices), 2-D or planar.
    :param hsiz: a constant target edge size.
    :param size: the target edge size at the boundary's points: a
        ``point_data`` name or an array.
    :param holes: one ``(x, y)`` point inside each hole.
    :param options: :func:`adapt`'s MMG parameters.
    :returns: the triangles and the boundary lines; each enclosed region is a
        cell region ``subdomain:<k>`` (``k`` MMG's reference), and the lines'
        cell regions are also side regions of the triangles. With
        ``return_report``, also :func:`adapt`'s report plus
        ``num_subdomains`` and ``num_holes``.
    :raises NotImplementedError: without the MMG backend.
    :raises ValueError: for other cells than lines, non-planar points, or a
        hole point outside the domain.
    """
    _core = _core_or_raise("mesh_domain")
    if hsiz is not None and size is not None:
        raise ValueError(f"{_PREFIX}mesh_domain: give one of hsiz and size")
    work, size_name = _field_name(
        boundary, size, "mesh_domain:size", "mesh_domain", "size"
    )
    params = _mmg_keywords(boundary, "mesh_domain", options)
    flat = [float(c) for point in holes or [] for c in point[:2]]
    out = _core.mesh_domain(
        work,
        hsiz=-1.0 if hsiz is None else float(hsiz),
        size=size_name,
        holes=flat,
        params=params,
    )
    result = out.pop("mesh")
    return (result, out) if return_report else result


def move(mesh, displacement, moving, *, lagrangian=1, return_report=False, **options):
    """Move ``mesh`` by a displacement of part of its boundary and remesh as it
    deforms, with MMG's Lagrangian mode (roadmap 3.2.4).

    ``MMG3D_mmg3dmov`` for tetrahedra, ``MMG2D_mmg2dmov`` for a planar
    triangle mesh. The boundary facets of the ``moving`` regions follow the
    displacement; MMG extends it into a layer of elements around them by
    linear elasticity (the ELAS library) and keeps the rest of the boundary
    in place, so the displacement given elsewhere is not read.

    :param displacement: the displacement per point: a ``point_data`` name or
        an ``(n, dim)`` array (``dim`` 3 for mmg3d, 2 or 3 for mmg2d).
    :param moving: the name of a cell or side region, or a list of them,
        whose boundary facets move.
    :param lagrangian: MMG's ``lag``: 0 moves only, 1 (the default) also
        swaps edges and faces, 2 also inserts and collapses points.
    :param options: :func:`adapt`'s MMG parameters.
    :returns: the moved mesh, regions as for :func:`adapt`. With
        ``return_report``, also :func:`adapt`'s report.
    :raises NotImplementedError: without the MMG backend, or with an MMG
        built without ELAS (:func:`has_mmg_elas`).
    """
    _core = _core_or_raise("move")
    work, name = _field_name(
        mesh, displacement, "move:displacement", "move", "displacement"
    )
    if not name:
        raise ValueError(f"{_PREFIX}move: give the displacement")
    names = [moving] if isinstance(moving, str) else list(moving)
    params = _mmg_keywords(mesh, "move", options)
    out = _core.move(
        work, name, moving=names, lagrangian=int(lagrangian), params=params
    )
    result = out.pop("mesh")
    return (result, out) if return_report else result


def metric_from_hessian(
    mesh, hessian, error=1e-2, hmin=1e-3, hmax=1.0, output="metric"
):
    """The anisotropic metric of a Hessian, for :func:`adapt` (roadmap 3.2.1.1).

    Per point, ``M = c / error * |H|``: the Hessian symmetrised, its
    eigenvalues taken in absolute value, scaled by ``c = (1/2)(d/(d+1))^2``
    for P1 interpolation in ``d`` dimensions and clamped so the sizes it asks
    for lie in ``[hmin, hmax]``.

    :param hessian: the ``point_data`` Hessian, ``(n, 9)`` (``(n, 4)`` in 2-D)
        row-major, such as :func:`hessian` writes.
    :param error: the target interpolation error.
    :param hmin, hmax: the smallest and largest edge size the metric asks for.
    :param output: the ``point_data`` name of the result.
    :returns: a copy of the mesh with ``point_data[output]``, six components
        (three for a 2-D Hessian) in the lower-triangle order of a Medit
        ``.sol``. Needs no MMG.
    """
    try:
        from . import _core
    except ImportError as exc:
        raise NotImplementedError(
            f"{_PREFIX}metric_from_hessian: needs the compiled C++ core "
            "(meshioplusplus._core); it has no numpy twin"
        ) from exc
    return _core.metric_from_hessian(
        mesh,
        hessian,
        error=float(error),
        hmin=float(hmin),
        hmax=float(hmax),
        output=output,
    )
