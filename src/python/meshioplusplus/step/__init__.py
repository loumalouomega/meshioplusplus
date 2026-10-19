"""STEP (ISO 10303-21, AP214): read through the optional OpenCASCADE backend,
written as a faceted B-rep solid without it.

Reading needs a CAD kernel: a core built with ``-DMESHIOPLUSPLUS_WITH_OCCT=ON``
(OpenCASCADE, LGPL-2.1 with the OCCT exception, bring your own). Without it,
:func:`read` raises ``NotImplementedError`` naming the flag; there is no
Python twin of a kernel.
"""

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._step import write as _py_write

_NEEDS_OCCT = (
    "reading needs the OpenCASCADE backend, which the PyPI wheels carry; rebuild with "
    "-DMESHIOPLUSPLUS_WITH_OCCT=ON (OpenCASCADE 7.6 or later, build/build-occt.sh "
    "builds one)"
)


def _read_cad(fmt, filename, deflection, angular_deflection):
    if is_buffer(filename, "r"):
        raise NotImplementedError(f"meshio++: {fmt}: reading needs a file path")
    reader = getattr(_core, f"{fmt}_read", None) if _core is not None else None
    if reader is None or not getattr(_core, "__has_occt__", False):
        raise NotImplementedError(f"meshio++: {fmt}: {_NEEDS_OCCT}")
    return reader(str(filename), float(deflection), float(angular_deflection))


def read(filename, deflection=0.0, angular_deflection=0.5):
    """Read a STEP file as its triangulation, through OpenCASCADE.

    The shape is triangulated by ``BRepMesh`` at ``deflection`` (a length; at
    most 0 takes 1e-3 of the bounding-box diagonal) and ``angular_deflection``
    (radians) into one ``triangle`` block with the ``cad:*`` cell data; the
    shape's faces, trimming loops and edges, as NURBS, are in
    ``mesh.info["brep"]``, which :func:`meshioplusplus.cad_json.write` writes as
    Kratos CAD JSON.

    :raises NotImplementedError: in a build without OpenCASCADE.
    """
    return _read_cad("step", filename, deflection, angular_deflection)


_FACES = ("auto", "planar", "fitted", "facet")


def write(filename, mesh, faces="auto", tolerance=0.0):
    """Write a closed surface (or a volume mesh's skin) as a STEP solid.

    :param faces: how facets become faces. ``"planar"`` merges edge-connected
        coplanar facets into one planar face (holes as inner bounds);
        ``"fitted"`` writes the patches :func:`meshioplusplus.segment` and
        :func:`meshioplusplus.fit_primitives` annotated as plane, sphere,
        cylinder, cone and torus faces, the rest as with ``"planar"``;
        ``"auto"`` (the default) is ``"fitted"`` when the mesh carries those
        patches, else ``"planar"``; ``"facet"`` writes one face per facet.
    :param tolerance: ``"fitted"`` only: how far, as a length, a patch's
        points may lie from its surface; at most 0 is automatic.

    The C++ core writes every mode for real file paths. The byte-identical
    Python twin writes ``"facet"`` only: ``"auto"`` falls back to it with a
    warning, and ``"planar"`` or ``"fitted"`` raise ``NotImplementedError``.
    """
    faces = str(faces)
    if faces not in _FACES:
        raise ValueError(
            "meshio++: step: faces must be 'auto', 'planar', 'fitted' or 'facet'"
        )
    if not is_buffer(filename, "w"):
        try:
            _core.step_write(str(filename), mesh, faces, float(tolerance))
            return
        except Exception as exc:
            if not core_declined(exc, "step", "write", filename):
                raise
    if faces in ("planar", "fitted"):
        raise NotImplementedError(
            f"meshio++: step: faces='{faces}' needs the C++ core "
            "(meshioplusplus._core); its Python twin writes faces='facet' only"
        )
    if faces == "auto":
        from .._common import warn

        warn(
            "step: one face per facet: merging faces needs the C++ core, "
            "so the Python twin writes faces='facet'"
        )
    return _py_write(filename, mesh)


register_format("step", [".step", ".stp"], read, {"step": write})

__all__ = ["read", "write"]
