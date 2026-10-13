import os
import pathlib

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._mdpa import read as _py_read
from ._mdpa import write as _py_write


def _has_mdpa_extras(mesh):
    """Whether ``mesh`` carries MDPA data only the Python reference can write.

    ``misc_data`` (SubModelPart hierarchy, original Kratos ids, ``Mesh``
    blocks), ``geometries_block`` and the ``field_data`` conventions
    (``properties_<id>``/``table_<id>`` dicts, gmsh physical-group tuples) have
    no counterpart in the C++ ``Mesh``, so a mesh holding any of them keeps the
    pure-Python writer.
    """
    if getattr(mesh, "misc_data", None):
        return True
    if getattr(mesh, "geometries_block", None):
        return True
    return bool(getattr(mesh, "field_data", None))


def cad_json_path(filename):
    """The Kratos CAD JSON geometry that goes with an MDPA deck, or ``None``.

    Kratos keeps IGA geometry in a ``*.cad.json`` file its ``CadIoModeler``
    reads, not in the deck: ``<stem>.cad.json`` beside the deck, else the
    ``geometry_file_name`` of the ``CadIoModeler`` in a ``ProjectParameters.json``
    beside it (``.cad.json`` appended when the name lacks it, as Kratos does).
    The twin of C++ ``mdpa_cad_json_path``.
    """
    from .. import _json_text
    from .._exceptions import ReadError

    path = pathlib.Path(filename)
    stem = (
        path.name[:-5]
        if path.name.endswith(".mdpa") and len(path.name) > 5
        else path.name
    )
    own = path.parent / (stem + ".cad.json")
    if own.is_file():
        return str(own)
    params = path.parent / "ProjectParameters.json"
    if not params.is_file():
        return None
    try:
        root = _json_text.parse(params.read_bytes(), "ProjectParameters.json")
    except ReadError:
        return None  # not ours to judge; there is simply no sidecar named
    modelers = (
        root.find("modelers") if isinstance(root, _json_text.JsonObject) else None
    )
    if not isinstance(modelers, list) or isinstance(modelers, _json_text.JsonObject):
        return None
    for m in modelers:
        if (
            not isinstance(m, _json_text.JsonObject)
            or m.find("modeler_name") != "CadIoModeler"
        ):
            continue
        parameters = m.find("Parameters")
        if not isinstance(parameters, _json_text.JsonObject):
            continue
        name = parameters.find("geometry_file_name")
        if not isinstance(name, str) or not name:
            continue
        if not name.endswith(".cad.json"):
            name += ".cad.json"
        candidate = path.parent / name
        if candidate.is_file():
            return str(candidate)
    return None


def _attach_cad_json(mesh, filename, cad_json):
    from .. import cad_json as _cad_json
    from .._exceptions import ReadError

    if cad_json is True:
        sidecar = cad_json_path(filename)
    else:
        sidecar = os.fspath(cad_json)
        if not os.path.isfile(sidecar):
            raise ReadError(
                f"meshio++: mdpa: the CAD JSON geometry '{sidecar}' does not exist"
            )
    if sidecar is None:
        return mesh
    try:
        brep = _core.cad_json_read_info(str(sidecar))
    except Exception as exc:
        if not core_declined(exc, "cad_json", "read", sidecar):
            raise
        brep = _cad_json.read_brep(sidecar)
    info = mesh.info if isinstance(getattr(mesh, "info", None), dict) else {}
    info["brep"] = brep
    mesh.info = info
    return mesh


def read(filename, cad_json=None):
    """Read a Kratos MDPA mesh file.

    ``cad_json`` also loads the Kratos IGA geometry that goes with the deck
    into ``mesh.info["brep"]`` (the CAD JSON side channel; the mesh itself is
    the deck's): ``True`` finds it as :func:`cad_json_path` does (none found
    is no error), a path names it. The default reads the deck alone.

    This deliberately stays on the pure-Python reference reader
    (``meshioplusplus.mdpa._mdpa``) rather than trying ``_core.mdpa_read``
    first, unlike the other format shims. The C++ core *does* implement the
    format (``meshioplusplus::read_mdpa``, reachable from the C API, Fortran,
    Julia, R, WebAssembly and the native CLI through the shared registry), but
    it necessarily returns a repo-standard ``Mesh``: ``cell_data`` keyed by
    array name with one array per cell block, no ``misc_data`` and no
    ``geometries_block``. The reference reader instead keys ``cell_data`` by
    cell type and carries MDPA's ``Tables``/``Properties``/``Geometries``/
    ``Mesh`` blocks and the SubModelPart hierarchy in ``mesh.misc_data`` — the
    documented behaviour of this module. Preferring the C++ reader here would
    silently change all of that, so it is reached explicitly
    (``_core.mdpa_read``) or through the flat bindings instead.
    """
    mesh = _py_read(filename)
    if cad_json is not None and cad_json is not False and not is_buffer(filename, "r"):
        _attach_cad_json(mesh, filename, cad_json)
    return mesh


def write(filename, mesh, float_fmt=".16e", binary=False, cad_json=None):
    """Write a Kratos MDPA mesh file (C++ core when possible, Python fallback).

    The C++ writer is used for real file paths and meshes carrying none of the
    MDPA-specific extras listed in :func:`_has_mdpa_extras`; anything else (and
    any failure in the C++ path) falls back to the pure-Python reference
    writer, which is the only one that can reproduce them.

    ``cad_json`` also writes ``mesh.info["brep"]`` (the IGA geometry
    :func:`read` loaded) as Kratos CAD JSON: ``True`` to ``<stem>.cad.json``
    beside the deck, a path to that path.
    """
    if cad_json is not None and cad_json is not False:
        from .. import cad_json as _cad_json
        from .._exceptions import WriteError

        brep = (
            (mesh.info or {}).get("brep")
            if isinstance(getattr(mesh, "info", None), dict)
            else None
        )
        if not brep:
            raise WriteError(
                "meshio++: mdpa: cad_json= needs a B-rep in mesh.info['brep']"
            )
        if cad_json is True:
            path = pathlib.Path(filename)
            stem = (
                path.name[:-5]
                if path.name.endswith(".mdpa") and len(path.name) > 5
                else path.name
            )
            sidecar = path.parent / (stem + ".cad.json")
        else:
            sidecar = cad_json
        _cad_json.write(sidecar, mesh, brep=brep)
    if not binary and not is_buffer(filename, "w") and not _has_mdpa_extras(mesh):
        try:
            _core.mdpa_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "mdpa", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt=float_fmt, binary=binary)


register_format("mdpa", [".mdpa"], read, {"mdpa": write})

__all__ = ["read", "write", "cad_json_path"]
