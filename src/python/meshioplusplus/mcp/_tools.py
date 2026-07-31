"""The pure payload layer of the MCP server: one function per tool.

Everything an MCP tool does lives here, and nothing here imports the ``mcp``
SDK — this module runs on every Python the package supports (>= 3.8) and is
tested in the default CI matrix with the ``[mcp]`` extra absent, exactly the
way :mod:`meshioplusplus._interop`'s payload layer is tested without PyVista.
The FastMCP registration in :mod:`meshioplusplus.mcp._server` is a thin,
typed wrapper around the functions below.

Design rules (the contracts the tests pin):

* **Stateless and file-path based.** Every tool takes input path(s), operates
  through the public :mod:`meshioplusplus` API, writes output path(s), and
  returns a JSON-safe ``dict`` report. No mesh ever lives between calls.
* **Strict JSON.** Every report survives ``json.dumps(..., allow_nan=False)``:
  :func:`_json_safe` converts numpy scalars/arrays, replaces non-finite
  values with ``None`` (counted in ``non_finite_replaced``), and truncates
  arrays longer than ``_MAX_ARRAY`` to a ``{"truncated": ...}`` wrapper.
* **Optional path sandbox.** With a root configured (``--root`` /
  ``MESHIOPLUSPLUS_MCP_ROOT``), every path is realpath-resolved and must stay
  inside it; violations raise a clean ``ValueError``, never a traceback the
  client cannot act on. Without a root, paths resolve against the CWD.
* **``TOOL_REGISTRY`` is the single source of truth.** ``_server`` registers
  from it and ``tests/python/test_mcp.py``'s parity guard asserts every public
  operation in ``meshioplusplus.__all__`` is claimed by some tool's ``wraps``.
"""

from __future__ import annotations

import os
import pathlib
import re
from collections import OrderedDict

import numpy as np

from .. import (
    attach_quality,
    cell_data_to_point_data,
    clean,
    compute_bandwidth,
    compute_quality,
    compute_stats,
    convert_cells,
    crop,
    data_calc,
    data_condition,
    data_info,
    data_manage,
    decimate,
    diff,
    extract_skin,
    extract_surface,
    gradient,
    interpolate,
    isosurface,
    merge,
    partition,
    point_data_to_cell_data,
    read,
    read_metadata,
    refine,
    reorder,
)
from .. import screenshot as _screenshot_fn
from .. import slice as _slice_op
from .. import smooth, sniff_format, split, transform, write
from .. import write_parquet as _write_parquet_fn
from .._helpers import (
    _filetypes_from_path,
    _writer_map,
    extension_to_filetypes,
    reader_map,
)

# --------------------------------------------------------------------------- #
# Path sandbox                                                                #
# --------------------------------------------------------------------------- #
_ROOT = os.environ.get("MESHIOPLUSPLUS_MCP_ROOT") or None

# Report-size caps: arrays longer than _MAX_ARRAY become a truncation wrapper
# carrying the first _PREVIEW values; region entry previews reuse _PREVIEW.
_MAX_ARRAY = 1000
_PREVIEW = 50


def set_root(root):
    """Confine every tool's paths to ``root`` (``None`` lifts the sandbox)."""
    global _ROOT
    _ROOT = os.path.realpath(str(root)) if root else None


def get_root():
    return _ROOT


def _resolve(path, must_exist=False, for_write=False):
    """Resolve a client-supplied path, enforcing the sandbox when configured."""
    raw = os.path.expanduser(str(path))
    base = _ROOT if _ROOT is not None else os.getcwd()
    resolved = os.path.realpath(os.path.join(base, raw))
    if _ROOT is not None:
        root = os.path.realpath(_ROOT)
        if resolved != root and not resolved.startswith(root + os.sep):
            raise ValueError(
                f"meshio++: mcp: path '{path}' resolves outside the configured "
                f"root '{root}'"
            )
    if must_exist and not os.path.isfile(resolved):
        raise ValueError(f"meshio++: mcp: input file not found: '{resolved}'")
    if for_write:
        parent = os.path.dirname(resolved)
        if parent:
            os.makedirs(parent, exist_ok=True)
    return resolved


# --------------------------------------------------------------------------- #
# Strict-JSON sanitizer                                                       #
# --------------------------------------------------------------------------- #
def _sanitize(obj, max_array, preview, stats):
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        if obj != obj or obj in (float("inf"), float("-inf")):
            stats["replaced"] += 1
            return None
        return obj
    if isinstance(obj, np.generic):
        return _sanitize(obj.item(), max_array, preview, stats)
    if isinstance(obj, np.ndarray):
        if obj.size > max_array:
            head = obj.reshape(-1)[:preview]
            return {
                "truncated": True,
                "size": int(obj.size),
                "shape": [int(s) for s in obj.shape],
                "dtype": str(obj.dtype),
                "preview": _sanitize(head.tolist(), max_array, preview, stats),
            }
        return _sanitize(obj.tolist(), max_array, preview, stats)
    if isinstance(obj, dict):
        return {str(k): _sanitize(v, max_array, preview, stats) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v, max_array, preview, stats) for v in obj]
    return str(obj)


def _json_safe(obj, max_array=_MAX_ARRAY, preview=_PREVIEW):
    """Return a strict-JSON copy of ``obj`` (numpy converted, NaN/Inf -> None).

    A ``non_finite_replaced`` count is added to dict results when any value
    had to be replaced, so lossy sanitization is visible rather than silent.
    """
    stats = {"replaced": 0}
    out = _sanitize(obj, int(max_array), int(preview), stats)
    if stats["replaced"] and isinstance(out, dict):
        out["non_finite_replaced"] = stats["replaced"]
    return out


# --------------------------------------------------------------------------- #
# Shared helpers                                                              #
# --------------------------------------------------------------------------- #
def _load(input_path, input_format=None, **read_kwargs):
    return read(
        _resolve(input_path, must_exist=True), file_format=input_format, **read_kwargs
    )


def _store(mesh, output_path, output_format=None, **write_kwargs):
    resolved = _resolve(output_path, for_write=True)
    write(resolved, mesh, file_format=output_format, **write_kwargs)
    return resolved


def _mesh_summary(mesh):
    return {
        "num_points": int(len(mesh.points)),
        "num_cells": int(sum(len(b) for b in mesh.cells)),
        "cell_blocks": [{"type": b.type, "num_cells": int(len(b))} for b in mesh.cells],
        "point_data": sorted(mesh.point_data),
        "cell_data": sorted(mesh.cell_data),
        "field_data": sorted(mesh.field_data),
        "num_regions": len(mesh.regions),
    }


def _result(output_path, mesh, **extra):
    report = {"output_path": output_path}
    report.update(_mesh_summary(mesh))
    report.update(extra)
    return _json_safe(report)


def _safe_key(key):
    return re.sub(r"[^\w.+-]", "_", str(key)) or "_"


# --------------------------------------------------------------------------- #
# Inspection tools (read-only)                                                #
# --------------------------------------------------------------------------- #
def tool_formats():
    """List readable/writable formats and the extension -> formats map."""
    return _json_safe(formats_payload())


def formats_payload():
    return {
        "readable": sorted(reader_map),
        "writable": sorted(_writer_map),
        "extensions": {
            ext: list(fmts) for ext, fmts in sorted(extension_to_filetypes.items())
        },
    }


def tool_sniff(input_path):
    """Identify a mesh file's format from its leading bytes and extension."""
    resolved = _resolve(input_path, must_exist=True)
    sniffed = sniff_format(resolved)
    try:
        from_extension = _filetypes_from_path(pathlib.Path(resolved))
    except Exception:
        from_extension = []
    return _json_safe(
        {"format": sniffed or None, "from_extension": list(from_extension)}
    )


def tool_info(input_path, file_format=None):
    """Summarize a mesh file without loading its heavy arrays (fast)."""
    meta = read_metadata(_resolve(input_path, must_exist=True), file_format=file_format)
    return _json_safe(meta)


def tool_stats(input_path, file_format=None, time_step=0):
    """Geometric statistics: bbox, centroid, areas/volumes, inverted cells."""
    mesh = _load(input_path, file_format, time_step=time_step)
    return _json_safe(compute_stats(mesh))


def tool_quality(input_path, file_format=None, output_path=None, output_format=None):
    """Per-metric mesh-quality summaries; optionally write the annotated mesh.

    The per-cell metric arrays are omitted from the report (they scale with the
    mesh); pass ``output_path`` to write a copy of the mesh carrying them as
    ``quality:<metric>`` cell_data instead.
    """
    mesh = _load(input_path, file_format)
    report = compute_quality(mesh)
    report.pop("cell_arrays", None)
    report["cell_arrays"] = (
        "omitted; pass output_path to write a mesh with quality:<metric> cell_data"
    )
    if output_path is not None:
        annotated = attach_quality(mesh)
        report["output_path"] = _store(annotated, output_path, output_format)
    return _json_safe(report)


def tool_data_info(input_path, file_format=None):
    """Describe every data array: dtype, shape, ranges, NaN/Inf counts."""
    mesh = _load(input_path, file_format)
    return _json_safe({"arrays": data_info(mesh)})


def tool_regions(input_path, file_format=None):
    """List the mesh's named regions (point/cell/side groups)."""
    mesh = _load(input_path, file_format)
    regions = [
        {
            "name": r.name,
            "kind": r.kind,
            "dim": r.dim,
            "tag": r.tag,
            "num_entries": len(r),
            "entries_preview": np.asarray(r.entries).reshape(len(r), -1)[:_PREVIEW],
        }
        for r in mesh.regions
    ]
    return _json_safe({"regions": regions})


def tool_bandwidth(input_path, file_format=None):
    """Node-numbering bandwidth (max per-cell node-index spread)."""
    mesh = _load(input_path, file_format)
    return _json_safe({"bandwidth": compute_bandwidth(mesh)})


def tool_data_preview(
    input_path, array, location="point", offset=0, limit=100, file_format=None
):
    """Return a bounded window of one data array's values."""
    mesh = _load(input_path, file_format)
    if location == "point":
        table = mesh.point_data
    elif location == "cell":
        table = mesh.cell_data
    elif location == "field":
        table = mesh.field_data
    else:
        raise ValueError(
            f"meshio++: mcp: unknown location '{location}' "
            "(expected 'point', 'cell' or 'field')"
        )
    if array not in table:
        raise ValueError(
            f"meshio++: mcp: no {location} data array named '{array}' "
            f"(available: {sorted(table)})"
        )
    if location == "cell":
        blocks = [np.asarray(a) for a in table[array]]
        values = (
            np.concatenate([b.reshape(len(b), -1) for b in blocks], axis=0)
            if blocks
            else np.empty((0,))
        )
    else:
        values = np.asarray(table[array])
    offset = max(int(offset), 0)
    limit = max(int(limit), 0)
    window = values[offset : offset + limit]
    return _json_safe(
        {
            "location": location,
            "name": array,
            "dtype": str(values.dtype),
            "shape": [int(s) for s in values.shape],
            "num_rows": int(len(values)),
            "offset": offset,
            "values": _json_safe(window, max_array=window.size + 1),
        }
    )


def tool_diff(
    path_a,
    path_b,
    format_a=None,
    format_b=None,
    atol=1e-12,
    rtol=1e-9,
    unordered=False,
    max_reported=10,
):
    """Compare two mesh files; the report carries a verdict and per-section detail."""
    mesh_a = _load(path_a, format_a)
    mesh_b = _load(path_b, format_b)
    report = diff(
        mesh_a,
        mesh_b,
        atol=atol,
        rtol=rtol,
        unordered=unordered,
        max_reported=max_reported,
    )
    report["equal"] = report.get("verdict") != "different"
    return _json_safe(report)


# --------------------------------------------------------------------------- #
# Conversion                                                                  #
# --------------------------------------------------------------------------- #
# The ascii/binary/compression writer-kwargs mirror the CLI's `ascii`,
# `binary`, `compress` and `decompress` verbs (src/python/meshioplusplus/_cli/)
# so both surfaces stay in step on which formats have which variant.
_ASCII_VARIANT = {
    "ansys": {"binary": False},
    "flac3d": {"binary": False},
    "gmsh": {"binary": False},
    "mdpa": {"binary": False},
    "ply": {"binary": False},
    "stl": {"binary": False},
    "vtk": {"binary": False},
    "vtu": {"binary": False},
    "xdmf": {"data_format": "XML"},
}
_BINARY_VARIANT = {
    "ansys": {"binary": True},
    "flac3d": {"binary": True},
    "gmsh": {"binary": True},
    "mdpa": {"binary": True},
    "ply": {"binary": True},
    "stl": {"binary": True},
    "vtk": {"binary": True},
    "vtu": {"binary": True},
    "xdmf": {"data_format": "HDF"},
}
_BLOCK_CODECS = ("zlib", "lz4", "zstd", "lzma")


def _variant_kwargs(out_fmt, mode, compression):
    kwargs = {}
    if mode not in ("auto", "ascii", "binary"):
        raise ValueError(
            f"meshio++: mcp: unknown mode '{mode}' "
            "(expected 'auto', 'ascii' or 'binary')"
        )
    if mode == "ascii" and compression is not None and compression != "none":
        raise ValueError(
            "meshio++: mcp: mode='ascii' cannot be combined with compression"
        )
    if mode != "auto":
        table = _ASCII_VARIANT if mode == "ascii" else _BINARY_VARIANT
        if out_fmt not in table:
            raise ValueError(
                f"meshio++: mcp: format '{out_fmt}' has no {mode} variant "
                f"(supported: {sorted(table)})"
            )
        kwargs.update(table[out_fmt])
    if compression is None:
        return kwargs
    if compression in _BLOCK_CODECS:
        if out_fmt not in ("vtu", "vtp"):
            raise ValueError(
                f"meshio++: mcp: compression '{compression}' selects the VTK XML "
                "block codec and only vtu/vtp have one"
            )
        kwargs.update({"binary": True, "compression": compression})
    elif compression == "gzip":
        if out_fmt in ("cgns", "h5m"):
            kwargs.update({"compression": "gzip", "compression_opts": 4})
        elif out_fmt == "xdmf":
            kwargs.update(
                {"data_format": "HDF", "compression": "gzip", "compression_opts": 4}
            )
        else:
            raise ValueError(
                f"meshio++: mcp: gzip compression applies to cgns/h5m/xdmf, "
                f"not '{out_fmt}'"
            )
    elif compression == "none":
        if out_fmt in ("cgns", "h5m"):
            kwargs.update({"compression": None})
        elif out_fmt == "vtu":
            kwargs.update({"binary": True, "compression": None})
        elif out_fmt == "xdmf":
            kwargs.update({"data_format": "HDF", "compression": None})
        else:
            raise ValueError(
                f"meshio++: mcp: compression='none' applies to cgns/h5m/vtu/xdmf, "
                f"not '{out_fmt}'"
            )
    else:
        raise ValueError(
            f"meshio++: mcp: unknown compression '{compression}' (expected one of "
            f"{list(_BLOCK_CODECS)}, 'gzip' or 'none')"
        )
    return kwargs


def tool_convert(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    points_only=False,
    arrays=None,
    time_step=0,
    mode="auto",
    compression=None,
):
    """Convert between mesh formats, optionally selecting variant/compression."""
    mesh = _load(
        input_path,
        input_format,
        points_only=points_only,
        arrays=arrays,
        time_step=time_step,
    )
    out_fmt = output_format
    if out_fmt is None:
        try:
            candidates = _filetypes_from_path(pathlib.Path(str(output_path)))
            out_fmt = candidates[0] if candidates else None
        except Exception:
            out_fmt = None
    write_kwargs = {}
    if mode != "auto" or compression is not None:
        if out_fmt is None:
            raise ValueError(
                "meshio++: mcp: cannot infer the output format needed to apply "
                "mode/compression; pass output_format explicitly"
            )
        write_kwargs = _variant_kwargs(out_fmt, mode, compression)
    resolved = _store(mesh, output_path, output_format, **write_kwargs)
    return _result(resolved, mesh, output_format=out_fmt)


# --------------------------------------------------------------------------- #
# Mesh operations (path in -> path out + report)                              #
# --------------------------------------------------------------------------- #
def tool_extract_surface(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    record_parent_ids=False,
):
    """Extract the boundary surface (3D -> faces, 2D -> edges)."""
    mesh = _load(input_path, input_format)
    out = extract_surface(mesh, record_parent_ids=record_parent_ids)
    return _result(_store(out, output_path, output_format), out)


def tool_extract_skin(
    input_path, output_path, input_format=None, output_format=None, linearize=False
):
    """Extract the outer skin of a volume mesh."""
    mesh = _load(input_path, input_format)
    out = extract_skin(mesh, linearize=linearize)
    return _result(_store(out, output_path, output_format), out)


def tool_reorder(
    input_path, output_path, input_format=None, output_format=None, method="rcm"
):
    """Renumber nodes/cells (rcm | morton | hilbert); reports bandwidth change."""
    mesh = _load(input_path, input_format)
    before = compute_bandwidth(mesh)
    out = reorder(mesh, method=method)
    after = compute_bandwidth(out)
    return _result(
        _store(out, output_path, output_format),
        out,
        method=method,
        bandwidth_before=before,
        bandwidth_after=after,
    )


def tool_clean(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    weld=False,
    atol=1e-8,
    remove_orphans=True,
    drop_degenerate=True,
    drop_duplicate_cells=True,
):
    """Weld/prune/de-duplicate a mesh; reports what was removed."""
    mesh = _load(input_path, input_format)
    out, report = clean(
        mesh,
        weld=weld,
        atol=atol,
        remove_orphans=remove_orphans,
        drop_degenerate=drop_degenerate,
        drop_duplicate_cells=drop_duplicate_cells,
        return_report=True,
    )
    return _result(_store(out, output_path, output_format), out, **report)


def tool_crop(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    bbox=None,
    plane_origin=None,
    plane_normal=None,
    mode="all",
    record_ids=False,
):
    """Crop to a bounding box (6 numbers) or a half-space (origin + normal)."""
    mesh = _load(input_path, input_format)
    plane = None
    if plane_origin is not None or plane_normal is not None:
        if plane_origin is None or plane_normal is None:
            raise ValueError(
                "meshio++: mcp: crop needs both plane_origin and plane_normal"
            )
        plane = (plane_origin, plane_normal)
    out = crop(mesh, bbox=bbox, plane=plane, mode=mode, record_ids=record_ids)
    return _result(_store(out, output_path, output_format), out)


def tool_slice(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    origin=(0.0, 0.0, 0.0),
    normal=(0.0, 0.0, 1.0),
    record_parent_ids=False,
):
    """Planar cross-section (one topological dimension below the input)."""
    mesh = _load(input_path, input_format)
    out = _slice_op(
        mesh, origin=origin, normal=normal, record_parent_ids=record_parent_ids
    )
    return _result(_store(out, output_path, output_format), out)


def tool_isosurface(
    input_path,
    output_path,
    array,
    isovalues,
    input_format=None,
    output_format=None,
    component=None,
    record_parent_ids=False,
):
    """Contour a point_data scalar at one or more isovalues."""
    mesh = _load(input_path, input_format)
    out = isosurface(
        mesh,
        array,
        isovalues,
        component=component,
        record_parent_ids=record_parent_ids,
    )
    return _result(_store(out, output_path, output_format), out)


def tool_gradient(
    input_path,
    output_path,
    array,
    input_format=None,
    output_format=None,
    operator="gradient",
    method="green-gauss",
    location="cell",
    output=None,
    component=None,
    overwrite=False,
):
    """Gradient, divergence or curl of a point_data field."""
    mesh = _load(input_path, input_format)
    out, report = gradient(
        mesh,
        array,
        operator=operator,
        method=method,
        location=location,
        output=output,
        component=component,
        overwrite=overwrite,
        return_report=True,
    )
    return _result(
        _store(out, output_path, output_format),
        out,
        num_skipped=int(report["num_skipped"]),
        num_fallback=int(report["num_fallback"]),
    )


def tool_transform(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    translate=None,
    scale=None,
    rotate_axis=None,
    rotate_degrees=None,
    matrix=None,
    scale_units=None,
    rotate_vector_data=False,
):
    """Affine-transform point coordinates (translate/scale/rotate/matrix)."""
    mesh = _load(input_path, input_format)
    rotate = None
    if rotate_axis is not None or rotate_degrees is not None:
        if rotate_axis is None or rotate_degrees is None:
            raise ValueError(
                "meshio++: mcp: transform needs both rotate_axis and rotate_degrees"
            )
        rotate = (rotate_axis, rotate_degrees)
    out = transform(
        mesh,
        translate=translate,
        scale=scale,
        rotate=rotate,
        matrix=matrix,
        scale_units=scale_units,
        rotate_vector_data=rotate_vector_data,
    )
    return _result(_store(out, output_path, output_format), out)


def tool_convert_cells(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    mode="linearize",
    record_parent_ids=False,
):
    """Convert the element representation (linearize | simplexify | elevate)."""
    mesh = _load(input_path, input_format)
    out = convert_cells(mesh, mode=mode, record_parent_ids=record_parent_ids)
    return _result(_store(out, output_path, output_format), out)


def tool_refine(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    levels=1,
    record_parent_ids=False,
    cells=None,
    region=None,
    where=None,
    closure="redgreen",
    record_levels=False,
):
    """Subdivide cells into same-type children, all or a selected subset."""
    mesh = _load(input_path, input_format)
    out = refine(
        mesh,
        levels=levels,
        record_parent_ids=record_parent_ids,
        cells=cells,
        region=region,
        where=where,
        closure=closure,
        record_levels=record_levels,
    )
    return _result(_store(out, output_path, output_format), out)


def tool_decimate(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    ratio=None,
    target_faces=None,
    max_error=None,
    placement="optimal",
    preserve_boundary=True,
    preserve_features=True,
    feature_angle=30.0,
):
    """Reduce a surface mesh's face count (exactly one stopping criterion)."""
    mesh = _load(input_path, input_format)
    out, report = decimate(
        mesh,
        ratio=ratio,
        target_faces=target_faces,
        max_error=max_error,
        placement=placement,
        preserve_boundary=preserve_boundary,
        preserve_features=preserve_features,
        feature_angle=feature_angle,
        return_report=True,
    )
    return _result(_store(out, output_path, output_format), out, **report)


def tool_smooth(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    method="taubin",
    iterations=10,
    lambda_factor=-1.0,
    mu=-0.34,
    fix_boundary=True,
    preserve_features=True,
    feature_angle=30.0,
    guard_inversion=True,
):
    """Smooth point coordinates (taubin | laplacian); geometry-only move."""
    mesh = _load(input_path, input_format)
    out, report = smooth(
        mesh,
        method=method,
        iterations=iterations,
        lambda_=lambda_factor,
        mu=mu,
        fix_boundary=fix_boundary,
        preserve_features=preserve_features,
        feature_angle=feature_angle,
        guard_inversion=guard_inversion,
        return_report=True,
    )
    return _result(_store(out, output_path, output_format), out, **report)


def tool_merge(
    input_paths,
    output_path,
    output_format=None,
    weld=False,
    atol=1e-8,
    source_tag=True,
    data_policy="intersection",
    drop_duplicate_cells=False,
):
    """Merge several mesh files into one (optionally welding coincident points)."""
    if not input_paths or len(input_paths) < 2:
        raise ValueError("meshio++: mcp: merge needs at least two input_paths")
    meshes = [_load(p) for p in input_paths]
    out = merge(
        meshes,
        weld=weld,
        atol=atol,
        source_tag=source_tag,
        data_policy=data_policy,
        drop_duplicate_cells=drop_duplicate_cells,
    )
    return _result(_store(out, output_path, output_format), out)


def tool_split(
    input_path,
    input_format=None,
    by="type",
    tag=None,
    output_dir=None,
    name_template="{stem}_{key}.vtu",
):
    """Split into submeshes by type/region/component/tag; writes one file per piece."""
    resolved_in = _resolve(input_path, must_exist=True)
    mesh = read(resolved_in, file_format=input_format)
    pieces = split(mesh, by=by, tag=tag)
    stem = os.path.splitext(os.path.basename(resolved_in))[0]
    out_dir = output_dir if output_dir is not None else os.path.dirname(resolved_in)
    written = OrderedDict()
    summaries = {}
    for key, piece in pieces.items():
        name = name_template.format(stem=stem, key=_safe_key(key))
        path = _resolve(os.path.join(out_dir, name), for_write=True)
        write(path, piece)
        written[str(key)] = path
        summaries[str(key)] = _mesh_summary(piece)
    return _json_safe({"by": by, "pieces": written, "summaries": summaries})


def tool_partition(
    input_path,
    nparts,
    input_format=None,
    method="auto",
    imbalance=0.03,
    mode="eco",
    seed=0,
    record_ids=False,
    ghost_layers=0,
    weights=None,
    output_dir=None,
    name_template="{stem}_part{part}.vtu",
):
    """Partition into exactly nparts balanced pieces; writes one file per part."""
    resolved_in = _resolve(input_path, must_exist=True)
    mesh = read(resolved_in, file_format=input_format)
    pieces = partition(
        mesh,
        nparts,
        method=method,
        imbalance=imbalance,
        mode=mode,
        seed=seed,
        record_ids=record_ids,
        ghost_layers=ghost_layers,
        weights=weights,
    )
    stem = os.path.splitext(os.path.basename(resolved_in))[0]
    out_dir = output_dir if output_dir is not None else os.path.dirname(resolved_in)
    written = []
    summaries = []
    for i, piece in enumerate(pieces):
        name = name_template.format(stem=stem, part=i)
        path = _resolve(os.path.join(out_dir, name), for_write=True)
        write(path, piece)
        written.append(path)
        summaries.append(_mesh_summary(piece))
    return _json_safe(
        {
            "nparts": int(nparts),
            "method": method,
            "parts": written,
            "summaries": summaries,
        }
    )


def tool_interpolate(
    source_path,
    target_path,
    output_path,
    source_format=None,
    target_format=None,
    output_format=None,
    method="nearest",
    arrays=None,
    extrapolate=False,
    default_value=0.0,
    on_conflict="error",
):
    """Sample the source mesh's data arrays onto the target mesh's geometry."""
    source = _load(source_path, source_format)
    target = _load(target_path, target_format)
    out = interpolate(
        source,
        target,
        method=method,
        arrays=arrays,
        extrapolate=extrapolate,
        default_value=default_value,
        on_conflict=on_conflict,
    )
    return _result(_store(out, output_path, output_format), out)


# --------------------------------------------------------------------------- #
# Data operations                                                             #
# --------------------------------------------------------------------------- #
def tool_data_manage(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    keep=None,
    drop=None,
    rename=None,
    ignore_missing=False,
):
    """Keep/drop/rename data arrays.

    ``keep``/``drop`` are lists of ``[location, name]`` pairs, ``rename`` a
    list of ``[location, old, new]`` triples (location: point/cell/field).
    """
    if not (keep or drop or rename):
        raise ValueError(
            "meshio++: mcp: data_manage needs at least one of keep/drop/rename"
        )
    mesh = _load(input_path, input_format)
    result = data_manage(
        mesh, keep=keep, drop=drop, rename=rename, ignore_missing=ignore_missing
    )
    out = result["mesh"]
    return _result(
        _store(out, output_path, output_format),
        out,
        dropped=result["dropped"],
        renamed=[list(pair) for pair in result["renamed"]],
    )


def tool_data_convert(
    input_path,
    output_path,
    direction,
    input_format=None,
    output_format=None,
    arrays=None,
    weighted=False,
    nan_policy="ignore",
    nan_replacement=0.0,
):
    """Average data between locations: point_to_cell or cell_to_point."""
    mesh = _load(input_path, input_format)
    if direction == "point_to_cell":
        out = point_data_to_cell_data(
            mesh, keys=arrays, nan_policy=nan_policy, nan_replacement=nan_replacement
        )
    elif direction == "cell_to_point":
        out = cell_data_to_point_data(
            mesh,
            keys=arrays,
            weighted=weighted,
            nan_policy=nan_policy,
            nan_replacement=nan_replacement,
        )
    else:
        raise ValueError(
            f"meshio++: mcp: unknown direction '{direction}' "
            "(expected 'point_to_cell' or 'cell_to_point')"
        )
    return _result(_store(out, output_path, output_format), out)


def tool_data_calc(
    input_path,
    output_path,
    expression,
    input_format=None,
    output_format=None,
    location="point",
    output_name="",
    overwrite=False,
):
    """Evaluate an arithmetic expression over data arrays into a new array.

    Accepts the CLI's ``"NAME = EXPR"`` spelling (split on the first ``=``,
    the documented colon rule) when ``output_name`` is not given separately.
    """
    if not output_name and "=" in expression:
        output_name, expression = (s.strip() for s in expression.split("=", 1))
    mesh = _load(input_path, input_format)
    out = data_calc(
        mesh, expression, location=location, output=output_name, overwrite=overwrite
    )
    return _result(_store(out, output_path, output_format), out)


def tool_data_condition(
    input_path,
    output_path,
    input_format=None,
    output_format=None,
    operation="clamp",
    location="point",
    arrays=None,
    scope="component",
    lo=0.0,
    hi=1.0,
    nan_policy="ignore",
    nan_replacement=0.0,
    suffix="",
    preserve_dtype=True,
):
    """Condition data values: clamp, normalize or standardize."""
    mesh = _load(input_path, input_format)
    out = data_condition(
        mesh,
        location=location,
        keys=arrays,
        mode=operation,
        scope=scope,
        lo=lo,
        hi=hi,
        nan_policy=nan_policy,
        nan_replacement=nan_replacement,
        suffix=suffix,
        preserve_dtype=preserve_dtype,
    )
    return _result(_store(out, output_path, output_format), out)


# --------------------------------------------------------------------------- #
# Gated tools (optional extras; the wrapped functions raise the named error)  #
# --------------------------------------------------------------------------- #
def tool_data_export(input_path, output_path, input_format=None, location="point"):
    """Export data arrays to a Parquet table (needs the [arrow] extra)."""
    mesh = _load(input_path, input_format)
    resolved = _resolve(output_path, for_write=True)
    _write_parquet_fn(mesh, resolved, location=location)
    return _json_safe(
        {"output_path": resolved, "location": location, **_mesh_summary(mesh)}
    )


def tool_screenshot(
    input_path,
    output_path,
    input_format=None,
    color_by=None,
    width=1280,
    height=960,
    transparent=False,
):
    """Render an off-screen screenshot PNG (needs the [viewer] extra)."""
    mesh = _load(input_path, input_format)
    resolved = _resolve(output_path, for_write=True)
    _screenshot_fn(
        mesh,
        resolved,
        color_by=color_by,
        size=(int(width), int(height)),
        transparent=transparent,
    )
    return _json_safe({"output_path": resolved, "size": [int(width), int(height)]})


# --------------------------------------------------------------------------- #
# Registry — the single source of truth                                       #
# --------------------------------------------------------------------------- #
# name -> {"fn", "wraps" (public meshioplusplus API names this tool covers,
# consumed by the parity-guard test), "gated" (extra name, or None)}.
TOOL_REGISTRY = OrderedDict(
    [
        (
            "formats",
            {"fn": tool_formats, "wraps": ("extension_to_filetypes",), "gated": None},
        ),
        ("sniff", {"fn": tool_sniff, "wraps": ("sniff_format",), "gated": None}),
        ("info", {"fn": tool_info, "wraps": ("read_metadata",), "gated": None}),
        ("stats", {"fn": tool_stats, "wraps": ("compute_stats",), "gated": None}),
        (
            "quality",
            {
                "fn": tool_quality,
                "wraps": ("compute_quality", "attach_quality"),
                "gated": None,
            },
        ),
        ("data_info", {"fn": tool_data_info, "wraps": ("data_info",), "gated": None}),
        ("regions", {"fn": tool_regions, "wraps": (), "gated": None}),
        (
            "bandwidth",
            {"fn": tool_bandwidth, "wraps": ("compute_bandwidth",), "gated": None},
        ),
        ("data_preview", {"fn": tool_data_preview, "wraps": (), "gated": None}),
        ("diff", {"fn": tool_diff, "wraps": ("diff", "meshes_equal"), "gated": None}),
        ("convert", {"fn": tool_convert, "wraps": ("read", "write"), "gated": None}),
        (
            "extract_surface",
            {"fn": tool_extract_surface, "wraps": ("extract_surface",), "gated": None},
        ),
        (
            "extract_skin",
            {"fn": tool_extract_skin, "wraps": ("extract_skin",), "gated": None},
        ),
        ("reorder", {"fn": tool_reorder, "wraps": ("reorder",), "gated": None}),
        ("clean", {"fn": tool_clean, "wraps": ("clean",), "gated": None}),
        ("crop", {"fn": tool_crop, "wraps": ("crop",), "gated": None}),
        ("slice", {"fn": tool_slice, "wraps": ("slice",), "gated": None}),
        (
            "isosurface",
            {"fn": tool_isosurface, "wraps": ("isosurface",), "gated": None},
        ),
        ("gradient", {"fn": tool_gradient, "wraps": ("gradient",), "gated": None}),
        ("transform", {"fn": tool_transform, "wraps": ("transform",), "gated": None}),
        (
            "convert_cells",
            {"fn": tool_convert_cells, "wraps": ("convert_cells",), "gated": None},
        ),
        ("refine", {"fn": tool_refine, "wraps": ("refine",), "gated": None}),
        ("decimate", {"fn": tool_decimate, "wraps": ("decimate",), "gated": None}),
        ("smooth", {"fn": tool_smooth, "wraps": ("smooth",), "gated": None}),
        ("merge", {"fn": tool_merge, "wraps": ("merge",), "gated": None}),
        ("split", {"fn": tool_split, "wraps": ("split",), "gated": None}),
        ("partition", {"fn": tool_partition, "wraps": ("partition",), "gated": None}),
        (
            "interpolate",
            {"fn": tool_interpolate, "wraps": ("interpolate",), "gated": None},
        ),
        (
            "data_manage",
            {
                "fn": tool_data_manage,
                "wraps": ("data_manage", "data_drop", "data_keep", "data_rename"),
                "gated": None,
            },
        ),
        (
            "data_convert",
            {
                "fn": tool_data_convert,
                "wraps": ("point_data_to_cell_data", "cell_data_to_point_data"),
                "gated": None,
            },
        ),
        ("data_calc", {"fn": tool_data_calc, "wraps": ("data_calc",), "gated": None}),
        (
            "data_condition",
            {"fn": tool_data_condition, "wraps": ("data_condition",), "gated": None},
        ),
        (
            "data_export",
            {"fn": tool_data_export, "wraps": ("write_parquet",), "gated": "arrow"},
        ),
        (
            "screenshot",
            {"fn": tool_screenshot, "wraps": ("screenshot",), "gated": "viewer"},
        ),
    ]
)
