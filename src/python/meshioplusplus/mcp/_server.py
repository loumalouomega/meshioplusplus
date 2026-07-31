"""FastMCP registration layer: typed wrappers over :mod:`._tools`.

This is the only module in the package that imports the ``mcp`` SDK, so it is
only importable with the ``[mcp]`` extra installed (which needs Python >= 3.10
— the SDK's own floor, not meshio++'s). All behaviour lives in the pure layer;
each wrapper here contributes exactly two things: the typed signature FastMCP
derives the tool's JSON schema from, and the docstring the client shows the
model. Wrappers never leak tracebacks — every call is guarded and failures
come back as ``{"error", "error_type"}`` payloads the model can act on.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import List, Optional, Union

from mcp.server.fastmcp import FastMCP, Image

from ..__about__ import __version__
from . import _tools

_INSTRUCTIONS = """\
meshio++ mesh I/O and processing. All tools are stateless and file-path based:
they read mesh files (40+ formats: VTK/VTU/VTP, Gmsh, Exodus, MED, Abaqus,
STL, OBJ, PLY, XDMF, CGNS, ...), operate, and write output files, returning a
JSON report. Start with `formats` to see supported formats, `info` for a fast
file summary, and `convert` to translate between formats. Inspection tools
(`stats`, `quality`, `data_info`, `regions`, `diff`) return reports without
writing anything. Mesh operations (clean, refine, decimate, smooth, slice,
partition, ...) take input_path/output_path. Paths may be confined to a root
directory by the server's --root option.
"""


def _guard(fn, /, **kwargs):
    try:
        return fn(**kwargs)
    except Exception as e:  # noqa: BLE001 — the client needs a payload, not a traceback
        return {"error": str(e), "error_type": type(e).__name__}


def create_server(root: Optional[str] = None) -> FastMCP:
    """Build the FastMCP server with every tool and resource registered."""
    if root is not None:
        _tools.set_root(root)
    server = FastMCP("meshioplusplus", instructions=_INSTRUCTIONS)
    _register_inspection(server)
    _register_conversion(server)
    _register_operations(server)
    _register_data(server)
    _register_gated(server)
    _register_resources(server)
    return server


# --------------------------------------------------------------------------- #
# Inspection                                                                  #
# --------------------------------------------------------------------------- #
def _register_inspection(server: FastMCP) -> None:
    @server.tool()
    def formats() -> dict:
        """List every readable and writable mesh format and the extension map."""
        return _guard(_tools.tool_formats)

    @server.tool()
    def sniff(input_path: str) -> dict:
        """Identify a mesh file's format from its leading bytes and extension."""
        return _guard(_tools.tool_sniff, input_path=input_path)

    @server.tool()
    def info(input_path: str, file_format: Optional[str] = None) -> dict:
        """Fast file summary without loading heavy arrays: point/cell counts,
        cell blocks, data-array names, named regions, time steps, format."""
        return _guard(_tools.tool_info, input_path=input_path, file_format=file_format)

    @server.tool()
    def stats(
        input_path: str, file_format: Optional[str] = None, time_step: int = 0
    ) -> dict:
        """Geometric statistics: bounding box, centroid, total area, signed and
        unsigned volume, per-type cell counts, inverted-cell count."""
        return _guard(
            _tools.tool_stats,
            input_path=input_path,
            file_format=file_format,
            time_step=time_step,
        )

    @server.tool()
    def quality(
        input_path: str,
        file_format: Optional[str] = None,
        output_path: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Mesh-quality metrics (scaled jacobian, aspect ratio, skewness,
        angles, ...) as per-metric min/max/mean + histograms. Pass output_path
        to also write a mesh copy carrying quality:<metric> cell data."""
        return _guard(
            _tools.tool_quality,
            input_path=input_path,
            file_format=file_format,
            output_path=output_path,
            output_format=output_format,
        )

    @server.tool()
    def data_info(input_path: str, file_format: Optional[str] = None) -> dict:
        """Describe every point/cell/field data array: dtype, shape, components,
        min/max/mean (whole-array and per component), NaN/Inf counts."""
        return _guard(
            _tools.tool_data_info, input_path=input_path, file_format=file_format
        )

    @server.tool()
    def regions(input_path: str, file_format: Optional[str] = None) -> dict:
        """List named regions (point/cell/side groups: gmsh physical groups,
        Exodus sets, Abaqus NSET/ELSET, ...) with kind, dim, tag and size."""
        return _guard(
            _tools.tool_regions, input_path=input_path, file_format=file_format
        )

    @server.tool()
    def bandwidth(input_path: str, file_format: Optional[str] = None) -> dict:
        """Node-numbering bandwidth (max node-index spread over any cell)."""
        return _guard(
            _tools.tool_bandwidth, input_path=input_path, file_format=file_format
        )

    @server.tool()
    def data_preview(
        input_path: str,
        array: str,
        location: str = "point",
        offset: int = 0,
        limit: int = 100,
        file_format: Optional[str] = None,
    ) -> dict:
        """Return a bounded window of one data array's values (location:
        point | cell | field; cell arrays are concatenated block-major)."""
        return _guard(
            _tools.tool_data_preview,
            input_path=input_path,
            array=array,
            location=location,
            offset=offset,
            limit=limit,
            file_format=file_format,
        )

    @server.tool()
    def diff(
        path_a: str,
        path_b: str,
        format_a: Optional[str] = None,
        format_b: Optional[str] = None,
        atol: float = 1e-12,
        rtol: float = 1e-9,
        unordered: bool = False,
        max_reported: int = 10,
    ) -> dict:
        """Compare two mesh files. Returns a verdict (identical / equal within
        tolerance / different), an `equal` boolean, and per-section detail."""
        return _guard(
            _tools.tool_diff,
            path_a=path_a,
            path_b=path_b,
            format_a=format_a,
            format_b=format_b,
            atol=atol,
            rtol=rtol,
            unordered=unordered,
            max_reported=max_reported,
        )


# --------------------------------------------------------------------------- #
# Conversion                                                                  #
# --------------------------------------------------------------------------- #
def _register_conversion(server: FastMCP) -> None:
    @server.tool()
    def convert(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        points_only: bool = False,
        arrays: Optional[List[str]] = None,
        time_step: int = 0,
        mode: str = "auto",
        compression: Optional[str] = None,
    ) -> dict:
        """Convert a mesh between formats (formats inferred from extensions
        unless given). points_only/arrays/time_step narrow the read. mode
        selects ascii|binary output where the format supports it; compression
        selects zlib|lz4|zstd|lzma (VTU/VTP block codecs), gzip (CGNS/H5M/XDMF)
        or 'none' to decompress."""
        return _guard(
            _tools.tool_convert,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            points_only=points_only,
            arrays=arrays,
            time_step=time_step,
            mode=mode,
            compression=compression,
        )


# --------------------------------------------------------------------------- #
# Mesh operations                                                             #
# --------------------------------------------------------------------------- #
def _register_operations(server: FastMCP) -> None:
    @server.tool()
    def extract_surface(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Extract the boundary surface (3D volume -> boundary faces, 2D
        surface -> boundary edges)."""
        return _guard(
            _tools.tool_extract_surface,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            record_parent_ids=record_parent_ids,
        )

    @server.tool()
    def extract_skin(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        linearize: bool = False,
    ) -> dict:
        """Extract the outer skin of a volume mesh (optionally linearized)."""
        return _guard(
            _tools.tool_extract_skin,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            linearize=linearize,
        )

    @server.tool()
    def reorder(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "rcm",
    ) -> dict:
        """Renumber nodes and cells (method: rcm | morton | hilbert). The
        report includes the bandwidth before and after."""
        return _guard(
            _tools.tool_reorder,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            method=method,
        )

    @server.tool()
    def clean(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        weld: bool = False,
        atol: float = 1e-8,
        remove_orphans: bool = True,
        drop_degenerate: bool = True,
        drop_duplicate_cells: bool = True,
    ) -> dict:
        """Clean a mesh: optionally weld coincident points, then drop
        degenerate cells, duplicate cells and orphan points. Reports counts."""
        return _guard(
            _tools.tool_clean,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            weld=weld,
            atol=atol,
            remove_orphans=remove_orphans,
            drop_degenerate=drop_degenerate,
            drop_duplicate_cells=drop_duplicate_cells,
        )

    @server.tool()
    def crop(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        bbox: Optional[List[float]] = None,
        plane_origin: Optional[List[float]] = None,
        plane_normal: Optional[List[float]] = None,
        mode: str = "all",
        record_ids: bool = False,
    ) -> dict:
        """Crop to a bounding box (bbox: [xmin,ymin,zmin,xmax,ymax,zmax]) or a
        half-space (plane_origin + plane_normal, keeps the +normal side). mode
        'all' keeps cells with every node inside, 'any' with any node inside."""
        return _guard(
            _tools.tool_crop,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            bbox=bbox,
            plane_origin=plane_origin,
            plane_normal=plane_normal,
            mode=mode,
            record_ids=record_ids,
        )

    @server.tool(name="slice")
    def slice_mesh(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        origin: Optional[List[float]] = None,
        normal: Optional[List[float]] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Planar cross-section: the true intersection with a plane, one
        topological dimension below the input (default plane: z = 0)."""
        return _guard(
            _tools.tool_slice,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            origin=origin if origin is not None else [0.0, 0.0, 0.0],
            normal=normal if normal is not None else [0.0, 0.0, 1.0],
            record_parent_ids=record_parent_ids,
        )

    @server.tool()
    def isosurface(
        input_path: str,
        output_path: str,
        array: str,
        isovalues: Union[List[float], float],
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        component: Optional[int] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Contour a point_data scalar at one or more isovalues; the output
        carries iso:value / iso:index cell data per contour."""
        return _guard(
            _tools.tool_isosurface,
            input_path=input_path,
            output_path=output_path,
            array=array,
            isovalues=isovalues,
            input_format=input_format,
            output_format=output_format,
            component=component,
            record_parent_ids=record_parent_ids,
        )

    @server.tool()
    def gradient(
        input_path: str,
        output_path: str,
        array: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        operator: str = "gradient",
        method: str = "green-gauss",
        location: str = "cell",
        output: Optional[str] = None,
        component: Optional[int] = None,
        overwrite: bool = False,
    ) -> dict:
        """Gradient, divergence or curl of a point_data field.

        operator: gradient|divergence|curl. method: green-gauss|least-squares.
        location: cell|point. The result is named <array>:<operator> unless
        `output` overrides it; a gradient of an nc-component field has 3*nc
        components laid out [component][derivative]. Divergence and curl need a
        2- or 3-component field. Cells that cannot be differentiated (below the
        mesh dimension, ragged, or a 3-D Lagrange type) yield NaN and are
        reported in num_skipped; least-squares cells with a degenerate
        neighbourhood fall back to Green-Gauss and are reported in
        num_fallback."""
        return _guard(
            _tools.tool_gradient,
            input_path=input_path,
            output_path=output_path,
            array=array,
            input_format=input_format,
            output_format=output_format,
            operator=operator,
            method=method,
            location=location,
            output=output,
            component=component,
            overwrite=overwrite,
        )

    @server.tool()
    def transform(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        translate: Optional[List[float]] = None,
        scale: Optional[Union[List[float], float]] = None,
        rotate_axis: Optional[List[float]] = None,
        rotate_degrees: Optional[float] = None,
        matrix: Optional[List[List[float]]] = None,
        scale_units: Optional[float] = None,
        rotate_vector_data: bool = False,
    ) -> dict:
        """Affine-transform point coordinates: translate [dx,dy,dz], scale
        (scalar or [sx,sy,sz]), rotate (axis + degrees), an explicit 4x4
        matrix, or a uniform unit-scale factor. rotate_vector_data also
        rotates 3-vector/9-tensor point data."""
        return _guard(
            _tools.tool_transform,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            translate=translate,
            scale=scale,
            rotate_axis=rotate_axis,
            rotate_degrees=rotate_degrees,
            matrix=matrix,
            scale_units=scale_units,
            rotate_vector_data=rotate_vector_data,
        )

    @server.tool()
    def convert_cells(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        mode: str = "linearize",
        record_parent_ids: bool = False,
    ) -> dict:
        """Convert the element representation: linearize (drop higher-order
        nodes), simplexify (split into triangles/tets) or elevate (linear ->
        quadratic)."""
        return _guard(
            _tools.tool_convert_cells,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            mode=mode,
            record_parent_ids=record_parent_ids,
        )

    @server.tool()
    def refine(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        levels: int = 1,
        record_parent_ids: bool = False,
        cells: Optional[List[int]] = None,
        region: Optional[str] = None,
        where: Optional[str] = None,
        closure: str = "redgreen",
        record_levels: bool = False,
    ) -> dict:
        """Refine: subdivide cells into congruent same-type children, `levels`
        times. With no selector every cell is refined; give at most one of
        `cells` (global block-major indices), `region` (a cell region selects its
        cells, a point region every cell touching it) or `where` (a threshold on
        a scalar cell_data array, e.g. "quality:scaled_jacobian < 0.3") and only
        those are, with the resulting hanging nodes resolved by `closure` —
        "redgreen" keeps that local, "propagate" reaches the whole connected
        component, and "balanced" keeps the hanging nodes and only enforces 2:1
        balance (the output is then NOT conforming; the constrained nodes are
        reported in refine:hanging). `record_levels` attaches refine:level."""
        return _guard(
            _tools.tool_refine,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            levels=levels,
            record_parent_ids=record_parent_ids,
            cells=cells,
            region=region,
            where=where,
            closure=closure,
            record_levels=record_levels,
        )

    @server.tool()
    def decimate(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        ratio: Optional[float] = None,
        target_faces: Optional[int] = None,
        max_error: Optional[float] = None,
        placement: str = "optimal",
        preserve_boundary: bool = True,
        preserve_features: bool = True,
        feature_angle: float = 30.0,
    ) -> dict:
        """Decimate a surface mesh by quadric edge collapse. Give exactly one
        stopping criterion: ratio (fraction of faces to KEEP), target_faces,
        or max_error. Reports faces/points removed and rejections."""
        return _guard(
            _tools.tool_decimate,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            ratio=ratio,
            target_faces=target_faces,
            max_error=max_error,
            placement=placement,
            preserve_boundary=preserve_boundary,
            preserve_features=preserve_features,
            feature_angle=feature_angle,
        )

    @server.tool()
    def smooth(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "taubin",
        iterations: int = 10,
        lambda_factor: float = -1.0,
        mu: float = -0.34,
        fix_boundary: bool = True,
        preserve_features: bool = True,
        feature_angle: float = 30.0,
        guard_inversion: bool = True,
    ) -> dict:
        """Smooth point coordinates (taubin, the feature-preserving default, or
        laplacian). Pure coordinate move; boundary and feature nodes pinned by
        default. lambda_factor < 0 means the method's own default."""
        return _guard(
            _tools.tool_smooth,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            method=method,
            iterations=iterations,
            lambda_factor=lambda_factor,
            mu=mu,
            fix_boundary=fix_boundary,
            preserve_features=preserve_features,
            feature_angle=feature_angle,
            guard_inversion=guard_inversion,
        )

    @server.tool()
    def merge(
        input_paths: List[str],
        output_path: str,
        output_format: Optional[str] = None,
        weld: bool = False,
        atol: float = 1e-8,
        source_tag: bool = True,
        data_policy: str = "intersection",
        drop_duplicate_cells: bool = False,
    ) -> dict:
        """Merge two or more mesh files into one; weld=true fuses coincident
        points within atol. data_policy: intersection | fill."""
        return _guard(
            _tools.tool_merge,
            input_paths=input_paths,
            output_path=output_path,
            output_format=output_format,
            weld=weld,
            atol=atol,
            source_tag=source_tag,
            data_policy=data_policy,
            drop_duplicate_cells=drop_duplicate_cells,
        )

    @server.tool()
    def split(
        input_path: str,
        input_format: Optional[str] = None,
        by: str = "type",
        tag: Optional[str] = None,
        output_dir: Optional[str] = None,
        name_template: str = "{stem}_{key}.vtu",
    ) -> dict:
        """Split into submeshes (by: type | component | tag | region | regions;
        tag names an integer cell-data array). Writes one file per piece into
        output_dir (default: the input's directory) and returns their paths."""
        return _guard(
            _tools.tool_split,
            input_path=input_path,
            input_format=input_format,
            by=by,
            tag=tag,
            output_dir=output_dir,
            name_template=name_template,
        )

    @server.tool()
    def partition(
        input_path: str,
        nparts: int,
        input_format: Optional[str] = None,
        method: str = "auto",
        imbalance: float = 0.03,
        mode: str = "eco",
        seed: int = 0,
        record_ids: bool = False,
        ghost_layers: int = 0,
        weights: Optional[str] = None,
        output_dir: Optional[str] = None,
        name_template: str = "{stem}_part{part}.vtu",
    ) -> dict:
        """Partition into exactly nparts balanced pieces (method: auto | sfc |
        kahip). Writes one file per part and returns their paths."""
        return _guard(
            _tools.tool_partition,
            input_path=input_path,
            nparts=nparts,
            input_format=input_format,
            method=method,
            imbalance=imbalance,
            mode=mode,
            seed=seed,
            record_ids=record_ids,
            ghost_layers=ghost_layers,
            weights=weights,
            output_dir=output_dir,
            name_template=name_template,
        )

    @server.tool()
    def interpolate(
        source_path: str,
        target_path: str,
        output_path: str,
        source_format: Optional[str] = None,
        target_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "nearest",
        arrays: Optional[List[str]] = None,
        extrapolate: bool = False,
        default_value: float = 0.0,
        on_conflict: str = "error",
    ) -> dict:
        """Sample the source mesh's data arrays onto the target mesh's
        geometry (method: nearest | barycentric); writes the target copy with
        the sampled arrays attached."""
        return _guard(
            _tools.tool_interpolate,
            source_path=source_path,
            target_path=target_path,
            output_path=output_path,
            source_format=source_format,
            target_format=target_format,
            output_format=output_format,
            method=method,
            arrays=arrays,
            extrapolate=extrapolate,
            default_value=default_value,
            on_conflict=on_conflict,
        )


# --------------------------------------------------------------------------- #
# Data operations                                                             #
# --------------------------------------------------------------------------- #
def _register_data(server: FastMCP) -> None:
    @server.tool()
    def data_manage(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        keep: Optional[List[List[str]]] = None,
        drop: Optional[List[List[str]]] = None,
        rename: Optional[List[List[str]]] = None,
        ignore_missing: bool = False,
    ) -> dict:
        """Keep/drop/rename data arrays. keep/drop: [[location, name], ...];
        rename: [[location, old, new], ...] (location: point | cell | field).
        Applied keep -> drop -> rename."""
        return _guard(
            _tools.tool_data_manage,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            keep=keep,
            drop=drop,
            rename=rename,
            ignore_missing=ignore_missing,
        )

    @server.tool()
    def data_convert(
        input_path: str,
        output_path: str,
        direction: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        arrays: Optional[List[str]] = None,
        weighted: bool = False,
        nan_policy: str = "ignore",
        nan_replacement: float = 0.0,
    ) -> dict:
        """Average data arrays between locations (direction: point_to_cell |
        cell_to_point; weighted uses cell measures for cell_to_point)."""
        return _guard(
            _tools.tool_data_convert,
            input_path=input_path,
            output_path=output_path,
            direction=direction,
            input_format=input_format,
            output_format=output_format,
            arrays=arrays,
            weighted=weighted,
            nan_policy=nan_policy,
            nan_replacement=nan_replacement,
        )

    @server.tool()
    def data_calc(
        input_path: str,
        output_path: str,
        expression: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        location: str = "point",
        output_name: str = "",
        overwrite: bool = False,
    ) -> dict:
        """Evaluate an arithmetic expression over data arrays into a new array
        (operators + - * /, functions abs/sqrt/min/max/norm; back-tick-quote
        names with spaces, e.g. "speed = norm(velocity)")."""
        return _guard(
            _tools.tool_data_calc,
            input_path=input_path,
            output_path=output_path,
            expression=expression,
            input_format=input_format,
            output_format=output_format,
            location=location,
            output_name=output_name,
            overwrite=overwrite,
        )

    @server.tool()
    def data_condition(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        operation: str = "clamp",
        location: str = "point",
        arrays: Optional[List[str]] = None,
        scope: str = "component",
        lo: float = 0.0,
        hi: float = 1.0,
        nan_policy: str = "ignore",
        nan_replacement: float = 0.0,
        suffix: str = "",
        preserve_dtype: bool = True,
    ) -> dict:
        """Condition data values: clamp to [lo, hi], normalize into [lo, hi],
        or standardize (zero mean, unit std). scope: component | magnitude."""
        return _guard(
            _tools.tool_data_condition,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            operation=operation,
            location=location,
            arrays=arrays,
            scope=scope,
            lo=lo,
            hi=hi,
            nan_policy=nan_policy,
            nan_replacement=nan_replacement,
            suffix=suffix,
            preserve_dtype=preserve_dtype,
        )


# --------------------------------------------------------------------------- #
# Gated tools (optional extras)                                               #
# --------------------------------------------------------------------------- #
def _register_gated(server: FastMCP) -> None:
    @server.tool()
    def data_export(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        location: str = "point",
    ) -> dict:
        """Export data arrays to a Parquet table (location: point | cell).
        Needs the [arrow] extra; a missing install returns a named error."""
        return _guard(
            _tools.tool_data_export,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            location=location,
        )

    @server.tool()
    def screenshot(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        color_by: Optional[str] = None,
        width: int = 1280,
        height: int = 960,
        transparent: bool = False,
    ):
        """Render an off-screen PNG screenshot of the mesh (optionally colored
        by a data array). Needs the [viewer] extra (polyscope); a missing
        install returns a named error."""
        try:
            report = _tools.tool_screenshot(
                input_path=input_path,
                output_path=output_path,
                input_format=input_format,
                color_by=color_by,
                width=width,
                height=height,
                transparent=transparent,
            )
        except Exception as e:  # noqa: BLE001
            return {"error": str(e), "error_type": type(e).__name__}
        return [Image(path=report["output_path"]), report]


# --------------------------------------------------------------------------- #
# Resources                                                                   #
# --------------------------------------------------------------------------- #
def _register_resources(server: FastMCP) -> None:
    @server.resource("meshioplusplus://formats")
    def formats_resource() -> str:
        """The supported-formats registry as JSON."""
        return json.dumps(_tools.formats_payload(), indent=2)

    @server.resource("meshioplusplus://version")
    def version_resource() -> str:
        """The installed meshio++ version."""
        return __version__


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="meshioplusplus-mcp",
        description="meshio++ MCP server (stdio transport).",
    )
    parser.add_argument(
        "--root",
        default=os.environ.get("MESHIOPLUSPLUS_MCP_ROOT") or None,
        help=(
            "confine every tool's input/output paths to this directory "
            "(default: the MESHIOPLUSPLUS_MCP_ROOT env var, else unrestricted)"
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"meshioplusplus {__version__}"
    )
    args = parser.parse_args(argv)
    create_server(root=args.root).run()
    return 0
