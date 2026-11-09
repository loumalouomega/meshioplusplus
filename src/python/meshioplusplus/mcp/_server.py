"""MCP registration layer: typed wrappers over :mod:`._tools`.

This is the only module in the package that imports the ``mcp`` SDK, so it is
only importable with the ``[mcp]`` extra installed (which needs Python >= 3.10
— the SDK's own floor, not meshio++'s). All behaviour lives in the pure layer;
each wrapper here contributes exactly two things: the typed signature the SDK
derives the tool's JSON schema from, and the docstring the client shows the
model. Wrappers never leak tracebacks — every call is guarded and failures
come back as ``{"error", "error_type"}`` payloads the model can act on.
"""

from __future__ import annotations

import argparse
import functools
import inspect
import json
import os
from typing import Dict, List, Optional, Union

try:  # mcp >= 2 renamed FastMCP to MCPServer (and moved it)
    from mcp.server.mcpserver import Image
    from mcp.server.mcpserver import MCPServer as FastMCP

    _HOST_IN_CTOR = False
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP, Image

    _HOST_IN_CTOR = True

from ..__about__ import __version__
from . import _tools

_INSTRUCTIONS = """\
meshio++ mesh I/O and processing. All tools are file-path based:
they read mesh files (40+ formats: VTK/VTU/VTP, Gmsh, Exodus, MED, Abaqus,
STL, OBJ, PLY, XDMF, CGNS, ...), operate, and write output files, returning a
JSON report. Start with `formats` to see supported formats, `info` for a fast
file summary, and `convert` to translate between formats. Inspection tools
(`stats`, `quality`, `data_info`, `regions`, `diff`) return reports without
writing anything. Mesh operations (clean, refine, decimate, smooth, slice,
partition, ...) take input_path/output_path. Paths may be confined to a root
directory by the server's --root option. Parsed inputs are cached between
calls and re-read whenever the file changes, so asking several questions of
one file is cheap.
"""


# The one guard, shared with the HTTP front-end's dispatch: a tool failure is
# a {"error", "error_type"} payload, never a traceback (or a 500).
_guard = _tools.guard


def _forward(target):
    """Make a typed tool wrapper forward its arguments to ``target``, guarded.

    The decorated function is the tool's declaration: its annotated signature
    is what the SDK turns into the tool's input schema, and its docstring the
    description, so both stay written out; only the body, ``return
    _guard(target, a=a, b=b, ...)`` with every parameter by name, is generated
    (roadmap §7). ``functools.wraps`` keeps the signature visible through
    ``__wrapped__`` (mcp 2.x) and this module's globals for the string
    annotations (mcp 1.x evaluates them against the wrapper's).
    """

    def decorate(fn):
        signature = inspect.signature(fn)

        @functools.wraps(fn)
        def forward(*args, **kwargs):
            bound = signature.bind(*args, **kwargs)
            bound.apply_defaults()
            return _guard(target, **bound.arguments)

        return forward

    return decorate


def create_server(root: Optional[str] = None, host: Optional[str] = None) -> FastMCP:
    """Build the MCP server with every tool and resource registered.

    ``host`` matters only for the HTTP front-end: the SDK enables its
    DNS-rebinding protection (Host/Origin allow-lists) when the host is a
    loopback name, so a ``--host 0.0.0.0`` server must be told, or the
    ``/mcp`` endpoint would reject every non-loopback ``Host`` header. mcp 1.x
    takes it in the constructor; 2.x takes it when the HTTP app is built, so
    it is remembered on the server for :func:`http_app`.
    """
    if root is not None:
        _tools.set_root(root)
    kwargs = {"host": host} if host is not None and _HOST_IN_CTOR else {}
    server = FastMCP("meshioplusplus", instructions=_INSTRUCTIONS, **kwargs)
    server._meshioplusplus_host = host
    _register_inspection(server)
    _register_conversion(server)
    _register_operations(server)
    _register_data(server)
    _register_dataset(server)
    _register_training(server)
    _register_gated(server)
    _register_resources(server)
    return server


def http_app(server: FastMCP):
    """The server's ASGI app and its transport name (``streamable-http``, or
    ``sse`` on an SDK older than streamable HTTP)."""
    host = None if _HOST_IN_CTOR else getattr(server, "_meshioplusplus_host", None)
    kwargs = {"host": host} if host is not None else {}
    if hasattr(server, "streamable_http_app"):
        return server.streamable_http_app(**kwargs), "streamable-http"
    return server.sse_app(**kwargs), "sse"  # pragma: no cover


# --------------------------------------------------------------------------- #
# Inspection                                                                  #
# --------------------------------------------------------------------------- #
def _register_inspection(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_formats)
    def formats() -> dict:
        """List every readable and writable mesh format and the extension map."""

    @server.tool()
    @_forward(_tools.tool_sniff)
    def sniff(input_path: str) -> dict:
        """Identify a mesh file's format from its leading bytes (a directory: the
        files it holds) and its extension."""

    @server.tool()
    @_forward(_tools.tool_info)
    def info(input_path: str, file_format: Optional[str] = None) -> dict:
        """Fast file summary without loading heavy arrays: point/cell counts,
        cell blocks, data-array names, named regions, time steps, format."""

    @server.tool()
    @_forward(_tools.tool_stats)
    def stats(
        input_path: str,
        file_format: Optional[str] = None,
        time_step: int = 0,
        piece: Optional[int] = None,
        ghosts: str = "keep",
    ) -> dict:
        """Geometric statistics: bounding box, centroid, total area, signed and
        unsigned volume, per-type cell counts, inverted-cell count. ghosts is
        keep|drop for the halo of a .pvtu/.pvtp/.pvd."""

    @server.tool()
    @_forward(_tools.tool_quality)
    def quality(
        input_path: str,
        file_format: Optional[str] = None,
        output_path: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Mesh-quality metrics (scaled jacobian, aspect ratio, skewness,
        angles, ...) as per-metric min/max/mean + histograms. Pass output_path
        to also write a mesh copy carrying quality:<metric> cell data."""

    @server.tool()
    @_forward(_tools.tool_data_info)
    def data_info(input_path: str, file_format: Optional[str] = None) -> dict:
        """Describe every point/cell/field data array: dtype, shape, components,
        min/max/mean (whole-array and per component), NaN/Inf counts."""

    @server.tool()
    @_forward(_tools.tool_data_integrate)
    def data_integrate(
        input_path: str,
        file_format: Optional[str] = None,
        arrays: Optional[List[str]] = None,
    ) -> dict:
        """Cell-measure-weighted total/mean of cell_data arrays (whole mesh,
        and independently for every named Cell region) -- gradient's
        integration companion. arrays defaults to every cell_data array."""

    @server.tool()
    @_forward(_tools.tool_regions)
    def regions(input_path: str, file_format: Optional[str] = None) -> dict:
        """List named regions (point/cell/side groups: gmsh physical groups,
        Exodus sets, Abaqus NSET/ELSET, ...) with kind, dim, tag and size."""

    @server.tool()
    @_forward(_tools.tool_bandwidth)
    def bandwidth(input_path: str, file_format: Optional[str] = None) -> dict:
        """Node-numbering bandwidth (max node-index spread over any cell)."""

    @server.tool()
    @_forward(_tools.tool_data_preview)
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

    @server.tool()
    @_forward(_tools.tool_diff)
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


# --------------------------------------------------------------------------- #
# Conversion                                                                  #
# --------------------------------------------------------------------------- #
def _register_conversion(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_med_multi)
    def med_multi(
        input_path: Optional[str] = None,
        output_path: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        mesh_name: Optional[str] = None,
        mesh_names: Optional[List[str]] = None,
        output_format: Optional[str] = None,
        med_version: str = "4.1.0",
    ) -> dict:
        """List named meshes in input_path, or extract mesh_name to output_path.
        Alternatively combine nonempty input_paths into a MED output_path with
        optional mesh_names and med_version (default 4.1.0). Input paths and
        output paths obey the server's root sandbox. Ordinary MED nodal/element
        profiles expand with NaN fill; enhanced fields may use Python/h5py."""

    @server.tool()
    @_forward(_tools.tool_convert)
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
        piece: Optional[int] = None,
        ghosts: str = "keep",
        grid_functions: Optional[Dict[str, str]] = None,
        write_grid_functions: bool = False,
        elmer_halo: bool = False,
        patran_results: Optional[Dict[str, str]] = None,
        write_nurbs: Optional[bool] = None,
        heal: bool = False,
        sew_tolerance: float = 0.0,
    ) -> dict:
        """Convert a mesh between formats (formats inferred from extensions
        unless given). grid_functions ({name: path}) reads MFEM .gf fields onto
        an MFEM .mesh input; patran_results ({name: path}) reads Patran 2.5
        result files (.nod/.dis nodal, .els element) onto a Patran neutral
        input; write_grid_functions writes the data of an MFEM
        output as <stem>.<name>.gf files beside it; elmer_halo adds ElmerGrid's
        -halo layer to the partitioning an Elmer output writes from
        partition:part. points_only/arrays/time_step narrow the read; piece
        keeps one partition/block of a partitioned file (VTKHDF, .pvtu/.pvtp,
        or one part of a .pvd step) instead of the merged mesh; ghosts=drop
        removes the ghost cells (halo) of a .pvtu/.pvtp/.pvd. mode selects ascii|binary output where the format
        supports it, or raw_appended (VTU: raw binary in one <AppendedData>
        section, no base64); compression selects zlib|lz4|zstd|lzma (VTU/VTP block
        codecs), gzip (CGNS/H5M/VTKHDF/XDMF), lzf (PCD binary_compressed) or
        'none' to decompress. MDPA-to-MDPA preserves nested geometry/constraint
        membership as raw file ids; constraints remain opaque, and other
        output formats may drop that side-channel content. Gmsh input versions
        2.2/4.0/4.1, including periodic links, use the native reader. Gmsh output
        is 4.1, or 2.2 via gmsh22; both preserve periodic links. Serial VTK XML
        inputs (VTU/VTP/VTS/VTR/VTI) accept appended raw/base64 arrays natively;
        Multiple pieces concatenate without welding; VTP/VTS/VTR/VTI retain
        inline output. Legacy VTK structured datasets also read natively.
        XDMF time-series reads preserve fixed shared point/cell/side regions.
        An MFEM NURBS input keeps its patches, so an MFEM output is written back
        as the NURBS mesh; write_nurbs=False writes the sampled cells instead,
        True refuses anything else. heal (STEP/IGES) heals the shape before it
        is read (ShapeFix_Shape, then sewing at sew_tolerance, a length; 0 is
        1e-5 of the bounding-box diagonal, and a positive value sews without
        heal)."""

    @server.tool()
    @_forward(_tools.tool_export_gltf)
    def export_gltf(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        split_angle: float = 30.0,
        normal_weight: str = "angle",
        normals: bool = True,
        fields: bool = True,
        color_by: Optional[str] = None,
        component: Optional[int] = None,
        cmap: str = "viridis",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        nan_color: str = "#808080",
        unlit: bool = True,
        up_axis: str = "auto",
        recenter: bool = True,
        scale: float = 1.0,
        by_region: bool = True,
        container: str = "auto",
    ) -> dict:
        """Export the surface of a mesh as glTF 2.0 -- a web-native format
        three.js, Blender, model viewers and dashboards load directly. A .glb
        output_path writes the binary container; .gltf writes JSON with a .bin
        beside it (container overrides the suffix).

        The skin of volume cells, 2-D cells, line cells (LINES) and vertex
        cells or a cell-less point cloud (POINTS) are exported, one named node
        per cell region. glTF normals are per vertex, so points are duplicated
        where the surface creases by more than split_angle degrees (0-180;
        normal_weight angle|area; normals=false writes none). Every one-to-four
        component point_data array is exported raw as a custom attribute
        (temperature -> _TEMPERATURE); fields=false skips them. color_by names
        a point_data or cell_data array to bake into COLOR_0 through cmap
        (viridis, coolwarm, turbo, magma, inferno, plasma, grey, synthwave,
        or a reversed *_r variant) over vmin..vmax (default: the finite range of
        what is exported), with an unlit material unless unlit=false;
        multi-component arrays reduce to component or their magnitude and
        non-finite values take nan_color (#rrggbb). The output is Y-up,
        right-handed, float32: up_axis (auto|x|y|z) names the source axis that
        points up (auto is y for a flat mesh, else z), and that rotation, scale
        (source unit to metres) and the bounding-box centre recenter subtracts
        all live on the root node, not in the coordinates, so the transform is
        exactly reversible, and glTF files read back (any glTF 2.0 or GLB
        file is a valid input to every other tool). A valid glTF-Validator pass with zero
        errors and warnings is part of the test suite."""

    @server.tool()
    @_forward(_tools.tool_pipeline)
    def pipeline(
        settings_path: str,
        input_path: Optional[str] = None,
        output_path: Optional[str] = None,
    ) -> dict:
        """Run a settings.json operation pipeline: read Input.Path, apply the
        Operations chain (Transform/Gradient/Refine/Clean/... -- PascalCase
        ops and keys, see doc/pipeline.md), write Output.Path. Version 2 adds
        Merge/Interpolate/UndoGreen Inputs and terminal Split/Partition outputs
        through Output.Pattern ({key}/{part}); a CAD JSON/STEP/IGES/MFEM NURBS
        Input carries its model to MeshCadSurface/IgaSample/IgaRefine steps
        until another step changes the mesh; all auxiliary inputs and every
        expanded output, including symlinks, are checked against the sandbox.
        Output.Codec='lzf' selects PCD compressed binary; ASCII is incompatible.
        input_path/output_path override the paths in the settings file; both
        the settings file and the paths inside it stay inside the sandbox
        root when one is configured."""

    @server.tool()
    def sequence(
        output_path: str,
        input_pattern: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        mode: Optional[str] = None,
        times: Optional[List[float]] = None,
        time_from: str = "auto",
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Treat a set of files (or the steps inside one multi-step file) as one
        transient dataset. Give exactly one of input_pattern (a glob -- '*' and
        '?' only) or input_paths. A '{step}'/'{index}' token in output_path
        writes one file per step (fan-out); a plain path writes one multi-step
        file (fan-in, only for xdmf, gid, usd, vtkhdf, pvd, femap and frd -- any
        other format fails by name rather than silently keeping step 0).
        XDMF fan-in stores the first mesh's fixed point/cell/side regions;
        later steps do not change their membership.
        Ordering is natural-numeric, so out_9 precedes
        out_10. mode optionally asserts 'sequence'/'fan-in'/'fan-out'.
        See doc/sequences.md."""
        return _guard(
            _tools.tool_sequence,
            input_pattern=input_pattern,
            input_paths=input_paths,
            output_path=output_path,
            mode=mode,
            times=times,
            time_from=time_from,
            input_format=input_format,
            output_format=output_format,
        )

    @server.tool()
    @_forward(_tools.tool_resample_sequence)
    def resample_sequence(
        output_path: str,
        input_pattern: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        times: Optional[List[float]] = None,
        times_from_pattern: Optional[str] = None,
        method: str = "linear",
        clamp: bool = False,
        blend_points: bool = False,
        time_from: str = "auto",
        input_format: Optional[str] = None,
    ) -> dict:
        """Resample a transient sequence onto new times -- align two solvers'
        timelines. Give exactly one of input_pattern (a glob) or input_paths,
        and exactly one of times or times_from_pattern (a second sequence whose
        step times are the targets). method is 'linear' (blend the bracketing
        steps' float data), 'nearest' or 'previous'; clamp takes the end steps
        for targets outside the source range instead of failing; blend_points
        also blends the coordinates. A '{step}'/'{index}' output_path writes
        one file per target time, a plain path one series file. At most two
        source meshes are held at once. See doc/sequences.md."""

    @server.tool()
    @_forward(_tools.tool_blend_steps)
    def blend_steps(
        path_a: str,
        path_b: str,
        output_path: str,
        weight: float,
        blend_points: bool = False,
        format_a: Optional[str] = None,
        format_b: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Blend two steps of one topology: every floating-point data array
        becomes (1 - weight) a + weight b, integer data comes from the nearer
        step, and points/cells stay a's unless blend_points. The steps must
        have the same point count and cell blocks. See doc/sequences.md."""


# --------------------------------------------------------------------------- #
# Mesh operations                                                             #
# --------------------------------------------------------------------------- #
def _register_operations(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_extract_surface)
    def extract_surface(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Extract the boundary surface (3D volume -> boundary faces, 2D
        surface -> boundary edges)."""

    @server.tool()
    @_forward(_tools.tool_extract_skin)
    def extract_skin(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        linearize: bool = False,
    ) -> dict:
        """Extract the outer skin of a volume mesh (optionally linearized)."""

    @server.tool()
    @_forward(_tools.tool_reorder)
    def reorder(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "rcm",
    ) -> dict:
        """Renumber nodes and cells (method: rcm | morton | hilbert). The
        report includes the bandwidth before and after."""

    @server.tool()
    @_forward(_tools.tool_clean)
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

    @server.tool()
    @_forward(_tools.tool_crop)
    def crop(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        bbox: Optional[List[float]] = None,
        plane_origin: Optional[List[float]] = None,
        plane_normal: Optional[List[float]] = None,
        where_array: Optional[str] = None,
        where_compare: str = "<",
        where_value: float = 0.0,
        mode: str = "all",
        record_ids: bool = False,
    ) -> dict:
        """Crop to a bounding box (bbox: [xmin,ymin,zmin,xmax,ymax,zmax]), a
        half-space (plane_origin + plane_normal, keeps the +normal side), or a
        cell_data predicate (where_array + where_compare in <, <=, >, >=, ==, !=
        + where_value). Give exactly one. A non-finite cell value never matches
        the predicate. mode 'all' keeps cells with every node inside, 'any' with
        any node inside -- it applies to bbox/plane only, since a cell_data value
        is already one per cell. Inside/outside a surface composes: run
        distance_to_surface with location='center', then crop on 'sdf:distance'
        < 0."""

    @server.tool(name="slice")
    def slice_mesh(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        origin: Optional[List[float]] = None,
        normal: Optional[List[float]] = None,
        record_parent_ids: bool = False,
        chain: bool = False,
        min_length: float = 0.0,
    ) -> dict:
        """Planar cross-section: the true intersection with a plane, one
        topological dimension below the input (default plane: z = 0). chain
        walks a surface's line section into polylines (section:chain and
        section:closed cell data), dropping chains shorter than min_length."""
        return _guard(
            _tools.tool_slice,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            origin=origin if origin is not None else [0.0, 0.0, 0.0],
            normal=normal if normal is not None else [0.0, 0.0, 1.0],
            record_parent_ids=record_parent_ids,
            chain=chain,
            min_length=min_length,
        )

    @server.tool()
    @_forward(_tools.tool_remove_islands)
    def remove_islands(
        input_path: str,
        output_path: str,
        connectivity: str = "edge",
        min_diagonal: float = 0.01,
        min_area: float = 0.0,
        min_cells: int = 0,
        relative: bool = True,
        dry_run: bool = False,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Drop a scan's crumbs: the connected components (by shared edges,
        or points) whose bounding-box diagonal, area or cell count is below a
        threshold, relative to the whole mesh by default (min_diagonal 0.01).
        dry_run labels islands:component and islands:removed instead. Reports
        num_components, num_removed, num_cells_removed, num_points_removed
        and per component its cells, diagonal, area and whether it goes."""

    @server.tool()
    @_forward(_tools.tool_flip_normals)
    def flip_normals(
        input_path: str,
        output_path: str,
        region: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Reverse the winding of every cell, or of a cell region's, with no
        point moved; Side regions move with their cells. Reports num_flipped."""

    @server.tool()
    @_forward(_tools.tool_clip)
    def clip(
        input_path: str,
        output_path: str,
        origin: List[float],
        normal: List[float],
        keep: str = "below",
        cap: bool = True,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Cut a surface by the plane through origin normal to normal, keep
        the side below it (or "above"), split the triangles across it at the
        crossings slice computes, and cap the cut's closed loops so a closed
        part stays closed. Reports num_kept, num_split, num_cap_triangles,
        num_loops, num_open_chains and warnings."""

    @server.tool()
    @_forward(_tools.tool_profile_solid)
    def profile_solid(
        input_path: str,
        output_path: str,
        features: Union[str, List[dict]],
        step_path: Optional[str] = None,
        scan_path: Optional[str] = None,
        provider: str = "",
        deflection: float = 0.0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Build a solid from a fit_curves sketch (input_path) by features
        applied in order, through a CAD provider (OpenCASCADE): first an
        extrusion ({"kind": "extrude", "distance": d}, 0 runs to scan_path's
        extent, "symmetric" half each way) or a revolution ({"kind":
        "revolve", "axis_origin", "axis_direction", "angle"}), then fillets
        ({"kind": "fillet", "radius"}), chamfers ({"kind": "chamfer",
        "distance"}) or a shell ({"kind": "shell", "distance": wall}) on
        every edge; features may also be text such as
        "extrude:10;fillet:0.5". Writes the tessellation (output_path) and
        STEP (step_path). Reports volume, area, num_faces, valid, provider,
        providers and warnings."""

    @server.tool()
    @_forward(_tools.tool_constrain_primitives)
    def constrain_primitives(
        input_path: str,
        output_path: str,
        face_tolerance: float = 0.0,
        angle_tolerance: float = -1.0,
        distance_tolerance: float = -1.0,
        min_fraction: float = -1.0,
        force: Optional[List[List[Union[int, str]]]] = None,
        forbid: Optional[List[List[Union[int, str]]]] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Fit related primitives of a fit_primitives skin together: coaxial
        groups, parallel and perpendicular planes, equal radii, as
        fit_cad's constrain_primitives does, with force and forbid lists of
        [patch_a, patch_b, relation] (coaxial, parallel, perpendicular,
        equal_radius; patches by segment:patch label) and angle (degrees)
        and distance tolerances (negative keeps the pass's own). Reports
        num_constrained, num_groups, face_tolerance, relations (each
        patch_a, patch_b, relation, angle_deg, distance, forced) and unmet."""

    @server.tool()
    @_forward(_tools.tool_face_support)
    def face_support(
        solid_path: str,
        scan_path: str,
        output_path: Optional[str] = None,
        tolerance: float = 0.0,
        measured_share: float = 0.9,
        inferred_share: float = 0.1,
        hole_data: str = "repair:hole",
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Class each face of a fit_cad re-mesh (solid_path) by how much of
        it the scan supports: measured (at least measured_share of its area
        within tolerance of the scan), inferred (below inferred_share: no
        scan near) or partial; scan cells marked by hole_data (repair's
        fills) count as unmeasured. Reports face_share, face_class,
        num_measured, num_partial, num_inferred, share_within and
        tolerance; output_path writes the solid with fit_cad:face_support."""

    @server.tool()
    @_forward(_tools.tool_fit_cad_iterate)
    def fit_cad_iterate(
        input_path: str,
        output_path: str,
        step_path: Optional[str] = None,
        part_class: str = "any",
        target_within: float = 0.0,
        target_explained: float = 0.0,
        segment_tolerance: float = 0.0,
        fit_tolerance: float = 0.0,
        max_rounds: int = 4,
        worst_faces: int = 3,
        refine_factor: float = 0.25,
        backend: str = "auto",
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Reverse-engineer a closed scan in rounds: segment, fit_primitives
        and fit_cad, stopping when the share of the scan within tolerance of
        the solid reaches target_within (and the explained area
        target_explained), at no improvement or after max_rounds, else
        re-segmenting the worst faces' patches finer. part_class (analytic,
        cast, freeform) sets the target (0.90, 0.85, 0.75) and tolerances.
        Writes the best solid (output_path, and step_path as STEP). Reports
        rounds, best_round, target_met, stop_reason and the face support."""

    @server.tool()
    @_forward(_tools.tool_offset_surface)
    def offset_surface(
        input_path: str,
        output_path: str,
        distance: float,
        smooth_iterations: int = 0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Offset a triangle and quad surface: each point moves along its
        offset direction (boundary_layer's) by distance over the least cosine
        to its faces, so every face's plane moves by distance (positive along
        the faces' normals). Keeps the connectivity, data and regions.
        Reports num_not_visible, num_clamped, num_folded, max_move and
        warnings."""

    @server.tool()
    @_forward(_tools.tool_thicken)
    def thicken(
        input_path: str,
        output_path: str,
        thickness: float,
        side: str = "inward",
        smooth_iterations: int = 0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Turn a surface into a closed shell thickness thick: side
        "inward" (the surface and its offset by -thickness), "outward" or
        "both" (+-thickness/2); the lower layer is reversed and the boundary
        edges are joined by quads (thicken:side cell data 0 upper, 1 lower, 2
        rim). Reports num_rim_quads, num_folded, max_move and warnings."""

    @server.tool()
    @_forward(_tools.tool_hollow)
    def hollow(
        input_path: str,
        output_path: str,
        wall: float,
        smooth_iterations: int = 0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Give a closed surface an inner wall: the surface plus its offset
        by -wall, reversed (an open surface is refused). Reports num_folded,
        max_move and warnings."""

    @server.tool()
    @_forward(_tools.tool_measure)
    def measure(
        input_path: str,
        spec: Union[dict, str],
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Measure a scan by a spec (a dict, its JSON text, or a JSON file's
        path): elements fitted ("fit": plane, sphere, cylinder, cone or auto,
        on a cell "region", from "seeds" [[x, y, z]], or a segment "patch"),
        constructed ("construct": intersection, axis, centre, apex,
        mid_plane, projection, plane_through; "of": [names]) or literal
        ("point"); dimensions (distance, angle between two elements, with
        nominal and tolerance) and tolerances (flatness, cylindricity,
        sphericity, circularity, parallelism, perpendicularity, angularity
        with basic_angle, coaxiality, concentricity; feature, datum, limit),
        and the form criterion. Reports elements, dimensions, tolerances
        (each value and pass), num_fail and warnings; output_path writes the
        mesh with its element table (measure:table) and measure:<name>
        regions."""

    @server.tool()
    @_forward(_tools.tool_align_datums)
    def align_datums(
        input_path: str,
        output_path: str,
        primary: str,
        secondary: str,
        origin: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """The 3-2-1 datum alignment of a measured mesh (measure): the
        primary element (a plane, or an axis) sets Z, the secondary X, and
        origin (a point-like element, or a third plane meeting the other two)
        the origin. Writes the moved mesh, its element table with it, and
        reports the 4x4 matrix transform applies."""

    @server.tool()
    @_forward(_tools.tool_write_step_elements)
    def write_step_elements(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        patches: Optional[List[int]] = None,
        segments: int = 64,
    ) -> dict:
        """Write the fitted primitives of a mesh (segment, then
        fit_primitives) to a STEP file as exact geometry a CAD user can snap
        to: each plane as a face bounded by its points' convex hull, a sphere
        as a ball, a cylinder as a solid with flat lids over its measured
        length, a cone as its frustum, a torus whole. patches selects patch
        labels (default every fitted one). Reports num_planes, num_spheres,
        num_cylinders, num_cones, num_tori, num_skipped, warnings and
        output_path."""

    @server.tool()
    @_forward(_tools.tool_silhouette)
    def silhouette(
        input_path: str,
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        normal: Optional[List[float]] = None,
        clearance: float = 0.0,
        round_angle: float = 180.0,
        arc_segments: int = 64,
        extrude: bool = False,
        region: Optional[str] = None,
        merge_tolerance: float = 1e-9,
    ) -> dict:
        """The outline of a surface (or a volume mesh's skin) seen along normal
        (default 0,0,1): projected and unioned exactly, as closed loops with
        its holes, in the plane at the part's lowest point. clearance > 0
        grows it (a loose fit), < 0 shrinks it (a press fit), corners turning
        by less than round_angle rounded with arcs; extrude returns the closed
        surface of the prism over it instead (a nest or a jig). Reports
        num_loops, num_outer, num_holes, num_triangles, area, perimeter,
        min_height, max_height and warnings; output_path writes the chained
        line mesh (ready for fit_curves) or the prism."""

    @server.tool()
    @_forward(_tools.tool_fit_curves)
    def fit_curves(
        input_path: str,
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        tolerance: float = 0.0,
        snap_angle: float = 0.0,
        normal: Optional[List[float]] = None,
        min_arc_points: int = 4,
        min_arc_sweep: float = 5.0,
    ) -> dict:
        """Lines, arcs and circles from a chained line mesh (slice with chain,
        or sections): each chain split into the runs a line or a circular arc
        explains within tolerance (0: 1e-3 of the diagonal), in its section
        plane's 2-D frame; snap_angle snaps lines onto the axes or diagonals.
        Reports pieces (type, start, end, centre, radius, sweep, residual,
        num_points, chain, frame), frames (4x4: e1, e2, normal, origin
        columns), num_lines, num_arcs, num_circles, max_residual, tolerance
        and warnings; output_path writes the mesh with sketch:piece and the
        sketch:* field data."""

    @server.tool()
    @_forward(_tools.tool_sections)
    def sections(
        input_path: str,
        output_path: str,
        offsets: List[float],
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        direction: Optional[List[float]] = None,
        min_length: float = 0.0,
        record_parent_ids: bool = False,
    ) -> dict:
        """Chained line sections of a surface at the parallel planes
        n . x = offsets[k] (n the unit direction, default +z), merged into one
        line mesh with section:offset, section:chain and section:closed cell
        data; the drawing profiles of a scanned part."""

    @server.tool()
    @_forward(_tools.tool_isosurface)
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

    @server.tool()
    @_forward(_tools.tool_grid)
    def grid(
        output_path: str,
        dims: List[int],
        output_format: Optional[str] = None,
        origin: Optional[List[float]] = None,
        spacing: Optional[List[float]] = None,
        max_cells: int = 20000000,
    ) -> dict:
        """Generate a regular hexahedron lattice of dims=(nx, ny, nz) cells. The
        only tool that reads no input mesh: it creates one. Useful as a background
        grid, a sampling domain, or a fixture that needs no file."""

    @server.tool()
    @_forward(_tools.tool_rectangle)
    def rectangle(
        output_path: str,
        output_format: Optional[str] = None,
        dims: Optional[List[int]] = None,
        size: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        cell_type: str = "quad",
        split: str = "zigzag",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a regular 2-D grid of dims=(nx, ny) quads or triangles at
        z = origin[2]. split (triangles) is zigzag, up, down or center; order=2
        gives quad8/triangle6 cells. The mesh carries side regions xmin, xmax,
        ymin, ymax. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_box)
    def box(
        output_path: str,
        output_format: Optional[str] = None,
        dims: Optional[List[int]] = None,
        size: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        cell_type: str = "hexahedron",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a regular 3-D grid of dims=(nx, ny, nz) hexahedra, tetrahedra
        (six Kuhn tetrahedra per cell) or wedges. order=2 gives hexahedron20,
        tetra10 or wedge15 cells. The mesh carries side regions xmin, xmax, ymin,
        ymax, zmin, zmax. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_disk)
    def disk(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        n: int = 8,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a disk of quads (an n x n grid pushed onto the disk) or
        triangles (n rings, 6 n^2 triangles). order=2 curves the rim. The mesh
        carries a side region rim. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_uv_sphere)
    def uv_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        nu: int = 16,
        nv: int = 8,
        center: Optional[List[float]] = None,
        cell_type: str = "triangle",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a latitude/longitude sphere surface with nu points per ring and
        nv bands, oriented outward. cell_type is triangle, or quad with triangle
        caps. order=2 puts the mid-edge nodes on the sphere. Reads no input mesh:
        it creates one."""

    @server.tool()
    @_forward(_tools.tool_icosa_sphere)
    def icosa_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        subdivisions: int = 2,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a sphere surface from an icosahedron split in four
        `subdivisions` times (20 * 4^subdivisions triangles), oriented outward.
        order=2 puts the mid-edge nodes on the sphere. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_cylinder)
    def cylinder(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        height: float = 1.0,
        nxy: int = 8,
        nz: int = 8,
        center: Optional[List[float]] = None,
        cell_type: str = "hexahedron",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a solid cylinder along z: hexahedra (the disk grid, extruded)
        or wedges. order=2 curves the side. The mesh carries side regions side,
        top and bottom. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_reference_triangle)
    def reference_triangle(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 4,
        size: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate the reference triangle (0,0) (1,0) (0,1), scaled by size and
        refined into n*n triangles in barycentric lattice order. The mesh carries
        side regions ymin, xmin and hypotenuse. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_ngon)
    def ngon(
        output_path: str,
        output_format: Optional[str] = None,
        sides: int = 6,
        radius: float = 1.0,
        rings: int = 1,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a regular polygon of `sides` sides fanned from its centre
        (rings=1), or refined into `rings` rings of triangles. The mesh carries a
        side region rim. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_annulus)
    def annulus(
        output_path: str,
        output_format: Optional[str] = None,
        inner_radius: float = 0.5,
        outer_radius: float = 1.0,
        ntheta: int = 16,
        nr: int = 4,
        angle: float = 6.283185307179586,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate an annulus of ntheta x nr quads or triangles, or a sector of
        `angle` radians whose two end rays stay separate node chains (the seam
        pair periodic links). order=2 curves the rims. The mesh carries side
        regions inner, outer and, for a sector, start and end. Reads no input
        mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_lshape)
    def lshape(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 4,
        size: float = 2.0,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate the L-shaped domain: the square of side `size` minus its
        fourth quadrant, with the re-entrant corner at the centre (a singular
        solution for refine and error estimation). The mesh carries side regions
        xmin, xmax, ymin, ymax, notch_x and notch_y. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_plate_with_hole)
    def plate_with_hole(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 0.25,
        size: float = 1.0,
        nt: int = 8,
        nr: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a square plate of side `size` around a circular hole of
        `radius`, in four quad blocks (the stress-concentration benchmark).
        order=2 makes the hole circular. The mesh carries side regions hole, xmin,
        xmax, ymin and ymax. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_honeycomb)
    def honeycomb(
        output_path: str,
        output_format: Optional[str] = None,
        dims: Optional[List[int]] = None,
        edge: float = 1.0,
        origin: Optional[List[float]] = None,
        max_cells: int = 20000000,
    ) -> dict:
        """Generate a tiling of pointy-top hexagons (dims = nx, ny) as one
        polygon block, a parametric fixture for the ragged polygon paths. The
        mesh carries a side region boundary. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_line)
    def line(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 8,
        start: Optional[List[float]] = None,
        end: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a straight segment of n line cells from start to end (default
        the origin to (1, 0, 0)), a 1-D and beam-format fixture. The mesh
        carries point regions start and end. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_circle)
    def circle(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        n: int = 16,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a closed polyline of n line cells on a circle in the z =
        center[2] plane, the seam welded. order=2 puts the mid-edge nodes on the
        arc. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_helix)
    def helix(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        pitch: float = 1.0,
        turns: float = 2.0,
        n: int = 32,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a helix (a spring) of n line cells around the z axis: turns
        turns rising pitch per turn. The mesh carries point regions start and
        end. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_sphere_points)
    def sphere_points(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 100,
        radius: float = 1.0,
        center: Optional[List[float]] = None,
        max_cells: int = 20000000,
    ) -> dict:
        """Generate n points on a sphere on a Fibonacci lattice, as vertex cells: a
        point-cloud fixture without poles or seams. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_box_points)
    def box_points(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 100,
        size: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        seed: int = 0,
        max_cells: int = 20000000,
    ) -> dict:
        """Generate n points drawn uniformly in a box from splitmix64 seeded with
        seed, as vertex cells: the same cloud on every platform and binding.
        Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_tetra_sphere)
    def tetra_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        subdivisions: int = 3,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a sphere surface from a regular tetrahedron split in four
        subdivisions times (4 * 4**subdivisions triangles), the most uneven
        seed: a mesh-quality probe. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_octa_sphere)
    def octa_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        subdivisions: int = 2,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a sphere surface from a regular octahedron split in four
        subdivisions times (8 * 4**subdivisions triangles), with the six axis
        points exact. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_cube_sphere)
    def cube_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        n: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a sphere surface from a cube's six faces of n x n quads (or
        triangles) pushed radially onto the sphere. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_geo_sphere)
    def geo_sphere(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        n: int = 4,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a geodesic sphere of frequency n: every icosahedron face cut
        into an n x n triangle lattice (20 n**2 triangles) on the sphere. Reads
        no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_torus)
    def torus(
        output_path: str,
        output_format: Optional[str] = None,
        major_radius: float = 1.0,
        minor_radius: float = 0.25,
        nu: int = 32,
        nv: int = 16,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a torus around the z axis, periodic both ways and welded: genus
        one (Euler characteristic zero), with negative Gaussian curvature on the
        inner rim. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_moebius)
    def moebius(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        width: float = 0.5,
        nu: int = 32,
        nv: int = 4,
        num_twists: int = 1,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a Moebius band with num_twists half-twists: an odd count has no
        consistent orientation (the compute_normals inconsistent_pairs probe).
        The mesh carries a side region boundary. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_tube)
    def tube(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        height: float = 1.0,
        ntheta: int = 16,
        nz: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate an open cylinder surface around the z axis, oriented outward.
        The mesh carries side regions bottom and top. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_cone)
    def cone(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        height: float = 1.0,
        ntheta: int = 16,
        nz: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a closed cone surface, the base disk at center and the apex
        height above it, the rim a sharp edge (a feature_edges probe). The mesh
        carries cell regions side and base. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_capsule)
    def capsule(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 0.5,
        height: float = 1.0,
        ntheta: int = 16,
        nz: int = 4,
        ncap: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a capsule: a cylinder of length height closed by two
        hemispheres, smooth where they meet. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_box_surface)
    def box_surface(
        output_path: str,
        output_format: Optional[str] = None,
        dims: Optional[List[int]] = None,
        size: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        cell_type: str = "quad",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate the boundary of a box as a closed surface, each face a grid,
        oriented outward. The mesh carries cell regions xmin ... zmax, one per
        face. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_platonic)
    def platonic(
        output_path: str,
        output_format: Optional[str] = None,
        name: str = "icosahedron",
        radius: float = 1.0,
        center: Optional[List[float]] = None,
        max_cells: int = 20000000,
    ) -> dict:
        """Generate a platonic solid's surface with one cell per face: tetrahedron,
        octahedron or icosahedron (triangles), cube (quads) or dodecahedron (a
        polygon block of pentagons). Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_ball)
    def ball(
        output_path: str,
        output_format: Optional[str] = None,
        radius: float = 1.0,
        n: int = 4,
        center: Optional[List[float]] = None,
        cell_type: str = "hexahedron",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a solid ball of hexahedra (or tetra) by the seven-block
        butterfly, so no cell degenerates at the centre; n cells from the centre
        to the surface. The mesh carries a side region surface. Reads no input
        mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_spherical_shell)
    def spherical_shell(
        output_path: str,
        output_format: Optional[str] = None,
        inner_radius: float = 0.5,
        outer_radius: float = 1.0,
        n: int = 4,
        nr: int = 2,
        ratio: float = 1.0,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a thick spherical shell of hexahedra between two radii, nr
        layers graded by ratio (each layer's thickness over the one inside it).
        The mesh carries side regions inner and outer. Reads no input mesh: it
        creates one."""

    @server.tool()
    @_forward(_tools.tool_pipe)
    def pipe(
        output_path: str,
        output_format: Optional[str] = None,
        inner_radius: float = 0.5,
        outer_radius: float = 1.0,
        height: float = 1.0,
        ntheta: int = 16,
        nr: int = 2,
        nz: int = 4,
        ratio: float = 1.0,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a thick-walled pipe of hexahedra around the z axis, nr layers
        through the wall graded by ratio. The mesh carries side regions inner,
        outer, bottom and top. Reads no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_fichera)
    def fichera(
        output_path: str,
        output_format: Optional[str] = None,
        n: int = 2,
        size: float = 2.0,
        center: Optional[List[float]] = None,
        cell_type: str = "hexahedron",
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate the Fichera corner: the cube of side size minus its positive
        octant, the 3-D L-shape with a re-entrant vertex at the centre. The mesh
        carries side regions xmin ... zmax, notch_x, notch_y and notch_z. Reads
        no input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_solid_torus)
    def solid_torus(
        output_path: str,
        output_format: Optional[str] = None,
        major_radius: float = 1.0,
        minor_radius: float = 0.25,
        n: int = 4,
        nu: int = 32,
        center: Optional[List[float]] = None,
        order: int = 1,
        max_cells: int = 20000000,
        complete: bool = False,
    ) -> dict:
        """Generate a solid torus of hexahedra: the squircle disk swept around the
        z axis and welded. The mesh carries a side region surface. Reads no
        input mesh: it creates one."""

    @server.tool()
    @_forward(_tools.tool_reference_cells)
    def reference_cells(
        output_path: str,
        output_format: Optional[str] = None,
        types: Optional[List[str]] = None,
        pairs: bool = False,
    ) -> dict:
        """Generate one cell of each requested type in reference coordinates, side
        by side (every supported type when types is omitted). pairs=True gives two
        cells sharing a facet, for linear types. Reads no input mesh: it creates
        one."""

    @server.tool()
    @_forward(_tools.tool_octree_mesh)
    def octree_mesh(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        output: str = "tet",
        min_size: float = 0.0,
        max_size: float = 0.0,
        max_level: int = -1,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.1,
        curvature_angle: float = 15.0,
        feature_angle: float = 30.0,
        corner_angle: float = 30.0,
        feature_size: float = 0.0,
        gap_cells: float = 2.0,
        gap_outside: bool = False,
        size_array: Optional[str] = None,
        surface_region: Optional[str] = None,
        grading_layers: int = 1,
        skin_grading: str = "graded",
        surface_size: float = 0.0,
        winding_threshold: float = 0.5,
        winding_beta: float = 2.0,
        conformity: str = "stuffing",
        alpha_long: float = 0.24799,
        alpha_short: float = 0.40173,
        preserve_features: bool = True,
        fit_iterations: int = 5,
        dihedral_floor_guard: bool = True,
        fill: str = "inside",
        keep_outside: bool = False,
        attach_distance: bool = True,
        attach_quality: bool = False,
        pillow: bool = True,
        hex_grading: str = "uniform",
        untangle_iterations: int = 50,
        max_cells: int = 20000000,
        max_tets: int = 100000000,
        volume_fix: bool = True,
        zero_volume_factor: float = 1.0e-6,
        cap_flatness: float = 0.01,
        crossing_margin: float = 0.01,
        crossing_gap: float = 0.05,
    ) -> dict:
        """Mesh the inside of a skin (any closed or slightly dirty surface:
        STL/OBJ/PLY/...) on a 2:1-balanced octree. output='castellated' gives
        the inside leaves as hexahedra with the 'inside' and 'boundary_cells'
        cell regions, the 'surrogate_boundary' side region (faces between
        inside and outside leaves) and sdf:distance / sdf:closest_point at the
        nodes, for embedded and Shifted Boundary Method solvers; output='tet'
        (default) conforming graded tetrahedra fitted to the skin; output='hex'
        all-hexahedra. Leaves are at most max_size, and near the skin the
        smallest of the curvature size (curvature_angle degrees per leaf), the
        feature size along edges sharper than feature_angle, the thickness over
        gap_cells and size_array, never below min_size (0: derived from the
        skin's diagonal). Inside is the fast winding number above
        winding_threshold, so holes and self-intersections are tolerated.
        skin_grading='uniform' keeps the skin at one level (the stuffing angle
        guarantee). fill='all'/'inside'/'cover' chooses the castellated leaves.
        hex_grading='graded' grades the hexahedra with the octree (dual
        transition schemes) instead of one grid at the skin's size.
        The tetrahedral output ends with a volume fix (volume_fix, default on):
        degenerate tetrahedra (volume under zero_volume_factor of the longest
        edge cubed) are collapsed, flat caps (cap_flatness) dropped or split,
        slivers (crossing_margin, crossing_gap) split at the crossing of their
        edges, inverted ones untangled or reordered, and overlaps counted.
        Reports the frame, leaf and cell counts, the quality (dihedral angles,
        scaled Jacobian, histograms), feature recovery, the volume-fix counts
        (num_collapsed, num_caps_removed, num_split, num_slivers,
        num_untangled, num_reordered, num_overlaps), timings and warnings."""

    @server.tool()
    @_forward(_tools.tool_mesh_cad_surface)
    def mesh_cad_surface(
        input_path: str,
        output_path: str,
        target_size: float,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        max_chordal_deviation: float = 0.0,
        min_points_per_edge: int = 0,
        tolerance: float = 0.0,
        record_face_ids: bool = True,
        edges: bool = True,
        match_unreferenced_trims: bool = True,
        strict: bool = False,
        max_points: int = 50000000,
        min_angle: float = 0.0,
        smooth_iterations: int = 0,
        smooth_method: str = "laplacian",
        order: int = 1,
        curvature_deviation: float = 0.0,
        cells_across_gaps: float = 0.0,
        grading: float = 0.0,
        pole_rings: int = 0,
        face_sizes: Optional[Dict[str, float]] = None,
        edge_sizes: Optional[Dict[str, float]] = None,
        background_path: Optional[str] = None,
        background_array: str = "size",
        metric: str = "mean",
        virtual_edges: bool = False,
        reparam_condition_number: float = 0.0,
        model_tolerances: bool = True,
        virtual_faces: bool = True,
        lop_corners: bool = False,
        metric_background_path: Optional[str] = None,
        metric_array: str = "",
        metric_error: float = 0.0,
        untangle: bool = False,
        quads: str = "none",
        quad_pairs: str = "matching",
        heal: bool = False,
        sew_tolerance: float = 0.0,
    ) -> dict:
        """Triangulate the faces of a CAD model (Kratos CAD JSON; STEP and
        IGES where OpenCASCADE is linked) into one conforming surface mesh with
        edges of about target_size, without a CAD kernel: every B-rep edge is
        sampled once and shared by the faces that use it (so closed solids come
        out watertight), and each face is a constrained Delaunay triangulation
        of its parameter domain lifted onto its surface. max_chordal_deviation
        caps the sagitta of the boundary segments (0: off); min_points_per_edge
        puts at least that many points inside every edge; tolerance (0: 1e-7 of
        the model's size) decides degenerate edges and pairs trims no edge
        names when match_unreferenced_trims. The output has one triangle block
        with cad:face cell data and a face:<k> cell region per face
        (record_face_ids) and a line block of the edges with cad:edge (edges).
        A face that cannot be meshed is left out and reported (status, reason)
        rather than meshed wrong; strict=True makes that an error. Reports the
        counts (failed faces, open, non-manifold and misoriented edges), the
        extreme angles and edge lengths, max_edge_gap, one record per face and
        the warnings. min_angle (degrees, under 60; 0: off) refines every face
        by Delaunay refinement until its triangles' smallest angles on the
        surface reach it (guaranteed up to 20.7 away from corners under 60
        degrees, a target with a reported shortfall above); the corners under
        60 degrees and the poles that keep triangles below it are listed per
        face (corners, below_bound, bound_deg), with num_refinement_points,
        num_faces_below_bound and an angle_histogram of 10-degree bins.
        smooth_iterations passes of guarded smoothing (smooth_method
        'laplacian' or 'odt') move interior points only where no triangle
        inverts and no star's smallest angle drops. order 2 or 3 writes curved
        triangle6 or triangle10 elements (line3, line4 for the edges) whose
        nodes lie on the surface and the edges; an element whose Jacobian is
        not positive everywhere is written straight-sided and counted in
        num_curved_straightened. curvature_deviation (0: off) sizes curved
        regions so a triangle sits at most that far from the surface, and
        cells_across_gaps (0: off) puts at least that many elements across a
        thin wall or a gap: both through one size field over the whole model,
        measured on a first uniform mesh, never above target_size, whose
        growth grading limits (a ratio above 1; 0: 1.5). pole_rings (0 to 8;
        0: the fan) places that many rings of points around each pole so its
        triangles are graded and meet min_angle. face_sizes and edge_sizes ({"index":
        size}, 0-based into the B-rep's faces and edges) and a background mesh
        (background_path, its point_data background_array the size at its
        points) give sizes that join the same field and only refine; each
        face's limited_by_edges names the edges finer than its given size
        allows. metric "varying" follows each face's metric as it varies
        instead of its mean; virtual_edges samples chains of short edges
        meeting smoothly as one, listed in virtual_edges.
        reparam_condition_number (0: off; else above 1; any order) meshes a
        face whose metric's condition number over its unrefined triangulation
        passes it in a flattening of that triangulation (a conformal map, or
        mean-value coordinates on a circle), every point still on the
        surface, listed in reparametrised (condition_number,
        condition_number_disc, taken, method 'conformal', 'mean_value' or
        'annulus' (a band between two circles), reason); model_tolerances
        accepts the gaps the model declares in its edges' and vertices'
        extras (a STEP or IGES read writes OpenCASCADE's) and the file's
        tolerances.model_tolerance, which also pairs trims no edge names,
        listed in tolerance_uses; virtual_faces (any order) meshes a face
        narrower than target_size with the neighbour it meets tangentially
        along the most edge, a chain of slivers joining the group the first
        forms, each triangle keeping its own face, listed in merged_faces
        (model_tolerances and virtual_faces are on by default since v19.0.0);
        a metric_background_path mesh whose 9-component point_data
        metric_array is a field's Hessian (hessian writes "<array>:hessian"),
        or whose 1-component one is the field itself (its Hessian then
        taken), with the target interpolation error metric_error, sizes the
        triangles by the anisotropic metric R |L| R^T / (8 metric_error)
        (Loseille and Alauzet 2011), triangles stretching where the field
        varies in one direction only (any order; min_metric_angle reports
        the least corner angle in the metric); lop_corners
        protects every acute corner by corner lopping (Pav and Walkington
        2005) in place of concentric shells, only with refinement, listed in
        lopped_corners (angle_deg, lfs0, lop_radius, turning_deg,
        stated_bound_deg, num_forced, min_forced_deg, within_two_radii,
        guaranteed, reason): the angle bound the paper proves on a planar
        face, and whether the run meets it. untangle (order 2 or 3) moves a
        face's own high-order nodes on the surface to untangle a folded curved
        element before it is straightened (Toulorge et al. 2013), then slides
        the nodes on the model's edges (quads' nodes too), counted in
        num_untangled. quads 'dominant' recombines each face's triangles in
        pairs into a quad block beside the triangles left; 'all' splits the
        quad-dominant mesh of twice the size (its size sources, a metric
        included, asked for twice their size) through its edges' midpoints
        into quads only, and takes no virtual_edges; with order 2 or 3 either
        writes curved quad9 or quad16 beside triangle6 or triangle10;
        quad_pairs 'matching' (the default, a maximum matching, Blossom-quad)
        or 'greedy' (best pair first) pairs them; reported as num_quads,
        num_leftover_triangles and min_quad_scaled_jacobian. heal (STEP/IGES)
        heals the shape before it is read (ShapeFix_Shape, then sewing at
        sew_tolerance; 0 is automatic, a positive value sews without heal)."""

    @server.tool()
    @_forward(_tools.tool_octree)
    def octree(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        root_cells: List[int] = [1, 1, 1],
        min_level: int = 0,
        max_level: int = 4,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.0,
        dimension: int = 0,
        sample_distance: float = 0.0,
        balance: str = "face",
        fill: str = "all",
        attach_distance: bool = False,
        even_odd: bool = False,
        max_cells: int = 20000000,
        sign: str = "pseudonormal",
        transition: str = "hanging",
        snap: bool = False,
        surface_level: Optional[int] = None,
        curvature_threshold: float = 0.0,
        feature_angle: float = 0.0,
        size_array: Optional[str] = None,
        plane: str = "auto",
    ) -> dict:
        """Build an octree (3-D, hexahedra) or quadtree (2-D, quads) around a mesh
        skin: leaves are fine only where the skin is. A closed polyline of line
        cells (or the boundary of a 2-D surface mesh) gives a quadtree, a surface
        an octree. Leaves are split while the skin cuts them (or lies within
        sample_distance) down to max_level; balance 'none'/'face'/'full' enforces
        2:1 size ratios. fill='all'/'surface'/'inside' keeps as for voxelize.
        Writes octree:level and octree:cut cell data and octree:hanging point
        data, and reports leaf counts, hanging nodes and the root lattice.
        transition='polyhedra' or 'templates' makes the output conforming (no
        hanging nodes); snap=True moves the outer nodes onto the skin (fill must
        be 'surface' or 'inside') and reports num_snapped/num_refused;
        surface_level, curvature_threshold, feature_angle and size_array refine
        beyond the cut leaves; plane chooses the plane of a 2-D skin."""

    @server.tool()
    @_forward(_tools.tool_voxelize)
    def voxelize(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        resolution: Optional[List[int]] = None,
        cell_size: Optional[float] = None,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.0,
        fill: str = "all",
        attach_occupancy: bool = False,
        max_cells: int = 20000000,
        sign: str = "pseudonormal",
        plane: str = "auto",
    ) -> dict:
        """Build a regular hexahedron grid around a mesh. Give exactly one of
        resolution (nx, ny, nz; nx, ny selects 2-D quad pixels around a planar
        polyline, in plane 'xy', 'xz' or 'yz', 'auto' taking the axis of zero
        extent) or cell_size. fill='all' keeps the whole bounding
        box, 'surface' keeps only cells a triangle passes through, and 'inside'
        keeps only cells whose centre is inside the surface. Reports the lattice
        origin/spacing/dims and how many cells were kept."""

    @server.tool()
    @_forward(_tools.tool_grid_sample)
    def grid_sample(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        resolution: Optional[List[int]] = None,
        cell_size: Optional[float] = None,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.0,
        fields: Optional[List[str]] = None,
        extrapolate: bool = False,
        fill_value: float = 0.0,
        max_cells: int = 20000000,
    ) -> dict:
        """Sample a mesh's point data onto a regular grid - the step a CNN or a
        superresolution model needs before it can see a mesh at all. Give exactly
        one of resolution (nx, ny, nz cell counts) or cell_size. Writes the grid
        as a lattice mesh; use .vti so the lattice survives the round trip.
        Reports the channel names (a k-component array expands to name_0..name_k)
        and the coverage, the fraction of grid points inside the mesh - a low
        coverage means the grid is mostly fill and a model would learn it."""

    @server.tool()
    @_forward(_tools.tool_grid_scatter)
    def grid_scatter(
        grid_path: str,
        target_path: str,
        output_path: str,
        grid_format: Optional[str] = None,
        target_format: Optional[str] = None,
        output_format: Optional[str] = None,
        fields: Optional[List[str]] = None,
        on_conflict: str = "error",
    ) -> dict:
        """Write a grid's fields back onto a mesh's points by trilinear
        interpolation - grid_sample's inverse, and the step that turns a model's
        output back into something every format can hold. Consecutive
        name_0/name_1/name_2 channels rebuild one multi-component array.
        on_conflict is 'error', 'overwrite' or 'suffix'."""

    @server.tool()
    @_forward(_tools.tool_grid_resample)
    def grid_resample(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        factor: Optional[int] = None,
        resolution: Optional[List[int]] = None,
        fields: Optional[List[str]] = None,
    ) -> dict:
        """Resample a grid onto a finer or coarser lattice covering the same box,
        by trilinear interpolation. This is the baseline a superresolution model
        has to beat. Give exactly one of factor (an integer upscale) or
        resolution (nx, ny, nz cell counts)."""

    @server.tool()
    @_forward(_tools.tool_grid_power_spectrum)
    def grid_power_spectrum(
        input_path: str,
        field: str,
        input_format: Optional[str] = None,
        max_bins: int = 256,
    ) -> dict:
        """The azimuthally averaged power spectrum of one field on a regular
        grid - whether a super-resolved or generated field carries the right
        small-scale content, which a pointwise error cannot see. Needs an
        isotropic grid (equal spacing on every axis) and reports power per
        wavenumber bin plus the mode count per bin. The power sums to
        mean(field**2) exactly."""

    @server.tool()
    @_forward(_tools.tool_subsample)
    def subsample(
        input_path: str,
        output_path: str,
        count: int,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "farthest",
        seed: int = 0,
        start: Optional[int] = 0,
        bounds: Optional[List[float]] = None,
        record_ids: bool = False,
    ) -> dict:
        """Reduce a mesh to a point cloud of exactly count points, so a large
        surface fits a transformer's token budget. method is 'farthest' (exact
        farthest-point sampling, O(N*count), the most uniform coverage), 'grid'
        (lattice representatives then farthest-point sampling, O(N + count^2),
        the scalable choice above a few hundred thousand points) or 'random'.
        start is the index of the first selected point (None draws it from
        seed); bounds (xlo, ylo, zlo, xhi, yhi, zhi) restricts the candidates to
        a box. The output keeps the selected points, their point_data and Point
        regions, and gets a vertex block so every format can hold it; cells and
        cell_data are dropped. record_ids attaches budget:original_point_id."""

    @server.tool()
    @_forward(_tools.tool_proximity_graph)
    def proximity_graph(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "radius",
        radius: Optional[float] = None,
        max_neighbors: Optional[int] = None,
        box_size: Optional[List[float]] = None,
        kind: str = "node",
    ) -> dict:
        """Build a graph from geometry rather than connectivity, for a particle
        state that has positions and no cells, or for contact edges between
        surfaces that are near but not connected. method is 'radius' (every pair
        closer than radius -- the physical interaction cutoff, degree varies with
        density) or 'knn' (each point's max_neighbors nearest, then symmetrized,
        so degree is near-constant). box_size (one value or one per axis) turns
        on the periodic minimum-image convention, so a pair either side of a
        boundary is linked by its short image; a radius past half the smallest
        box side is refused, the minimum image being ambiguous there. kind is
        'node' (mesh points) or 'cell' (cell centroids). The output is the graph
        as line cells over those positions, with a 'degree' point array."""

    @server.tool()
    @_forward(_tools.tool_compute_sdf)
    def compute_sdf(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        structure: str = "voxel",
        resolution: Optional[List[int]] = None,
        cell_size: Optional[float] = None,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.1,
        root_resolution: int = 8,
        max_depth: int = 4,
        band_cells: float = 1.0,
        max_cells: int = 20000000,
        sign: str = "pseudonormal",
        location: str = "corner",
        band: float = 0.0,
    ) -> dict:
        """Generate a grid over a surface and fill it with signed distances --
        the one call that turns a surface into a field. structure='voxel' is a
        dense lattice sized by resolution (nx, ny, nz) or cell_size, exactly one
        of which must be given; structure='octree' refines near the surface
        instead and sizes itself from root_resolution and max_depth, so passing
        resolution or cell_size with it is an error. An octree's output is
        1-irregular (it has hanging nodes). Reports the root dims, the finest
        spacing, the octree depth and the surface's watertightness. Write the
        result as .vti to keep the grid header; no other format carries it."""

    @server.tool()
    @_forward(_tools.tool_distance_to_surface)
    def distance_to_surface(
        input_path: str,
        surface_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        surface_format: Optional[str] = None,
        output_format: Optional[str] = None,
        sign: str = "pseudonormal",
        location: str = "corner",
        band: float = 0.0,
        record_inside: bool = False,
        record_closest_cell: bool = False,
    ) -> dict:
        """Attach the signed distance from a mesh's points (location='corner') or
        cell centres (location='center') to a surface, as sdf:distance. Negative
        is inside. sign='winding-number' is robust to holes and inconsistent
        winding but costs O(triangles) per query."""

    @server.tool()
    @_forward(_tools.tool_sample_distance)
    def sample_distance(
        input_path: str,
        points: List[List[float]],
        input_format: Optional[str] = None,
        sign: str = "pseudonormal",
        band: float = 0.0,
    ) -> dict:
        """Signed distances from a list of [x, y, z] points to a surface, without
        writing a mesh. Negative is inside."""

    @server.tool()
    @_forward(_tools.tool_surface_watertight_check)
    def surface_watertight_check(
        input_path: str, input_format: Optional[str] = None
    ) -> dict:
        """Report what is wrong with a surface, in numbers: boundary edges,
        non-manifold edges, inconsistently wound pairs and degenerate triangles.
        A signed distance is only meaningful where these are zero."""

    @server.tool()
    @_forward(_tools.tool_gradient)
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

    @server.tool()
    @_forward(_tools.tool_hessian)
    def hessian(
        input_path: str,
        output_path: str,
        array: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "green-gauss",
        location: str = "cell",
        output: Optional[str] = None,
        overwrite: bool = False,
    ) -> dict:
        """Hessian (second derivative) of a scalar point_data field --
        gradient's companion one order further, for curvature-based adaptive
        refinement.

        A composition of TWO gradient calls, not a new numerical kernel: the
        field is differentiated once (point location), then that (n,3)
        gradient is differentiated again with the default gradient operator,
        producing (n,9) -- the flattened row-major 3x3 Hessian. method is
        forwarded to BOTH internal passes. The result is named
        <array>:hessian unless `output` overrides it. Input must have exactly
        one component -- a vector field's Hessian is a separate quantity per
        component, call this once per component. A field that is at most
        linear has an exactly zero Hessian everywhere; for a genuinely
        quadratic field the composition is exact on a structured/symmetric
        mesh away from its own boundary and a good, standard, but genuinely
        approximate curvature estimate on an irregular mesh. Cells that
        cannot be evaluated yield NaN and are reported in num_skipped;
        least-squares cells with a degenerate neighbourhood in either
        internal pass fall back to Green-Gauss and are reported in
        num_fallback. A curvature-driven refinement indicator needs no new
        tool: `norm(...)` in data_calc on the 9-component output is exactly
        its Frobenius norm, ready for refine's where selector."""

    @server.tool()
    @_forward(_tools.tool_curvature)
    def curvature(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        mean: bool = True,
        gaussian: bool = True,
        dual_area: str = "mixed-voronoi",
        include_boundary: bool = False,
        record_area: bool = False,
        record_principal: bool = False,
        region: str = "",
    ) -> dict:
        """Per-vertex mean (H) and Gaussian (K) curvature of a surface mesh --
        the signed distance's natural companion as a node feature: sdf says
        how far a point is from the surface, this says how the surface bends
        there.

        K is the angle defect and H the cotangent Laplace-Beltrami operator,
        the estimators the discrete-differential-geometry convergence results
        are about. Writes curvature:mean and curvature:gaussian, optionally
        curvature:area (the dual area each was divided by) and
        curvature:principal ((n,2), k1 >= k2). Triangles come from the same
        fan convert_cells(simplexify) uses, so a quad mesh's curvature is the
        curvature of its canonical triangulation; a volume or polyhedron block
        is refused by name pointing at extract_surface and a higher-order one
        pointing at linearize. dual_area is "mixed-voronoi" (default, better
        on an irregular tessellation) or "barycentric" (cruder but
        branch-free). Boundary vertices are NaN unless include_boundary, and
        isolated ones always; both are counted. total_angle_defect is the
        oracle: on a CLOSED surface it is 2*pi*chi exactly -- 4*pi for
        anything sphere-like -- whatever the tessellation and whichever dual
        area, so a value that is not that means the input is not closed or the
        result is not sane. H's sign is orientation-dependent and K's is not,
        so check quality.inconsistent_pairs before trusting a sign: a nonzero
        count means facets disagree about which side is out and H is
        sign-flipped in patches. This never repairs its input."""

    @server.tool()
    @_forward(_tools.tool_normals)
    def normals(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        point_normals: bool = True,
        cell_normals: bool = False,
        weight: str = "angle",
        split_angle: Optional[float] = None,
        record_parent_ids: bool = False,
        region: str = "",
    ) -> dict:
        """Point and cell normals of a surface mesh, optionally splitting
        vertices at creases so every point carries exactly one normal.

        A vertex normal is a property of a smooth patch, not of a position: at
        the edge of a cube the one position has three normals. With
        split_angle unset the result is one smooth normal per point -- the
        angle- (weight "angle", default) or area- (weight "area") weighted
        mean of the incident faces. With split_angle set to a number of degrees
        in [0, 180] the corners around a vertex are grouped into smooth fans
        and every fan beyond the first gets its own copy of the point,
        appended after the original points; cells keep their numbering,
        point_data is gathered by row and a copy joins its source's point
        regions. Writes `normals` as point data (n, 3) and, with cell_normals,
        as cell data; record_parent_ids adds normals:parent_point. A `.pcd` or
        `.xyz` output then carries the normal columns. Triangles come from the
        same fan convert_cells(simplexify) uses; a volume or polyhedron block
        is refused by name pointing at extract_surface and a higher-order one
        pointing at linearize. Never reorients: a nonzero
        quality.inconsistent_pairs means some normals average faces that
        disagree about which side is out, and a split always cuts at such an
        edge; run repair first. Reports num_added_points, num_split_points,
        num_isolated / num_undefined (points whose normal is NaN) and
        num_degenerate."""

    @server.tool()
    @_forward(_tools.tool_check_quality)
    def check_quality(
        input_path: str,
        require: Optional[Union[str, List[str]]] = None,
        max_inverted: int = 0,
        max_degenerate: int = 0,
        file_format: Optional[str] = None,
    ) -> dict:
        """Quality gate: pass/fail thresholds over every cell's quality metrics
        -- what a CI job over meshes asserts on.

        require is the threshold text (or a list of texts): clauses separated
        by ";", "," or newlines, each "METRIC >= VALUE" or "METRIC <= VALUE",
        optionally "@ FRACTION" -- the fraction of evaluated cells allowed to
        violate it, a number in [0, 1] or a percentage like "@1%". METRIC is a
        compute_quality metric (scaled_jacobian, aspect_ratio, skewness,
        min_angle, max_angle, warpage, min_dihedral, max_dihedral, volume,
        inverted, degenerate), with or without "quality:". Values are tested
        per cell, not against the histograms; cells where a metric does not
        apply are not evaluated. max_inverted / max_degenerate bound those
        counts (default 0; negative disables). Reports passed, num_cells,
        num_inverted, num_degenerate and one entry per check: name, metric,
        min, max, max_fraction, evaluated, violations, fraction, worst,
        worst_cell and passed."""

    @server.tool()
    @_forward(_tools.tool_check_conformity)
    def check_conformity(
        input_path: str,
        output_path: Optional[str] = None,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        tolerance: float = 1e-9,
    ) -> dict:
        """Whether a mesh is conforming: every two cells meet in a whole
        shared face, edge or corner with the same point ids.

        Examines the volume cells (else the 2-D cells). Reports
        num_hanging_nodes (a point on a volume cell's face or edge that the
        cell does not use), num_t_junctions (the same on a 2-D cell's edge),
        num_partial_overlaps (unmatched facets in one plane overlapping over a
        positive area), num_unwelded_facets (coincident facets with different
        point ids), num_boundary_facets, num_non_manifold_facets, conforming,
        warnings and interface_mode: "conforming", "weld" (merge the
        coincident points first) or "proximity" -- the find_interface mode
        that suits the mesh. tolerance is a fraction of the examined facets'
        bounding-box diagonal. Pass output_path to write the mesh with
        conformity:hanging point_data (1 hanging node, 2 T-junction)."""

    @server.tool()
    @_forward(_tools.tool_check_tetrahedra)
    def check_tetrahedra(
        input_path: str,
        file_format: Optional[str] = None,
        zero_volume_factor: float = 1e-6,
    ) -> dict:
        """Whether the tetrahedra of a mesh are usable. Examines the tetra and
        tetra10 cells (a quadratic one by its four corners) and reports
        num_tetrahedra, num_inverted (signed volume below minus the tolerance),
        num_degenerate (volume within it), num_overlaps (pairs of tetrahedra
        sharing a face, used by no other, whose opposite nodes lie on the same
        side of it), min_volume, max_volume, valid and warnings. The tolerance is
        zero_volume_factor times the longest edge cubed plus the effect of six
        significant digits of written coordinates. The mesh is not changed."""

    @server.tool()
    @_forward(_tools.tool_find_intersections)
    def find_intersections(
        input_path: str,
        output_path: Optional[str] = None,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        other_path: Optional[str] = None,
        other_format: Optional[str] = None,
        segments_path: Optional[str] = None,
        region: str = "",
        region_other: str = "",
        max_pairs: int = 0,
    ) -> dict:
        """Where a surface intersects itself, or where it and a second
        surface (other_path) interpenetrate or touch.

        Quads and polygons are fanned into triangles and every candidate
        pair is decided by an exact triangle-triangle test, so every hit is
        exact. On one surface, triangles sharing corners count only when they
        meet beyond them; on two, any contact counts. Reports pairs (global
        cell ids, the cell of input_path first, capped by max_pairs, 0 = all),
        num_pairs, num_triangle_pairs, num_coplanar (overlaps in one plane,
        which give no segment), num_segments, num_cells_a, num_cells_b,
        truncated, clean and warnings. Surfaces only: run extract_surface on a
        volume mesh first. Pass output_path to write the mesh with the
        intersections:count cell data and segments_path for the intersection
        segments as a line mesh."""

    @server.tool()
    @_forward(_tools.tool_thickness)
    def thickness(
        input_path: str,
        output_path: Optional[str] = None,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        other_path: Optional[str] = None,
        other_format: Optional[str] = None,
        direction: str = "auto",
        max_distance: float = 0.0,
        min_thickness: Optional[float] = None,
        region: str = "",
        region_other: str = "",
        method: str = "ray",
        cone_angle: float = 15.0,
        cone_rays: int = 16,
        two_sided: bool = False,
        max_normal_deviation: float = 90.0,
    ) -> dict:
        """Wall thickness or clearance per facet of a surface, by a ray from
        each cell's centroid along its normal.

        direction "inward" stops at the far side of the wall (thickness),
        "outward" at the next wall facing the ray (clearance); "auto" is
        inward for one surface and outward against other_path. The surface
        must be consistently oriented with outward normals. max_distance
        bounds the ray (0 = the bounding-box diagonal). method "sphere" reads
        twice the distance from the ray's midpoint to the nearest other
        surface (the perpendicular distance where walls diverge), "cone" the
        shortest of cone_rays rays within cone_angle degrees (a chamfer's
        narrow spot); two_sided casts both ways with either face counting (a
        scan's unreliable winding; the ray is then a fifth of the diagonal by
        default) and max_normal_deviation steps over faces more edge-on than
        that. Reports direction, method, num_cells, num_hit, num_missed,
        num_below (below min_thickness), min, max, mean, p2, p95, min_cell,
        min_point, max_distance and warnings; with
        min_thickness also threshold and passed. Pass output_path to write the
        mesh with thickness:distance (NaN for a miss) and thickness:hit_cell
        cell data."""

    @server.tool()
    @_forward(_tools.tool_deviation)
    def deviation(
        input_path: str,
        reference_path: str,
        output_path: Optional[str] = None,
        file_format: Optional[str] = None,
        reference_format: Optional[str] = None,
        output_format: Optional[str] = None,
        patch: Optional[int] = None,
        tolerance: float = 0.0,
        min_share: Optional[float] = None,
        bins: int = 20,
        range_min: Optional[float] = None,
        range_max: Optional[float] = None,
        max_normal_angle: float = 90.0,
        reference_region: str = "",
        register: bool = False,
        element: Optional[str] = None,
    ) -> dict:
        """Signed deviation of every point of a scan (input_path) from a
        reference surface (reference_path): positive outside it (excess
        material), negative inside (missing material). With patch, measure
        against that patch label's fitted primitive of a mesh fit_primitives
        wrote instead. max_normal_angle is the far-side rule: a reference
        triangle whose normal turns away from the scan point's by more is
        skipped, so a thin wall is not measured to its back face (180 turns
        it off). Reports num_points, num_measured, min, max, mean, rms, sigma,
        percentiles (p1, p5, p50, p95, p99), histogram_edges and
        histogram_counts (bins over range_min..range_max, default +-max|d|),
        num_outside_range, share_within (within +-tolerance), suggested_range,
        num_without_normal and warnings; with min_share also passed. Pass
        output_path to write the scan with deviation:distance point data.
        register registers the scan to the reference first (an automatic
        start, point-to-plane ICP, adaptive trimming): the written scan is
        then the registered one, and the report adds registered, matrix and
        registration_rms. element measures against an element of the
        reference's element table (measure) instead: a plane, sphere,
        cylinder or cone."""

    @server.tool()
    @_forward(_tools.tool_deviation_regions)
    def deviation_regions(
        input_path: str,
        output_path: str,
        tolerance: float,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        field: str = "deviation:distance",
        prefix: str = "deviation:",
        min_area: float = 0.0,
    ) -> dict:
        """Connected regions (cells joined across shared edges) of a measured
        surface whose deviation, the mean of each cell's corners' field
        values, lies beyond +tolerance (excess material) or -tolerance
        (missing material). Each is written as the cell region
        <prefix>excess_<k> or <prefix>missing_<k>, with the cell data
        deviation:region. Reports regions (name, sign, num_cells, area,
        extreme, mean, volume, centroid, bounds; largest |volume| first),
        num_excess, num_missing, num_dropped (below min_area) and warnings."""

    @server.tool()
    @_forward(_tools.tool_exaggerate)
    def exaggerate(
        input_path: str,
        output_path: str,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        field: str = "deviation:distance",
        scale: float = 10.0,
        clamp: float = 0.0,
        smooth_fraction: float = 0.01,
    ) -> dict:
        """A copy of a measured surface with each point displaced along its
        angle-weighted normal by scale times its deviation (field), clamped
        to +-clamp (0: the 95th percentile of |d| rounded up to two
        significant figures) and averaged over the points within
        smooth_fraction of the bounding-box diagonal: a mesh for a figure or
        an animation. Reports clamp, radius, num_moved and warnings."""

    @server.tool()
    @_forward(_tools.tool_colorize)
    def colorize(
        input_path: str,
        output_path: str,
        field: str,
        file_format: Optional[str] = None,
        output_format: Optional[str] = None,
        cmap: str = "viridis",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        symmetric: bool = False,
        steps: int = 0,
        nan_color: str = "#808080",
        component: Optional[int] = None,
    ) -> dict:
        """Bake a point or cell data field into UInt8 red, green and blue data
        through a colormap (viridis, coolwarm, turbo, ...): the vertex or face
        colours a PLY file carries and glTF takes with color_by "rgb". vmin
        and vmax set the range (default the field's min and max; symmetric
        +-max|v|), steps of 2 or more bands it at fixed steps, nan_color
        colours non-finite values. Reports vmin, vmax and location."""

    @server.tool()
    def render_mesh(
        input_path: str,
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        width: int = 640,
        height: int = 480,
        cols: int = 100,
        rows: int = 40,
        encoding: str = "halfblock",
        text_format: str = "plain",
        view: Optional[str] = None,
        azimuth: Optional[float] = None,
        elevation: Optional[float] = None,
        perspective: bool = False,
        zoom: Optional[float] = None,
        shading: str = "flat",
        edges: str = "none",
        color_by: Optional[str] = None,
        component: Optional[int] = None,
        cmap: Optional[str] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        colorbar: bool = False,
        supersample: int = 2,
        background: Optional[str] = None,
        axes: bool = False,
        expr: Optional[str] = None,
        reduce: Optional[str] = None,
        clip: Optional[list] = None,
        symmetric: bool = False,
        scale: Optional[str] = None,
        categorical: bool = False,
        color_regions: bool = False,
        category_edges: bool = False,
        isolines: int = 0,
        vectors: Optional[str] = None,
        warp: Optional[str] = None,
        warp_scale: float = 1.0,
        warp_outline: bool = False,
        diagnostic: Optional[str] = None,
        quality_metric: Optional[str] = None,
        cutaway: Optional[list] = None,
        cutaway_tint: Optional[str] = None,
        streamlines: Optional[str] = None,
        stream_seeds: int = 40,
        stream_length: float = 0.5,
        stream_color: Optional[str] = None,
    ):
        """Look at a mesh: draw it with the software rasterizer, which needs
        no display, GPU or optional extra (unlike screenshot).

        Without output_path the picture comes back as text, cols x rows
        terminal cells: text_format "plain" (glyphs only, the default),
        "ansi" (24-bit colour escapes) or "html"; encoding is halfblock,
        quadrant, sextant, braille or ascii. With output_path the extension
        picks the form: .png (width x height pixels, returned as an image),
        .txt, .ansi, .html, or .cast (an asciinema orbit). The camera is a
        named view (iso, +x, -x, +y, -y, +z, -z: the side the camera sits
        on) or azimuth/elevation in degrees, orthographic unless
        perspective; zoom magnifies the fitted view. shading is none, flat
        or smooth; edges none, all or feature. color_by maps a point or cell
        data array through cmap (viridis, coolwarm, turbo, ...) over
        vmin..vmax, with an optional colorbar; component picks one component
        of a vector array (else its magnitude). background is #rrggbb,
        #rrggbbaa or none (transparent, the default); axes draws the world
        axes. Volume meshes are drawn through their boundary skin.

        Field rendering: expr colours by a data_calc expression ("mag(u) /
        max(p)") instead of color_by, and reduce ("mises", "hydrostatic" or
        "principal", with component) by a tensor invariant of color_by; clip
        is [low, high] percentiles bounding the range (null for either end),
        symmetric centres it on zero, and scale is "linear", "log" or
        "symlog". categorical colours integer data (a material id) and
        color_regions the named cell regions from a fixed palette, with a key
        in the notes under the picture; category_edges draws where they meet.
        isolines draws that many contour lines of the point array color_by,
        vectors draws arrows for a vector point array, and warp moves the
        points by a displacement point array (times warp_scale), with
        warp_outline drawing the undeformed outline. diagnostic is one of
        quality (with quality_metric, e.g. scaled_jacobian), inverted,
        degenerate, orientation (front and back faces), free_edges (open,
        non-manifold and inconsistent edges) or edge_length. cutaway clips
        the geometry away to look inside: a list of at most two planes, each
        a string like "+x:0.5" (keep x >= 0.5; "-x:0.5" keeps x <= 0.5) or
        six numbers [px, py, pz, nx, ny, nz] (a point and the normal of the
        side kept); the back faces then in view, the inside of the cut, are
        drawn in cutaway_tint (#rrggbb). streamlines draws lines that follow
        a vector point array over the mesh's cells from about stream_seeds
        seeds (default 40), each at most stream_length diagonals of the model
        long (default 0.5) in each direction, in stream_color (#rrggbb); on a
        volume they run inside it, so combine them with cutaway. Returned
        text ends with the colour range, ticks and keys."""
        report = _guard(
            _tools.tool_render_mesh,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            width=width,
            height=height,
            cols=cols,
            rows=rows,
            encoding=encoding,
            text_format=text_format,
            view=view,
            azimuth=azimuth,
            elevation=elevation,
            perspective=perspective,
            zoom=zoom,
            shading=shading,
            edges=edges,
            color_by=color_by,
            component=component,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            colorbar=colorbar,
            supersample=supersample,
            background=background,
            axes=axes,
            expr=expr,
            reduce=reduce,
            clip=clip,
            symmetric=symmetric,
            scale=scale,
            categorical=categorical,
            color_regions=color_regions,
            category_edges=category_edges,
            isolines=isolines,
            vectors=vectors,
            warp=warp,
            warp_scale=warp_scale,
            warp_outline=warp_outline,
            diagnostic=diagnostic,
            quality_metric=quality_metric,
            cutaway=cutaway,
            cutaway_tint=cutaway_tint,
            streamlines=streamlines,
            stream_seeds=stream_seeds,
            stream_length=stream_length,
            stream_color=stream_color,
        )
        if isinstance(report, dict) and str(
            report.get("output_path", "")
        ).lower().endswith(".png"):
            return [Image(path=report["output_path"]), report]
        return report

    @server.tool()
    @_forward(_tools.tool_feature_edges)
    def feature_edges(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        feature_angle: float = 30.0,
        feature: bool = True,
        boundary: bool = True,
        non_manifold: bool = True,
        inconsistent: bool = True,
        region: str = "",
    ) -> dict:
        """The sharp, open, non-manifold and inconsistently wound edges of a
        surface -- or, for a volume mesh, of its boundary skin -- written as a
        mesh of `line` cells over the input's points.

        Each edge shared by two faces is compared by the dihedral angle
        between their normals: above feature_angle degrees (0-180, default 30)
        it is a feature edge. An edge used by one face is a boundary edge, by
        three or more non-manifold, and a pair that walks the edge the same
        way is inconsistent (its angle is measured after reorienting one
        face). Each category can be switched off. The output carries cell data
        feature:kind (1 feature, 2 boundary, 3 non-manifold, 4 inconsistent)
        and feature:angle (degrees, NaN where undefined); region restricts to
        a named cell region. Reports num_edges written and num_feature /
        num_boundary / num_non_manifold / num_inconsistent over the whole
        surface. This is the crease test decimate and smooth pin nodes with."""

    @server.tool()
    @_forward(_tools.tool_segment)
    def segment(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        feature_angle: float = 30.0,
        normal_tolerance: float = 0.0,
        curvature_tolerance: float = 0.0,
        min_patch_cells: int = 1,
        region: str = "",
        prefix: str = "patch",
        method: str = "angle",
        tolerance: float = 1e-3,
    ) -> dict:
        """Split a surface into smooth patches -- the first step of turning a
        triangle skin (an STL, a scan) into a CAD model.

        A patch grows from its lowest-numbered cell across every edge two
        cells share whose dihedral angle is at most feature_angle degrees
        (default 30); open, non-manifold and sharp edges are never crossed.
        normal_tolerance (degrees) also bounds a cell's normal against its
        seed's, and curvature_tolerance (1/length) its curvature estimate
        against the patch's mean so far -- both off at 0. Patches with fewer than
        min_patch_cells cells merge into the neighbour sharing the most edges.
        The mesh is written unchanged plus cell data segment:patch, field data
        segment:adjacency (rows a, b, shared edges) and one cell region per
        patch named <prefix>:<k>. Reports num_patches, num_merged,
        num_feature_edges, num_boundary_edges and num_non_manifold_edges.
        method="primitive" ends a patch where the surface its facets lie on
        changes instead: a smooth area one plane, sphere, cylinder, cone or
        torus explains is one patch, any other is split into the regions
        primitives explain within tolerance (relative to the diagonal), so a
        fillet comes apart from the planes it is tangent to -- what a
        CAD-exported STL wants. Feed the output to fit_primitives."""

    @server.tool(name="select")
    @_forward(_tools.tool_select)
    def select_region(
        input_path: str,
        output_path: str,
        name: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        kind: str = "cell",
        bbox: Optional[List[float]] = None,
        plane_origin: Optional[List[float]] = None,
        plane_normal: Optional[List[float]] = None,
        sphere_center: Optional[List[float]] = None,
        sphere_radius: Optional[float] = None,
        where_array: Optional[str] = None,
        where_compare: str = "<",
        where_value: float = 0.0,
        quality_metric: Optional[str] = None,
        quality_compare: str = "<",
        quality_value: float = 0.0,
        component: Optional[int] = None,
        seed: Optional[List[float]] = None,
        feature_angle: float = 30.0,
        normal_tolerance: float = 0.0,
        mode: str = "all",
        boundary_only: bool = True,
        tag: int = -1,
    ) -> dict:
        """Write one named region -- a boundary-condition patch, a contact
        surface, a material zone -- so a solver deck or find_interface can
        name it. The mesh is otherwise written unchanged; a region of the same
        name and kind is replaced.

        kind is 'cell' (global cell ids), 'point', or 'side' ((cell, local
        facet) pairs: faces of the 3-D cells, or edges of the 2-D cells when
        there are none; boundary facets only unless boundary_only=false).
        Every criterion given must pass, and at least one is required: bbox
        [xmin,ymin,zmin,xmax,ymax,zmax], a half-space (plane_origin +
        plane_normal, the +normal side), a ball (sphere_center +
        sphere_radius), a data comparison (where_array on point_data for
        points, cell_data for cells and sides, where_compare in <, <=, >, >=,
        ==, !=, where_value), a quality comparison (quality_metric such as
        scaled_jacobian or min_dihedral, quality_compare, quality_value;
        cells and sides only), the k-th node-connected component (0-based, as
        split numbers its pieces), and seed [x,y,z]: a patch grown from the
        facet nearest it across edges no sharper than feature_angle degrees,
        within normal_tolerance of the seed's normal (0 = off). mode 'all',
        'any' or 'centroid' says how a cell or side reduces its nodes for the
        box, plane and ball. Reports num_candidates, num_selected, num_grown
        (-1 without a seed) and num_components (-1 without that criterion)."""

    @server.tool()
    @_forward(_tools.tool_iga_refine)
    def iga_refine(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        subdivisions: Union[int, List[int]] = 1,
        elevate: Union[int, List[int]] = 0,
        sample_subdivisions: int = 0,
        model: str = "auto",
        trim_samples: int = 16,
        edges: bool = True,
        trim: bool = True,
    ) -> dict:
        """Refine every patch of a NURBS model exactly -- the isogeometric
        counterpart of refine, acting on the control net and knot vectors, not
        on a sampling. Works on an MFEM NURBS mesh, Kratos CAD JSON, STEP or
        IGES (the model a reader keeps beside its cells).

        elevate raises the degree by that many (an int, or one per parametric
        direction), then subdivisions splits every nonzero knot span into that
        many equal parts by knot insertion (an int or one per direction):
        k-refinement, the highest continuity. The geometry, the
        parametrisation, the NURBS fields and a trimmed face's trimming loops
        do not change. An MFEM NURBS mesh must refine every direction alike
        (its patches share knot vectors along differently oriented directions).
        The refined model is written: an MFEM NURBS mesh or CAD JSON keeps its
        patches, any other output format gets the sampled cells
        (sample_subdivisions per knot span, 0 automatic). model picks "nurbs"
        or "brep" when the input carries both; trim_samples, edges and trim
        shape the B-rep sampling as in iga_sample."""

    @server.tool()
    @_forward(_tools.tool_iga_sample)
    def iga_sample(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        subdivisions: int = 0,
        model: str = "auto",
        trim_samples: int = 16,
        edges: bool = True,
        trim: bool = True,
    ) -> dict:
        """Sample a NURBS model (MFEM NURBS mesh, Kratos CAD JSON, STEP, IGES)
        into Lagrange cells: the exact rational geometry at the lattice of
        each knot-span cell, subdivisions cells per span and direction (0 is
        automatic: one per span, about 16 per direction on a trimmed face; an
        MFEM model with its patch topology samples as the MFEM reader does).
        A NURBS mesh gets one block of the highest degree's quad/hexahedron/
        line cells, the NURBS fields as point data and iga:patch / iga:span
        cell data; a B-rep is sampled as read_cad_json does (trim_samples per
        knot span for the clip, edges for the 3-D edge curves, trim=false to
        keep every cell)."""

    @server.tool()
    @_forward(_tools.tool_iga_evaluate)
    def iga_evaluate(
        input_path: str,
        params: List[List[float]],
        input_format: Optional[str] = None,
        patch: int = 0,
        derivatives: int = 0,
        model: str = "auto",
        trim_samples: int = 16,
    ) -> dict:
        """Evaluate a patch of a NURBS model (MFEM NURBS mesh) or a face of a
        B-rep (Kratos CAD JSON, STEP, IGES) at parametric coordinates.

        params is a list of [u], [u, v] or [u, v, w] points, each inside the
        patch's knot range (at most 10000 per call). The answer carries the
        points (N x space dimension) and the NURBS fields; derivatives=1 adds
        the first derivatives (N x d x space dimension) and, for a surface in
        3-D, the unit normals; a B-rep face also reports inside: whether each
        point lies within its trimming loops (trim_samples per knot span)."""

    @server.tool()
    @_forward(_tools.tool_fit_primitives)
    def fit_primitives(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        patch_data: str = "segment:patch",
        tolerance: float = 1e-3,
        relative: bool = True,
        min_inlier_fraction: float = 0.9,
        types: Optional[List[str]] = None,
        max_iterations: int = 50,
        region: Optional[str] = None,
        seeds: Optional[List[List[float]]] = None,
        feature_angle: float = 30.0,
        rim_angle: float = 15.0,
        output_region: str = "fit:feature",
        used_points: str = "auto",
        criterion: str = "least_squares",
    ) -> dict:
        """Fit a plane, sphere, cylinder, cone or torus to each patch of a
        surface and report how well it fits -- a recovered radius, axis or
        plane normal is a measurement.

        Patches are the integer cell-data label patch_data (segment:patch,
        written by the segment tool; "" fits the whole surface as one patch).
        A type is accepted when its RMS orthogonal distance is at most
        tolerance (a fraction of the patch's bounding-box diagonal unless
        relative is false) and at least min_inlier_fraction of the points lie
        within it; the first accepted of plane, sphere, cylinder, cone, torus
        wins (types narrows the list), else the patch is reported unfitted
        with its best residual. The mesh is written with cell data fit:type /
        fit:rms / fit:max and field data fit:patch_*; the report lists each
        patch's type, params (x, y, z, ax, ay, az, r1, r2, angle_deg, 0), rms,
        max and fraction, the form error and extremes, sigma, num_used and
        num_rejected, plus counts per type.

        region fits one named cell region as one patch; seeds (a list of
        [x, y, z]) grows a feature from the facet nearest each point across
        edges below feature_angle while one primitive of types explains it,
        peels rim faces turning more than rim_angle from the fit, fits it as
        one patch and writes its cells as the cell region output_region
        (num_feature_cells). used_points "3sigma"/"2sigma"/"1sigma" refits on
        the points within k RMS distances until the set stops changing;
        "auto" is "3sigma" with a region or seeds, "all" otherwise.
        criterion "minimum_zone" (the narrowest zone holding every point, as
        GD&T form tolerances are defined), "min_circumscribed" or
        "max_inscribed" (a shaft's or a hole's envelope) refits planes,
        spheres and cylinders under that reference."""

    @server.tool()
    @_forward(_tools.tool_refine_patches)
    def refine_patches(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        segment_method: str = "angle",
        segment_tolerance: float = 1e-3,
        feature_angle: float = 30.0,
        tolerance: float = 1e-3,
        relative: bool = True,
        normal_weight: float = 1.0,
        epsilon: float = 1e-3,
        max_iterations: int = 10,
        min_patch_cells: int = 3,
        ransac_seed: int = 0x6D657368696F2B2B,
        types: Optional[List[str]] = None,
        allow_split: bool = True,
        allow_merge: bool = True,
        prefix: str = "patch",
    ) -> dict:
        """Refine a segmentation until one plane, sphere, cylinder, cone or
        torus explains each patch -- between segment and fit_primitives /
        fit_cad, so more of a skin becomes real CAD faces (a fillet strip a
        plane claimed, a cylinder cut in two, a cylinder inside a freeform
        patch).

        Hybrid variational shape approximation: each patch is explained by a
        primitive holding every point within tolerance (a fraction of the
        bounding-box diagonal unless relative is false; types narrows which),
        else a freeform quadric; patches grow back from their best-explained
        cells, cheapest first (cost: squared distance over tolerance squared
        plus normal_weight times the normal deviation), never across a
        feature edge above feature_angle degrees; with allow_split a primitive
        is found inside a freeform patch by RANSAC seeded by ransac_seed (the
        output is byte-identical for a seed), with allow_merge neighbours one
        primitive explains merge, and a freeform patch under min_patch_cells
        cells joins a neighbour. Stops when nothing changes, an iteration adds
        less than epsilon to the explained area share, or after max_iterations.
        A surface without segment:patch is segmented first (segment_method
        "angle" by default, feature_angle, segment_tolerance for "primitive").
        The mesh is written with segment:patch, the <prefix>:<k> regions and
        the fit:* data rewritten (fit_cad reads it unchanged) and refine:*
        per-patch field data. The report gives num_patches, num_iterations,
        num_splits, num_merges, num_absorbed, num_primitive, num_freeform,
        num_facets, the absolute tolerance, explained_area_fraction and
        initial_explained_area_fraction, hausdorff_max / hausdorff_rms,
        initial_energy / energy and stop (unchanged, small_gain,
        iteration_cap, energy_rose). Needs the compiled core."""

    @server.tool()
    @_forward(_tools.tool_region_adjacency)
    def region_adjacency(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        regions: Optional[list[str]] = None,
    ) -> dict:
        """Write conforming facets shared by Cell regions, or by cell blocks
        when fewer than two Cell regions exist. Reports total facets/measure
        and writes parent-cell/facet ids with per-facet length or area."""

    @server.tool()
    @_forward(_tools.tool_find_interface)
    def find_interface(
        input_path: str,
        output_path: str,
        region_a: str,
        region_b: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        mesh_b_path: Optional[str] = None,
        format_b: Optional[str] = None,
        mode: str = "conforming",
        master: str = "a",
        gap_tolerance: float = 0.0,
        angle_tolerance: float = 30.0,
        overlap_tolerance: float = 0.0,
    ) -> dict:
        """Write conforming or proximity-matched interface facets between
        two Cell regions. Optionally read part B from a second mesh file. The
        report includes side-region `(cell, facet)` pairs for each input."""

    @server.tool()
    @_forward(_tools.tool_contact_pairs)
    def contact_pairs(
        slave_path: str,
        slave_region: str,
        master_region: str,
        master_path: Optional[str] = None,
        slave_format: Optional[str] = None,
        master_format: Optional[str] = None,
        tolerance: float = 0.0,
        require_complete: bool = False,
    ) -> dict:
        """Project every node of a Point region onto the closest facets of a
        Cell region; return source ids, fan-subfacet ordinals, local coordinates,
        signed gaps and master normals."""

    @server.tool()
    @_forward(_tools.tool_split_interface)
    def split_interface(
        input_path: str,
        output_path: str,
        side_region: Optional[str] = None,
        side_entries: Optional[list[list[int]]] = None,
        side_name: str = "interface:side_a",
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        add_cohesive: bool = False,
    ) -> dict:
        """Duplicate point fans along a Side region, optionally adding
        cohesive cells. `side_entries` can be the `(cell, facet)` list returned
        by `find_interface`; polyhedron inputs are rejected."""

    @server.tool()
    @_forward(_tools.tool_hausdorff)
    def hausdorff(
        path_a: str,
        path_b: str,
        format_a: Optional[str] = None,
        format_b: Optional[str] = None,
        face_samples: int = 0,
        region_a: str = "",
        region_b: str = "",
        max_distance: Optional[float] = None,
    ) -> dict:
        """The Hausdorff distance between the surfaces of two mesh files (a
        volume mesh contributes its skin): how far apart they are at their
        worst. Each surface is sampled -- its vertices, plus with
        face_samples = s > 0 the centroids of the s*s sub-triangles of every
        triangle -- and each sample's unsigned distance to the other surface
        measured. Vertex sampling alone is exact when the farthest point is a
        vertex and a lower bound otherwise. Reports distance, the one-sided
        a_to_b / b_to_a maxima with their means and RMS, the sample counts and
        the worst sample of each side; with max_distance, also passed."""

    @server.tool()
    @_forward(_tools.tool_edit_regions)
    def edit_regions(
        input_path: str,
        output_path: str,
        edits: List[Dict[str, object]],
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Edit a mesh's named regions and write the result; points, cells
        and data are untouched.

        edits is a list applied in order, each {"op", "inputs", ...}: op is
        union / intersection / difference (two or more inputs of one kind; the
        result is named "output", keeps its inputs unless "keep_inputs" is
        false, and may set "dim"/"tag"), rename (one input, to "output"),
        retag (one input; new "tag" and/or "dim") or delete. An input is a
        region name or {"name", "kind", "dim", "tag"} and must match exactly
        one region; side regions combine their (cell, facet) pairs, and an
        empty result is kept. A result that would overwrite an unrelated
        region is an error. Reports the resulting region list."""

    @server.tool()
    @_forward(_tools.tool_graph)
    def graph(
        input_path: str,
        output_path: str,
        kind: str = "node",
        input_format: Optional[str] = None,
        file_format: Optional[str] = None,
    ) -> dict:
        """Write a mesh's graph as a METIS or Scotch graph file, for an
        external partitioner or a sparse-matrix preallocation.

        kind "node" joins nodes that share a cell edge; "clique" joins nodes
        that share any cell (the sparsity pattern of an assembled FEM matrix);
        "cell" joins cells that share a facet (a face in 3-D, an edge in 2-D),
        numbered block-major: the dual graph partitioners cut. file_format is
        "metis" (.graph: a header "n m", then each vertex's 1-based
        neighbours) or "scotch" (.grf, 0-based), inferred from the extension
        (.graph/.metis, .grf) when omitted. Reports output_path, kind, format,
        num_vertices, num_edges (undirected) and max_degree."""

    @server.tool()
    @_forward(_tools.tool_periodic)
    def periodic(
        input_path: str,
        slave: str,
        master: str,
        input_format: Optional[str] = None,
        translate: Optional[List[float]] = None,
        rotate_axis: Optional[Union[str, List[float]]] = None,
        rotate_degrees: Optional[float] = None,
        origin: Optional[List[float]] = None,
        matrix: Optional[List[float]] = None,
        atol: float = 1e-8,
        require_complete: bool = True,
        pairs_path: Optional[str] = None,
    ) -> dict:
        """Periodic node pairs: which master-region node each slave-region
        node maps onto under an affine transform, within atol -- what a
        periodic boundary condition (a Kratos periodic condition, a Gmsh
        $Periodic section) ties together.

        The transform is translate [dx, dy, dz], and/or a rotation of
        rotate_degrees about rotate_axis ("x", "y", "z" or a 3-vector) turning
        about origin (applied first), or a row-major 4x4 matrix of 16 numbers
        overriding both. Point regions contribute their nodes, cell regions
        their cells' nodes, side regions their facets' nodes. The nearest
        master within atol wins; a master claimed twice is an error, as is an
        unmatched slave node unless require_complete is false. Slave nodes the
        transform leaves in place that are also master nodes (on a rotation
        axis) are counted as num_fixed. Reports num_pairs, pairs
        ([slave, master], 0-based, ascending in the slave id), unmatched and
        max_residual; pairs_path also writes them as 'slave,master' CSV."""

    @server.tool()
    @_forward(_tools.tool_repair)
    def repair(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        fix_orientation: bool = True,
        orient_outward: bool = True,
        fill_holes: bool = True,
        split_non_manifold: bool = True,
        max_hole_edges: int = 10,
        weld_tolerance: float = 0.0,
        record_provenance: bool = False,
        fill_style: str = "fan",
    ) -> dict:
        """Repair a surface mesh's orientation, holes and pinched vertices --
        the three defects `clean` does not touch (it welds, drops degenerate
        and duplicate cells, and prunes orphans; it never rewinds a triangle
        or closes a hole).

        The passes run in order: weld (opt-in, `weld_tolerance`) ->
        triangulate (quads and polygons fan exactly as convert_cells
        simplexify does, blocks staying 1:1) -> split bowties (a vertex whose
        triangle star is edge-disconnected is duplicated, geometry unchanged)
        -> orient (a BFS per connected component by the TOPOLOGICAL half-edge
        rule -- two triangles sharing an edge agree iff they traverse it in
        opposite directions, which is exact where a normal-angle test
        mis-orients across a sharp crease; fewest flips wins ties) -> fill
        holes (every traceable boundary loop of at most max_hole_edges edges,
        `<= 0` for no limit, gets one centroid point and one triangle per loop
        edge, wound to AGREE with the surrounding surface) -> orient outward
        (every closed component whose signed volume is negative is flipped
        whole). Lower-dimensional blocks (boundary lines) ride along; fill
        triangles land in one trailing triangle block, added only when there
        is one, so a closed input keeps its block count.

        Reports the input's and the output's defect counts, so what was fixed
        and what remains are both visible. Non-manifold EDGES (used by three or
        more triangles) are neither split nor crossed -- they are counted;
        nested cavities are not detected, so every closed component is
        oriented outward on its own. Point, Cell and Side regions survive (a
        split copy joins its source's; an edge a flip renumbers is found again
        by its nodes).

        fill_style picks how a hole is filled: "fan" (one centroid point,
        the default), "flat" (minimum-area, refined to the rim's edge
        density, membrane-faired) or "curved" (thin-plate-faired, continuing
        the surface's curvature; roadmap 7.8.6)."""

    @server.tool()
    @_forward(_tools.tool_shrinkwrap)
    def shrinkwrap(
        input_path: str,
        target_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        target_format: Optional[str] = None,
        output_format: Optional[str] = None,
        offset: float = 0.0,
        max_distance: float = 0.0,
        weights: str = "",
        target_region: str = "",
        normal_weight: str = "angle",
        record_distance: bool = False,
        record_closest_cell: bool = False,
    ) -> dict:
        """Project a mesh's points onto a target triangle surface:
        x' = x + w (p + offset*n - x), with p the closest point on the target
        and n the unit pseudonormal there -- the fitting step a scanned skin, a
        CAD shell or a coarse solve needs before it can be used as a template.

        It is ONE projection, not an iteration: there is no self-intersection
        guard and no inversion guard, because a wrap is a fit and not a
        smoothing. Every point of the input moves whatever cells it carries --
        a volume mesh's interior points are projected too; only the TARGET must
        be a surface (quads and polygons are fanned, a volume or higher-order
        block is refused by name). Use `weights` to select or blend: an
        integer/bool point_data array selects (nonzero moves), a float one
        blends, applied unclamped so a caller can overshoot on purpose. A
        point farther than max_distance (`<= 0` means unlimited) is left where
        it is and counted, as is one whose hit feature has no direction to
        offset along.

        The offset goes along the pseudonormal of the hit FEATURE (face, edge
        or vertex), not of the selected triangle: at a crease the offset
        surface's normal is the bisector, so a face normal would land the
        point off that surface and make the result depend on which of two
        equidistant faces won the tie-break. Connectivity, data, regions and
        property sets pass through verbatim, and the points keep their
        dtype."""

    @server.tool()
    def convert_coordinates(
        input_path: str,
        output_path: str,
        system: str = "cylindrical",
        axis: Optional[List[float]] = None,
        origin: Optional[List[float]] = None,
        reference: Optional[List[float]] = None,
        arrays: Optional[List[str]] = None,
        inverse: bool = False,
        record_coordinates: bool = True,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Cylindrical or spherical components about an axis, for rotating
        machinery and pressure vessels. The points stay put;
        record_coordinates attaches coords:cylindrical (r, theta, z) or
        coords:spherical (r, theta polar from the axis, phi) point data, and
        each array in `arrays` -- point or cell data, a vector (n, 3) or a
        tensor (n, 9)/(n, 3, 3) -- is rewritten in the local basis of its point
        (a cell's at its centroid); inverse=true converts local components back
        to Cartesian. The axis (default +z) passes through origin; angles are
        measured from reference projected off the axis. Reports num_arrays and
        num_singular (points on the axis)."""
        return _guard(
            _tools.tool_convert_coordinates,
            input_path=input_path,
            output_path=output_path,
            system=system,
            axis=axis or (0.0, 0.0, 1.0),
            origin=origin or (0.0, 0.0, 0.0),
            reference=reference or (1.0, 0.0, 0.0),
            arrays=arrays,
            inverse=inverse,
            record_coordinates=record_coordinates,
            input_format=input_format,
            output_format=output_format,
        )

    @server.tool()
    @_forward(_tools.tool_prepare_metric)
    def prepare_metric(
        input_path: str,
        output_path: str,
        arrays: List[str],
        output: str = "metric",
        target_vertices: float = 0.0,
        norm: float = 2.0,
        hmin: float = 0.0,
        hmax: float = 0.0,
        gradation: float = 0.0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """The metric toolbox, no MMG needed: intersect several anisotropic
        metric fields (point data, 6 components in 3-D or 3 in 2-D, as
        metric_from_hessian writes) so one mesh resolves every field,
        normalise to target_vertices (L^p, norm = p, 0 = infinity), clamp
        sizes to [hmin, hmax] and grade (gradation = largest size growth per
        unit length, > 1). Writes point data `output` for adapt. Reports
        complexity_before/after (about the vertex count), num_graded,
        sweeps."""

    @server.tool()
    @_forward(_tools.tool_untangle)
    def untangle(
        input_path: str,
        output_path: str,
        fix_boundary: bool = True,
        frozen: Optional[str] = None,
        rings: int = 2,
        smooth: bool = False,
        max_sweeps: int = 100,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Repair inverted tetrahedra and hexahedra by moving nodes (Escobar
        et al.'s untangling functional), keeping the connectivity -- after a
        large morph, or in a result file from an explicit run. By default
        only nodes within `rings` edge rings of an inverted cell move and the
        boundary stays fixed; smooth=true untangles and smooths everything.
        `frozen` names a point_set or point-data mask of nodes to keep.
        Reports num_inverted_before/after, min_quality_before/after (least
        corner scaled Jacobian), sweeps, moves."""

    @server.tool()
    @_forward(_tools.tool_morph)
    def morph(
        input_path: str,
        output_path: str,
        array: str,
        fix_boundary: bool = True,
        constrained: Optional[str] = None,
        stiffening: float = 1.0,
        steps: int = 1,
        allow_inverted: bool = False,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Boundary-driven mesh morphing (shape optimisation, fluid-structure
        mesh motion): the constrained points -- the boundary unless
        fix_boundary is false, plus the point_set or nonzero point-data array
        named by `constrained` -- move by the point-data displacement `array`
        exactly; every other point moves by its harmonic extension with
        Jacobian-based stiffening (small cells move nearly rigidly). Give
        boundary points that stay put a zero displacement. Larger motions:
        more `steps`. A result with inverted cells fails unless
        allow_inverted. Reports num_constrained, num_free, num_inverted,
        max_displacement, converged."""

    @server.tool()
    @_forward(_tools.tool_laplacian_eigenmodes)
    def laplacian_eigenmodes(
        input_path: str,
        output_path: str,
        num_modes: int = 10,
        tolerance: float = 1e-8,
        max_iterations: int = 0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """The lowest eigenpairs of a triangle surface's Laplace-Beltrami
        operator (cotangent stiffness against lumped mass, Neumann on a
        boundary): writes point data laplacian:modes (n, num_modes),
        M-orthonormal, and field data laplacian:eigenvalues, ascending -- a
        spectral positional encoding for graph networks, or a shape
        descriptor. Quads and polygons are refused naming convert_cells
        simplexify, volumes naming extract_surface. Reports eigenvalues,
        iterations, max_residual, converged."""

    @server.tool()
    @_forward(_tools.tool_register_rigid)
    def register_rigid(
        source_path: str,
        target_path: str,
        output_path: Optional[str] = None,
        method: str = "point_to_point",
        max_iterations: int = 50,
        tolerance: float = 1e-10,
        max_distance: float = 0.0,
        rotate_data: bool = False,
        source_format: Optional[str] = None,
        target_format: Optional[str] = None,
        output_format: Optional[str] = None,
        initial: str = "identity",
        symmetric_check: bool = False,
        source_region: str = "",
        target_region: str = "",
        adaptive_trim: float = 0.0,
    ) -> dict:
        """Rigid registration by iterative closest points: the 4x4 matrix
        (rows) laying the source's points onto the target -- its surface when
        it has surface cells, else its points -- e.g. a scan onto its design.
        method "point_to_point" or "point_to_plane" (tangent planes, fewer
        iterations; a point-cloud target needs normals point data, see
        estimate_normals). max_distance > 0 drops pairs farther apart
        (outliers). ICP finds the nearest optimum, so the two should start
        roughly aligned -- or initial="auto": the best of the 24 candidates
        laying the two principal frames onto each other, coarse ICP from each
        (any start; symmetric_check warns when the part looks symmetric).
        source_region / target_region restrict the pairs to cell regions (the
        part, not its fixture); adaptive_trim drops pairs beyond that multiple
        of the median pair distance. With output_path the moved source is
        written too. Reports matrix, initial_rms, rms, iterations, num_pairs,
        converged, num_rejected, candidate_rms, chosen_candidate, warnings."""

    @server.tool()
    @_forward(_tools.tool_register_landmarks)
    def register_landmarks(
        source_points: List[List[float]],
        target_points: List[List[float]],
        weights: Optional[List[float]] = None,
        input_path: Optional[str] = None,
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        rotate_data: bool = False,
    ) -> dict:
        """The rigid motion laying picked points ([[x, y, z], ...], three or
        more, the sources not collinear) onto their targets: the weighted Horn
        solve in closed form, a start for register_rigid. With input_path and
        output_path the moved mesh is written too. Reports matrix, residuals
        and rms."""

    @server.tool()
    @_forward(_tools.tool_find_symmetry_plane)
    def find_symmetry_plane(
        input_path: str,
        input_format: Optional[str] = None,
        seed_origin: Optional[List[float]] = None,
        seed_normal: Optional[List[float]] = None,
        tolerance: float = 0.0,
        min_share: float = 0.9,
        max_samples: int = 4000,
    ) -> dict:
        """The plane a mesh or point cloud is most nearly symmetric about:
        each principal plane (and the seed plane) is mirrored, registered back
        by ICP and read off the reflection-motion product. Refuses a plane
        whose share of samples within tolerance (0: 0.005 of the diagonal) is
        below min_share. Reports origin, normal, rms, share, tolerance,
        candidate_rms, chosen_candidate and warnings."""

    @server.tool()
    @_forward(_tools.tool_auto_align)
    def auto_align(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        angle_tolerance: float = 10.0,
        subdivisions: int = 3,
        region: str = "",
        rotate_vector_data: bool = False,
    ) -> dict:
        """Pose a scan from itself alone: on its largest flat face, or a
        cylindrical part on the end of its axis, from an area-weighted
        histogram of its facet normals (an open scan's opening counts as a
        face). The result sits on z = 0, centred on the Z axis. Reports
        matrix, mode ("face" or "axis"), primary, secondary, share,
        face_share, axis_share, opening_area and warnings alongside the mesh
        summary."""

    @server.tool()
    def estimate_normals(
        input_path: str,
        output_path: str,
        k: int = 16,
        orientation: str = "consistent",
        viewpoint: Optional[List[float]] = None,
        record_variation: bool = False,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Normals of a point cloud (a PLY, PCD or XYZ scan, or any mesh's
        points): a plane fitted to each point and its k nearest neighbours,
        written as `normals` point data (the columns the PCD and XYZ writers
        emit). orientation: "consistent" (propagated along a spanning tree of
        the neighbour graph, then each piece turned outward; a flat patch
        faces +z), "viewpoint" (towards `viewpoint`, the scanner) or "none".
        record_variation adds normals:variation (0 on a plane). Points with
        no plane in their neighbourhood get NaN. Reports num_degenerate and
        num_components."""
        return _guard(
            _tools.tool_estimate_normals,
            input_path=input_path,
            output_path=output_path,
            k=k,
            orientation=orientation,
            viewpoint=viewpoint or (0.0, 0.0, 0.0),
            record_variation=record_variation,
            input_format=input_format,
            output_format=output_format,
        )

    @server.tool()
    @_forward(_tools.tool_reconstruct_surface)
    def reconstruct_surface(
        input_path: str,
        output_path: str,
        normals: str = "normals",
        k: int = 8,
        cell_size: float = 0.0,
        resolution: int = 0,
        max_distance: float = 0.0,
        estimate_normals: bool = True,
        normals_k: int = 16,
        max_cells: int = 20_000_000,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """A triangle surface through an oriented point cloud (a PLY, PCD or
        XYZ scan): the zero level set of the signed distance to the samples'
        tangent planes (Hoppe et al. 1992), cut by marching tetrahedra on a
        band of lattice cells about the samples. `normals` names the point
        data (estimated when absent unless estimate_normals is false, with
        normals_k neighbours); k tangent planes are blended per lattice node
        (1 is Hoppe's field). cell_size 0 takes the extent over `resolution`
        when given, else the sample spacing; max_distance 0 trims at twice
        the larger of spacing and cell size, which keeps an open scan open.
        Reports num_samples, spacing, cell_size, max_distance,
        num_lattice_points, num_undefined, num_triangles, num_components and
        estimated_normals."""

    @server.tool()
    @_forward(_tools.tool_convex_hull)
    def convex_hull(
        input_path: str,
        output_path: str,
        dimension: int = 0,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """The convex hull of a file's points (a scan, a point cloud, or any
        mesh: its cells are ignored). 3-D: a closed triangle surface wound
        outward; 2-D (the x, y footprint): a counter-clockwise line loop.
        dimension is 2, 3, or 0 to decide from the points. Only extreme
        points are vertices; each carries hull:point_id (its 0-based input
        point). Exact predicates, so coplanar points are handled without a
        tolerance. Reports dimension, num_vertices, num_facets, area (enclosed
        in 2-D, surface in 3-D), volume and perimeter."""

    @server.tool()
    @_forward(_tools.tool_delaunay_2d)
    def delaunay_2d(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """The Delaunay triangulation of a file's points in the x, y plane:
        scattered xyz/pcd points turned into a terrain (each point keeps its
        z and point data). One counter-clockwise triangle block covering the
        convex hull; a point coinciding in x, y with an earlier one is left
        out. Reports num_triangles, num_duplicates and num_boundary."""

    @server.tool()
    @_forward(_tools.tool_geodesic_distance)
    def geodesic_distance(
        input_path: str,
        output_path: str,
        sources: Optional[List[int]] = None,
        source_region: str = "",
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        time_factor: float = 1.0,
        tolerance: float = 1e-10,
        max_iterations: int = 0,
    ) -> dict:
        """Geodesic distance on a triangle surface: the distance along the
        surface from the source points to every point, written as point_data
        geodesic:distance, by the heat method (Crane et al. 2013) -- a short
        heat flow from the sources (time time_factor * h^2, h the mean edge
        length), its normalised gradient, and a Poisson solve.

        Sources are point ids (0-based) and/or the points of a region
        (source_region). The surface is the mesh's triangle blocks; a volume,
        quad, polygon or higher-order cell is refused naming the operation that
        makes triangles. A point in no triangle, or in a piece with no source,
        gets NaN (num_unreached). The method is first-order accurate (about 1%
        on a sphere at h = 0.03); a larger time_factor smooths it. Reports
        num_sources, num_unreached, num_degenerate, max_distance, time_step,
        heat_iterations, poisson_iterations and converged."""

    @server.tool()
    @_forward(_tools.tool_parameterize)
    def parameterize(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        padding: float = 0.02,
        max_depth: int = 8,
        max_stretch: float = 2.0,
    ) -> dict:
        """A UV atlas of a triangle surface: point_data uv, two texture
        coordinates per point in the unit square, which the OBJ (vt), glTF
        (TEXCOORD_0) and USD (primvars:st) writers turn into texture
        coordinates.

        Each edge-connected piece that is a disc is one chart; a closed or holed
        piece is cut into discs, and one that folds or whose L2 stretch exceeds
        max_stretch is halved (up to max_depth times). Charts are flattened by
        least-squares conformal maps (else mean-value coordinates) with no
        flipped triangle, keep their area, and are shelf-packed padding *
        sqrt(area) apart. Points on seams are split (point_data
        parameterize:parent_point names the input point; the first copy keeps
        its index); cell_data parameterize:chart gives each triangle's chart.
        The surface is the mesh's triangle blocks; a volume, quad, polygon or
        higher-order cell is refused naming the operation that makes triangles.
        Reports num_charts, num_flipped, stretch (1 keeps lengths),
        packing_efficiency and num_seam_points."""

    @server.tool()
    @_forward(_tools.tool_sobolev_deform)
    def sobolev_deform(
        input_path: str,
        output_path: str,
        array: str,
        length_scale: float,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        fixed_points_array: str = "",
        fix_boundary: bool = False,
        record_filtered: bool = False,
        max_iterations: int = 128,
        tolerance: float = 1e-10,
    ) -> dict:
        """Sobolev (Helmholtz-filtered) deformation: smooth a raw per-point
        displacement field through the mesh's own P1 finite-element operators,
        then move the points by the smoothed field.

        Solves (M + l^2 K) u = M d per ambient component and sets x' = x + u,
        with K the P1 stiffness matrix assembled from the simplex edge Gram
        matrix, M a uniform vertex mass and l = length_scale -- a
        screened-Poisson low-pass filter whose cutoff wavelength is l, which
        is what turns a jagged per-node displacement (a shape gradient, a
        scattered measurement, a model's raw output) into one a mesh can
        follow without tangling. length_scale = 0 applies the raw field at the
        free points.

        `array` names the (n, dim) point_data displacement. Every cell block
        at the mesh's top topological dimension must be a linear simplex --
        line, triangle or tetra -- since that is what the assembly is defined
        on; a quadratic block is refused naming linearize, anything else
        naming convert_cells(simplexify). Lower-dimensional blocks ride along,
        and a point in no top-dimensional cell receives its raw displacement.

        Nothing is pinned by default: an unfixed boundary carries the natural
        Neumann condition, so a constant displacement is preserved exactly.
        fixed_points_array (an integer/bool point_data mask) and fix_boundary
        impose zero-Dirichlet rows instead. Non-convergence within
        max_iterations is reported through `converged` with the last iterate
        returned -- a partially smoothed field is still usable, and `residual`
        says how far it got. This is a pure coordinate move: connectivity,
        every data array, regions and property sets pass through."""

    @server.tool()
    @_forward(_tools.tool_estimate_error)
    def estimate_error(
        input_path: str,
        output_path: str,
        array: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        method: str = "zz",
        marking: str = "none",
        marking_value: float = 0.0,
        output: Optional[str] = None,
        marked: Optional[str] = None,
        overwrite: bool = False,
    ) -> dict:
        """ZZ recovery-based error indicator of a point_data field, plus
        optional marking of cells for refinement.

        A composition of `gradient` (Green-Gauss, cell location) with the
        measure-weighted point<->cell averaging round trip: the indicator is
        sqrt(|measure| * sum((recovered - raw)^2)) per cell, attached as
        `output` (default "error:zz"). marking: none|absolute|fraction|dorfler;
        when not "none" a second Int64 0/1 array `marked` (default
        "error:marked") is attached too, so refine's own `where` selector
        needs no change at all -- the intended use is
        `refine(..., where="error:marked > 0.5")`. marking_value's meaning
        depends on marking: an absolute indicator threshold, a fraction in
        (0, 1] of cells, or the Doerfler bulk fraction theta in (0, 1].
        Cells that cannot be evaluated read NaN in the indicator array and 0
        (never NaN) in the marking array, and are reported in num_skipped
        (excluded from global_error and from num_marked)."""

    @server.tool()
    @_forward(_tools.tool_error_size)
    def error_size(
        input_path: str,
        output_path: str,
        target_error: float,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        error_array: str = "error:zz",
        order: float = 1.0,
        min_size: float = 0.0,
        max_size: float = 0.0,
        output: str = "size",
        overwrite: bool = False,
    ) -> dict:
        """A per-point size field from a per-cell error indicator (roadmap
        6.10.2), by the adaptive procedure of Zienkiewicz and Zhu (1987).

        With N cells evaluated and the target global error target_error, a
        cell of diameter h whose one-component cell_data indicator
        error_array (estimate_error writes "error:zz") is eta asks
        h * (target_error / sqrt(N) / eta)**(1 / order), order being the
        indicator's rate in the cell size (1 for linear elements), clamped to
        [min_size, max_size] (0: no bound); a point takes the least of its
        cells' (NaN where none is evaluated). The sizes are written as
        point_data output (default "size", what mesh_cad_surface's
        background_array reads), so estimate_error -> error_size ->
        mesh_cad_surface(background_path=...) closes the adaptive loop.
        Reports num_cells, num_skipped, num_clamped, min_size and max_size."""

    @server.tool()
    @_forward(_tools.tool_transform)
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

    @server.tool()
    @_forward(_tools.tool_convert_cells)
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

    @server.tool()
    @_forward(_tools.tool_tessellate)
    def tessellate(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        levels: int = 2,
        curved: bool = True,
        fields: bool = True,
        record_stencil: bool = False,
    ) -> dict:
        """Isoparametric subdivision of a mesh's curved cells (quad9, quad8,
        triangle6, tetra10, hexahedron27) onto a levels-divisions-per-axis
        reference lattice mapped through each cell's own shape functions;
        every other cell (linear types, hexahedron20, wedge/pyramid,
        VTK-Lagrange, ragged/polyhedron) passes through unchanged. curved=
        False disables curved handling entirely (a full no-op). fields=True
        interpolates point_data/broadcasts cell_data onto the output.
        Attaches tessellate:source_point/source_cell/sub_index provenance
        (a synthetic point's source_point is -1); record_stencil also
        attaches tessellate:stencil/weights so the tessellation can be
        reconstructed (Tessellation.from_mesh) after a file round trip."""

    @server.tool()
    @_forward(_tools.tool_subdivide)
    def subdivide(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Polyhedrally refine: split every eligible 3D cell into one
        polyhedral child per face, connected to a new interior point. No
        per-type template table is needed -- it handles tabulated types and
        existing polyhedron blocks uniformly. Automatically conforming."""

    @server.tool()
    @_forward(_tools.tool_agglomerate)
    def agglomerate(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        target_group_size: int = 8,
        merge_coplanar_faces: bool = False,
        coplanar_angle: float = 1.0,
        min_sphericity: float = 0.0,
    ) -> dict:
        """Polyhedrally coarsen: merge groups of cells into single larger
        polyhedral cells via greedy seed-and-grow over the shared-face dual.
        Non-volume blocks pass through unchanged; points are never pruned or
        renumbered (clean with remove_orphans=True is the follow-up for
        that). merge_coplanar_faces fuses the coplanar faces two groups (or a
        group and the boundary) share into one polygon, within coplanar_angle
        degrees; min_sphericity (0 = off) refuses an absorption that would make
        a group less round than that. Reports num_faces_merged and
        num_rejected. See doc/agglomerate.md."""

    @server.tool()
    @_forward(_tools.tool_tet_to_poly)
    def tet_to_poly(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        feature_angle: float = 1.0,
        point_data_to_cell_data: bool = False,
        record_vertex_ids: bool = False,
    ) -> dict:
        """Convert a tetrahedral mesh into its polyhedral (median) dual: one
        polyhedron per vertex, bounded by the ring of tetrahedron and triangle
        centroids around each edge and, on the boundary, by the boundary
        triangles' pieces, fused while their normals differ by at most
        feature_angle degrees (so feature edges and corners survive). A vertex
        on an interface between tetra blocks gets one cell per block. Each
        tetra block becomes a polyhedron block, each triangle block a polygon
        block of its dual faces; other blocks are dropped. cell_data is
        volume-weighted (floats) or the largest tet's value (integers);
        point_data_to_cell_data copies each vertex's point data onto its cells;
        record_vertex_ids adds dual:vertex. Reports num_new_points,
        num_boundary_faces and num_split_vertices. See doc/tet_poly.md."""

    @server.tool()
    @_forward(_tools.tool_poly_to_tet)
    def poly_to_tet(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        record_parent_ids: bool = False,
    ) -> dict:
        """Decompose a polyhedral (or mixed) mesh into tetrahedra:
        convert_cells(mode="simplexify"), except that a polyhedron that already
        is a tetrahedron becomes that one tetra with no new points. Other
        polyhedra use the volume-conserving corner-average fan; hexahedra,
        wedges and pyramids the fixed templates. Reports num_exact_tets and
        num_fanned_cells. See doc/tet_poly.md."""

    @server.tool()
    @_forward(_tools.tool_refine)
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
        record_hierarchy: bool = False,
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
        reported in refine:hanging). `record_levels` attaches refine:level.
        `record_hierarchy` attaches refine:cell_id/refine:parent_id -- the
        persistent parent/child hierarchy a multigrid caller resolves across
        the sequence of meshes it keeps; also forces refine:entity to be
        attached even when the closure leaves no hanging node."""

    @server.tool()
    @_forward(_tools.tool_boundary_layer)
    def boundary_layer(
        input_path: str,
        output_path: str,
        first_height: float,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        num_layers: int = 3,
        growth: float = 1.2,
        reverse: bool = False,
        collision_fraction: float = 0.45,
        smooth_iterations: int = 5,
        max_shrink: int = 6,
        collapse: str = "points",
        normal_smooth_iterations: int = 0,
        front_check: bool = True,
    ) -> dict:
        """Grow prismatic boundary layers from a triangle and quad surface
        (roadmap 6.9): wedges from triangles and hexahedra from quads along
        each point's most normal normal (Aubry and Löhner 2009), the first
        first_height thick and each next growth times the one before (against
        the faces' normals with reverse), stopped at collision_fraction of the
        distance to the surface ahead, smoothed smooth_iterations times; the
        points of a cell that would fold are halved max_shrink times and then
        collapsed into pyramids and tetrahedra (collapse="points") or their
        top layer dropped (collapse="layers"), so no cell is written folded.
        normal_smooth_iterations Laplacian passes smooth the growth
        directions; front_check tests the last layer's faces against each
        other and the surface. The wall and the cap come out as faces with
        boundary_layer:layer cell data. Reports num_cells, num_wedges,
        num_hexahedra, num_pyramids, num_tetrahedra, num_shrunk,
        num_collapsed, num_not_visible, num_fewer_layers,
        num_front_crossings, min_scaled_jacobian, min_thickness,
        max_thickness and cap_closed."""

    @server.tool()
    def extrude(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        direction: Optional[List[float]] = None,
        num_layers: int = 1,
        ratio: float = 1.0,
        offsets: Optional[List[float]] = None,
        record_layers: bool = False,
        max_cells: int = 20_000_000,
    ) -> dict:
        """Sweep a surface into a volume, or a curve into a surface, along a
        vector (roadmap 3.1.2): triangles to wedges, quads to hexahedra, lines
        (when there are no 2-D cells) to quads, in num_layers layers spaced
        uniformly, by a geometric ratio, or at explicit offsets (fractions of
        direction, increasing, ending at 1). direction (default [0, 0, 1]) is
        the whole displacement; faces are turned to point along it. The
        output holds the swept cells only, with Side regions bottom, top and
        side, and the input's regions and data carried; record_layers adds
        extrude:layer cell data. Reports num_cells, num_hexahedra, num_wedges,
        num_pyramids, num_tetrahedra, num_quads, num_triangles, num_flipped,
        min_scaled_jacobian and warnings."""
        return _guard(
            _tools.tool_extrude,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            direction=direction if direction is not None else (0.0, 0.0, 1.0),
            num_layers=num_layers,
            ratio=ratio,
            offsets=offsets,
            record_layers=record_layers,
            max_cells=max_cells,
        )

    @server.tool()
    def revolve(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        axis: Optional[List[float]] = None,
        axis_point: Optional[List[float]] = None,
        angle: float = 360.0,
        num_layers: int = 16,
        tolerance: float = 1e-9,
        record_layers: bool = False,
        max_cells: int = 20_000_000,
    ) -> dict:
        """Sweep a surface (or a curve) around an axis (roadmap 3.1.3):
        axis (default [0, 1, 0]) through axis_point (default the origin),
        angle degrees (right-hand rule; 360 closes the turn) in num_layers
        steps. Points within tolerance (a fraction of the bounding-box
        diagonal) of the axis stay put, so the cells around them collapse
        into pyramids, tetrahedra and wedges. Side regions start, end (a
        partial turn) and side. Reports the extrude counts, num_on_axis,
        full, min_scaled_jacobian and warnings."""
        return _guard(
            _tools.tool_revolve,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            axis=axis if axis is not None else (0.0, 1.0, 0.0),
            axis_point=axis_point if axis_point is not None else (0.0, 0.0, 0.0),
            angle=angle,
            num_layers=num_layers,
            tolerance=tolerance,
            record_layers=record_layers,
            max_cells=max_cells,
        )

    @server.tool()
    def mirror(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        origin: Optional[List[float]] = None,
        normal: Optional[List[float]] = None,
        keep_original: bool = False,
        weld: bool = True,
        tolerance: float = 1e-9,
        reflect_vector_data: bool = False,
        plane: Optional[str] = None,
    ) -> dict:
        """Reflect a mesh across the plane through origin (default the
        origin) with normal (default [1, 0, 0]), every cell right side out
        (roadmap 3.2.5; a reflecting transform leaves them inverted). With
        keep_original the input followed by its image, the points on the
        plane shared when weld (default); reflect_vector_data reflects float
        vector and tensor data. plane="auto" mirrors across the mesh's own
        symmetry plane (find_symmetry_plane), reported as plane. Reports
        num_flipped_cells, num_points_on_plane, num_dropped_cells,
        num_dropped_sides and warnings."""
        return _guard(
            _tools.tool_mirror,
            input_path=input_path,
            output_path=output_path,
            input_format=input_format,
            output_format=output_format,
            origin=origin if origin is not None else (0.0, 0.0, 0.0),
            normal=normal if normal is not None else (1.0, 0.0, 0.0),
            keep_original=keep_original,
            weld=weld,
            tolerance=tolerance,
            reflect_vector_data=reflect_vector_data,
            plane=plane,
        )

    @server.tool()
    @_forward(_tools.tool_tetrahedralize)
    def tetrahedralize(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        backend: str = "auto",
        max_size: float = 0.0,
        max_volume: float = 0.0,
        grading: float = 0.3,
        optimize_steps: int = 3,
        max_cells: int = 50_000_000,
        face_array: str = "cad:face",
    ) -> dict:
        """Fill a closed, manifold, consistently oriented triangle surface
        (a mesh_cad_surface result, a boundary layer's cap) with tetrahedra
        (roadmap 6.7). backend "netgen" keeps its triangles and points
        exactly; "octree", the native backend always built, fits the boundary
        within the reported hausdorff_to_surface; "auto" (default) takes
        Netgen when built in, else octree with a warning; "tetgen" and "gmsh"
        run those installed programs on the surface as separate processes,
        its triangles kept (checked). Each shell bounds a
        region (a box inside a box gives two), meshed on its own. max_size and max_volume (0: none) bound the
        tetrahedra, grading (0 to 1) is how fast the size grows inward,
        optimize_steps the backend's volume optimisation passes; face_array
        (cad:face) is carried to the triangles. backend "netgen" in a build
        without it (the meshioplusplus-netgen package, or
        -DMESHIOPLUSPLUS_WITH_NETGEN in a source build) returns an error naming
        both. Reports num_tetrahedra, num_points,
        num_new_points, num_regions, num_boundary_triangles,
        boundary_preserved, hausdorff_to_surface, min_dihedral_deg,
        max_dihedral_deg, volume, backend, warnings, netgen_available and
        backends_available (every backend name besides auto)."""

    @server.tool()
    @_forward(_tools.tool_cut_level_set)
    def cut_level_set(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        level_sets: Optional[List[str]] = None,
        planes: Optional[List[List[float]]] = None,
        selection: str = "",
        include_uncut: bool = True,
        backend: str = "auto",
        degree: int = 2,
        triangulate: bool = False,
        record_parent_ids: bool = True,
    ) -> dict:
        """Cut a mesh by level sets through CutCells (doc/cut.md) and write
        the part the selection picks: the pieces of the cut cells (every cell
        family, quadratic level sets on triangle6/tetra10) and, with
        include_uncut, the whole cells inside. level_sets are scalar
        point_data arrays; planes are six numbers each (origin, then normal),
        named plane0, plane1, ...; selection is an expression over the names
        such as "phi < 0", "phi = 0" (the interface) or "a < 0 and b > 0"
        (empty: "<first> < 0"). backend "lut" gives straight pieces,
        "quadrays" curved Lagrange pieces of degree, "auto" picks quadrays for
        a quadratic level set. Point data is interpolated, cell data and cell
        regions follow each piece's parent (cut:parent_cell). Without
        CutCells it returns an error naming MESHIOPLUSPLUS_WITH_CUTCELLS.
        Reports num_cells, num_cut_cells, num_uncut_cells, num_pieces,
        dimension, measure, selection, backend, warnings and
        cutcells_available."""

    @server.tool()
    @_forward(_tools.tool_cut_quadrature)
    def cut_quadrature(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        level_sets: Optional[List[str]] = None,
        planes: Optional[List[List[float]]] = None,
        selection: str = "",
        order: int = 2,
        include_uncut: bool = True,
        backend: str = "auto",
        triangulate: bool = False,
        record_parent_ids: bool = True,
    ) -> dict:
        """Quadrature rules (cut-cell quadrature for unfitted finite elements)
        of degree order on the part of a mesh that level sets select, through
        CutCells (doc/cut.md). The level sets, selection and backend are
        cut_level_set's. Writes one vertex per quadrature point at its
        physical position with point data quadrature:weight (physical; they
        sum to the part's measure), quadrature:reference (in the parent
        cell), quadrature:rule and quadrature:parent_cell. Reports num_rules,
        num_points, sum_weights, dimension, selection, backend, warnings and
        cutcells_available."""

    @server.tool()
    @_forward(_tools.tool_remesh_quads)
    def remesh_quads(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        target_faces: int = 0,
        backend: str = "auto",
        sharp: bool = False,
        boundary: bool = True,
        adaptive: bool = False,
        dominant: bool = False,
        seed: int = 0,
        max_cells: int = 50_000_000,
    ) -> dict:
        """Remesh a triangle surface into field-aligned quads (roadmap
        7.7.3) through an installed QuadriFlow or Instant Meshes, run as a
        separate process: a sphere comes out all quads, with no pole of poor
        quads as recombining a mesher's triangles leaves. target_faces (0:
        half the triangles) is the number asked for; backend "auto"
        (QuadriFlow, else Instant Meshes), "quadriflow" or "instant_meshes";
        sharp keeps sharp edges, boundary (default true) an open surface's
        boundary, adaptive lets the size follow curvature, dominant allows
        triangles (Instant Meshes), seed repeats a run. No remesher installed
        returns an error naming the programs. Reports num_quads,
        num_triangles, num_points, min_scaled_jacobian,
        mean_scaled_jacobian, hausdorff_to_surface (to the input), backend,
        warnings and backends_available."""

    @server.tool()
    @_forward(_tools.tool_project_to_cad)
    def project_to_cad(
        input_path: str,
        model_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        face_array: str = "cad:face",
        edge_array: str = "cad:edge",
        max_distance: float = 0.0,
        trim_samples: int = 16,
    ) -> dict:
        """Move a mesh's points back onto the CAD model it was made from
        (roadmap 6.10.1): a point on one edge's lines (cad:edge) onto that
        edge, on two or more (a corner) onto the lowest, a point whose cells
        name one face (cad:face) onto its surface by Newton point inversion;
        a point on two faces and no edge stays. model_path is a CAD JSON,
        STEP or IGES file. A point farther than max_distance (0: no limit)
        stays. Reports num_projected, num_on_edges, num_on_curves (of those
        on edges, onto the edge's own 3-D curve), num_corners,
        num_ambiguous, num_too_far, num_failed, num_outside_trim (projections
        outside their face's trims, kept, never clamped), max_move and
        outside_trim."""

    @server.tool()
    @_forward(_tools.tool_undo_green)
    def undo_green(
        coarse_path: str,
        fine_path: str,
        output_path: str,
        coarse_format: Optional[str] = None,
        fine_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Restore fine's transitional (green) cells to their original
        parent, read verbatim from coarse (the mesh a prior refine() call
        with record_hierarchy=True, record_levels=True was run on) --
        undoes refine's known quality-degradation issue with repeated
        selective passes over the same region. Reports
        num_groups_undone/num_cells_removed."""

    @server.tool()
    @_forward(_tools.tool_decimate)
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
        max_deviation: Optional[float] = None,
        coverage: float = 0.999,
        max_faces: Optional[int] = None,
    ) -> dict:
        """Decimate a surface mesh by quadric edge collapse. Give exactly one
        stopping criterion: ratio (fraction of faces to KEEP), target_faces,
        max_error, or max_deviation (the largest two-sided distance from the
        input on a coverage share of samples: the lightest mesh within it);
        max_faces caps any of them (alone, it is the target). Reports
        faces/points removed and rejections, with max_deviation or max_faces
        also deviation, deviation_met and rounds."""

    @server.tool()
    @_forward(_tools.tool_decimate_volume)
    def decimate_volume(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        ratio: Optional[float] = None,
        target_cells: Optional[int] = None,
        max_error: Optional[float] = None,
        placement: str = "optimal",
        preserve_boundary: bool = False,
        preserve_features: bool = True,
        feature_angle: float = 30.0,
    ) -> dict:
        """Decimate a tetrahedral mesh by quadric-error tet-edge collapse. Give
        exactly one stopping criterion: ratio (fraction of tets to KEEP),
        target_cells, or max_error. Boundary vertices participate by default
        (preserve_boundary=False, unlike surface decimate). Reports
        tets/points removed and rejections."""

    @server.tool()
    @_forward(_tools.tool_adapt)
    def adapt(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        hsiz: Optional[float] = None,
        size: Optional[str] = None,
        metric: Optional[str] = None,
        hmin: Optional[float] = None,
        hmax: Optional[float] = None,
        hausd: Optional[float] = None,
        hgrad: Optional[float] = None,
        hgradreq: Optional[float] = None,
        angle_detection: Optional[float] = 45.0,
        optim: bool = False,
        noinsert: bool = False,
        noswap: bool = False,
        nomove: bool = False,
        nosurf: bool = False,
        nreg: bool = False,
        opnbdy: bool = False,
        frozen: Optional[str] = None,
        required: Optional[List[str]] = None,
        local: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> dict:
        """Remesh a simplicial mesh with MMG (needs a build with
        MESHIOPLUSPLUS_WITH_MMG; otherwise returns error_type
        NotImplementedError and mmg_available False). mmg3d for
        tetrahedra, mmg2d for a planar triangle mesh, mmgs for another
        triangle surface. Give at most one of hsiz (constant edge size),
        size (point_data name of a size field) and metric (point_data name
        of a six-component, three in 2-D, metric in Medit .sol order); none
        lets MMG size from the mesh. hmin/hmax bound the edge sizes, hausd
        bounds how far the boundary moves, hgrad/hgradreq set the
        gradation, angle_detection (45 degrees; None turns it off) finds
        ridges. optim/noinsert/noswap/nomove/nosurf/nreg/opnbdy are MMG's
        switches. frozen names a point region MMG must keep; required
        names regions whose entities it must keep; local maps a region to
        {"hmin", "hmax", "hausd"}. Cell and side regions and a medit:ref
        array come back under their names; other point/cell data is
        dropped (use interpolate), field_data is carried. Reports status
        ("success" or "low_failure"), library, num_points_in/out,
        num_elements_in/out, quality (aspect ratio min/mean/max/histogram)
        and warnings."""

    @server.tool()
    @_forward(_tools.tool_discretize_levelset)
    def discretize_levelset(
        input_path: str,
        output_path: str,
        level_set: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        isovalue: float = 0.0,
        materials: Optional[List[List[int]]] = None,
        remove_small_components: Optional[float] = None,
        hmin: Optional[float] = None,
        hmax: Optional[float] = None,
        hausd: Optional[float] = None,
        hgrad: Optional[float] = None,
        hgradreq: Optional[float] = None,
        angle_detection: Optional[float] = 45.0,
        optim: bool = False,
        noinsert: bool = False,
        noswap: bool = False,
        nomove: bool = False,
        nosurf: bool = False,
        nreg: bool = False,
        opnbdy: bool = False,
        frozen: Optional[str] = None,
        required: Optional[List[str]] = None,
        local: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> dict:
        """Embed the isovalue of the point_data level set `level_set` (a
        signed distance, say) conformally in a simplicial mesh and remesh
        both sides with MMG (needs a build with MESHIOPLUSPLUS_WITH_MMG;
        otherwise returns error_type NotImplementedError and mmg_available
        False): mmg3d, mmg2d or mmgs by the mesh, as adapt. The elements
        carry cell_data medit:ref (3 inside, below the isovalue; 2 outside)
        and the cell regions inside/outside; the interface is the side region
        interface, on the inside elements. materials lists [ref, split,
        ref_inside, ref_outside] over the input's medit:ref (MMG's
        multi-material mode: split 0 keeps the cells whole);
        remove_small_components is MMG's rmc volume fraction. The other MMG
        parameters are adapt's. The input's side regions come back, its cell
        regions do not. Reports as adapt."""

    @server.tool()
    @_forward(_tools.tool_mesh_domain)
    def mesh_domain(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        hsiz: Optional[float] = None,
        size: Optional[str] = None,
        holes: Optional[List[List[float]]] = None,
        hmin: Optional[float] = None,
        hmax: Optional[float] = None,
        hausd: Optional[float] = None,
        hgrad: Optional[float] = None,
        hgradreq: Optional[float] = None,
        angle_detection: Optional[float] = 45.0,
        optim: bool = False,
        noinsert: bool = False,
        noswap: bool = False,
        nomove: bool = False,
        nosurf: bool = False,
        nreg: bool = False,
        opnbdy: bool = False,
        frozen: Optional[str] = None,
        required: Optional[List[str]] = None,
        local: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> dict:
        """Triangulate the planar domain the line cells of a mesh bound with
        MMG's MMG2D_mmg2dmesh (needs a build with MESHIOPLUSPLUS_WITH_MMG;
        otherwise returns error_type NotImplementedError and mmg_available
        False). hsiz is a constant edge size, size a point_data size field on
        the boundary points. MMG meshes every region the lines enclose, each
        a cell region subdomain:<k>; it has no notion of a hole, so holes
        lists one [x, y] point inside each region to leave out. The lines'
        cell regions also come back as side regions of the triangles. The
        other MMG parameters are adapt's. Reports as adapt plus
        num_subdomains and num_holes."""

    @server.tool()
    @_forward(_tools.tool_move)
    def move(
        input_path: str,
        output_path: str,
        displacement: str,
        moving: List[str],
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        lagrangian: int = 1,
        hmin: Optional[float] = None,
        hmax: Optional[float] = None,
        hausd: Optional[float] = None,
        hgrad: Optional[float] = None,
        hgradreq: Optional[float] = None,
        angle_detection: Optional[float] = 45.0,
        optim: bool = False,
        noinsert: bool = False,
        noswap: bool = False,
        nomove: bool = False,
        nosurf: bool = False,
        nreg: bool = False,
        opnbdy: bool = False,
        frozen: Optional[str] = None,
        required: Optional[List[str]] = None,
        local: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> dict:
        """Move a tetrahedral or planar triangle mesh by the point_data
        vector `displacement` of the boundary facets of the regions `moving`
        (cell or side regions), remeshing as it deforms with MMG's
        Lagrangian mode (MMG3D_mmg3dmov, MMG2D_mmg2dmov). MMG extends the
        displacement inward by linear elasticity (the ELAS library) and
        keeps the rest of the boundary in place, so the displacement given
        elsewhere is not read. lagrangian is MMG's lag: 0 moves only, 1
        also swaps, 2 also inserts and collapses. Needs a build with
        MESHIOPLUSPLUS_WITH_MMG and an MMG built with ELAS; otherwise returns
        error_type NotImplementedError with mmg_available and
        mmg_elas_available. The other MMG parameters are adapt's. Reports as
        adapt."""

    @server.tool()
    @_forward(_tools.tool_metric_from_hessian)
    def metric_from_hessian(
        input_path: str,
        output_path: str,
        hessian: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        error: float = 1e-2,
        hmin: float = 1e-3,
        hmax: float = 1.0,
        output: str = "metric",
    ) -> dict:
        """The anisotropic metric of the point_data Hessian `hessian`
        ((n, 9), or (n, 4) in 2-D, row-major as the hessian tool writes)
        for the adapt tool: M = c/error |H| with c = (1/2)(d/(d+1))^2, its
        eigenvalues clamped so the sizes lie in [hmin, hmax], written as
        point_data `output` (six components, three in 2-D, in the
        lower-triangle order of a Medit .sol). Needs no MMG."""

    @server.tool()
    def remesh(
        input_path: str,
        output_path: str,
        num_clusters: int,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        subdivide: Optional[int] = None,
        subsample_ratio: float = 10.0,
        max_subdivide: int = 4,
        max_iterations: int = 100,
        max_repair_passes: int = 10,
        metric: str = "isotropic",
        gradation: float = 0.0,
        preserve_boundary: bool = True,
        max_anisotropy: float = 4.0,
    ) -> dict:
        """Replace a surface mesh's triangulation with a new,
        near-uniformly-sized, well-shaped one at num_clusters vertices
        (approximated centroidal Voronoi diagram clustering). Unlike every
        other resolution-changing tool, the output has NO correspondence to
        the input -- point_data/cell_data/named regions are dropped,
        field_data is carried; use interpolate/conservative_interpolate to
        transfer a field onto the result. metric is "isotropic" (default,
        fast, rounds sharp features), "quadric" (Garland-Heckbert error,
        preserves sharp edges/corners), or "anisotropic" (clusters shaped by
        a local curvature tensor, elongated along low-curvature directions --
        see max_anisotropy). gradation is the curvature-gradation exponent
        gamma in the item weight area * kappa**gamma (0.0 default disables
        gradation). preserve_boundary (default True) seeds and pins the
        input's open boundary and emits a line dual cell along it, a no-op
        on a closed mesh. max_anisotropy (default 4.0) is, under
        metric="anisotropic", the maximum ratio between the two in-plane
        target edge lengths a curvature tensor may request; 1.0 recovers
        the isotropic shape; an error to set away from the default under
        any other metric."""
        return _guard(
            _tools.tool_remesh,
            input_path=input_path,
            output_path=output_path,
            num_clusters=num_clusters,
            input_format=input_format,
            output_format=output_format,
            subdivide=subdivide,
            subsample_ratio=subsample_ratio,
            max_subdivide=max_subdivide,
            max_iterations=max_iterations,
            max_repair_passes=max_repair_passes,
            metric=metric,
            max_anisotropy=max_anisotropy,
            gradation=gradation,
            preserve_boundary=preserve_boundary,
        )

    @server.tool()
    @_forward(_tools.tool_remesh_volume)
    def remesh_volume(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        resolution: Optional[List[int]] = None,
        cell_size: Optional[float] = None,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.1,
        max_cells: int = 20000000,
        max_tets: int = 20000000,
        warp_fraction: float = 0.35,
        watertight_check: str = "warn",
    ) -> dict:
        """Retetrahedralize a volume mesh (or closed surface) at a chosen
        resolution by isosurface stuffing -- the volumetric sibling of remesh.
        Give exactly one of resolution (nx, ny, nz) or cell_size. Like remesh,
        the output has NO correspondence to the input -- point_data/cell_data/
        named regions are dropped, field_data is carried. warp_fraction (0-1,
        default 0.35) trades boundary tet quality for a small, measured chance
        of non-manifold boundary edges (reported as num_non_manifold_edges);
        0 gives an exactly watertight but lower-quality boundary."""

    @server.tool()
    @_forward(_tools.tool_optimize_volume)
    def optimize_volume(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        max_iterations: int = 10,
        relocate: bool = True,
        flip: bool = True,
        preserve_boundary: bool = True,
        min_improvement: float = 1e-6,
        objective: str = "worst",
        max_ring_size: int = 7,
    ) -> dict:
        """ODT-remesh a tetrahedral mesh: raise its worst element quality by
        relocating vertices AND changing connectivity, predicate-free.
        objective="worst" (the default) raises the least dihedral angle and
        never lowers it -- guarded ODT relocation, 2-3 flips, edge removal and
        multi-face removal over rings of up to max_ring_size tets, and
        worst-cell smoothing; "energy" is the ODT energy with scaled-Jacobian
        flips, which can lower it. Reports the flip counts, num_edge_removals,
        num_multi_face_removals, num_worst_smoothed, and the least scaled
        Jacobian and dihedral angle (min_dihedral_before/after) in and out.
        The genuine "ODT remeshing" sibling of remesh_volume (which generates a
        fresh lattice mesh) and of smooth method="odt" (which only moves points
        on fixed connectivity). Tet-only and C++-core-only. The point set is
        invariant, so point_data and named Point regions carry; cell_data and
        Cell/Side regions are dropped. With preserve_boundary the boundary
        surface is exactly preserved."""

    @server.tool()
    @_forward(_tools.tool_smooth)
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
        """Smooth point coordinates (taubin, the feature-preserving default;
        laplacian; or odt, optimal-Delaunay-triangulation smoothing, tet-only
        and C++-core-only with no pure-Python fallback). Pure coordinate move;
        boundary and feature nodes pinned by default. lambda_factor < 0 means
        the method's own default."""

    @server.tool()
    @_forward(_tools.tool_merge)
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

    @server.tool()
    @_forward(_tools.tool_split)
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

    @server.tool()
    @_forward(_tools.tool_partition)
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
        kahip). Writes one file per part and returns their paths; a
        name_template ending in .pvtu/.pvtp with no {part} writes one
        parallel index over every part instead (halo layers kept as
        vtkGhostType) and also returns its path as `index`."""

    @server.tool()
    @_forward(_tools.tool_interpolate)
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

    @server.tool()
    @_forward(_tools.tool_conservative_interpolate)
    def conservative_interpolate(
        source_path: str,
        target_path: str,
        output_path: str,
        source_format: Optional[str] = None,
        target_format: Optional[str] = None,
        output_format: Optional[str] = None,
        arrays: Optional[List[str]] = None,
        default_value: float = 0.0,
        on_conflict: str = "error",
    ) -> dict:
        """Mass-preservingly (overlap-measure weighted) sample the source
        mesh's data arrays onto the target mesh's geometry; unlike
        interpolate, conserves sum(value * measure) over the shared region."""


# --------------------------------------------------------------------------- #
# Data operations                                                             #
# --------------------------------------------------------------------------- #
def _register_data(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_sets_data)
    def sets_data(
        input_path: str,
        output_path: str,
        direction: str = "sets_to_data",
        location: str = "cell",
        key: Optional[str] = None,
        data_name: Optional[str] = None,
        join_char: str = "-",
        order: Optional[List[str]] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
    ) -> dict:
        """Convert point/cell sets and scalar integer labels. direction is
        sets_to_data|data_to_sets; location is point|cell. Labels are zero-based,
        later overlaps win, uncovered rows are -1. data_to_sets requires key;
        sets_to_data can choose data_name, join_char and explicit set order.
        The converted sets or source field are removed; side regions survive."""

    @server.tool()
    @_forward(_tools.tool_data_manage)
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

    @server.tool()
    @_forward(_tools.tool_data_convert)
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

    @server.tool()
    @_forward(_tools.tool_data_calc)
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

    @server.tool()
    @_forward(_tools.tool_data_condition)
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

    @server.tool()
    @_forward(_tools.tool_tensor_invariants)
    def tensor_invariants(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        location: str = "point",
        arrays: Optional[List[str]] = None,
        outputs: Optional[List[str]] = None,
        prefix: str = "",
        suffix: str = "",
        overwrite: bool = True,
    ) -> dict:
        """von Mises / principal / hydrostatic / deviatoric of a symmetric
        (6-component) or general 3x3 (9-component) tensor array. outputs: any
        of mises | principal | hydrostatic | deviatoric (default: all four)."""


# --------------------------------------------------------------------------- #
# Dataset manifests (doc/datasets.md)                                         #
# --------------------------------------------------------------------------- #
def _register_dataset(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_dataset_add)
    def dataset_add(
        manifest_path: str,
        input_pattern: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        entry_id: Optional[str] = None,
        input_format: Optional[str] = None,
        time_from: Optional[str] = None,
        times: Optional[List[float]] = None,
        sort: bool = False,
        target_pattern: Optional[str] = None,
        target_paths: Optional[List[str]] = None,
        target_format: Optional[str] = None,
        target_time_from: Optional[str] = None,
        target_times: Optional[List[float]] = None,
        target_sort: bool = False,
        split: Optional[str] = None,
        tags: Optional[List[str]] = None,
        group: Optional[str] = None,
        notes: Optional[str] = None,
        metadata: Optional[dict] = None,
        thumbnail: Optional[str] = None,
        render_thumbnail: bool = False,
        thumbnail_width: int = 256,
        thumbnail_height: int = 192,
    ) -> dict:
        """Add a case to a dataset manifest JSON (created if absent). Give
        exactly one of input_pattern (a glob) or input_paths; the source is
        validated now and stored relative to the manifest's directory.
        The optional target_* family records a paired coarse/fine series for a
        superresolution dataset: the two must have the same number of steps at
        the same instants, checked here. Leave it out for the ordinary case --
        an entry without a target is self-supervised, one mesh supplying both
        sides. Optional curation: split, tags, group, notes, metadata, and a
        picture of the case: thumbnail (an existing PNG under the manifest's
        directory) or render_thumbnail=true, which draws the first step with
        the software rasterizer to thumbnails/<id>.png beside the manifest."""

    @server.tool()
    @_forward(_tools.tool_dataset_list)
    def dataset_list(
        manifest_path: str,
        split: Optional[str] = None,
        tags: Optional[List[str]] = None,
        group: Optional[str] = None,
        resolve: bool = False,
    ) -> dict:
        """List a dataset manifest's entries, optionally filtered by split /
        tags (must carry all) / group (path or descendant). resolve=true also
        expands each entry's file/step/time plan (no mesh is read)."""

    @server.tool()
    @_forward(_tools.tool_dataset_update)
    def dataset_update(
        manifest_path: str,
        entry_ids: Optional[List[str]] = None,
        all_entries: bool = False,
        split: Optional[str] = None,
        assign_splits: Optional[dict] = None,
        seed: int = 0,
        by_group: bool = False,
        add_tags: Optional[List[str]] = None,
        remove_tags: Optional[List[str]] = None,
        group: Optional[str] = None,
        notes: Optional[str] = None,
        metadata: Optional[dict] = None,
        drop_metadata: Optional[List[str]] = None,
        thumbnail: Optional[str] = None,
        render_thumbnail: bool = False,
        thumbnail_step: int = 0,
        thumbnail_width: int = 256,
        thumbnail_height: int = 192,
    ) -> dict:
        """Curate a dataset manifest: set a split on selected entries, assign
        splits by fractions (assign_splits={"train": 0.8, ...}, deterministic
        via seed; by_group keeps groups together), add/remove tags, or set one
        entry's group/notes/metadata. thumbnail records an existing PNG under
        the manifest's directory as one entry's picture ("none" clears it);
        render_thumbnail=true draws each selected entry's step
        thumbnail_step to thumbnails/<id>.png beside the manifest."""

    @server.tool()
    @_forward(_tools.tool_dataset_find)
    def dataset_find(root_dir: str = ".", max_depth: int = 2) -> dict:
        """Find dataset manifests: every *.json at most max_depth levels
        below root_dir that parses as a DatasetManifest, with its name,
        entry count, splits, modification time and SHA-256 content hash."""

    @server.tool()
    @_forward(_tools.tool_dataset_health)
    def dataset_health(
        manifest_path: str,
        split: Optional[str] = None,
        entry_ids: Optional[List[str]] = None,
        quality: bool = True,
        all_steps: bool = False,
    ) -> dict:
        """Scan a dataset manifest's entries (optionally one split / given
        ids) and report their health: per entry the step count, NaN/Inf
        counts over data arrays, inverted/degenerate cells, worst scaled
        Jacobian and arrays present; per manifest the split balance, totals,
        fields missing across entries and the bad entries. Reads one mesh at
        a time (step 0, or every step with all_steps)."""


# --------------------------------------------------------------------------- #
# Training jobs (doc/dashboard.md)                                            #
# --------------------------------------------------------------------------- #
def _register_training(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_train_defaults)
    def train_defaults(
        manifest_path: str,
        fields: Optional[List[str]] = None,
        target_fields: Optional[List[str]] = None,
    ) -> dict:
        """What a training launch needs: the point/cell data arrays the
        manifest's first entry carries, its splits, whether the frameworks
        are installed, and a complete default training spec."""

    @server.tool()
    @_forward(_tools.tool_train_start)
    def train_start(
        manifest_path: str,
        fields: List[str],
        target_fields: List[str],
        train_split: str = "train",
        valid_split: str = "valid",
        epochs: int = 100,
        batch_size: int = 8,
        learning_rate: float = 1e-3,
        seed: int = 0,
        model_name: str = "meshgraphnet",
        processor_size: int = 8,
        hidden_dim: int = 64,
        aggregation: str = "sum",
        scaling_factor: int = 2,
        conv_layer_size: int = 32,
        resid_blocks: int = 8,
        resolution: Optional[List[int]] = None,
        cell_size: Optional[float] = None,
        bounds: Optional[List[float]] = None,
        padding: float = 0.0,
        padding_relative: float = 0.0,
        extrapolate: bool = False,
        fill_value: float = 0.0,
        squeeze: Optional[int] = None,
        squeeze_index: Optional[int] = None,
        latent_channels: int = 32,
        num_fno_layers: int = 4,
        num_fno_modes: int = 16,
        spectral_padding: int = 8,
        patch_size: Optional[List[int]] = None,
        embed_dim: int = 256,
        depth: int = 4,
        num_blocks: int = 16,
        parameters: Optional[List[str]] = None,
        trunk: str = "points",
        trunk_count: Optional[int] = None,
        trunk_method: str = "farthest",
        trunk_seed: int = 0,
        branch_layers: int = 4,
        branch_layer_size: int = 128,
        trunk_layers: int = 4,
        trunk_layer_size: int = 128,
        width: int = 64,
        regions: bool = False,
        kind: str = "node",
        undirected: bool = True,
        edge_features: bool = True,
        float32: bool = True,
        target_offset: int = 0,
        target_delta: bool = False,
        checkpoint_every: int = 10,
        device: str = "auto",
        notes: Optional[str] = None,
        tags: Optional[List[str]] = None,
    ) -> dict:
        """Start a PhysicsNeMo training run on a manifest split (fields ->
        target_fields) as a background job under the runs directory; returns
        the job id and initial status. model_name picks the family --
        meshgraphnet (graph), srresnet / fno / afno (grid; fno is 2-D when
        squeeze names a world axis, afno requires it) or deeponet (parameters
        in, field out: no fields, the per-entry Metadata keys in parameters)
        -- and each reads only its own hyperparameters. Needs
        nvidia-physicsnemo, plus torch_geometric for meshgraphnet (no pip
        extra); a missing framework is a named error."""

    @server.tool()
    @_forward(_tools.tool_train_status)
    def train_status(job_id: str) -> dict:
        """A training job's status (running/finished/failed/stopped), epoch
        progress, best validation loss, ETA and last metrics row."""

    @server.tool()
    @_forward(_tools.tool_train_list)
    def train_list(
        status: Optional[str] = None, manifest_path: Optional[str] = None
    ) -> dict:
        """Every training job (newest first) with its hyperparameters and
        final/best losses; optionally filtered by status or manifest."""

    @server.tool()
    @_forward(_tools.tool_train_stop)
    def train_stop(job_id: str, grace_seconds: float = 10.0) -> dict:
        """Stop a job: SIGTERM lets the trainer finish its epoch and write
        final.mdlus; SIGKILL after grace_seconds."""

    @server.tool()
    @_forward(_tools.tool_train_log)
    def train_log(job_id: str, offset: int = 0, max_bytes: int = 65536) -> dict:
        """A window of a job's stdout/stderr from a byte offset; poll with
        next_offset to tail it."""

    @server.tool()
    @_forward(_tools.tool_train_metrics)
    def train_metrics(job_id: str, since_epoch: int = 0) -> dict:
        """A job's per-epoch metrics rows (train/valid loss, lr, timing)."""

    @server.tool()
    @_forward(_tools.tool_train_checkpoints)
    def train_checkpoints(job_id: str) -> dict:
        """A job's .mdlus checkpoints (periodic/best/final) with epoch,
        validation loss, size and which is marked best."""

    @server.tool()
    @_forward(_tools.tool_train_mark_best)
    def train_mark_best(job_id: str, checkpoint: str) -> dict:
        """Mark one of a job's checkpoints as best (copied to best.mdlus)."""

    @server.tool()
    @_forward(_tools.tool_train_predict)
    def train_predict(
        manifest_path: str,
        job_id: Optional[str] = None,
        checkpoint: Optional[str] = None,
        entry_ids: Optional[List[str]] = None,
        split: Optional[str] = "test",
        step: int = 0,
        output_dir: Optional[str] = None,
    ) -> dict:
        """Predict over a manifest split with a job's (best or named)
        checkpoint or an explicit .mdlus, writing <column>_pred /
        <column>_error data arrays into output_dir/<entry_id>.vtu; returns
        per-entry RMSE. Needs the frameworks."""

    @server.tool()
    @_forward(_tools.tool_guard_fit)
    def guard_fit(
        manifest_path: str,
        output_path: str,
        split: Optional[str] = "train",
        margin: float = 1.5,
        quality: bool = True,
    ) -> dict:
        """Fit a geometry guardrail over a manifest split and write it as JSON.

        A trained surrogate answers any mesh it is given, and the answer for a
        part unlike anything it saw is finite, plausible and wrong. This
        describes the shapes in a split — extents, centroid, area, volume,
        counts, surface area and quality summaries — so a new mesh can be
        scored against them. The descriptors are deliberately NOT invariant:
        a scaled part is a different part, and a model trained on brackets 10
        cm across has learnt physics at that scale. The threshold is margin
        times the worst training score."""

    @server.tool()
    @_forward(_tools.tool_guard_check)
    def guard_check(
        input_path: str,
        guard_path: Optional[str] = None,
        input_format: Optional[str] = None,
        top: int = 3,
    ) -> dict:
        """Describe one mesh's shape, and score it against a guardrail if given.

        Returns the raw descriptors always; with guard_path (a fitted guard, or
        a model card carrying one) it adds the score, the threshold, a verdict
        of 'in' or 'out', and the top descriptors that put it there — so a flag
        is actionable rather than a bare number. Advisory: nothing is refused."""

    @server.tool()
    @_forward(_tools.tool_predict_file)
    def predict_file(
        checkpoint: str,
        input_path: str,
        output_path: str,
        time_step: Optional[int] = None,
        target_path: Optional[str] = None,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        device: str = "auto",
        parameters: Optional[dict] = None,
    ) -> dict:
        """Predict with a trained .mdlus checkpoint on ONE mesh file — no
        manifest, no split, no entry, for a mesh that was never catalogued.
        Everything the prediction needs comes from the checkpoint's own model
        card: which model family wrote it, the sample options, the read
        options, the column contract and the normalization. time_step picks a
        step of a multi-step input; target_path supplies the paired mesh a
        t->t+n or coarse/fine checkpoint compares against; parameters supplies
        a deeponet checkpoint's per-case inputs (an object of the Metadata
        keys it was trained on). A file carrying no truth predicts anyway,
        with rmse/max_error reported as null rather than measured against
        itself. Needs the frameworks."""


# --------------------------------------------------------------------------- #
# Gated tools (optional extras)                                               #
# --------------------------------------------------------------------------- #
def _register_gated(server: FastMCP) -> None:
    @server.tool()
    @_forward(_tools.tool_fit_cad)
    def fit_cad(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        mesh_output_path: Optional[str] = None,
        curvature_tolerance: float = 0.0,
        min_patch_cells: int = 1,
        fit_tolerance: Optional[float] = None,
        tolerance: Optional[float] = None,
        exact_edges: bool = True,
        types: Optional[List[str]] = None,
        mesh_deflection: Optional[float] = None,
        max_demotion_rounds: int = 8,
        freeform: bool = True,
        freeform_degree: int = 3,
        freeform_max_control: int = 64,
        cad_json_output_path: Optional[str] = None,
        freeform_split: bool = True,
        fillets: bool = True,
        threads: bool = True,
        segment_method: str = "primitive",
        segment_tolerance: float = 1e-3,
        merge_faces: bool = True,
        backend: str = "auto",
        refine_segmentation: bool = False,
        max_refine_iterations: int = 10,
        refine_epsilon: float = 1e-3,
        normal_weight: float = 1.0,
        allow_merge: bool = True,
        allow_split: bool = True,
        ransac_seed: int = 0x6D657368696F2B2B,
        refine_tolerance: float = 1e-3,
        relative_tolerance: bool = True,
        freeform_general: bool = False,
        freeform_tolerance: float = 0.0,
        freeform_smoothing: float = 1e-12,
        max_freeform_poles: int = 40,
        max_freeform_iterations: int = 30,
        max_freeform_splits: int = 2,
        sew_tolerance: float = 0.0,
        extend_surfaces: bool = False,
        snap_junctions: bool = False,
        face_deviation: bool = False,
        freeform_mode: str = "legacy",
        freeform_parameterization: str = "auto",
        extraordinary_vertex: str = "g1x",
        local_refine: str = "knots",
        freeform_g1_tolerance_deg: float = 1.0,
        freeform_max_control_points: int = 400,
        layout_target_patch_size: float = 0.1,
        layout_max_valence_deviation: int = 2,
        resegment_rounds: int = 2,
        layout_snap_to_features: bool = True,
        validate_by_surface_sampling: bool = True,
        constrain_primitives: bool = False,
    ) -> dict:
        """Turn a closed triangle skin into a CAD solid with analytic faces,
        written to output_path as AP214 STEP -- the last step of skin to CAD
        after segment and fit_primitives (run here first when the skin has no
        segment:patch: segment_method "primitive" (the default) ends a patch
        where the plane, sphere, cylinder, cone or torus its facets lie on
        changes, at segment_tolerance of the diagonal, and fits at that
        tolerance unless fit_tolerance (relative) is given; "angle" segments
        by angles with curvature_tolerance and min_patch_cells, which a noisy
        scan wants in the tens).

        Each accepted fit becomes a plane, sphere, cylinder, cone or torus face
        (types narrows which); a patch no primitive explains becomes, with
        fillets, a rolling-ball fillet face when it runs between two analytic
        faces (of one radius or a radius varying along it), with threads, a
        helical thread face (about a neighbouring cylinder's or cone's axis or
        the screw the patch fits, straight flank or rounded profile), and
        otherwise, with freeform, a B-spline face (degree up to
        freeform_degree, fitted on at most freeform_max_control grid points
        per direction; with freeform_split, a patch folding over every
        flattening is split into pieces first; with fillets, a patch between
        two fitted faces, one of them freeform, then becomes a rolling-ball
        blend); with fillets, a patch no rolling ball or thread explains
        between two analytic faces becomes a chamfer or a conic blend, and a
        corner where three or more blends meet a vertex blend; with threads, a
        flank of varying pitch or cut on a taper is a helical face too; the
        rest coplanar polygons;
        OpenCASCADE sews and checks the solid, a face that fails or lies
        further than tolerance from the skin is demoted to facets, and with
        merge_faces faces on one surface merge. The report
        gives backend, valid, the volume ratio against the skin, max_deviation
        (Hausdorff), faces by type (num_bspline the freeform ones, num_fillet,
        num_helical, num_chamfer, num_vertex_blend; of those num_variable_fillet,
        num_blend, num_conic_blend, num_rounded_thread, num_tapered_thread and
        num_variable_pitch_thread), num_split_patches, num_demoted and exact/approximated
        edges; mesh_output_path also writes the
        re-meshed solid and cad_json_output_path the solid as Kratos CAD JSON
        (NURBS faces, trimming loops, edges). backend "occt" needs a core
        built with OpenCASCADE (-DMESHIOPLUSPLUS_WITH_OCCT=ON; without it the
        error names the flag); "native" needs no CAD kernel: the fitted
        patches become analytic STEP faces bounded by line, circle and
        polyline edges, the rest coplanar polygons, no freeform or feature
        faces, NaN volumes; "auto" (the default) is "occt" when the core
        carries it (occt_available), else "native". refine_segmentation (off
        by default) first refines the patches as the refine_patches tool does
        (max_refine_iterations, refine_epsilon, normal_weight, allow_merge,
        allow_split, ransac_seed, refine_tolerance relative to the diagonal
        unless relative_tolerance is false), so more of the skin becomes
        analytic faces; the report then adds refined, refinement_kept (whether
        the solid was built from the refined patches; OpenCASCADE also builds
        it from the patches as given and keeps the better), refine_iterations,
        num_refine_splits, num_refine_merges, initial_explained_area_fraction
        and explained_area_fraction. freeform_general (off by default) tries a
        patch the plane and unrolled-cylinder freeform fit leaves unexplained
        as a general B-spline face (a disc parametrisation and an iterative
        fit to freeform_tolerance, 0 a quarter of the tolerance, with
        thin-plate weight freeform_smoothing, at most max_freeform_poles per
        direction and max_freeform_iterations, halved up to
        max_freeform_splits times) and adds num_freeform_general,
        num_freeform_refinements and avg_freeform_poles to the report.
        OpenCASCADE only: sew_tolerance sews at that tolerance (0, the
        default, at tolerance); extend_surfaces extends freeform surfaces past
        their patches before the edges are built (adds
        num_extended_edges); snap_junctions places the vertices where
        analytic faces meet by a least-squares solve on their tangent planes
        (adds num_snapped_junctions); face_deviation measures each face's
        two-sided distance to the skin as fit_cad:face_hausdorff and
        fit_cad:face_rms field data and adds hausdorff_max and hausdorff_rms
        (null when not measured). The report always holds closed (every shell
        edge borders two faces). freeform_mode picks how the patches no
        primitive or feature fit explains are treated: "legacy" (the default,
        every output as before), "adaptive" (the general fit with the
        freeform_parameterization ladder "auto" | "projection" | "lscm" |
        "harmonic", surface sampling when validate_by_surface_sampling, at most
        freeform_max_control_points poles per face, and an escalation ladder
        before demotion with resegment_rounds re-segmentation rounds),
        "network" (a quad-layout B-spline patch network per residual region:
        layout_target_patch_size of the diagonal, layout_max_valence_deviation,
        layout_snap_to_features, G1 within freeform_g1_tolerance_deg, and
        extraordinary_vertex "g0" | "g1" | "g1x" | "g1xx") or "auto" (adaptive,
        then a network for what still fails); local_refine "knots" or "thb"
        (experimental). Outside "legacy" the report adds freeform_mode,
        brep_checked, facet_area_fraction (the share of the area left on facet
        faces), num_freeform_adaptive, num_network_patches, num_network_edges,
        num_extraordinary_vertices, num_resegmented, num_escalations,
        max_g1_angle_deg, max_freeform_distance, rms_freeform_distance and
        faces (per patch tried: param_method, num_ctrl, iterations,
        refinements, distances, escalation_steps, demotion_reason).
        constrain_primitives (off) fits related primitives together on either
        backend: coaxial cylinders, cones, tori and spheres share one axis,
        parallel and perpendicular planes and nearly equal radii are made
        exact, each change kept only where its patch still passes its fit;
        adds num_constrained_primitives and num_constraint_groups."""

    @server.tool()
    @_forward(_tools.tool_data_export)
    def data_export(
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        location: str = "point",
    ) -> dict:
        """Export data arrays to a Parquet table (location: point | cell).
        Needs the [arrow] extra; a missing install returns a named error."""

    @server.tool()
    @_forward(_tools.tool_export_dataset)
    def export_dataset(
        input_pattern: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        output_path: Optional[str] = None,
        input_format: Optional[str] = None,
        location: str = "point",
        dataset_format: str = "parquet",
        mesh_id: str = "stem",
    ) -> dict:
        """Export a SET of meshes as one dataset keyed by mesh_id — a
        hive-partitioned Parquet directory (dataset_format: parquet, needs the
        [arrow] extra) or chunked zarr/hdf5 groups ([zarr] / h5py). Give
        exactly one of input_pattern (a glob) or input_paths; mesh_id: stem |
        index. Returns the dataset manifest."""

    @server.tool()
    @_forward(_tools.tool_export_cae)
    def export_cae(
        output_dir: str,
        input_pattern: Optional[str] = None,
        input_paths: Optional[List[str]] = None,
        input_format: Optional[str] = None,
        surface_fields: Optional[List[str]] = None,
        volume_fields: Optional[List[str]] = None,
        global_params: Optional[Dict[str, float]] = None,
        global_params_reference: Optional[Dict[str, float]] = None,
        global_params_order: Optional[List[str]] = None,
        name_template: str = "case_{index}.npz",
    ) -> dict:
        """Export a SET of meshes as one .npz per case in the CAE sample
        layout PhysicsNeMo's DoMINO/Transolver datapipes read: the
        triangulated skin with normals and areas, the volume's nodes, the
        named field blocks and the case's global parameters. Give exactly one
        of input_pattern (a glob) or input_paths. Returns the files written
        and the keys of the first one."""

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
        by a data array). Uses polyscope when the [viewer] extra is installed
        (size and transparent apply to it); otherwise the software rasterizer
        draws the picture (see render_mesh for its options)."""
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
    @server.prompt(
        name="measure_scan", description="Measure a scan against its drawing"
    )
    def measure_scan_prompt(
        scan_path: str = "scan.stl",
        spec_path: str = "inspection.json",
        part_class: str = "analytic",
    ) -> str:
        """The documented order of the measuring tools (roadmap 7.9.2)."""
        return _tools.prompt_measure_scan(scan_path, spec_path, part_class)

    @server.prompt(
        name="reverse_engineer_scan",
        description="Reverse-engineer a scan into a CAD solid",
    )
    def reverse_engineer_scan_prompt(
        scan_path: str = "scan.stl",
        step_path: str = "part.step",
        part_class: str = "analytic",
    ) -> str:
        """The documented order of the reverse-engineering tools (roadmap 7.9.2)."""
        return _tools.prompt_reverse_engineer_scan(scan_path, step_path, part_class)

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
    http = parser.add_argument_group(
        "HTTP front-end",
        "serve the browser dataset manager's JSON API and MCP over HTTP from "
        "this process instead of stdio (needs the [dashboard] extra; "
        "doc/dashboard.md)",
    )
    http.add_argument("--http", action="store_true", help="serve over HTTP")
    http.add_argument("--host", default=None, help="bind address (default 127.0.0.1)")
    http.add_argument("--port", type=int, default=None, help="port (default 8765)")
    http.add_argument(
        "--token",
        default=None,
        help="the bearer token to require (default: a fresh random one, printed)",
    )
    http.add_argument(
        "--no-token",
        action="store_true",
        help="require no token (only on a machine you alone use)",
    )
    http.add_argument(
        "--allow-origin",
        action="append",
        default=[],
        metavar="ORIGIN",
        help="an extra browser origin to admit (loopback origins and the hosted "
        "docs site are always admitted); repeatable",
    )
    http.add_argument(
        "--runs-dir",
        default=None,
        help="where training runs land (default: <root or cwd>/runs)",
    )
    http.add_argument(
        "--webhook",
        default=None,
        metavar="URL",
        help="POST a JSON payload to this URL when a training job finishes, "
        "fails or is stopped (server-side by design: a client-supplied URL "
        "would be server-side request forgery)",
    )
    args = parser.parse_args(argv)
    if args.runs_dir:
        _tools.set_runs_dir(args.runs_dir)
    if args.webhook:
        _tools.set_webhook(args.webhook)
    if not args.http:
        create_server(root=args.root).run()
        return 0
    from . import _require_http

    try:
        _http = _require_http()
    except ImportError as e:
        import sys

        print(str(e), file=sys.stderr)
        return 1
    host = args.host or _http.DEFAULT_HOST
    port = args.port if args.port is not None else _http.DEFAULT_PORT
    token = None if args.no_token else (args.token or _http.new_token())
    app = _http.build_app(
        create_server(root=args.root, host=host),
        token=token,
        allowed_origins=[*_http.DEFAULT_ALLOWED_ORIGINS, *args.allow_origin],
        root=_tools.get_root(),
        runs_dir=args.runs_dir,
        watch_jobs=bool(args.webhook),
    )
    _http.serve(app, host=host, port=port)
    return 0
