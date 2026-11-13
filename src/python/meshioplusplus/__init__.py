import importlib
import os as _os
import sys as _sys
from typing import TYPE_CHECKING

from .__about__ import __version__
from ._exceptions import ReadError, WriteError
from ._helpers import (
    deregister_format,
    extension_to_filetypes,
    formats,
    read,
    read_metadata,
    register_format,
    write,
    write_points_cells,
)
from ._mesh import CellBlock, Mesh, topological_dimension
from ._regions import Region

if TYPE_CHECKING:
    from . import (
        abaqus,
        abaqus_fil,
        ansys,
        ansys_rst,
        ansysInp,
        avsucd,
        cad_json,
        cae,
        cgns,
        code_aster,
        dex,
        dolfin,
        dxf,
        elmer,
        ensight,
        exodus,
        febio,
        femap,
        flac3d,
        flux,
        frd,
        freefem,
        gid,
        gltf,
        gmsh,
        h5m,
        hmf,
        iges,
        ip,
        libmesh,
        lsdyna,
        lsdyna_binout,
        lsdyna_d3plot,
        marc,
        mdpa,
        med,
        medit,
        medit_sol,
        mfem,
        mff,
        mfm,
        mphbin,
        mphtxt,
        nastran,
        nastran_h5,
        nastran_op2,
        netgen,
        neuroglancer,
        obj,
        off,
        openfoam,
        patran,
        pcd,
        permas,
        ply,
        pmsh,
        pvd,
        pvti,
        pvtp,
        pvtr,
        pvts,
        pvtu,
        radioss,
        radioss_anim,
        radioss_th,
        step,
        stl,
        su2,
        svg,
        szplt,
        tecplot,
        tetgen,
        tikz,
        triangle,
        ugrid,
        unv,
        usd,
        vti,
        vtk,
        vtkhdf,
        vtm,
        vtp,
        vtr,
        vts,
        vtu,
        vtx,
        wkt,
        xdmf,
        xplt,
        xyz,
        z88,
        zarr,
    )
    from ._adapt_volume import adapt_volume
    from ._adjacency import adjacency, write_graph
    from ._agglomerate import agglomerate
    from ._blender import from_blender, has_blender, to_blender
    from ._boundary_layer import boundary_layer
    from ._check_conformity import check_conformity
    from ._check_tetrahedra import check_tetrahedra
    from ._clean import clean
    from ._clip import clip
    from ._colorize import colorize
    from ._conservative_interpolate import conservative_interpolate
    from ._convert_cells import convert_cells
    from ._coordinates import convert_coordinates
    from ._crop import crop
    from ._curvature import compute_curvature
    from ._cut import cut_has_cutcells, cut_level_set, cut_quadrature
    from ._data_average import cell_data_to_point_data, point_data_to_cell_data
    from ._data_calc import data_calc
    from ._data_condition import data_condition
    from ._data_info import data_info
    from ._data_integrate import data_integrate
    from ._data_manage import data_drop, data_keep, data_manage, data_rename
    from ._dataset import DatasetEntry, DatasetManifest
    from ._decimate import decimate
    from ._decimate_volume import decimate_volume
    from ._deviation import deviation, deviation_regions, exaggerate
    from ._diff import diff, meshes_equal
    from ._eigenmodes import laplacian_eigenmodes
    from ._error import estimate_error
    from ._error_size import error_size
    from ._extrude import extrude, revolve
    from ._feature_edges import feature_edges
    from ._find_intersections import find_intersections
    from ._fit_cad import fit_cad
    from ._fit_cad_loop import constrain_primitives, face_support, fit_cad_iterate
    from ._fit_curves import fit_curves
    from ._fit_primitives import fit_primitives
    from ._flip_normals import flip_normals
    from ._geodesic import geodesic_distance
    from ._gpu import (
        from_cupy,
        has_cuda_device,
        has_cupy,
        has_jax,
        has_torch,
        to_cupy,
        to_dlpack,
        to_jax,
        to_torch,
    )
    from ._gradient import gradient
    from ._grid import grid
    from ._grid_transfer import (
        GridArray,
        GridSpec,
        PowerSpectrum,
        expand_grid,
        interpolate_grid,
        power_spectrum,
        resample_grid,
        sample_grid,
        scatter_grid,
        squeeze_grid,
    )
    from ._guard import GeometryGuard, geometry_descriptors
    from ._hausdorff import hausdorff_distance
    from ._hessian import hessian
    from ._hull import convex_hull, delaunay_2d
    from ._iga import iga_evaluate, iga_refine, iga_sample
    from ._interfaces import (
        contact_pairs,
        find_interface,
        region_adjacency,
        split_interface,
    )
    from ._interop import (
        from_arrow,
        from_pyvista,
        from_trimesh,
        has_arrow,
        has_dolfinx,
        has_open3d,
        has_pandas,
        has_polars,
        has_pyvista,
        has_trimesh,
        read_parquet,
        to_arrow,
        to_pandas,
        to_polars,
        to_pyvista,
        to_trimesh,
        write_parquet,
    )
    from ._interpolate import interpolate
    from ._isosurface import isosurface
    from ._measure import align_datums, measure, measure_table
    from ._merge import merge
    from ._mesh_cad import mesh_cad_surface, triangulate_polygon_2d
    from ._metric import prepare_metric
    from ._mirror import mirror
    from ._ml import FeatureMatrix, edge_index, feature_matrix, has_zarr, write_dataset
    from ._mmg import (
        adapt,
        discretize_levelset,
        has_mmg,
        has_mmg_elas,
        mesh_domain,
        metric_from_hessian,
        move,
    )
    from ._morph import morph
    from ._normals import compute_normals
    from ._octree import octree
    from ._octree_mesh import octree_mesh
    from ._offset_surface import hollow, offset_surface, thicken
    from ._optimize_volume import optimize_volume
    from ._parameterize import parameterize
    from ._partition import partition, partition_labels
    from ._periodic import match_periodic_nodes
    from ._pipeline import apply_operations, run_pipeline
    from ._point_budget import PointBudget, select_points, subsample_points
    from ._point_normals import estimate_normals
    from ._primitives import (
        annulus,
        ball,
        box,
        box_points,
        box_surface,
        capsule,
        circle,
        cone,
        cube_sphere,
        cylinder,
        disk,
        fichera,
        geo_sphere,
        helix,
        honeycomb,
        icosa_sphere,
        line,
        lshape,
        moebius,
        ngon,
        octa_sphere,
        pipe,
        plate_with_hole,
        platonic,
        rectangle,
        reference_cells,
        reference_triangle,
        solid_torus,
        sphere_points,
        spherical_shell,
        tetra_sphere,
        torus,
        tube,
        uv_sphere,
    )
    from ._profile_solid import (
        profile_features_from_text,
        profile_solid,
        profile_solid_providers,
    )
    from ._project_to_cad import project_to_cad
    from ._proximity import (
        BistrideHierarchy,
        bistride_hierarchy,
        edge_vectors,
        proximity_graph,
    )
    from ._quality import attach_quality, compute_quality
    from ._quality_gate import check_quality
    from ._reconstruct_surface import reconstruct_surface
    from ._refine import refine
    from ._refine_patches import refine_patches
    from ._region_ops import edit_regions
    from ._registration import (
        auto_align,
        find_symmetry_plane,
        register_landmarks,
        register_rigid,
    )
    from ._remesh import remesh
    from ._remesh_quads import remesh_quads, remesh_quads_backends
    from ._remesh_volume import remesh_volume
    from ._remove_islands import remove_islands
    from ._render import render_image, render_text, snapshot
    from ._reorder import compute_bandwidth, reorder
    from ._repair import repair
    from ._sdf import (
        compute_sdf,
        distance_to_surface,
        sample_distance,
        surface_watertight_check,
    )
    from ._segment import segment
    from ._select import select
    from ._sequence import (
        TimeSeries,
        blend_steps,
        read_sequence,
        resample_sequence,
        run_sequence_pipeline,
        sequence_entries,
        write_sequence,
    )
    from ._shrinkwrap import shrinkwrap
    from ._silhouette import silhouette
    from ._skin import extract_skin
    from ._slice import sections, slice
    from ._smooth import smooth
    from ._sniff import detect, sniff_format
    from ._sobolev_deform import sobolev_deform
    from ._split import split
    from ._stats import compute_stats
    from ._step_elements import write_step_elements
    from ._subdivide import subdivide
    from ._surface import extract_surface
    from ._tensor_invariants import tensor_invariants
    from ._tessellation import Tessellation, tessellate
    from ._tet_poly import poly_to_tet, tet_to_poly
    from ._tetrahedralize import has_netgen, tetrahedralize, tetrahedralize_backends
    from ._thickness import thickness
    from ._transform import transform
    from ._tui import tui
    from ._undo_green import undo_green
    from ._untangle import untangle
    from ._viewer import has_viewer, screenshot, view
    from ._voxelize import voxelize


# Operations, classes and helpers are loaded on first attribute access (PEP 562)
# rather than with the package: importing all of them costs more than the three
# registries and ``read``/``write`` that every user needs. ``name -> module``
# for every public name that is not imported above; format subpackages and
# private modules need no entry, ``__getattr__`` imports a submodule of that name.
# A split wheel ships the shared core in lib/ beside this package; Linux and
# macOS find it through the extension's RPATH, Windows only through a DLL
# directory registered before `_core` is first imported.
_LIB_DIR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "lib")
if _sys.platform == "win32" and _os.path.isdir(_LIB_DIR):  # pragma: no cover
    _os.add_dll_directory(_LIB_DIR)


def get_cmake_dir() -> str:
    """The directory of the installed CMake package (``meshioplusplusConfig.cmake``).

    Only a wheel built with the shared core ships it; an add-on wheel finds the
    core through ``find_package(meshioplusplus CONFIG)`` with
    ``meshioplusplus_DIR`` set to this, the way pybind11's ``get_cmake_dir()``
    works.

    :raises FileNotFoundError: in a build without the shared core.
    """
    path = _os.path.join(_LIB_DIR, "cmake", "meshioplusplus")
    if not _os.path.isdir(path):
        raise FileNotFoundError(
            "meshioplusplus: this build has no shared core (the CMake package is only in the "
            "wheels built with MESHIOPLUSPLUS_PYTHON_SHARED_CORE=ON)"
        )
    return path


_LAZY_ATTRS = {
    "agglomerate": "_agglomerate",
    "from_blender": "_blender",
    "has_blender": "_blender",
    "to_blender": "_blender",
    "adjacency": "_adjacency",
    "geodesic_distance": "_geodesic",
    "parameterize": "_parameterize",
    "laplacian_eigenmodes": "_eigenmodes",
    "convex_hull": "_hull",
    "delaunay_2d": "_hull",
    "estimate_normals": "_point_normals",
    "reconstruct_surface": "_reconstruct_surface",
    "morph": "_morph",
    "prepare_metric": "_metric",
    "untangle": "_untangle",
    "register_rigid": "_registration",
    "register_landmarks": "_registration",
    "find_symmetry_plane": "_registration",
    "auto_align": "_registration",
    "convert_coordinates": "_coordinates",
    "write_graph": "_adjacency",
    "clean": "_clean",
    "conservative_interpolate": "_conservative_interpolate",
    "convert_cells": "_convert_cells",
    "crop": "_crop",
    "compute_curvature": "_curvature",
    "cell_data_to_point_data": "_data_average",
    "point_data_to_cell_data": "_data_average",
    "data_calc": "_data_calc",
    "data_condition": "_data_condition",
    "data_info": "_data_info",
    "data_integrate": "_data_integrate",
    "data_drop": "_data_manage",
    "data_keep": "_data_manage",
    "data_manage": "_data_manage",
    "data_rename": "_data_manage",
    "DatasetEntry": "_dataset",
    "DatasetManifest": "_dataset",
    "check_conformity": "_check_conformity",
    "check_tetrahedra": "_check_tetrahedra",
    "find_intersections": "_find_intersections",
    "thickness": "_thickness",
    "deviation": "_deviation",
    "deviation_regions": "_deviation",
    "exaggerate": "_deviation",
    "colorize": "_colorize",
    "decimate": "_decimate",
    "decimate_volume": "_decimate_volume",
    "diff": "_diff",
    "meshes_equal": "_diff",
    "estimate_error": "_error",
    "error_size": "_error_size",
    "feature_edges": "_feature_edges",
    "fit_cad": "_fit_cad",
    "fit_primitives": "_fit_primitives",
    "iga_evaluate": "_iga",
    "iga_refine": "_iga",
    "iga_sample": "_iga",
    "segment": "_segment",
    "select": "_select",
    "refine_patches": "_refine_patches",
    "from_cupy": "_gpu",
    "has_cuda_device": "_gpu",
    "has_cupy": "_gpu",
    "has_jax": "_gpu",
    "has_torch": "_gpu",
    "to_cupy": "_gpu",
    "to_dlpack": "_gpu",
    "to_jax": "_gpu",
    "to_torch": "_gpu",
    "gradient": "_gradient",
    "grid": "_grid",
    "rectangle": "_primitives",
    "box": "_primitives",
    "disk": "_primitives",
    "uv_sphere": "_primitives",
    "icosa_sphere": "_primitives",
    "cylinder": "_primitives",
    "reference_triangle": "_primitives",
    "ngon": "_primitives",
    "annulus": "_primitives",
    "lshape": "_primitives",
    "plate_with_hole": "_primitives",
    "honeycomb": "_primitives",
    "reference_cells": "_primitives",
    "line": "_primitives",
    "circle": "_primitives",
    "helix": "_primitives",
    "sphere_points": "_primitives",
    "box_points": "_primitives",
    "tetra_sphere": "_primitives",
    "octa_sphere": "_primitives",
    "cube_sphere": "_primitives",
    "geo_sphere": "_primitives",
    "torus": "_primitives",
    "moebius": "_primitives",
    "tube": "_primitives",
    "cone": "_primitives",
    "capsule": "_primitives",
    "box_surface": "_primitives",
    "platonic": "_primitives",
    "ball": "_primitives",
    "spherical_shell": "_primitives",
    "pipe": "_primitives",
    "fichera": "_primitives",
    "solid_torus": "_primitives",
    "GridArray": "_grid_transfer",
    "GridSpec": "_grid_transfer",
    "PowerSpectrum": "_grid_transfer",
    "expand_grid": "_grid_transfer",
    "interpolate_grid": "_grid_transfer",
    "power_spectrum": "_grid_transfer",
    "resample_grid": "_grid_transfer",
    "sample_grid": "_grid_transfer",
    "scatter_grid": "_grid_transfer",
    "squeeze_grid": "_grid_transfer",
    "GeometryGuard": "_guard",
    "geometry_descriptors": "_guard",
    "hausdorff_distance": "_hausdorff",
    "hessian": "_hessian",
    "contact_pairs": "_interfaces",
    "find_interface": "_interfaces",
    "region_adjacency": "_interfaces",
    "split_interface": "_interfaces",
    "from_arrow": "_interop",
    "from_pyvista": "_interop",
    "from_trimesh": "_interop",
    "has_arrow": "_interop",
    "has_dolfinx": "_interop",
    "has_open3d": "_interop",
    "has_pandas": "_interop",
    "has_polars": "_interop",
    "has_pyvista": "_interop",
    "has_trimesh": "_interop",
    "read_parquet": "_interop",
    "to_arrow": "_interop",
    "to_pandas": "_interop",
    "to_polars": "_interop",
    "to_pyvista": "_interop",
    "to_trimesh": "_interop",
    "write_parquet": "_interop",
    "interpolate": "_interpolate",
    "isosurface": "_isosurface",
    "merge": "_merge",
    "FeatureMatrix": "_ml",
    "edge_index": "_ml",
    "feature_matrix": "_ml",
    "has_zarr": "_ml",
    "write_dataset": "_ml",
    "compute_normals": "_normals",
    "optimize_volume": "_optimize_volume",
    "adapt_volume": "_adapt_volume",
    "partition": "_partition",
    "partition_labels": "_partition",
    "match_periodic_nodes": "_periodic",
    "run_pipeline": "_pipeline",
    "apply_operations": "_pipeline",
    "PointBudget": "_point_budget",
    "select_points": "_point_budget",
    "subsample_points": "_point_budget",
    "BistrideHierarchy": "_proximity",
    "bistride_hierarchy": "_proximity",
    "edge_vectors": "_proximity",
    "proximity_graph": "_proximity",
    "attach_quality": "_quality",
    "compute_quality": "_quality",
    "check_quality": "_quality_gate",
    "refine": "_refine",
    "edit_regions": "_region_ops",
    "remesh": "_remesh",
    "adapt": "_mmg",
    "metric_from_hessian": "_mmg",
    "has_mmg": "_mmg",
    "has_mmg_elas": "_mmg",
    "discretize_levelset": "_mmg",
    "mesh_domain": "_mmg",
    "move": "_mmg",
    "remesh_volume": "_remesh_volume",
    "compute_bandwidth": "_reorder",
    "render_image": "_render",
    "render_text": "_render",
    "snapshot": "_render",
    "tui": "_tui",
    "reorder": "_reorder",
    "repair": "_repair",
    "compute_sdf": "_sdf",
    "distance_to_surface": "_sdf",
    "sample_distance": "_sdf",
    "surface_watertight_check": "_sdf",
    "TimeSeries": "_sequence",
    "blend_steps": "_sequence",
    "read_sequence": "_sequence",
    "resample_sequence": "_sequence",
    "run_sequence_pipeline": "_sequence",
    "sequence_entries": "_sequence",
    "write_sequence": "_sequence",
    "shrinkwrap": "_shrinkwrap",
    "extract_skin": "_skin",
    "slice": "_slice",
    "sections": "_slice",
    "fit_curves": "_fit_curves",
    "silhouette": "_silhouette",
    "measure": "_measure",
    "remove_islands": "_remove_islands",
    "offset_surface": "_offset_surface",
    "profile_solid": "_profile_solid",
    "profile_solid_providers": "_profile_solid",
    "profile_features_from_text": "_profile_solid",
    "constrain_primitives": "_fit_cad_loop",
    "face_support": "_fit_cad_loop",
    "fit_cad_iterate": "_fit_cad_loop",
    "thicken": "_offset_surface",
    "hollow": "_offset_surface",
    "flip_normals": "_flip_normals",
    "clip": "_clip",
    "measure_table": "_measure",
    "align_datums": "_measure",
    "write_step_elements": "_step_elements",
    "smooth": "_smooth",
    "sniff_format": "_sniff",
    "detect": "_sniff",
    "sobolev_deform": "_sobolev_deform",
    "split": "_split",
    "compute_stats": "_stats",
    "subdivide": "_subdivide",
    "extract_surface": "_surface",
    "tensor_invariants": "_tensor_invariants",
    "Tessellation": "_tessellation",
    "tessellate": "_tessellation",
    "transform": "_transform",
    "undo_green": "_undo_green",
    "has_viewer": "_viewer",
    "screenshot": "_viewer",
    "view": "_viewer",
    "voxelize": "_voxelize",
    "octree": "_octree",
    "octree_mesh": "_octree_mesh",
    "mesh_cad_surface": "_mesh_cad",
    "project_to_cad": "_project_to_cad",
    "boundary_layer": "_boundary_layer",
    "extrude": "_extrude",
    "revolve": "_extrude",
    "mirror": "_mirror",
    "tet_to_poly": "_tet_poly",
    "poly_to_tet": "_tet_poly",
    "tetrahedralize": "_tetrahedralize",
    "has_netgen": "_tetrahedralize",
    "tetrahedralize_backends": "_tetrahedralize",
    "remesh_quads": "_remesh_quads",
    "cut_has_cutcells": "_cut",
    "cut_level_set": "_cut",
    "cut_quadrature": "_cut",
    "remesh_quads_backends": "_remesh_quads",
    "triangulate_polygon_2d": "_mesh_cad",
}

__all__ = [
    "abaqus",
    "abaqus_fil",
    "ansys",
    "ansys_rst",
    "ansysInp",
    "avsucd",
    "cad_json",
    "cae",
    "cgns",
    "code_aster",
    "dex",
    "dolfin",
    "dxf",
    "elmer",
    "ensight",
    "exodus",
    "febio",
    "femap",
    "flac3d",
    "flux",
    "frd",
    "freefem",
    "gid",
    "gltf",
    "gmsh",
    "h5m",
    "hmf",
    "iges",
    "ip",
    "libmesh",
    "lsdyna",
    "lsdyna_binout",
    "lsdyna_d3plot",
    "marc",
    "mdpa",
    "med",
    "medit",
    "medit_sol",
    "mfem",
    "mff",
    "mfm",
    "mphbin",
    "mphtxt",
    "nastran",
    "nastran_h5",
    "nastran_op2",
    "netgen",
    "neuroglancer",
    "obj",
    "off",
    "openfoam",
    "patran",
    "pcd",
    "permas",
    "ply",
    "pmsh",
    "pvd",
    "pvtp",
    "pvti",
    "pvtr",
    "pvts",
    "pvtu",
    "radioss",
    "radioss_anim",
    "radioss_th",
    "step",
    "stl",
    "su2",
    "svg",
    "szplt",
    "tecplot",
    "tetgen",
    "tikz",
    "triangle",
    "ugrid",
    "unv",
    "usd",
    "vti",
    "vtk",
    "vtkhdf",
    "vtx",
    "vtm",
    "vtp",
    "vtr",
    "vts",
    "vtu",
    "wkt",
    "xdmf",
    "xplt",
    "xyz",
    "z88",
    "zarr",
    "_cli",
    "read",
    "read_metadata",
    "write",
    "register_format",
    "deregister_format",
    "write_points_cells",
    "extension_to_filetypes",
    "formats",
    "extract_skin",
    "extract_surface",
    "compute_quality",
    "compute_curvature",
    "compute_normals",
    "feature_edges",
    "fit_cad",
    "fit_primitives",
    "iga_evaluate",
    "iga_refine",
    "iga_sample",
    "segment",
    "select",
    "refine_patches",
    "hausdorff_distance",
    "match_periodic_nodes",
    "edit_regions",
    "attach_quality",
    "check_quality",
    "check_conformity",
    "check_tetrahedra",
    "find_intersections",
    "thickness",
    "deviation",
    "deviation_regions",
    "exaggerate",
    "colorize",
    "sniff_format",
    "detect",
    "reorder",
    "compute_bandwidth",
    "diff",
    "meshes_equal",
    "merge",
    "interpolate",
    "conservative_interpolate",
    "gradient",
    "hessian",
    "region_adjacency",
    "find_interface",
    "contact_pairs",
    "split_interface",
    "data_integrate",
    "estimate_error",
    "error_size",
    "slice",
    "sections",
    "fit_curves",
    "silhouette",
    "write_step_elements",
    "measure",
    "measure_table",
    "align_datums",
    "remove_islands",
    "offset_surface",
    "profile_solid",
    "profile_solid_providers",
    "profile_features_from_text",
    "constrain_primitives",
    "face_support",
    "fit_cad_iterate",
    "thicken",
    "hollow",
    "flip_normals",
    "clip",
    "isosurface",
    "transform",
    "clean",
    "repair",
    "crop",
    "split",
    "convert_cells",
    "tessellate",
    "Tessellation",
    "subdivide",
    "agglomerate",
    "adjacency",
    "write_graph",
    "geodesic_distance",
    "parameterize",
    "laplacian_eigenmodes",
    "convex_hull",
    "delaunay_2d",
    "estimate_normals",
    "reconstruct_surface",
    "morph",
    "prepare_metric",
    "untangle",
    "register_rigid",
    "register_landmarks",
    "find_symmetry_plane",
    "auto_align",
    "convert_coordinates",
    "refine",
    "undo_green",
    "decimate",
    "decimate_volume",
    "remesh",
    "remesh_volume",
    "adapt_volume",
    "adapt",
    "metric_from_hessian",
    "has_mmg",
    "has_mmg_elas",
    "discretize_levelset",
    "mesh_domain",
    "move",
    "optimize_volume",
    "grid",
    "rectangle",
    "box",
    "disk",
    "uv_sphere",
    "icosa_sphere",
    "cylinder",
    "reference_triangle",
    "ngon",
    "annulus",
    "lshape",
    "plate_with_hole",
    "honeycomb",
    "reference_cells",
    "line",
    "circle",
    "helix",
    "sphere_points",
    "box_points",
    "tetra_sphere",
    "octa_sphere",
    "cube_sphere",
    "geo_sphere",
    "torus",
    "moebius",
    "tube",
    "cone",
    "capsule",
    "box_surface",
    "platonic",
    "ball",
    "spherical_shell",
    "pipe",
    "fichera",
    "solid_torus",
    "voxelize",
    "octree",
    "octree_mesh",
    "mesh_cad_surface",
    "project_to_cad",
    "boundary_layer",
    "extrude",
    "revolve",
    "mirror",
    "tet_to_poly",
    "poly_to_tet",
    "tetrahedralize",
    "has_netgen",
    "tetrahedralize_backends",
    "remesh_quads",
    "remesh_quads_backends",
    "cut_has_cutcells",
    "cut_level_set",
    "cut_quadrature",
    "triangulate_polygon_2d",
    "sample_distance",
    "distance_to_surface",
    "surface_watertight_check",
    "compute_sdf",
    "smooth",
    "shrinkwrap",
    "sobolev_deform",
    "partition",
    "partition_labels",
    "run_pipeline",
    "apply_operations",
    "read_sequence",
    "write_sequence",
    "sequence_entries",
    "run_sequence_pipeline",
    "resample_sequence",
    "blend_steps",
    "TimeSeries",
    "compute_stats",
    "data_manage",
    "data_drop",
    "data_keep",
    "data_rename",
    "point_data_to_cell_data",
    "cell_data_to_point_data",
    "data_calc",
    "data_condition",
    "data_info",
    "tensor_invariants",
    "view",
    "screenshot",
    "has_viewer",
    "render_image",
    "render_text",
    "snapshot",
    "tui",
    "to_pyvista",
    "from_pyvista",
    "to_trimesh",
    "from_trimesh",
    "to_blender",
    "from_blender",
    "has_blender",
    "to_arrow",
    "from_arrow",
    "write_parquet",
    "read_parquet",
    "to_pandas",
    "to_polars",
    "has_pyvista",
    "has_trimesh",
    "has_arrow",
    "has_pandas",
    "has_polars",
    "has_open3d",
    "has_dolfinx",
    "edge_index",
    "feature_matrix",
    "FeatureMatrix",
    "write_dataset",
    "has_zarr",
    "GridSpec",
    "GridArray",
    "PowerSpectrum",
    "sample_grid",
    "scatter_grid",
    "interpolate_grid",
    "resample_grid",
    "squeeze_grid",
    "expand_grid",
    "power_spectrum",
    "PointBudget",
    "select_points",
    "subsample_points",
    "GeometryGuard",
    "geometry_descriptors",
    "BistrideHierarchy",
    "bistride_hierarchy",
    "edge_vectors",
    "proximity_graph",
    "DatasetManifest",
    "DatasetEntry",
    "to_dlpack",
    "to_cupy",
    "from_cupy",
    "has_cupy",
    "has_cuda_device",
    "to_torch",
    "to_jax",
    "has_torch",
    "has_jax",
    "Mesh",
    "CellBlock",
    "Region",
    "ReadError",
    "WriteError",
    "topological_dimension",
    "__version__",
]


def __getattr__(name):
    # ``from meshioplusplus import *`` reaches every name in ``__all__`` here,
    # so it still loads everything; the pinned order of formats that share an
    # extension (``_helpers._EXTENSION_PRIORITY``) does not depend on it.
    module = _LAZY_ATTRS.get(name)
    if module is not None:
        value = getattr(importlib.import_module(f".{module}", __name__), name)
        globals()[name] = value
        return value
    if not name.startswith("__"):
        # A format subpackage (``meshioplusplus.stl``), ``_cli``, or any other
        # submodule that was an attribute while the package imported them all.
        try:
            return importlib.import_module(f".{name}", __name__)
        except ModuleNotFoundError as exc:
            if exc.name != f"{__name__}.{name}":
                raise
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
