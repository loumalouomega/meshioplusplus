"""The declarative operation pipeline: read a mesh, apply a chain of
operations, write the result — described by a ``settings.json`` document (or
the equivalent dict) and exposed as ``meshioplusplus.run_pipeline`` and the
``pipeline`` CLI verb.

This is the **pure-Python twin** of the C++ engine in
``operations/pipeline.cpp`` (reachable as ``_core.run_pipeline_file`` /
``_core.run_pipeline_json`` when the build carries the vendored nlohmann/json).
It deliberately dispatches over the *public Python API* rather than calling the
C++ engine, for the same reason every shim in this repo has a Python path: the
Python ops inherit each operation's C++/numpy parity contract **and** the
per-format Python fallbacks, so a pipeline whose input needs one (a gmsh
``$Periodic`` file, say) still runs — which ``_core.run_pipeline_file``, being
C++ end to end, could never do. It also makes ``run_pipeline`` fully
functional on an sdist install without the submodule.

Vocabulary (PascalCase ops and keys, lowercase enum values, ``Compare`` — not
``Op`` — for refine's comparator) and the report shape
``{"steps": [{"op", ...counters}], "warnings": [...]}`` are the C++ engine's,
transcribed; ``tests/python/test_pipeline.py`` pins the op/key table against
``_core.pipeline_op_table()`` and the two engines' outputs against each other,
so the twins cannot drift silently. See ``doc/pipeline.md`` for the schema.
"""

import json
import os
import pathlib
from contextvars import ContextVar

import numpy as np

from . import _primitives
from ._agglomerate import agglomerate
from ._check_conformity import check_conformity
from ._clean import clean
from ._convert_cells import convert_cells
from ._crop import crop
from ._curvature import compute_curvature
from ._data_average import cell_data_to_point_data, point_data_to_cell_data
from ._data_calc import data_calc
from ._data_condition import data_condition
from ._data_manage import data_drop, data_keep, data_rename
from ._decimate import decimate
from ._decimate_volume import decimate_volume
from ._error import estimate_error
from ._error_size import error_size
from ._feature_edges import feature_edges
from ._fit_cad import DEFAULT_RANSAC_SEED, fit_cad
from ._fit_primitives import TYPES, fit_primitives, types_mask
from ._gradient import gradient
from ._helpers import _filetypes_from_path, _write_format_for_path, read, write
from ._hessian import hessian
from ._interfaces import find_interface, region_adjacency, split_interface
from ._interpolate import interpolate
from ._isosurface import isosurface
from ._merge import merge
from ._mesh import Mesh
from ._mmg import adapt, discretize_levelset, mesh_domain, move
from ._normals import compute_normals
from ._octree import octree
from ._octree_mesh import octree_mesh
from ._optimize_volume import optimize_volume
from ._partition import partition, partition_labels
from ._primitives import (
    annulus,
    box,
    cylinder,
    disk,
    honeycomb,
    icosa_sphere,
    lshape,
    ngon,
    plate_with_hole,
    rectangle,
    reference_cells,
    reference_triangle,
    uv_sphere,
)
from ._provenance import add_operation as _prov_add_operation
from ._provenance import set_source as _prov_set_source
from ._provenance import set_target as _prov_set_target
from ._quality import attach_quality
from ._quality_gate import check_quality, format_quality_gate
from ._refine import refine
from ._refine_patches import STOPS as REFINE_STOPS
from ._refine_patches import refine_patches
from ._region_ops import edit_regions
from ._remesh import remesh
from ._remesh_volume import remesh_volume
from ._render import snapshot
from ._reorder import reorder
from ._repair import repair
from ._sdf import compute_sdf
from ._segment import segment
from ._skin import extract_skin
from ._slice import slice as _slice
from ._smooth import smooth
from ._sobolev_deform import sobolev_deform
from ._split import split
from ._subdivide import subdivide
from ._surface import extract_surface
from ._tensor_invariants import tensor_invariants
from ._tessellation import tessellate
from ._transform import transform
from ._undo_green import undo_green
from ._voxelize import voxelize

__all__ = ["run_pipeline"]

# MCP resolves every expanded output before any writer runs, including existing
# file symlinks. Context-local, so concurrent non-MCP runs inherit no sandbox.
_OUTPUT_PATH_GUARD = ContextVar("pipeline_output_path_guard", default=None)

# The step vocabulary: op name -> allowed parameter keys. A transcription of
# pipe_op_table() in src/cpp/src/operations/pipeline.cpp, pinned against
# _core.pipeline_op_table() by tests/python/test_pipeline.py -- edit both.
_OP_TABLE = {
    "Quality": (),
    "Clean": ("Weld", "Atol", "RemoveOrphans", "DropDegenerate", "DropDuplicateCells"),
    "Smooth": (
        "Method",
        "Iterations",
        "Lambda",
        "Mu",
        "FixBoundary",
        "PreserveFeatures",
        "FeatureAngle",
        "GuardInversion",
    ),
    "Refine": (
        "Levels",
        "Cells",
        "Region",
        "Array",
        "Compare",
        "Value",
        "Closure",
        "RecordLevels",
        "RecordHierarchy",
    ),
    "Decimate": (
        "Ratio",
        "TargetFaces",
        "MaxError",
        "Placement",
        "PreserveBoundary",
        "PreserveFeatures",
        "FeatureAngle",
    ),
    "DecimateVolume": (
        "Ratio",
        "TargetCells",
        "MaxError",
        "Placement",
        "PreserveBoundary",
        "PreserveFeatures",
        "FeatureAngle",
    ),
    "Partition": ("Nparts", "Method", "Imbalance", "Mode", "Seed", "WeightsKey"),
    "Slice": ("Point", "Normal", "RecordParentIds"),
    "Section": ("Point", "Normal", "RecordParentIds"),  # alias of Slice
    "Gradient": ("Array", "Operator", "Method", "Location", "Output", "Component"),
    "Hessian": ("Array", "Method", "Location", "Output"),
    "Curvature": (
        "Mean",
        "Gaussian",
        "DualArea",
        "IncludeBoundary",
        "RecordArea",
        "RecordPrincipal",
        "Region",
    ),
    "Normals": (
        "PointNormals",
        "CellNormals",
        "Weight",
        "SplitAngle",
        "RecordParentIds",
        "Region",
    ),
    "FeatureEdges": (
        "FeatureAngle",
        "Feature",
        "Boundary",
        "NonManifold",
        "Inconsistent",
        "Region",
    ),
    "Segment": (
        "FeatureAngle",
        "NormalTolerance",
        "CurvatureTolerance",
        "MinPatchCells",
        "Region",
        "Prefix",
        "Method",
        "Tolerance",
    ),
    "FitPrimitives": (
        "PatchData",
        "Tolerance",
        "Relative",
        "MinInlierFraction",
        "Types",
        "MaxIterations",
    ),
    "FitCad": (
        "Path",
        "Tolerance",
        "ExactEdges",
        "Types",
        "MeshDeflection",
        "MaxDemotionRounds",
        "Freeform",
        "FreeformDegree",
        "FreeformMaxControl",
        "FreeformSplit",
        "Fillets",
        "Threads",
        "MergeFaces",
        "Backend",
        "RefineSegmentation",
        "MaxRefineIterations",
        "RefineEpsilon",
        "NormalWeight",
        "AllowMerge",
        "AllowSplit",
        "RansacSeed",
        "RefineTolerance",
        "RelativeTolerance",
        "FreeformGeneral",
        "FreeformTolerance",
        "FreeformSmoothing",
        "MaxFreeformPoles",
        "MaxFreeformIterations",
        "MaxFreeformSplits",
        "SewTolerance",
        "ExtendSurfaces",
        "SnapJunctions",
        "FaceDeviation",
        "FreeformMode",
        "FreeformParameterization",
        "ExtraordinaryVertex",
        "LocalRefine",
        "FreeformG1ToleranceDeg",
        "FreeformMaxControlPoints",
        "LayoutTargetPatchSize",
        "LayoutMaxValenceDeviation",
        "ResegmentRounds",
        "LayoutSnapToFeatures",
        "ValidateBySurfaceSampling",
        "ConstrainPrimitives",
    ),
    "RefinePatches": (
        "Tolerance",
        "Relative",
        "NormalWeight",
        "Epsilon",
        "FeatureAngle",
        "MaxIterations",
        "MinPatchCells",
        "RansacSeed",
        "Types",
        "AllowSplit",
        "AllowMerge",
        "Prefix",
    ),
    "EditRegions": ("Edit", "Inputs", "Output", "Kind", "Dim", "Tag", "KeepInputs"),
    "RegionAdjacency": ("Regions",),
    "FindInterface": (
        "RegionA",
        "RegionB",
        "Mode",
        "Master",
        "GapTolerance",
        "AngleTolerance",
        "OverlapTolerance",
    ),
    "SplitInterface": ("Region", "AddCohesive"),
    "QualityGate": ("Require", "MaxInverted", "MaxDegenerate"),
    "CheckConformity": ("Tolerance", "Require", "AttachFlags"),
    "Snapshot": (
        "Path",
        "Width",
        "Height",
        "Supersample",
        "Cols",
        "Rows",
        "Encoding",
        "ColorDepth",
        "PngCompress",
        "View",
        "Azimuth",
        "Elevation",
        "Roll",
        "Perspective",
        "Fov",
        "Zoom",
        "Shading",
        "Edges",
        "FeatureAngle",
        "Background",
        "FillColor",
        "ColorBy",
        "Component",
        "Cmap",
        "VMin",
        "VMax",
        "Colorbar",
        "Axes",
        "ScaleBar",
        "Reduce",
        "Expr",
        "ClipLow",
        "ClipHigh",
        "Symmetric",
        "Scale",
        "ScaleThreshold",
        "Categorical",
        "ColorRegions",
        "CategoryEdges",
        "Isolines",
        "IsoLevels",
        "Vectors",
        "VectorCount",
        "Warp",
        "WarpScale",
        "WarpOutline",
        "Diagnostic",
        "QualityMetric",
        "Cutaway",
        "CutawayTint",
        "Streamlines",
        "StreamSeeds",
        "StreamLength",
        "StreamColor",
        "Theme",
        "Scanlines",
        "Bloom",
        "Fringe",
        "GridFloor",
    ),
    "Repair": (
        "FixOrientation",
        "OrientOutward",
        "FillHoles",
        "SplitNonManifold",
        "MaxHoleEdges",
        "WeldTolerance",
        "RecordProvenance",
    ),
    "SobolevDeform": (
        "Array",
        "LengthScale",
        "FixedPointsArray",
        "FixBoundary",
        "RecordFiltered",
        "MaxIterations",
        "Tolerance",
    ),
    "EstimateError": ("Array", "Method", "Marking", "MarkingValue", "Output", "Marked"),
    "ErrorSize": ("Array", "TargetError", "Order", "MinSize", "MaxSize", "Output"),
    "Adapt": (
        "Hsiz",
        "Size",
        "Metric",
        "Hmin",
        "Hmax",
        "Hausd",
        "Hgrad",
        "Hgradreq",
        "AngleDetection",
        "Optim",
        "NoInsert",
        "NoSwap",
        "NoMove",
        "NoSurf",
        "Nreg",
        "Opnbdy",
        "Frozen",
        "Required",
    ),
    "DiscretizeLevelset": (
        "LevelSet",
        "Isovalue",
        "RemoveSmallComponents",
        "Materials",
        "Hmin",
        "Hmax",
        "Hausd",
        "Hgrad",
        "Hgradreq",
        "AngleDetection",
        "Optim",
        "NoInsert",
        "NoSwap",
        "NoMove",
        "NoSurf",
        "Nreg",
        "Opnbdy",
        "Frozen",
        "Required",
    ),
    "MeshDomain": (
        "Hsiz",
        "Size",
        "Holes",
        "Hmin",
        "Hmax",
        "Hausd",
        "Hgrad",
        "Hgradreq",
        "AngleDetection",
        "Optim",
        "NoInsert",
        "NoSwap",
        "NoMove",
        "NoSurf",
        "Nreg",
        "Opnbdy",
        "Frozen",
        "Required",
    ),
    "Move": (
        "Displacement",
        "Moving",
        "Lagrangian",
        "Hmin",
        "Hmax",
        "Hausd",
        "Hgrad",
        "Hgradreq",
        "AngleDetection",
        "Optim",
        "NoInsert",
        "NoSwap",
        "NoMove",
        "NoSurf",
        "Nreg",
        "Opnbdy",
        "Frozen",
        "Required",
    ),
    "Remesh": (
        "NumClusters",
        "Subdivide",
        "SubsampleRatio",
        "MaxSubdivide",
        "MaxIterations",
        "MaxRepairPasses",
        "Metric",
        "Gradation",
        "PreserveBoundary",
        "MaxAnisotropy",
    ),
    "RemeshVolume": (
        "Resolution",
        "CellSize",
        "Bounds",
        "Padding",
        "PaddingRelative",
        "MaxCells",
        "MaxTets",
        "WarpFraction",
        "Sign",
        "WatertightCheck",
    ),
    "OptimizeVolume": (
        "MaxIterations",
        "Relocate",
        "Flip",
        "PreserveBoundary",
        "MinImprovement",
    ),
    "Isosurface": ("Array", "Isovalue", "Isovalues", "Component", "RecordParentIds"),
    "Voxelize": (
        "Resolution",
        "CellSize",
        "Bounds",
        "Padding",
        "PaddingRelative",
        "Fill",
        "AttachOccupancy",
        "MaxCells",
        "Sign",
        "Plane",
    ),
    "Generate": (
        "Shape",
        "Dims",
        "Size",
        "Origin",
        "Center",
        "Radius",
        "Height",
        "N",
        "Nu",
        "Nv",
        "Nxy",
        "Nz",
        "Subdivisions",
        "Sides",
        "Rings",
        "InnerRadius",
        "OuterRadius",
        "Ntheta",
        "Nr",
        "Nt",
        "Angle",
        "Edge",
        "CellType",
        "Split",
        "Order",
        "MaxCells",
        "Types",
        "Pairs",
        "Complete",
        "Start",
        "End",
        "Pitch",
        "Turns",
        "Seed",
        "MajorRadius",
        "MinorRadius",
        "Width",
        "NumTwists",
        "Ncap",
        "Ratio",
        "Name",
    ),
    "OctreeMesh": (
        "Output",
        "MinSize",
        "MaxSize",
        "MaxLevel",
        "Bounds",
        "Padding",
        "PaddingRelative",
        "CurvatureAngle",
        "FeatureAngle",
        "CornerAngle",
        "FeatureSize",
        "GapCells",
        "GapOutside",
        "SizeArray",
        "SurfaceRegion",
        "GradingLayers",
        "SkinGrading",
        "SurfaceSize",
        "WindingThreshold",
        "WindingBeta",
        "Conformity",
        "AlphaLong",
        "AlphaShort",
        "PreserveFeatures",
        "FitIterations",
        "DihedralFloorGuard",
        "Fill",
        "KeepOutside",
        "AttachDistance",
        "AttachQuality",
        "Pillow",
        "HexGrading",
        "UntangleIterations",
        "MaxCells",
        "MaxTets",
    ),
    "Octree": (
        "RootCells",
        "MinLevel",
        "MaxLevel",
        "Bounds",
        "Padding",
        "PaddingRelative",
        "Dimension",
        "SampleDistance",
        "Balance",
        "Fill",
        "AttachDistance",
        "EvenOdd",
        "MaxCells",
        "Sign",
        "Transition",
        "Snap",
        "SurfaceLevel",
        "CurvatureThreshold",
        "FeatureAngle",
        "SizeArray",
        "Plane",
    ),
    "ComputeSdf": (
        "Structure",
        "Resolution",
        "CellSize",
        "Bounds",
        "Padding",
        "PaddingRelative",
        "RootResolution",
        "MaxDepth",
        "BandCells",
        "RecordLevels",
        "MaxCells",
        "Sign",
        "Location",
        "Band",
    ),
    "Transform": (
        "Translate",
        "Scale",
        "RotateAxis",
        "RotateDegrees",
        "Matrix",
        "ScaleUnits",
        "RotateData",
    ),
    "ConvertCells": ("Mode", "RecordParentIds"),
    "Subdivide": ("RecordParentIds",),
    "Tessellate": ("Levels", "Curved", "Fields", "RecordStencil"),
    "Agglomerate": (
        "TargetGroupSize",
        "MergeCoplanarFaces",
        "CoplanarAngle",
        "MinSphericity",
    ),
    "Crop": (
        "Bbox",
        "Point",
        "Normal",
        "Where",
        "Compare",
        "Value",
        "Mode",
        "RecordIds",
    ),
    "ExtractSurface": ("RecordParentIds",),
    "ExtractSkin": ("Linearize",),
    "Reorder": ("Method",),
    "DataDrop": ("Point", "Cell", "Field", "IgnoreMissing"),
    "DataKeep": ("Point", "Cell", "Field"),
    "DataRename": ("Point", "Cell", "Field"),
    "SetsToData": ("Location", "Name", "Join", "Order"),
    "DataToSets": ("Location", "Key"),
    "DataCalc": ("Expr", "Location", "Overwrite"),
    "DataCondition": (
        "Mode",
        "Location",
        "Names",
        "Scope",
        "Lo",
        "Hi",
        "NanPolicy",
        "NanReplacement",
        "Suffix",
    ),
    "TensorInvariants": (
        "Location",
        "Names",
        "Outputs",
        "Prefix",
        "Suffix",
        "Overwrite",
    ),
    "ToCell": ("Names",),
    "ToPoint": ("Names", "Weight"),
    # The model steps (roadmap 6.1.1): they act on the CAD model the Input
    # read carries beside its mesh (see _apply_model_step).
    "MeshCadSurface": (
        "TargetSize",
        "MaxChordalDeviation",
        "MinPointsPerEdge",
        "Tolerance",
        "RecordFaceIds",
        "Edges",
        "MatchUnreferencedTrims",
        "Strict",
        "MaxPoints",
        "MinAngle",
        "SmoothIterations",
        "SmoothMethod",
        "Order",
        "CurvatureDeviation",
        "CellsAcrossGaps",
        "Grading",
        "PoleRings",
        "FaceSizes",
        "EdgeSizes",
        "Background",
        "BackgroundArray",
        "Metric",
        "VirtualEdges",
        "ReparamConditionNumber",
        "ModelTolerances",
        "VirtualFaces",
        "LopCorners",
        "MetricBackground",
        "MetricArray",
        "MetricError",
        "Untangle",
        "Quads",
        "QuadPairs",
    ),
    "IgaSample": ("Subdivisions", "TrimSamples", "Trim", "Edges"),
    "IgaRefine": (
        "Subdivisions",
        "Elevate",
        "SampleSubdivisions",
        "TrimSamples",
        "Trim",
        "Edges",
    ),
    # Back onto the CAD (roadmap 6.10.1): the B-rep at Model, or the one a
    # CAD Input carries when no step has changed the mesh since.
    "ProjectToCad": ("Model", "FaceArray", "EdgeArray", "MaxDistance", "TrimSamples"),
    # Prismatic layers grown from the surface (roadmap 6.9).
    "BoundaryLayer": (
        "FirstHeight",
        "NumLayers",
        "Growth",
        "Reverse",
        "CollisionFraction",
        "SmoothIterations",
        "MaxShrink",
        "Collapse",
        "NormalSmoothIterations",
        "FrontCheck",
    ),
    # Sweeps (roadmap 3.1.2, 3.1.3) and the reflection (3.2.5).
    "Extrude": (
        "Direction",
        "NumLayers",
        "Ratio",
        "Offsets",
        "RecordLayers",
        "MaxCells",
    ),
    "Revolve": (
        "Axis",
        "AxisPoint",
        "Angle",
        "NumLayers",
        "Tolerance",
        "RecordLayers",
        "MaxCells",
    ),
    "Mirror": ("Point", "Normal", "KeepOriginal", "Weld", "Tolerance", "ReflectData"),
    # Tetrahedra filling the closed triangle surface (roadmap 6.7).
    "Tetrahedralize": (
        "Backend",
        "MaxSize",
        "MaxVolume",
        "Grading",
        "OptimizeSteps",
        "MaxCells",
        "FaceArray",
    ),
}

_V2_OP_TABLE = {
    **_OP_TABLE,
    "Merge": (
        "Inputs",
        "Weld",
        "Atol",
        "SourceTag",
        "DataPolicy",
        "DropDuplicateCells",
    ),
    "Interpolate": (
        "Inputs",
        "Method",
        "Arrays",
        "Extrapolate",
        "DefaultValue",
        "OnConflict",
    ),
    "UndoGreen": ("Inputs",),
    "Split": ("By", "Tag"),
    "Partition": _OP_TABLE["Partition"] + ("RecordIds", "GhostLayers"),
}

# Ops that exist in meshio++ but cannot be a single-mesh pipeline step; the
# error names where to go instead (the C++ engine's exact rule).
_EXCLUDED_OPS = {
    "Merge": "'Merge' needs several input meshes and is not a pipeline step; "
    "use the `merge` CLI verb",
    "Interpolate": "'Interpolate' needs a second (source) mesh and is not a "
    "pipeline step; use the `interpolate` CLI verb",
    "Split": "'Split' produces several output meshes and is not a pipeline "
    "step; use the `split` CLI verb",
    "Shrinkwrap": "'Shrinkwrap' needs a second (target) mesh and is not a "
    "pipeline step; use the `shrinkwrap` CLI verb",
    "Diff": "'Diff' compares two meshes and is not a pipeline step; use the "
    "`diff` CLI verb",
    "UndoGreen": "'UndoGreen' needs a second (coarse) mesh and is not a "
    "pipeline step; use the `undo-green` CLI verb",
    "ConservativeInterpolate": "'ConservativeInterpolate' needs a second "
    "(source) mesh and is not a pipeline step; use the "
    "`conservative-interpolate` CLI verb",
}

# The model steps, and the formats whose read carries a model (the C++
# engine's pipe_is_model_step and pipe_read_with_model).
_MODEL_OPS = ("MeshCadSurface", "IgaSample", "IgaRefine")
_MODEL_FORMATS = ("cad_json", "step", "iges", "mfem")


def _err(op, what):
    return ValueError(f"meshio++: pipeline: step '{op}': {what}")


def _mmg_kwargs(step):
    """The shared MMG parameters of a step (roadmap 3.2); ``Frozen`` names a point region."""

    def size(key):
        value = _number(step, key, -1.0)
        return None if value < 0 else value

    return dict(
        hmin=size("Hmin"),
        hmax=size("Hmax"),
        hausd=size("Hausd"),
        hgrad=size("Hgrad"),
        hgradreq=size("Hgradreq"),
        angle_detection=_number(step, "AngleDetection", 45.0),
        optim=_flag(step, "Optim", False),
        noinsert=_flag(step, "NoInsert", False),
        noswap=_flag(step, "NoSwap", False),
        nomove=_flag(step, "NoMove", False),
        nosurf=_flag(step, "NoSurf", False),
        nreg=_flag(step, "Nreg", False),
        opnbdy=_flag(step, "Opnbdy", False),
        frozen=_text(step, "Frozen", "") or None,
        required=_svec(step, "Required") or None,
    )


def _mmg_entry(entry, warnings, report):
    """An MMG step's counters, as the core's pipe_mmg_report writes them."""
    entry["NumPointsIn"] = float(report["num_points_in"])
    entry["NumPointsOut"] = float(report["num_points_out"])
    entry["NumElementsIn"] = float(report["num_elements_in"])
    entry["NumElementsOut"] = float(report["num_elements_out"])
    entry["LowFailure"] = 1.0 if report["status"] == "low_failure" else 0.0
    entry["QualityMin"] = float(report["quality"]["min"])
    entry["QualityMean"] = float(report["quality"]["mean"])
    warnings.extend(report["warnings"])


def _number(step, key, fallback):
    v = step.get(key)
    if v is None:
        return fallback
    # bool is an int subclass; a boolean here is a type error, not a number.
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise _err(step["Op"], f"parameter '{key}' must be a number")
    return float(v)


def _flag(step, key, fallback):
    v = step.get(key)
    if v is None:
        return fallback
    if not isinstance(v, bool):
        raise _err(step["Op"], f"parameter '{key}' must be a boolean")
    return v


def _segment_method(step):
    method = _text(step, "Method", "angle")
    if method not in ("angle", "primitive"):
        raise _err(step["Op"], "'Method' must be 'angle' or 'primitive'")
    return method


def _fit_cad_backend(step):
    backend = _text(step, "Backend", "auto")
    if backend not in ("auto", "occt", "native"):
        raise _err(step["Op"], "'Backend' must be 'auto', 'occt' or 'native'")
    return backend


def _seed(step, key, fallback):
    """A ``RansacSeed``: an integer in ``[0, 2**53]``, the largest a JSON
    number (a double, in the C++ engine) carries exactly."""
    v = step.get(key)
    if v is None:
        return fallback
    what = f"parameter '{key}' must be an integer in [0, 2^53] (9007199254740992)"
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise _err(step["Op"], what)
    if isinstance(v, float) and not v.is_integer():
        raise _err(step["Op"], what)
    if not 0 <= v <= 1 << 53:
        raise _err(step["Op"], what)
    return int(v)


def _text(step, key, fallback):
    v = step.get(key)
    if v is None:
        return fallback
    if not isinstance(v, str):
        raise _err(step["Op"], f"parameter '{key}' must be a string")
    return v


def _snap_color(step, key):
    if key not in step:
        return None
    from ._render import parse_color

    try:
        return parse_color(_text(step, key, ""))
    except ValueError:
        raise _err(
            step["Op"], f"parameter '{key}' must be #rrggbb, #rrggbbaa or none"
        ) from None


def _snap_choice(step, key, names):
    if key not in step:
        return None
    name = _text(step, key, "")
    if name not in names:
        raise _err(
            step["Op"],
            f"parameter '{key}' must be one of {', '.join(names)}, not '{name}'",
        )
    return name


def _snapshot_options(step):
    """The keyword arguments of `snapshot` from a Snapshot step's parameters,
    with the defaults of the native engine's `pipe_snap_render_options`."""
    iso_levels = _dvec(step, "IsoLevels")
    flat_cut = _dvec(step, "Cutaway") or []
    if len(flat_cut) % 6 != 0 or len(flat_cut) > 12:
        raise _err(
            step["Op"],
            "parameter 'Cutaway' must be one or two planes of six numbers each "
            "(a point, then the normal of the side kept)",
        )
    cutaway = [tuple(flat_cut[i : i + 6]) for i in range(0, len(flat_cut), 6)] or None
    options = {
        "width": int(_number(step, "Width", 320)),
        "height": int(_number(step, "Height", 240)),
        "supersample": int(_number(step, "Supersample", 1)),
        "cols": int(_number(step, "Cols", 100)),
        "rows": int(_number(step, "Rows", 40)),
        "encoding": _snap_choice(
            step, "Encoding", ("halfblock", "quadrant", "sextant", "braille", "ascii")
        )
        or "halfblock",
        "color_depth": _snap_choice(
            step, "ColorDepth", ("truecolor", "256", "16", "mono")
        )
        or "truecolor",
        "png_compress": int(_number(step, "PngCompress", 0)),
        "view": _text(step, "View", "") or None,
        "azimuth": _number(step, "Azimuth", 45.0),
        "elevation": _number(step, "Elevation", 35.264389682754654),
        "roll": _number(step, "Roll", 0.0),
        "projection": "perspective" if _flag(step, "Perspective", False) else None,
        "fov": _number(step, "Fov", 30.0),
        "zoom": _number(step, "Zoom", 1.0),
        "shading": _snap_choice(step, "Shading", ("none", "flat", "smooth")),
        "edges": _snap_choice(step, "Edges", ("none", "all", "feature")),
        "feature_angle": _number(step, "FeatureAngle", 30.0),
        "background": _snap_color(step, "Background"),
        "fill_color": _snap_color(step, "FillColor"),
        "color_by": _text(step, "ColorBy", "") or None,
        "component": (
            int(_number(step, "Component", 0)) if "Component" in step else None
        ),
        "theme": _snap_choice(step, "Theme", ("none", "synthwave")),
        "cmap": _text(
            step,
            "Cmap",
            "synthwave" if _text(step, "Theme", "none") == "synthwave" else "viridis",
        ),
        "scanlines": _flag(step, "Scanlines", False),
        "bloom": _flag(step, "Bloom", False),
        "fringe": _flag(step, "Fringe", False),
        "grid_floor": _flag(step, "GridFloor", False),
        "vmin": _number(step, "VMin", 0.0) if "VMin" in step else None,
        "vmax": _number(step, "VMax", 0.0) if "VMax" in step else None,
        "colorbar": _flag(step, "Colorbar", False),
        "axes": _flag(step, "Axes", False),
        "scale_bar": _flag(step, "ScaleBar", False),
        "reduce": _text(step, "Reduce", "") or None,
        "expr": _text(step, "Expr", "") or None,
        "clip": (
            (
                _number(step, "ClipLow", 0.0) if "ClipLow" in step else None,
                _number(step, "ClipHigh", 100.0) if "ClipHigh" in step else None,
            )
            if "ClipLow" in step or "ClipHigh" in step
            else None
        ),
        "symmetric": _flag(step, "Symmetric", False),
        "scale": _snap_choice(step, "Scale", ("linear", "log", "symlog")),
        "scale_threshold": _number(step, "ScaleThreshold", 1.0),
        "categorical": _flag(step, "Categorical", False),
        "color_regions": _flag(step, "ColorRegions", False),
        "category_edges": _flag(step, "CategoryEdges", False),
        "isolines": int(_number(step, "Isolines", 0)),
        "iso_levels": iso_levels if iso_levels else None,
        "vectors": _text(step, "Vectors", "") or None,
        "vector_count": int(_number(step, "VectorCount", 200)),
        "warp": _text(step, "Warp", "") or None,
        "warp_scale": _number(step, "WarpScale", 1.0),
        "warp_outline": _flag(step, "WarpOutline", False),
        "diagnostic": _snap_choice(
            step,
            "Diagnostic",
            (
                "none",
                "quality",
                "inverted",
                "degenerate",
                "orientation",
                "free_edges",
                "edge_length",
            ),
        ),
        "quality_metric": _text(step, "QualityMetric", "") or None,
        "cutaway": cutaway,
        "cutaway_tint": _snap_color(step, "CutawayTint"),
        "streamlines": _text(step, "Streamlines", "") or None,
        "stream_seeds": int(_number(step, "StreamSeeds", 40)),
        "stream_length": _number(step, "StreamLength", 0.5),
        "stream_color": _snap_color(step, "StreamColor"),
    }
    return {k: v for k, v in options.items() if v is not None}


def _dvec(step, key):
    v = step.get(key)
    if v is None:
        return None
    if not isinstance(v, (list, tuple)) or any(
        isinstance(e, bool) or not isinstance(e, (int, float)) for e in v
    ):
        raise _err(step["Op"], f"parameter '{key}' must be an array of numbers")
    return [float(e) for e in v]


def _svec(step, key):
    v = step.get(key)
    if v is None:
        return None
    if not isinstance(v, (list, tuple)) or any(not isinstance(e, str) for e in v):
        raise _err(step["Op"], f"parameter '{key}' must be an array of strings")
    return list(v)


def _vec3(step, key):
    v = _dvec(step, key)
    if v is None or len(v) != 3:
        raise _err(step["Op"], f"parameter '{key}' must be an array of 3 numbers")
    return v


def _apply_sets_data(mesh, step):
    """Python reference for the region-backed sets/data pipeline steps."""
    location = _text(step, "Location", "cell")
    if location not in ("point", "cell"):
        raise _err(step["Op"], "Location must be 'point' or 'cell'")
    mesh = mesh.copy()
    if step["Op"] == "DataToSets":
        if "Key" not in step:
            raise _err(step["Op"], "'DataToSets' requires 'Key'")
        getattr(mesh, f"{location}_data_to_sets")(_text(step, "Key", ""))
        return mesh
    sets = getattr(mesh, f"{location}_sets")
    names = _svec(step, "Order") or list(sets)
    if (
        len(names) != len(sets)
        or len(set(names)) != len(names)
        or set(names) != set(sets)
    ):
        raise _err(step["Op"], "Order must name every set exactly once")
    if not names:
        return mesh
    name = _text(step, "Name", _text(step, "Join", "-").join(names))
    if mesh._sets_data_core(
        "sets_to_data", location, name=name, join=_text(step, "Join", "-"), order=names
    ):
        return mesh
    if location == "point":
        labels = np.full(len(mesh.points), -1, dtype=np.int64)
        for label, set_name in enumerate(names):
            labels[sets[set_name]] = label
        mesh.point_data[name] = labels
        mesh.point_sets = {}
    else:
        blocks = [np.full(len(cb), -1, dtype=np.int64) for cb in mesh.cells]
        for label, set_name in enumerate(names):
            for b, entries in enumerate(sets[set_name]):
                if entries is not None:
                    blocks[b][entries] = label
        mesh.cell_data[name] = blocks
        mesh.cell_sets = {}
    return mesh


def _total_cells(mesh):
    return int(sum(len(block.data) for block in mesh.cells))


def _validate_step(step, version=1):
    if not isinstance(step, dict):
        raise ValueError("meshio++: pipeline: every operation must be an object")
    op = step.get("Op")
    if not isinstance(op, str) or not op:
        raise ValueError("meshio++: pipeline: every operation needs a string 'Op'")
    table = _V2_OP_TABLE if version == 2 else _OP_TABLE
    if op not in table:
        if op in _EXCLUDED_OPS:
            raise ValueError(f"meshio++: pipeline: {_EXCLUDED_OPS[op]}")
        known = ", ".join(table)
        raise ValueError(
            f"meshio++: pipeline: unknown operation '{op}' (known: {known})"
        )
    allowed = table[op]
    for key in step:
        if key != "Op" and key not in allowed:
            keys = ", ".join(allowed)
            hint = f" (known: {keys})" if keys else " (the op takes no parameters)"
            raise _err(op, f"unknown parameter '{key}'{hint}")
    if version == 2 and op in ("Merge", "Interpolate", "UndoGreen"):
        paths = _svec(step, "Inputs") or []
        if not paths or (op != "Merge" and len(paths) != 1):
            raise _err(
                op,
                "Inputs requires "
                + ("at least one path" if op == "Merge" else "exactly one path"),
            )
        if not all(paths):
            raise _err(op, "Inputs paths must not be empty")


def _apply_data_manage(mesh, op, step):
    locations = (("Point", "point"), ("Cell", "cell"), ("Field", "field"))
    if op == "DataDrop":
        ignore_missing = _flag(step, "IgnoreMissing", False)
        for key, location in locations:
            if key in step:
                mesh = data_drop(mesh, location, _svec(step, key) or [], ignore_missing)
    elif op == "DataKeep":
        for key, location in locations:
            if key in step:
                mesh = data_keep(mesh, location, _svec(step, key) or [])
    else:  # DataRename
        for key, location in locations:
            for entry in _svec(step, key) or []:
                # "OLD:NEW" splits on the LAST colon (names carry colons:
                # `gmsh:physical`) -- the CLI's documented rule.
                old, sep, new = entry.rpartition(":")
                if not sep or not old or not new:
                    raise _err(op, f"rename entry '{entry}' must be 'OLD:NEW'")
                mesh = data_rename(mesh, location, old, new)
    return mesh


def _apply_transform(mesh, step):
    kwargs = {}
    if "Translate" in step:
        kwargs["translate"] = _vec3(step, "Translate")
    if "Scale" in step:
        v = step["Scale"]
        kwargs["scale"] = (
            _dvec(step, "Scale")
            if isinstance(v, (list, tuple))
            else _number(step, "Scale", 1.0)
        )
    if "RotateAxis" in step or "RotateDegrees" in step:
        axis = _vec3(step, "RotateAxis")
        if "RotateDegrees" not in step:
            raise _err("Transform", "'RotateAxis' requires 'RotateDegrees'")
        kwargs["rotate"] = (axis, _number(step, "RotateDegrees", 0.0))
    if "Matrix" in step:
        matrix = _dvec(step, "Matrix")
        if matrix is None or len(matrix) != 16:
            raise _err("Transform", "'Matrix' expects 16 values (row-major 4x4)")
        kwargs["matrix"] = matrix
    if "ScaleUnits" in step:
        kwargs["scale_units"] = _number(step, "ScaleUnits", 1.0)
    if len(kwargs) != 1:
        raise _err(
            "Transform",
            "give exactly one of 'Translate'/'Scale'/'RotateAxis'+'RotateDegrees'/"
            "'Matrix'/'ScaleUnits'",
        )
    return transform(
        mesh, rotate_vector_data=_flag(step, "RotateData", False), **kwargs
    )


_GENERATE_SHAPES = (
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
)

_ORIGIN = ("origin", "triple", [0.0, 0.0, 0.0])
_CENTER = ("center", "triple", [0.0, 0.0, 0.0])
_QUAD = ("cell_type", "text", "quad")
_HEX = ("cell_type", "text", "hexahedron")
_SURFACE = {
    "Radius": ("radius", "number", 1.0),
    "Height": ("height", "number", 1.0),
    "Ntheta": ("ntheta", "count", 16),
    "Nz": ("nz", "count", 4),
    "Center": _CENTER,
    "CellType": _QUAD,
}
#: The second tier of shapes: pipeline key -> (keyword, kind, default). The
#: ``order``, ``complete`` and ``max_cells`` keywords are added by
#: :func:`_generate` unless the shape is in ``_GENERATE_LINEAR``.
_GENERATE_SPECS = {
    "line": {
        "N": ("n", "count", 8),
        "Start": ("start", "triple", [0.0, 0.0, 0.0]),
        "End": ("end", "triple", [1.0, 0.0, 0.0]),
    },
    "circle": {
        "Radius": ("radius", "number", 1.0),
        "N": ("n", "count", 16),
        "Center": _CENTER,
    },
    "helix": {
        "Radius": ("radius", "number", 1.0),
        "Pitch": ("pitch", "number", 1.0),
        "Turns": ("turns", "number", 2.0),
        "N": ("n", "count", 32),
        "Center": _CENTER,
    },
    "sphere_points": {
        "N": ("n", "count", 100),
        "Radius": ("radius", "number", 1.0),
        "Center": _CENTER,
    },
    "box_points": {
        "N": ("n", "count", 100),
        "Size": ("size", "triple", [1.0, 1.0, 1.0]),
        "Origin": _ORIGIN,
        "Seed": ("seed", "count", 0),
    },
    "tetra_sphere": {
        "Radius": ("radius", "number", 1.0),
        "Subdivisions": ("subdivisions", "count", 3),
        "Center": _CENTER,
    },
    "octa_sphere": {
        "Radius": ("radius", "number", 1.0),
        "Subdivisions": ("subdivisions", "count", 2),
        "Center": _CENTER,
    },
    "cube_sphere": {
        "Radius": ("radius", "number", 1.0),
        "N": ("n", "count", 4),
        "Center": _CENTER,
        "CellType": _QUAD,
    },
    "geo_sphere": {
        "Radius": ("radius", "number", 1.0),
        "N": ("n", "count", 4),
        "Center": _CENTER,
    },
    "torus": {
        "MajorRadius": ("major_radius", "number", 1.0),
        "MinorRadius": ("minor_radius", "number", 0.25),
        "Nu": ("nu", "count", 32),
        "Nv": ("nv", "count", 16),
        "Center": _CENTER,
        "CellType": _QUAD,
    },
    "moebius": {
        "Radius": ("radius", "number", 1.0),
        "Width": ("width", "number", 0.5),
        "Nu": ("nu", "count", 32),
        "Nv": ("nv", "count", 4),
        "NumTwists": ("num_twists", "count", 1),
        "Center": _CENTER,
        "CellType": _QUAD,
    },
    "tube": _SURFACE,
    "cone": _SURFACE,
    "capsule": dict(
        _SURFACE, Radius=("radius", "number", 0.5), Ncap=("ncap", "count", 4)
    ),
    "box_surface": {
        "Dims": ("dims", "dims3", 1),
        "Size": ("size", "triple", [1.0, 1.0, 1.0]),
        "Origin": _ORIGIN,
        "CellType": _QUAD,
    },
    "platonic": {
        "Name": ("name", "text", "icosahedron"),
        "Radius": ("radius", "number", 1.0),
        "Center": _CENTER,
    },
    "ball": {
        "Radius": ("radius", "number", 1.0),
        "N": ("n", "count", 4),
        "Center": _CENTER,
        "CellType": _HEX,
    },
    "spherical_shell": {
        "InnerRadius": ("inner_radius", "number", 0.5),
        "OuterRadius": ("outer_radius", "number", 1.0),
        "N": ("n", "count", 4),
        "Nr": ("nr", "count", 2),
        "Ratio": ("ratio", "number", 1.0),
        "Center": _CENTER,
    },
    "pipe": {
        "InnerRadius": ("inner_radius", "number", 0.5),
        "OuterRadius": ("outer_radius", "number", 1.0),
        "Height": ("height", "number", 1.0),
        "Ntheta": ("ntheta", "count", 16),
        "Nr": ("nr", "count", 2),
        "Nz": ("nz", "count", 4),
        "Ratio": ("ratio", "number", 1.0),
        "Center": _CENTER,
    },
    "fichera": {
        "N": ("n", "count", 2),
        "Size": ("size", "number", 2.0),
        "Center": _CENTER,
        "CellType": _HEX,
    },
    "solid_torus": {
        "MajorRadius": ("major_radius", "number", 1.0),
        "MinorRadius": ("minor_radius", "number", 0.25),
        "N": ("n", "count", 4),
        "Nu": ("nu", "count", 32),
        "Center": _CENTER,
    },
}
#: Shapes of the second tier without ``order`` (points and polygon faces).
_GENERATE_LINEAR = {"sphere_points", "box_points", "platonic"}


def _generate(step):
    """The mesh a ``Generate`` step builds (the engine twin of the C++ branch)."""
    shape = _text(step, "Shape", "")
    max_cells = int(_number(step, "MaxCells", 20000000.0))
    order = int(_number(step, "Order", 1.0))
    complete = _flag(step, "Complete", False)

    def triple(key, default, pad_zero=False):
        v = _dvec(step, key)
        if not v:
            return default
        v = list(v)
        if len(v) == 2 and pad_zero:
            v.append(0.0)
        if len(v) != 3:
            raise _err("Generate", f"'{key}' expects three values")
        return v

    def count(key, default):
        return int(_number(step, key, float(default)))

    def dims_of(n, default):
        v = _dvec(step, "Dims")
        if not v:
            return [default] * n
        if len(v) != n:
            raise _err("Generate", f"'Dims' expects {n} counts")
        return [int(x) for x in v]

    zero = [0.0, 0.0, 0.0]
    if shape == "rectangle":
        size = _dvec(step, "Size")
        if size and len(size) != 2:
            raise _err("Generate", "'Size' expects two values")
        return rectangle(
            dims=dims_of(2, 1),
            size=list(size) if size else [1.0, 1.0],
            origin=triple("Origin", zero, True),
            cell_type=_text(step, "CellType", "quad"),
            split=_text(step, "Split", "zigzag"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "box":
        return box(
            dims=dims_of(3, 1),
            size=triple("Size", [1.0, 1.0, 1.0]),
            origin=triple("Origin", zero),
            cell_type=_text(step, "CellType", "hexahedron"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "disk":
        return disk(
            radius=_number(step, "Radius", 1.0),
            n=count("N", 8),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "quad"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "uv_sphere":
        return uv_sphere(
            radius=_number(step, "Radius", 1.0),
            nu=count("Nu", 16),
            nv=count("Nv", 8),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "triangle"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "icosa_sphere":
        return icosa_sphere(
            radius=_number(step, "Radius", 1.0),
            subdivisions=count("Subdivisions", 2),
            center=triple("Center", zero),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "cylinder":
        return cylinder(
            radius=_number(step, "Radius", 1.0),
            height=_number(step, "Height", 1.0),
            nxy=count("Nxy", 8),
            nz=count("Nz", 8),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "hexahedron"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "reference_triangle":
        size = _dvec(step, "Size")
        if size and len(size) != 2:
            raise _err("Generate", "'Size' expects two values")
        return reference_triangle(
            n=count("N", 4),
            size=list(size) if size else [1.0, 1.0],
            origin=triple("Origin", zero, True),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "ngon":
        return ngon(
            sides=count("Sides", 6),
            radius=_number(step, "Radius", 1.0),
            rings=count("Rings", 1),
            center=triple("Center", zero),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "annulus":
        return annulus(
            inner_radius=_number(step, "InnerRadius", 0.5),
            outer_radius=_number(step, "OuterRadius", 1.0),
            ntheta=count("Ntheta", 16),
            nr=count("Nr", 4),
            angle=_number(step, "Angle", 6.283185307179586),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "quad"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "lshape":
        return lshape(
            n=count("N", 4),
            size=_number(step, "Size", 2.0),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "quad"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "plate_with_hole":
        return plate_with_hole(
            radius=_number(step, "Radius", 0.25),
            size=_number(step, "Size", 1.0),
            nt=count("Nt", 8),
            nr=count("Nr", 4),
            center=triple("Center", zero),
            cell_type=_text(step, "CellType", "quad"),
            order=order,
            max_cells=max_cells,
            complete=complete,
        )
    if shape == "honeycomb":
        return honeycomb(
            dims=dims_of(2, 4),
            edge=_number(step, "Edge", 1.0),
            origin=triple("Origin", zero, True),
            max_cells=max_cells,
        )
    if shape == "reference_cells":
        return reference_cells(
            types=_svec(step, "Types") or None, pairs=_flag(step, "Pairs", False)
        )
    if shape in _GENERATE_SPECS:
        kwargs = {}
        for key, (name, kind, default) in _GENERATE_SPECS[shape].items():
            if kind == "number":
                kwargs[name] = _number(step, key, default)
            elif kind == "count":
                kwargs[name] = count(key, default)
            elif kind == "text":
                kwargs[name] = _text(step, key, default)
            elif kind == "dims3":
                kwargs[name] = dims_of(3, default)
            else:
                kwargs[name] = triple(key, default)
        if shape not in _GENERATE_LINEAR:
            kwargs["order"] = order
            kwargs["complete"] = complete
        kwargs["max_cells"] = max_cells
        return getattr(_primitives, shape)(**kwargs)
    raise _err(
        "Generate",
        f"'Shape' must be one of {', '.join(_GENERATE_SHAPES)}, not '{shape}'",
    )


def _render_op(step):
    """Renders one step as ``Op(Key=value, Key=value, ...)`` for the
    provenance operation chain -- the Python twin of `pipe_render_op` in
    pipeline.cpp. Keys are sorted (the C++ side's std::map iteration order
    is already sorted, so this is what keeps the two textually identical)."""
    parts = []
    for key in sorted(k for k in step if k != "Op"):
        value = step[key]
        if isinstance(value, bool):
            rendered = "true" if value else "false"
        elif isinstance(value, list):
            rendered = "[" + ", ".join(str(v) for v in value) + "]"
        else:
            rendered = str(value)
        parts.append(f"{key}={rendered}")
    return f"{step['Op']}(" + ", ".join(parts) + ")"


def _bl_collapse(step):
    collapse = _text(step, "Collapse", "points")
    if collapse not in ("points", "layers"):
        raise _err(
            "BoundaryLayer",
            f"'Collapse' must be 'points' or 'layers', not '{collapse}'",
        )
    return collapse


def _tetrahedralize_step(mesh, step, entry, warnings):
    """``Tetrahedralize``: :func:`meshioplusplus.tetrahedralize`, its counters
    and warnings the C++ engine's."""
    import warnings as _warnings

    from ._tetrahedralize import tetrahedralize

    backend = _text(step, "Backend", "auto")
    if backend not in ("auto", "netgen", "octree"):
        raise _err(
            "Tetrahedralize",
            f"'Backend' must be 'auto', 'netgen' or 'octree', not '{backend}'",
        )
    with _warnings.catch_warnings():
        # The report carries them, as the C++ engine's does.
        _warnings.simplefilter("ignore", RuntimeWarning)
        mesh, report = tetrahedralize(
            mesh,
            backend=backend,
            max_size=_number(step, "MaxSize", 0.0),
            max_volume=_number(step, "MaxVolume", 0.0),
            grading=_number(step, "Grading", 0.3),
            optimize_steps=int(_number(step, "OptimizeSteps", 3)),
            max_cells=int(_number(step, "MaxCells", 50_000_000)),
            face_array=_text(step, "FaceArray", "cad:face"),
            return_report=True,
        )
    warnings.extend("Tetrahedralize: " + w for w in report["warnings"])
    entry.update(
        NumTetrahedra=float(report["num_tetrahedra"]),
        NumPoints=float(report["num_points"]),
        NumNewPoints=float(report["num_new_points"]),
        NumRegions=float(report["num_regions"]),
        NumBoundaryTriangles=float(report["num_boundary_triangles"]),
        BoundaryPreserved=1.0 if report["boundary_preserved"] else 0.0,
        HausdorffToSurface=float(report["hausdorff_to_surface"]),
        MinDihedralDeg=float(report["min_dihedral_deg"]),
        MaxDihedralDeg=float(report["max_dihedral_deg"]),
        Volume=float(report["volume"]),
    )
    return mesh


def _apply_step(mesh, step, steps, warnings):
    """One validated step. Counters and warning strings are the C++ engine's,
    transcribed -- keep the two in lockstep."""
    op = step["Op"]
    entry = {"op": op}

    if op in _MODEL_OPS:
        raise _err(
            op,
            "needs the CAD model a pipeline's Input reads, which only a whole run "
            "carries (run_pipeline), not a step applied on its own",
        )
    if op == "BoundaryLayer":
        from ._boundary_layer import boundary_layer

        if "FirstHeight" not in step:
            raise _err(op, "'FirstHeight' is required")
        mesh, report = boundary_layer(
            mesh,
            _number(step, "FirstHeight", 0.0),
            num_layers=int(_number(step, "NumLayers", 3)),
            growth=_number(step, "Growth", 1.2),
            reverse=_flag(step, "Reverse", False),
            collision_fraction=_number(step, "CollisionFraction", 0.45),
            smooth_iterations=int(_number(step, "SmoothIterations", 5)),
            max_shrink=int(_number(step, "MaxShrink", 6)),
            collapse=_bl_collapse(step),
            normal_smooth_iterations=int(_number(step, "NormalSmoothIterations", 0)),
            front_check=_flag(step, "FrontCheck", True),
            return_report=True,
        )
        entry.update(
            NumCells=float(report["num_cells"]),
            NumShrunk=float(report["num_shrunk"]),
            NumCollapsed=float(report["num_collapsed"]),
            NumFewerLayers=float(report["num_fewer_layers"]),
            NumFrontCrossings=float(report["num_front_crossings"]),
            MinScaledJacobian=float(report["min_scaled_jacobian"]),
        )
    elif op == "Extrude":
        from ._extrude import extrude

        mesh, report = extrude(
            mesh,
            _vec3(step, "Direction") if "Direction" in step else (0.0, 0.0, 1.0),
            num_layers=int(_number(step, "NumLayers", 1)),
            ratio=_number(step, "Ratio", 1.0),
            offsets=_dvec(step, "Offsets"),
            record_layers=_flag(step, "RecordLayers", False),
            max_cells=int(_number(step, "MaxCells", 20_000_000)),
            return_report=True,
        )
        entry.update(
            NumCells=float(report["num_cells"]),
            NumFlipped=float(report["num_flipped"]),
            MinScaledJacobian=float(report["min_scaled_jacobian"]),
        )
    elif op == "Revolve":
        from ._extrude import revolve

        mesh, report = revolve(
            mesh,
            _vec3(step, "Axis") if "Axis" in step else (0.0, 1.0, 0.0),
            axis_point=(
                _vec3(step, "AxisPoint") if "AxisPoint" in step else (0.0, 0.0, 0.0)
            ),
            angle=_number(step, "Angle", 360.0),
            num_layers=int(_number(step, "NumLayers", 16)),
            tolerance=_number(step, "Tolerance", 1e-9),
            record_layers=_flag(step, "RecordLayers", False),
            max_cells=int(_number(step, "MaxCells", 20_000_000)),
            return_report=True,
        )
        entry.update(
            NumCells=float(report["num_cells"]),
            NumOnAxis=float(report["num_on_axis"]),
            MinScaledJacobian=float(report["min_scaled_jacobian"]),
        )
    elif op == "Mirror":
        from ._mirror import mirror

        mesh, report = mirror(
            mesh,
            _vec3(step, "Point") if "Point" in step else (0.0, 0.0, 0.0),
            _vec3(step, "Normal") if "Normal" in step else (1.0, 0.0, 0.0),
            keep_original=_flag(step, "KeepOriginal", False),
            weld=_flag(step, "Weld", True),
            tolerance=_number(step, "Tolerance", 1e-9),
            reflect_vector_data=_flag(step, "ReflectData", False),
            return_report=True,
        )
        entry.update(
            NumFlippedCells=float(report["num_flipped_cells"]),
            NumPointsOnPlane=float(report["num_points_on_plane"]),
            NumDroppedCells=float(report["num_dropped_cells"]),
        )
    elif op == "Tetrahedralize":
        mesh = _tetrahedralize_step(mesh, step, entry, warnings)
    elif op == "ProjectToCad":
        mesh = _project_to_cad_step(mesh, step, entry, None)
    elif op == "Quality":
        mesh = attach_quality(mesh)
    elif op == "Clean":
        mesh, report = clean(
            mesh,
            weld=_flag(step, "Weld", False),
            atol=_number(step, "Atol", 1e-8),
            remove_orphans=_flag(step, "RemoveOrphans", True),
            drop_degenerate=_flag(step, "DropDegenerate", True),
            drop_duplicate_cells=_flag(step, "DropDuplicateCells", True),
            return_report=True,
        )
        entry["PointsWelded"] = report["points_welded"]
        entry["PointsRemovedOrphan"] = report["points_removed_orphan"]
        entry["CellsDroppedDegenerate"] = report["cells_dropped_degenerate"]
        entry["CellsDroppedDuplicate"] = report["cells_dropped_duplicate"]
    elif op == "Smooth":
        mesh, report = smooth(
            mesh,
            method=_text(step, "Method", "taubin"),
            iterations=int(_number(step, "Iterations", 10)),
            lambda_=_number(step, "Lambda", -1.0),
            mu=_number(step, "Mu", -0.34),
            fix_boundary=_flag(step, "FixBoundary", True),
            preserve_features=_flag(step, "PreserveFeatures", True),
            feature_angle=_number(step, "FeatureAngle", 30.0),
            guard_inversion=_flag(step, "GuardInversion", True),
            return_report=True,
        )
        entry["NumNodesMoved"] = report["num_nodes_moved"]
        entry["MaxDisplacement"] = report["max_displacement"]
        entry["NumSkippedInversion"] = report["num_skipped_inversion"]
    elif op == "Refine":
        where = None
        array = _text(step, "Array", "")
        if array:
            # `Compare`, not `Op`: `Op` is the step's own discriminant.
            compare = _text(step, "Compare", "<")
            where = f"{array} {compare} {_number(step, 'Value', 0.0)}"
        cells = _dvec(step, "Cells")
        mesh = refine(
            mesh,
            levels=int(_number(step, "Levels", 1)),
            cells=[int(c) for c in cells] if cells is not None else None,
            region=_text(step, "Region", "") or None,
            where=where,
            closure=_text(step, "Closure", "redgreen"),
            record_levels=_flag(step, "RecordLevels", False),
            record_hierarchy=_flag(step, "RecordHierarchy", False),
        )
    elif op == "Decimate":
        ratio = _number(step, "Ratio", -1.0)
        target_faces = _number(step, "TargetFaces", -1.0)
        max_error = _number(step, "MaxError", -1.0)
        if ratio < 0.0 and target_faces < 0.0 and max_error < 0.0:
            ratio = 0.5  # a usable pipeline default (the C++ engine's rule)
        mesh, report = decimate(
            mesh,
            ratio=ratio if ratio >= 0.0 else None,
            target_faces=int(target_faces) if target_faces >= 0.0 else None,
            max_error=max_error if max_error >= 0.0 else None,
            placement=_text(step, "Placement", "optimal"),
            preserve_boundary=_flag(step, "PreserveBoundary", True),
            preserve_features=_flag(step, "PreserveFeatures", True),
            feature_angle=_number(step, "FeatureAngle", 30.0),
            return_report=True,
        )
        entry["FacesRemoved"] = report["faces_removed"]
        entry["CollapsesRejected"] = report["collapses_rejected"]
    elif op == "DecimateVolume":
        ratio = _number(step, "Ratio", -1.0)
        target_cells = _number(step, "TargetCells", -1.0)
        max_error = _number(step, "MaxError", -1.0)
        if ratio < 0.0 and target_cells < 0.0 and max_error < 0.0:
            ratio = 0.5  # a usable pipeline default (the C++ engine's rule)
        mesh, report = decimate_volume(
            mesh,
            ratio=ratio if ratio >= 0.0 else None,
            target_cells=int(target_cells) if target_cells >= 0.0 else None,
            max_error=max_error if max_error >= 0.0 else None,
            placement=_text(step, "Placement", "optimal"),
            preserve_boundary=_flag(step, "PreserveBoundary", False),
            preserve_features=_flag(step, "PreserveFeatures", True),
            feature_angle=_number(step, "FeatureAngle", 30.0),
            return_report=True,
        )
        entry["TetsRemoved"] = report["tets_removed"]
        entry["CollapsesRejected"] = report["collapses_rejected"]
    elif op == "Partition":
        nparts = int(_number(step, "Nparts", 2))
        labels = partition_labels(
            mesh,
            nparts,
            method=_text(step, "Method", "auto"),
            imbalance=_number(step, "Imbalance", 0.03),
            mode=_text(step, "Mode", "eco"),
            seed=int(_number(step, "Seed", 0)),
            weights=_text(step, "WeightsKey", "") or None,
        )
        mesh.cell_data["partition:part"] = labels
        entry["Nparts"] = float(nparts)
    elif op in ("Slice", "Section"):
        mesh = _slice(
            mesh,
            origin=_vec3(step, "Point"),
            normal=_vec3(step, "Normal"),
            record_parent_ids=_flag(step, "RecordParentIds", False),
        )
        entry["SectionFaces"] = float(_total_cells(mesh))
        if _total_cells(mesh) == 0:
            warnings.append(
                "the section is empty; the plane misses the mesh -- move or flip it"
            )
    elif op == "Gradient":
        component = int(_number(step, "Component", -1.0))
        mesh, report = gradient(
            mesh,
            _text(step, "Array", ""),
            operator=_text(step, "Operator", "gradient"),
            method=_text(step, "Method", "green-gauss"),
            location=_text(step, "Location", "cell"),
            output=_text(step, "Output", "") or None,
            component=component if component >= 0 else None,
            overwrite=True,
            return_report=True,
        )
        entry["NumSkipped"] = report["num_skipped"]
        entry["NumFallback"] = report["num_fallback"]
        if report["num_skipped"] > 0:
            warnings.append(
                f"gradient: {report['num_skipped']} cell(s) could not be "
                "differentiated and are NaN"
            )
    elif op == "Hessian":
        mesh, report = hessian(
            mesh,
            _text(step, "Array", ""),
            method=_text(step, "Method", "green-gauss"),
            location=_text(step, "Location", "cell"),
            output=_text(step, "Output", "") or None,
            overwrite=True,
            return_report=True,
        )
        entry["NumSkipped"] = report["num_skipped"]
        entry["NumFallback"] = report["num_fallback"]
        if report["num_skipped"] > 0:
            warnings.append(
                f"hessian: {report['num_skipped']} cell(s) could not be "
                "evaluated and are NaN"
            )
    elif op == "Curvature":
        mesh, report = compute_curvature(
            mesh,
            mean=_flag(step, "Mean", True),
            gaussian=_flag(step, "Gaussian", True),
            dual_area=_text(step, "DualArea", "mixed-voronoi"),
            include_boundary=_flag(step, "IncludeBoundary", False),
            record_area=_flag(step, "RecordArea", False),
            record_principal=_flag(step, "RecordPrincipal", False),
            region=_text(step, "Region", ""),
            return_report=True,
        )
        entry["NumBoundary"] = report["num_boundary"]
        entry["NumIsolated"] = report["num_isolated"]
        entry["NumDegenerate"] = report["num_degenerate"]
        entry["TotalAngleDefect"] = report["total_angle_defect"]
        # H's sign comes from the surface's own winding, so a mesh whose facets
        # disagree about which side is out yields sign-flipped patches with no
        # error raised. Say so rather than letting it pass silently.
        if report["quality"]["inconsistent_pairs"] > 0:
            warnings.append(
                f"curvature: {report['quality']['inconsistent_pairs']} edge "
                "pair(s) wind the same way, so the sign of 'curvature:mean' is "
                "not trustworthy"
            )
    elif op == "FeatureEdges":
        mesh, report = feature_edges(
            mesh,
            feature_angle=_number(step, "FeatureAngle", 30.0),
            feature=_flag(step, "Feature", True),
            boundary=_flag(step, "Boundary", True),
            non_manifold=_flag(step, "NonManifold", True),
            inconsistent=_flag(step, "Inconsistent", True),
            region=_text(step, "Region", ""),
            return_report=True,
        )
        entry["NumFeature"] = report["num_feature"]
        entry["NumBoundary"] = report["num_boundary"]
        entry["NumNonManifold"] = report["num_non_manifold"]
        entry["NumInconsistent"] = report["num_inconsistent"]
    elif op == "Segment":
        mesh, report = segment(
            mesh,
            feature_angle=_number(step, "FeatureAngle", 30.0),
            normal_tolerance=_number(step, "NormalTolerance", 0.0),
            curvature_tolerance=_number(step, "CurvatureTolerance", 0.0),
            min_patch_cells=int(_number(step, "MinPatchCells", 1.0)),
            region=_text(step, "Region", ""),
            prefix=_text(step, "Prefix", "patch"),
            return_report=True,
            method=_segment_method(step),
            tolerance=_number(step, "Tolerance", 1e-3),
        )
        entry["NumPatches"] = report["num_patches"]
        entry["NumMerged"] = report["num_merged"]
        entry["NumFeatureEdges"] = report["num_feature_edges"]
        entry["NumBoundaryEdges"] = report["num_boundary_edges"]
    elif op == "FitPrimitives":
        types = _svec(step, "Types")
        try:
            mask_types = tuple(types) if types else TYPES[1:]
            if types:
                types_mask(mask_types)
        except ValueError:
            bad = next(t for t in types if t not in TYPES[1:])
            raise _err(op, f"unknown primitive type '{bad}' in 'Types'") from None
        mesh, report = fit_primitives(
            mesh,
            patch_data=_text(step, "PatchData", "segment:patch"),
            tolerance=_number(step, "Tolerance", 1e-3),
            relative=_flag(step, "Relative", True),
            min_inlier_fraction=_number(step, "MinInlierFraction", 0.9),
            types=mask_types,
            max_iterations=int(_number(step, "MaxIterations", 50.0)),
            return_report=True,
        )
        entry["NumPatches"] = report["num_patches"]
        for t in TYPES[1:]:
            entry["Num" + t.capitalize()] = report["num_" + t]
        entry["NumUnfitted"] = report["num_unfitted"]
    elif op == "FitCad":
        # The fitted skin becomes a CAD solid, written to `Path` when given; the
        # mesh continues as the solid re-meshed.
        types = _svec(step, "Types") if "Types" in step else None
        bad = [t for t in (types or ()) if t not in TYPES[1:]]
        if bad:
            raise _err(op, f"unknown primitive type '{bad[0]}' in 'Types'")
        tolerance = _number(step, "Tolerance", 0.0)
        deflection = _number(step, "MeshDeflection", 0.0)
        refine_seg = _flag(step, "RefineSegmentation", False)
        mesh, report = fit_cad(
            mesh,
            _text(step, "Path", "") or None,
            tolerance=tolerance if tolerance > 0.0 else None,
            exact_edges=_flag(step, "ExactEdges", True),
            types=types,
            mesh_deflection=deflection if deflection > 0.0 else None,
            max_demotion_rounds=int(_number(step, "MaxDemotionRounds", 8.0)),
            freeform=_flag(step, "Freeform", True),
            freeform_degree=int(_number(step, "FreeformDegree", 3.0)),
            freeform_max_control=int(_number(step, "FreeformMaxControl", 64.0)),
            freeform_split=_flag(step, "FreeformSplit", True),
            fillets=_flag(step, "Fillets", True),
            threads=_flag(step, "Threads", True),
            return_report=True,
            merge_faces=_flag(step, "MergeFaces", True),
            backend=_fit_cad_backend(step),
            refine_segmentation=refine_seg,
            max_refine_iterations=int(_number(step, "MaxRefineIterations", 10.0)),
            refine_epsilon=_number(step, "RefineEpsilon", 1e-3),
            normal_weight=_number(step, "NormalWeight", 1.0),
            allow_merge=_flag(step, "AllowMerge", True),
            allow_split=_flag(step, "AllowSplit", True),
            ransac_seed=_seed(step, "RansacSeed", DEFAULT_RANSAC_SEED),
            refine_tolerance=_number(step, "RefineTolerance", 1e-3),
            relative_tolerance=_flag(step, "RelativeTolerance", True),
            freeform_general=_flag(step, "FreeformGeneral", False),
            freeform_tolerance=_number(step, "FreeformTolerance", 0.0),
            freeform_smoothing=_number(step, "FreeformSmoothing", 1e-12),
            max_freeform_poles=int(_number(step, "MaxFreeformPoles", 40.0)),
            max_freeform_iterations=int(_number(step, "MaxFreeformIterations", 30.0)),
            max_freeform_splits=int(_number(step, "MaxFreeformSplits", 2.0)),
            sew_tolerance=_number(step, "SewTolerance", 0.0),
            extend_surfaces=_flag(step, "ExtendSurfaces", False),
            snap_junctions=_flag(step, "SnapJunctions", False),
            face_deviation=_flag(step, "FaceDeviation", False),
            freeform_mode=_text(step, "FreeformMode", "legacy"),
            freeform_parameterization=_text(step, "FreeformParameterization", "auto"),
            extraordinary_vertex=_text(step, "ExtraordinaryVertex", "g1x"),
            local_refine=_text(step, "LocalRefine", "knots"),
            freeform_g1_tolerance_deg=_number(step, "FreeformG1ToleranceDeg", 1.0),
            freeform_max_control_points=int(
                _number(step, "FreeformMaxControlPoints", 400.0)
            ),
            layout_target_patch_size=_number(step, "LayoutTargetPatchSize", 0.1),
            layout_max_valence_deviation=int(
                _number(step, "LayoutMaxValenceDeviation", 2.0)
            ),
            resegment_rounds=int(_number(step, "ResegmentRounds", 2.0)),
            layout_snap_to_features=_flag(step, "LayoutSnapToFeatures", True),
            validate_by_surface_sampling=_flag(step, "ValidateBySurfaceSampling", True),
            constrain_primitives=_flag(step, "ConstrainPrimitives", False),
        )
        entry["Valid"] = 1.0 if report["valid"] else 0.0
        entry["Native"] = 1.0 if report["backend"] == "native" else 0.0
        entry["VolumeRatio"] = report["volume_ratio"]
        entry["MaxDeviation"] = report["max_deviation"]
        entry["NumFaces"] = report["num_faces"]
        entry["NumFacetedFaces"] = report["num_faceted_faces"]
        entry["NumBspline"] = report["num_bspline"]
        entry["NumSplitPatches"] = report["num_split_patches"]
        entry["NumFillet"] = report["num_fillet"]
        entry["NumHelical"] = report["num_helical"]
        entry["NumVariableFillet"] = report["num_variable_fillet"]
        entry["NumBlend"] = report["num_blend"]
        entry["NumRoundedThread"] = report["num_rounded_thread"]
        entry["NumVertexBlend"] = report["num_vertex_blend"]
        entry["NumChamfer"] = report["num_chamfer"]
        entry["NumConicBlend"] = report["num_conic_blend"]
        entry["NumTaperedThread"] = report["num_tapered_thread"]
        entry["NumVariablePitchThread"] = report["num_variable_pitch_thread"]
        entry["NumDemoted"] = report["num_demoted"]
        entry["Closed"] = 1.0 if report["closed"] else 0.0
        # Only a run that refined reports the refinement, so a report without
        # `RefineSegmentation` stays what it was before the option existed.
        if report.get("refined"):
            entry["Refined"] = 1.0
            entry["NumRefineSplits"] = report["num_refine_splits"]
            entry["NumRefineMerges"] = report["num_refine_merges"]
            entry["ExplainedAreaFraction"] = report["explained_area_fraction"]
            entry["RefinementKept"] = 1.0 if report["refinement_kept"] else 0.0
        # Likewise each Phase 3 option reports only when it is set.
        if "num_freeform_general" in report:
            entry["NumFreeformGeneral"] = report["num_freeform_general"]
            entry["NumFreeformRefinements"] = report["num_freeform_refinements"]
            entry["AvgFreeformPoles"] = report["avg_freeform_poles"]
        if "num_extended_edges" in report:
            entry["NumExtendedEdges"] = report["num_extended_edges"]
        if "num_snapped_junctions" in report:
            entry["NumSnappedJunctions"] = report["num_snapped_junctions"]
        if "hausdorff_max" in report:
            entry["HausdorffMax"] = report["hausdorff_max"]
            entry["HausdorffRms"] = report["hausdorff_rms"]
        if "num_constrained_primitives" in report:
            entry["NumConstrainedPrimitives"] = report["num_constrained_primitives"]
            entry["NumConstraintGroups"] = report["num_constraint_groups"]
        # The freeform modes (v17.12.0) report only outside the legacy mode.
        if "freeform_mode" in report:
            entry["BrepChecked"] = 1.0 if report["brep_checked"] else 0.0
            entry["FacetAreaFraction"] = report["facet_area_fraction"]
            entry["NumFreeformAdaptive"] = report["num_freeform_adaptive"]
            entry["NumNetworkPatches"] = report["num_network_patches"]
            entry["NumNetworkEdges"] = report["num_network_edges"]
            entry["NumExtraordinaryVertices"] = report["num_extraordinary_vertices"]
            entry["NumResegmented"] = report["num_resegmented"]
            entry["NumEscalations"] = report["num_escalations"]
            entry["MaxG1AngleDeg"] = report["max_g1_angle_deg"]
    elif op == "RefinePatches":
        # A data step: the `segment:patch` labels are refined, the cells never move.
        types = _svec(step, "Types") if "Types" in step else None
        bad = [t for t in (types or ()) if t not in TYPES[1:]]
        if bad:
            raise _err(op, f"unknown primitive type '{bad[0]}' in 'Types'")
        mesh, report = refine_patches(
            mesh,
            tolerance=_number(step, "Tolerance", 1e-3),
            relative=_flag(step, "Relative", True),
            normal_weight=_number(step, "NormalWeight", 1.0),
            epsilon=_number(step, "Epsilon", 1e-3),
            feature_angle=_number(step, "FeatureAngle", 30.0),
            max_iterations=int(_number(step, "MaxIterations", 10.0)),
            min_patch_cells=int(_number(step, "MinPatchCells", 3.0)),
            ransac_seed=_seed(step, "RansacSeed", DEFAULT_RANSAC_SEED),
            types=types,
            allow_split=_flag(step, "AllowSplit", True),
            allow_merge=_flag(step, "AllowMerge", True),
            prefix=_text(step, "Prefix", "patch"),
            return_report=True,
        )
        entry["NumPatches"] = report["num_patches"]
        entry["NumIterations"] = report["num_iterations"]
        entry["NumSplits"] = report["num_splits"]
        entry["NumMerges"] = report["num_merges"]
        entry["NumPrimitive"] = report["num_primitive"]
        entry["NumFreeform"] = report["num_freeform"]
        entry["ExplainedAreaFraction"] = report["explained_area_fraction"]
        entry["HausdorffMax"] = report["hausdorff_max"]
        # Counters are numbers: the `RefineStop` code, 0 unchanged, 1 small_gain,
        # 2 iteration_cap, 3 energy_rose.
        entry["Stop"] = float(REFINE_STOPS.index(report["stop"]))
    elif op == "Snapshot":
        # A side output, not a transform: the mesh passes through untouched and
        # the frame is written to `Path` (.png, .txt, .ansi, .html or .cast).
        path = _text(step, "Path", "")
        if not path:
            raise _err(step["Op"], "'Path' is required")
        options = _snapshot_options(step)
        snapshot(mesh, path, **options)
        entry["Width"] = options["width"]
        entry["Height"] = options["height"]
    elif op == "QualityGate":
        # A gate, not a transform: the mesh passes through untouched, and a
        # failed check stops the pipeline with the summary as the error.
        result = check_quality(
            mesh,
            _svec(step, "Require"),
            max_inverted=int(_number(step, "MaxInverted", 0)),
            max_degenerate=int(_number(step, "MaxDegenerate", 0)),
        )
        if not result["passed"]:
            raise RuntimeError(
                "meshio++: pipeline: QualityGate failed\n" + format_quality_gate(result)
            )
        entry["NumChecks"] = len(result["checks"])
        entry["NumFailed"] = sum(not c["passed"] for c in result["checks"])
    elif op == "CheckConformity":
        # A gate by default: a mesh that is not conforming stops the pipeline;
        # with Require false it only reports. AttachFlags writes
        # conformity:hanging onto the mesh that flows on.
        attach = _flag(step, "AttachFlags", False)
        flagged, rep = check_conformity(
            mesh,
            tolerance=_number(step, "Tolerance", 1e-9),
            attach_flags=attach,
            return_report=True,
        )
        if _flag(step, "Require", True) and not rep["conforming"]:
            raise RuntimeError(
                "meshio++: pipeline: CheckConformity failed: "
                f"{rep['num_hanging_nodes']} hanging node(s), "
                f"{rep['num_t_junctions']} T-junction(s), "
                f"{rep['num_partial_overlaps']} partial overlap(s), "
                f"{rep['num_unwelded_facets']} unwelded facet pair(s)"
            )
        entry["NumHangingNodes"] = float(rep["num_hanging_nodes"])
        entry["NumTJunctions"] = float(rep["num_t_junctions"])
        entry["NumPartialOverlaps"] = float(rep["num_partial_overlaps"])
        entry["NumUnweldedFacets"] = float(rep["num_unwelded_facets"])
        if attach:
            mesh = flagged
    elif op == "EditRegions":
        kind = _text(step, "Kind", "") or None
        edit = {
            "op": _text(step, "Edit", ""),
            "inputs": [{"name": n, "kind": kind} for n in _svec(step, "Inputs")],
            "output": _text(step, "Output", ""),
            "keep_inputs": _flag(step, "KeepInputs", True),
        }
        if "Dim" in step:
            edit["dim"] = int(_number(step, "Dim", -1))
        if "Tag" in step:
            edit["tag"] = int(_number(step, "Tag", -1))
        mesh = edit_regions(mesh, [edit])
        entry["NumRegions"] = len(mesh.regions)
    elif op == "RegionAdjacency":
        mesh = region_adjacency(mesh, _svec(step, "Regions") or None)
        entry["NumFacets"] = _total_cells(mesh)
        entry["Measure"] = sum(
            float(values.sum()) for values in mesh.cell_data["interface:measure"]
        )
    elif op == "FindInterface":
        mesh, report = find_interface(
            mesh,
            _text(step, "RegionA", ""),
            _text(step, "RegionB", ""),
            mode=_text(step, "Mode", "conforming"),
            master=_text(step, "Master", "a"),
            gap_tolerance=_number(step, "GapTolerance", 0.0),
            angle_tolerance=_number(step, "AngleTolerance", 30.0),
            overlap_tolerance=_number(step, "OverlapTolerance", 0.0),
            return_report=True,
        )
        for key, target in (
            ("num_pairs", "NumPairs"),
            ("area", "Area"),
            ("max_gap", "MaxGap"),
            ("unmatched_a", "UnmatchedA"),
            ("unmatched_b", "UnmatchedB"),
        ):
            entry[target] = report[key]
    elif op == "SplitInterface":
        mesh, report = split_interface(
            mesh,
            _text(step, "Region", ""),
            add_cohesive=_flag(step, "AddCohesive", False),
            return_report=True,
        )
        entry["NumDuplicatedPoints"] = report["num_duplicated_points"]
        entry["NumCohesiveCells"] = report["num_cohesive_cells"]
    elif op == "Normals":
        # An absent SplitAngle means one smooth normal per point; a number is
        # the crease angle in degrees, and the split appends points.
        split = _number(step, "SplitAngle", 30.0) if "SplitAngle" in step else None
        mesh, report = compute_normals(
            mesh,
            point_normals=_flag(step, "PointNormals", True),
            cell_normals=_flag(step, "CellNormals", False),
            weight=_text(step, "Weight", "angle"),
            split_angle=split,
            region=_text(step, "Region", ""),
            record_parent_ids=_flag(step, "RecordParentIds", False),
            return_report=True,
        )
        entry["NumIsolated"] = report["num_isolated"]
        entry["NumUndefined"] = report["num_undefined"]
        entry["NumDegenerate"] = report["num_degenerate"]
        entry["NumSplitPoints"] = report["num_split_points"]
        entry["NumAddedPoints"] = report["num_added_points"]
        if report["quality"]["inconsistent_pairs"] > 0 and split is None:
            warnings.append(
                f"normals: {report['quality']['inconsistent_pairs']} edge pair(s) "
                "wind the same way, so the normals there average faces that "
                "disagree about which side is out"
            )
    elif op == "Repair":
        mesh, report = repair(
            mesh,
            fix_orientation=_flag(step, "FixOrientation", True),
            orient_outward=_flag(step, "OrientOutward", True),
            fill_holes=_flag(step, "FillHoles", True),
            split_non_manifold=_flag(step, "SplitNonManifold", True),
            max_hole_edges=int(_number(step, "MaxHoleEdges", 10)),
            weld_tolerance=_number(step, "WeldTolerance", 0.0),
            record_provenance=_flag(step, "RecordProvenance", False),
            return_report=True,
        )
        entry["NumFlipped"] = report["num_flipped"]
        entry["NumComponents"] = report["num_components"]
        entry["NumVerticesSplit"] = report["num_vertices_split"]
        entry["NumHolesFilled"] = report["num_holes_filled"]
        entry["NumHolesSkipped"] = report["num_holes_skipped"]
        entry["NumFacesAdded"] = report["num_faces_added"]
        entry["NumPointsAdded"] = report["num_points_added"]
        entry["PointsWelded"] = report["points_welded"]
        # What repair could NOT fix is worth saying.
        if report["num_unorientable"] > 0:
            warnings.append(
                f"repair: {report['num_unorientable']} component(s) are not "
                "orientable; their winding is a best effort"
            )
        if report["quality_after"]["non_manifold_edges"] > 0:
            warnings.append(
                f"repair: {report['quality_after']['non_manifold_edges']} "
                "non-manifold edge(s) remain; repair counts them, it does not "
                "split them"
            )
    elif op == "SobolevDeform":
        mesh, report = sobolev_deform(
            mesh,
            _text(step, "Array", ""),
            _number(step, "LengthScale", 0.0),
            fixed_points=_text(step, "FixedPointsArray", "") or None,
            fix_boundary=_flag(step, "FixBoundary", False),
            record_filtered=_flag(step, "RecordFiltered", False),
            max_iterations=int(_number(step, "MaxIterations", 128)),
            tolerance=_number(step, "Tolerance", 1e-10),
            return_report=True,
        )
        entry["NumIterations"] = report["num_iterations"]
        entry["Residual"] = report["residual"]
        entry["Converged"] = 1.0 if report["converged"] else 0.0
        entry["NumFixed"] = report["num_fixed"]
        entry["NumIsolated"] = report["num_isolated"]
        entry["MaxDisplacement"] = report["max_displacement"]
        if not report["converged"]:
            warnings.append(
                f"sobolev_deform: conjugate gradients did not converge in "
                f"{report['num_iterations']} iteration(s); the last iterate is "
                "returned"
            )
    elif op == "EstimateError":
        mesh, report = estimate_error(
            mesh,
            _text(step, "Array", ""),
            method=_text(step, "Method", "zz"),
            marking=_text(step, "Marking", "none"),
            marking_value=_number(step, "MarkingValue", 0.0),
            output=_text(step, "Output", "") or None,
            marked_name=_text(step, "Marked", "") or None,
            overwrite=True,
            return_report=True,
        )
        entry["GlobalError"] = report["global_error"]
        entry["NumSkipped"] = report["num_skipped"]
        entry["NumMarked"] = report["num_marked"]
        if report["num_skipped"] > 0:
            warnings.append(
                f"estimate_error: {report['num_skipped']} cell(s) could not be "
                "evaluated and are NaN"
            )
    elif op == "ErrorSize":
        # A pure data step (roadmap 6.10.2): EstimateError's per-cell
        # indicator becomes the per-point size a later mesh_cad_surface's
        # background reads. Geometry is untouched.
        mesh, report = error_size(
            mesh,
            _text(step, "Array", "error:zz"),
            target_error=_number(step, "TargetError", 0.0),
            order=_number(step, "Order", 1.0),
            min_size=_number(step, "MinSize", 0.0),
            max_size=_number(step, "MaxSize", 0.0),
            output=_text(step, "Output", "size"),
            overwrite=True,
            return_report=True,
        )
        entry["NumCells"] = report["num_cells"]
        entry["NumSkipped"] = report["num_skipped"]
        entry["NumClamped"] = report["num_clamped"]
        entry["MinSize"] = report["min_size"]
        entry["MaxSize"] = report["max_size"]
        if report["num_skipped"] > 0:
            warnings.append(
                f"error_size: {report['num_skipped']} cell(s) have no finite size "
                "and size no point"
            )
    elif op == "Adapt":
        # MMG (roadmap 3.2.1): new points and connectivity, regions by name.
        hsiz = _number(step, "Hsiz", -1.0)
        mesh, report = adapt(
            mesh,
            hsiz=hsiz if hsiz > 0 else None,
            size=_text(step, "Size", "") or None,
            metric=_text(step, "Metric", "") or None,
            return_report=True,
            **_mmg_kwargs(step),
        )
        _mmg_entry(entry, warnings, report)
    elif op == "DiscretizeLevelset":
        # MMG (roadmap 3.2.2): the level set's isovalue becomes a conformal interface.
        level_set = _text(step, "LevelSet", "")
        if not level_set:
            raise _err(op, "parameter 'LevelSet' is required")
        mats = _dvec(step, "Materials") or []
        if len(mats) % 4:
            raise _err(
                op,
                "parameter 'Materials' must hold Ref, Split, RefInside, RefOutside "
                "per material",
            )
        rmc = _number(step, "RemoveSmallComponents", -1.0)
        mesh, report = discretize_levelset(
            mesh,
            level_set,
            _number(step, "Isovalue", 0.0),
            materials=[
                (int(mats[k]), mats[k + 1] != 0, int(mats[k + 2]), int(mats[k + 3]))
                for k in range(0, len(mats), 4)
            ],
            remove_small_components=rmc if rmc > 0 else None,
            return_report=True,
            **_mmg_kwargs(step),
        )
        _mmg_entry(entry, warnings, report)
    elif op == "Move":
        # MMG (roadmap 3.2.4): the Lagrangian mode, with ELAS.
        displacement = _text(step, "Displacement", "")
        if not displacement:
            raise _err(op, "parameter 'Displacement' is required")
        moving = _svec(step, "Moving")
        if not moving:
            raise _err(op, "parameter 'Moving' is required")
        mesh, report = move(
            mesh,
            displacement,
            moving,
            lagrangian=int(_number(step, "Lagrangian", 1.0)),
            return_report=True,
            **_mmg_kwargs(step),
        )
        _mmg_entry(entry, warnings, report)
    elif op == "MeshDomain":
        # MMG (roadmap 3.2.3): the planar domain the line cells bound, triangulated.
        holes = _dvec(step, "Holes") or []
        if len(holes) % 2:
            raise _err(op, "parameter 'Holes' must hold x, y per hole")
        hsiz = _number(step, "Hsiz", -1.0)
        mesh, report = mesh_domain(
            mesh,
            hsiz=hsiz if hsiz > 0 else None,
            size=_text(step, "Size", "") or None,
            holes=[(holes[k], holes[k + 1]) for k in range(0, len(holes), 2)],
            return_report=True,
            **_mmg_kwargs(step),
        )
        _mmg_entry(entry, warnings, report)
        entry["NumSubdomains"] = float(report["num_subdomains"])
        entry["NumHoles"] = float(report["num_holes"])
    elif op == "Remesh":
        # Unlike every other step, the output has NO correspondence to the
        # input: new points, new connectivity. point_data/cell_data/regions
        # are dropped (remesh's own contract), field_data carries through.
        n_subdivide = _number(step, "Subdivide", -1)
        mesh, report = remesh(
            mesh,
            int(_number(step, "NumClusters", 0)),
            subdivide=None if n_subdivide < 0 else int(n_subdivide),
            subsample_ratio=_number(step, "SubsampleRatio", 10.0),
            max_subdivide=int(_number(step, "MaxSubdivide", 4)),
            max_iterations=int(_number(step, "MaxIterations", 100)),
            max_repair_passes=int(_number(step, "MaxRepairPasses", 10)),
            metric=_text(step, "Metric", "isotropic"),
            gradation=_number(step, "Gradation", 0.0),
            preserve_boundary=_flag(step, "PreserveBoundary", True),
            max_anisotropy=_number(step, "MaxAnisotropy", 4.0),
            return_report=True,
        )
        entry["NumClusters"] = float(report["num_clusters"])
        entry["NumIterations"] = float(report["num_iterations"])
        entry["SubdivideApplied"] = float(report["subdivide_applied"])
        entry["NumIsolatedClusters"] = float(report["num_isolated_clusters"])
        entry["NumNonManifoldVertices"] = float(report["num_non_manifold_vertices"])
        if (
            report["num_isolated_clusters"] > 0
            or report["num_non_manifold_vertices"] > 0
        ):
            warnings.append(
                f"remesh: {report['num_isolated_clusters']} isolated cluster(s), "
                f"{report['num_non_manifold_vertices']} non-manifold vertex/vertices could "
                "not be repaired; output may be non-manifold near them"
            )
    elif op == "RemeshVolume":
        # Remesh's volumetric sibling: same no-correspondence-with-the-input
        # contract, same lattice-sizing vocabulary as Voxelize/ComputeSdf.
        resolution = _dvec(step, "Resolution")
        bounds = _dvec(step, "Bounds")
        mesh, report = remesh_volume(
            mesh,
            resolution=None if not resolution else [int(v) for v in resolution],
            cell_size=_number(step, "CellSize", 0.0) if "CellSize" in step else None,
            bounds=None if not bounds else list(bounds),
            padding=_number(step, "Padding", 0.0),
            padding_relative=_number(step, "PaddingRelative", 0.1),
            max_cells=int(_number(step, "MaxCells", 20000000.0)),
            max_tets=int(_number(step, "MaxTets", 20000000.0)),
            warp_fraction=_number(step, "WarpFraction", 0.35),
            sign=_text(step, "Sign", "pseudonormal"),
            watertight_check=_text(step, "WatertightCheck", "warn"),
            return_report=True,
        )
        entry["NumTets"] = float(report["num_tets"])
        entry["NumVerticesWarped"] = float(report["num_vertices_warped"])
        entry["NumTetsRejected"] = float(report["num_tets_rejected"])
        entry["NumNonManifoldEdges"] = float(report["num_non_manifold_edges"])
        if report["num_non_manifold_edges"] > 0:
            warnings.append(
                f"remesh_volume: {report['num_non_manifold_edges']} non-manifold boundary "
                "edge(s), from warping (see doc/remesh_volume.md)"
            )
    elif op == "OptimizeVolume":
        # ODT remeshing: relocate vertices AND flip connectivity (2-3/3-2). The
        # point set is invariant, so point_data/field_data and Point regions
        # carry; cell_data + Cell/Side regions are dropped.
        mesh, report = optimize_volume(
            mesh,
            max_iterations=int(_number(step, "MaxIterations", 10.0)),
            relocate=_flag(step, "Relocate", True),
            flip=_flag(step, "Flip", True),
            preserve_boundary=_flag(step, "PreserveBoundary", True),
            min_improvement=_number(step, "MinImprovement", 1e-6),
            return_report=True,
        )
        entry["NumFlips"] = float(report["num_flips"])
        entry["Num23Flips"] = float(report["num_23_flips"])
        entry["Num32Flips"] = float(report["num_32_flips"])
        entry["NumVerticesMoved"] = float(report["num_vertices_moved"])
        entry["NumTets"] = float(report["num_tets"])
        entry["MinQualityBefore"] = float(report["min_quality_before"])
        entry["MinQualityAfter"] = float(report["min_quality_after"])
    elif op == "Generate":
        # A mesh from parameters; like Voxelize it REPLACES the geometry, and it
        # never reads the incoming mesh, so a document may omit Input when this
        # is its first step.
        mesh = _generate(step)
        entry["NumPoints"] = float(len(mesh.points))
        entry["NumCells"] = float(_total_cells(mesh))

    elif op == "Voxelize":
        resolution = _dvec(step, "Resolution")
        bounds = _dvec(step, "Bounds")
        mesh, rep = voxelize(
            mesh,
            resolution=None if not resolution else [int(v) for v in resolution],
            cell_size=_number(step, "CellSize", 0.0) if "CellSize" in step else None,
            bounds=None if not bounds else list(bounds),
            padding=_number(step, "Padding", 0.0),
            padding_relative=_number(step, "PaddingRelative", 0.0),
            fill=_text(step, "Fill", "all"),
            attach_occupancy=_flag(step, "AttachOccupancy", False),
            max_cells=int(_number(step, "MaxCells", 20000000.0)),
            sign=_text(step, "Sign", "pseudonormal"),
            watertight_check="off",
            plane=_text(step, "Plane", "auto"),
            return_report=True,
        )
        entry["NumOccupied"] = float(rep["num_occupied"])

    elif op == "OctreeMesh":
        # A volume mesh of the skin's inside: it REPLACES the input's geometry.
        bounds = _dvec(step, "Bounds")
        mesh, rep = octree_mesh(
            mesh,
            output=_text(step, "Output", "tet"),
            min_size=_number(step, "MinSize", 0.0),
            max_size=_number(step, "MaxSize", 0.0),
            max_level=int(_number(step, "MaxLevel", -1.0)),
            bounds=None if not bounds else list(bounds),
            padding=_number(step, "Padding", 0.0),
            padding_relative=_number(step, "PaddingRelative", 0.1),
            curvature_angle=_number(step, "CurvatureAngle", 15.0),
            feature_angle=_number(step, "FeatureAngle", 30.0),
            corner_angle=_number(step, "CornerAngle", 30.0),
            feature_size=_number(step, "FeatureSize", 0.0),
            gap_cells=_number(step, "GapCells", 2.0),
            gap_outside=_flag(step, "GapOutside", False),
            size_array=_text(step, "SizeArray", "") or None,
            surface_region=_text(step, "SurfaceRegion", "") or None,
            grading_layers=int(_number(step, "GradingLayers", 1.0)),
            skin_grading=_text(step, "SkinGrading", "graded"),
            surface_size=_number(step, "SurfaceSize", 0.0),
            winding_threshold=_number(step, "WindingThreshold", 0.5),
            winding_beta=_number(step, "WindingBeta", 2.0),
            conformity=_text(step, "Conformity", "stuffing"),
            alpha_long=_number(step, "AlphaLong", 0.24799),
            alpha_short=_number(step, "AlphaShort", 0.40173),
            preserve_features=_flag(step, "PreserveFeatures", True),
            fit_iterations=int(_number(step, "FitIterations", 5.0)),
            dihedral_floor_guard=_flag(step, "DihedralFloorGuard", True),
            fill=_text(step, "Fill", "inside"),
            keep_outside=_flag(step, "KeepOutside", False),
            attach_distance=_flag(step, "AttachDistance", True),
            attach_quality=_flag(step, "AttachQuality", False),
            pillow=_flag(step, "Pillow", True),
            hex_grading=_text(step, "HexGrading", "uniform"),
            untangle_iterations=int(_number(step, "UntangleIterations", 50.0)),
            max_cells=int(_number(step, "MaxCells", 20000000.0)),
            max_tets=int(_number(step, "MaxTets", 100000000.0)),
            return_report=True,
        )
        entry["NumLeaves"] = float(rep["num_leaves"])
        entry["NumCells"] = float(rep["num_cells"])
        entry["NumInside"] = float(rep["num_inside"])
        entry["NumSurrogateFaces"] = float(rep["num_surrogate_faces"])
        entry["NumHanging"] = float(rep["num_hanging"])
        entry["NumInverted"] = float(rep["num_inverted"])
        entry["MinDihedral"] = float(rep["min_dihedral"])
        entry["MinScaledJacobian"] = float(rep["min_scaled_jacobian"])

    elif op == "Octree":
        # Voxelize's adaptive sibling: it likewise REPLACES the input's geometry.
        roots = _dvec(step, "RootCells")
        bounds = _dvec(step, "Bounds")
        mesh, rep = octree(
            mesh,
            root_cells=[int(v) for v in roots] if roots else (1, 1, 1),
            min_level=int(_number(step, "MinLevel", 0.0)),
            max_level=int(_number(step, "MaxLevel", 4.0)),
            bounds=None if not bounds else list(bounds),
            padding=_number(step, "Padding", 0.0),
            padding_relative=_number(step, "PaddingRelative", 0.0),
            dimension=int(_number(step, "Dimension", 0.0)),
            sample_distance=_number(step, "SampleDistance", 0.0),
            balance=_text(step, "Balance", "face"),
            fill=_text(step, "Fill", "all"),
            attach_distance=_flag(step, "AttachDistance", False),
            even_odd=_flag(step, "EvenOdd", False),
            max_cells=int(_number(step, "MaxCells", 20000000.0)),
            sign=_text(step, "Sign", "pseudonormal"),
            watertight_check="off",
            return_report=True,
            transition=_text(step, "Transition", "hanging"),
            snap=_flag(step, "Snap", False),
            surface_level=(
                int(_number(step, "SurfaceLevel", -1.0))
                if "SurfaceLevel" in step
                else None
            ),
            curvature_threshold=_number(step, "CurvatureThreshold", 0.0),
            feature_angle=_number(step, "FeatureAngle", 0.0),
            size_array=_text(step, "SizeArray", ""),
            plane=_text(step, "Plane", "auto"),
        )
        entry["NumLeaves"] = float(rep["num_leaves"])
        entry["NumHanging"] = float(rep["num_hanging"])
        entry["NumSnapped"] = float(rep["num_snapped"])
        entry["NumRefused"] = float(rep["num_refused"])

    elif op == "ComputeSdf":
        # Voxelize's sibling: it likewise REPLACES the input's geometry, taking a
        # surface in and handing a filled grid back.
        resolution = _dvec(step, "Resolution")
        bounds = _dvec(step, "Bounds")
        mesh, rep = compute_sdf(
            mesh,
            structure=_text(step, "Structure", "voxel"),
            resolution=None if not resolution else [int(v) for v in resolution],
            cell_size=_number(step, "CellSize", 0.0) if "CellSize" in step else None,
            bounds=None if not bounds else list(bounds),
            padding=_number(step, "Padding", 0.0),
            padding_relative=_number(step, "PaddingRelative", 0.1),
            root_resolution=int(_number(step, "RootResolution", 8.0)),
            max_depth=int(_number(step, "MaxDepth", 4.0)),
            band_cells=_number(step, "BandCells", 1.0),
            record_levels=_flag(step, "RecordLevels", True),
            max_cells=int(_number(step, "MaxCells", 20000000.0)),
            sign=_text(step, "Sign", "pseudonormal"),
            location=_text(step, "Location", "corner"),
            band=_number(step, "Band", 0.0),
            watertight_check="off",
            return_report=True,
        )
        entry["MaxDepth"] = float(rep["max_depth"])
        entry["NumBanded"] = float(rep["num_banded"])

    elif op == "Isosurface":
        isovalues = _dvec(step, "Isovalues")
        if not isovalues:
            isovalues = [_number(step, "Isovalue", 0.0)]
        component = int(_number(step, "Component", -1.0))
        mesh = isosurface(
            mesh,
            _text(step, "Array", ""),
            isovalues,
            component=component if component >= 0 else None,
            record_parent_ids=_flag(step, "RecordParentIds", False),
        )
        entry["ContourCells"] = float(_total_cells(mesh))
        if _total_cells(mesh) == 0:
            warnings.append(
                "the contour is empty; the isovalue lies outside the field's range"
            )
    elif op == "Transform":
        mesh = _apply_transform(mesh, step)
    elif op == "ConvertCells":
        mesh = convert_cells(
            mesh,
            mode=_text(step, "Mode", "linearize"),
            record_parent_ids=_flag(step, "RecordParentIds", False),
        )
    elif op == "Subdivide":
        mesh = subdivide(
            mesh,
            record_parent_ids=_flag(step, "RecordParentIds", False),
        )
    elif op == "Tessellate":
        tess = tessellate(
            mesh,
            levels=int(_number(step, "Levels", 2.0)),
            curved=_flag(step, "Curved", True),
            fields=_flag(step, "Fields", True),
            record_stencil=_flag(step, "RecordStencil", False),
        )
        mesh = tess.mesh
        schema = tess.schema
        entry["NumCurvedSourceCells"] = float(schema["num_curved_source_cells"])
        entry["NumPassThroughSourceCells"] = float(
            schema["num_pass_through_source_cells"]
        )
        entry["NumPoints"] = float(schema["num_points"])
        entry["NumCells"] = float(schema["num_cells"])
        entry["NumFacetsWithBadCount"] = float(
            schema["watertight"]["num_facets_with_bad_count"]
        )
    elif op == "Agglomerate":
        mesh, report = agglomerate(
            mesh,
            target_group_size=int(_number(step, "TargetGroupSize", 8.0)),
            merge_coplanar_faces=_flag(step, "MergeCoplanarFaces", False),
            coplanar_angle=_number(step, "CoplanarAngle", 1.0),
            min_sphericity=_number(step, "MinSphericity", 0.0),
            return_report=True,
        )
        entry["NumFacesMerged"] = float(report["num_faces_merged"])
        entry["NumRejected"] = float(report["num_rejected"])
    elif op == "Crop":
        has_bbox = "Bbox" in step
        has_plane = "Point" in step or "Normal" in step
        has_where = "Where" in step
        if has_bbox + has_plane + has_where != 1:
            raise _err(op, "give one of 'Bbox', 'Point'+'Normal' or 'Where'")
        if has_where and "Mode" in step:
            raise _err(
                op,
                "'Mode' applies to 'Bbox' and 'Point'+'Normal', which test points; "
                "a 'Where' predicate is already one value per cell and has nothing "
                "to reduce",
            )
        if has_bbox:
            kwargs = {"bbox": _dvec(step, "Bbox")}
            if kwargs["bbox"] is None or len(kwargs["bbox"]) != 6:
                raise _err(op, "'Bbox' must be [xmin, ymin, zmin, xmax, ymax, zmax]")
        elif has_plane:
            kwargs = {"plane": (_vec3(step, "Point"), _vec3(step, "Normal"))}
        else:
            kwargs = {
                "where": (
                    _text(step, "Where", ""),
                    _text(step, "Compare", "<"),
                    _number(step, "Value", 0.0),
                )
            }
        mesh = crop(
            mesh,
            record_ids=_flag(step, "RecordIds", False),
            **({} if has_where else {"mode": _text(step, "Mode", "all")}),
            **kwargs,
        )
        entry["CellsKept"] = float(_total_cells(mesh))
    elif op == "ExtractSurface":
        mesh = extract_surface(
            mesh, record_parent_ids=_flag(step, "RecordParentIds", False)
        )
    elif op == "ExtractSkin":
        mesh = extract_skin(mesh, linearize=_flag(step, "Linearize", False))
    elif op == "Reorder":
        mesh = reorder(mesh, method=_text(step, "Method", "rcm"))
    elif op in ("DataDrop", "DataKeep", "DataRename"):
        mesh = _apply_data_manage(mesh, op, step)
    elif op in ("SetsToData", "DataToSets"):
        mesh = _apply_sets_data(mesh, step)
    elif op == "DataCalc":
        spec = _text(step, "Expr", "")
        # "NAME = EXPR" splits on the FIRST '=' (the CLI's documented rule).
        name, sep, expr = spec.partition("=")
        if not sep:
            raise _err(op, "'Expr' must be 'NAME = EXPRESSION'")
        name = name.strip()
        if not name:
            raise _err(op, "'Expr' has an empty output name")
        mesh = data_calc(
            mesh,
            expr.strip(),
            location=_text(step, "Location", "point"),
            output=name,
            overwrite=_flag(step, "Overwrite", False),
        )
    elif op == "DataCondition":
        mesh = data_condition(
            mesh,
            location=_text(step, "Location", "point"),
            keys=_svec(step, "Names"),
            mode=_text(step, "Mode", "clamp"),
            scope=_text(step, "Scope", "component"),
            lo=_number(step, "Lo", 0.0),
            hi=_number(step, "Hi", 1.0),
            nan_policy=_text(step, "NanPolicy", "ignore"),
            nan_replacement=_number(step, "NanReplacement", 0.0),
            suffix=_text(step, "Suffix", ""),
        )
    elif op == "TensorInvariants":
        outputs_str = _text(step, "Outputs", "")
        mesh = tensor_invariants(
            mesh,
            location=_text(step, "Location", "point"),
            keys=_svec(step, "Names"),
            outputs=outputs_str.split(",") if outputs_str else None,
            prefix=_text(step, "Prefix", ""),
            suffix=_text(step, "Suffix", ""),
            overwrite=_flag(step, "Overwrite", True),
        )
    elif op == "ToCell":
        mesh = point_data_to_cell_data(mesh, keys=_svec(step, "Names"))
    elif op == "ToPoint":
        mesh = cell_data_to_point_data(
            mesh,
            keys=_svec(step, "Names"),
            weighted=_text(step, "Weight", "uniform") == "measure",
        )
    else:  # unreachable after _validate_step
        raise ValueError(f"meshio++: pipeline: unknown operation '{op}'")

    steps.append(entry)
    # A no-op when no scope is open (the common case); mirrors what
    # pipe_push_step does on the C++ side, from the step's own INPUT
    # parameters (not the result entry above, which holds counters).
    _prov_add_operation(_render_op(step))
    return mesh


def _project_to_cad_step(mesh, step, entry, model):
    """``ProjectToCad`` on the B-rep at ``Model``, or else on the run's model
    (the C++ engine's branch of pipe_apply_step)."""
    from ._project_to_cad import project_to_cad

    op = step["Op"]
    if "Model" in step:
        path = _text(step, "Model", "")
        fmt = (_filetypes_from_path(pathlib.Path(path)) or [""])[0]
        if fmt not in ("cad_json", "step", "iges"):
            raise _err(
                op, f"'Model' must be a CAD JSON, STEP or IGES file, not '{fmt}'"
            )
        brep = read(path, file_format=fmt).info["brep"]
    elif model is not None and model.kind == "brep":
        brep = model.data
    else:
        why = model.why if model is not None and not model.has() and model.why else ""
        raise _err(
            op,
            "needs 'Model' (a CAD JSON, STEP or IGES file), or the B-rep a CAD Input "
            "carries while no step has changed the mesh" + (": " + why if why else ""),
        )
    mesh, report = project_to_cad(
        mesh,
        brep,
        face_array=_text(step, "FaceArray", "cad:face"),
        edge_array=_text(step, "EdgeArray", "cad:edge"),
        max_distance=_number(step, "MaxDistance", 0.0),
        trim_samples=int(_number(step, "TrimSamples", 16)),
        return_report=True,
    )
    entry.update(
        NumProjected=float(report["num_projected"]),
        NumOnEdges=float(report["num_on_edges"]),
        NumOnCurves=float(report["num_on_curves"]),
        NumOutsideTrim=float(report["num_outside_trim"]),
        MaxMove=float(report["max_move"]),
    )
    return mesh


def _check_keys(obj, where, allowed):
    for key in obj:
        if key not in allowed:
            raise ValueError(
                f"meshio++: pipeline: unknown key '{key}' in {where} "
                f"(known: {', '.join(allowed)})"
            )


def _read_kwargs_from(inp, warnings):
    """``Input.Options`` -> ``read()`` kwargs, validating strictly.

    Shared with ``_sequence.py`` rather than transcribed: two copies of this
    mapping would drift on the next option added, and a sequence document's
    ``Input.Options`` must mean exactly what a single-file one's does.
    """
    read_kwargs = {}
    options = inp.get("Options") or {}
    if not isinstance(options, dict):
        raise ValueError("meshio++: pipeline: Input.Options must be an object")
    _check_keys(
        options,
        "Input.Options",
        ("PointsOnly", "DataArrays", "TimeStep", "Lenient", "Mmap"),
    )
    if options.get("PointsOnly"):
        read_kwargs["points_only"] = True
    if "DataArrays" in options:
        read_kwargs["arrays"] = list(options["DataArrays"])
    if "TimeStep" in options:
        read_kwargs["time_step"] = int(options["TimeStep"])
    if options.get("Mmap") not in (None, "auto", "on", "off"):
        raise ValueError(
            "meshio++: pipeline: Input.Options.Mmap must be 'auto', 'on' or 'off'"
        )
    # Mmap is advisory by contract (a preference, never a demand), so the
    # Python runner -- whose readers do not memory-map -- may ignore it.
    if options.get("Lenient"):
        # Lenient exists to make the *C++* readers tolerant; this runner reads
        # through the Python API, whose per-format fallbacks already tolerate
        # those constructs. Say so rather than silently dropping the flag.
        warnings.append(
            "Input.Options.Lenient has no effect in the Python runner (its "
            "per-format fallbacks already read what lenient would skip)"
        )
    return read_kwargs


def _write_kwargs_from(out, out_path):
    """``Output`` -> ``write()`` kwargs. Shared with ``_sequence.py``.

    ``out_path`` is passed separately because a sequence's output is a pattern
    whose *expanded* name is what carries the extension.
    """
    encoding = out.get("Encoding", "default")
    if encoding not in ("default", "ascii", "binary", "raw_appended"):
        raise ValueError(
            "meshio++: pipeline: Output.Encoding must be 'ascii', 'binary' or "
            "'raw_appended'"
        )
    codec = out.get("Codec")
    if codec is not None and codec not in ("none", "zlib", "lz4", "zstd", "lzf"):
        raise ValueError(
            "meshio++: pipeline: Output.Codec must be 'none', 'zlib', 'lz4', 'zstd' or 'lzf'"
        )
    write_kwargs = {}
    if encoding != "default" or codec is not None:
        out_fmt = out.get("Format")
        if not out_fmt:
            out_fmt = _write_format_for_path(pathlib.Path(str(out_path)))
        if out_fmt is None:
            raise ValueError(
                "meshio++: pipeline: cannot infer the output format needed to "
                "apply Encoding/Codec; set Output.Format"
            )
        # The MCP convert tool's writer-kwargs mapping is the single Python
        # owner of "this format's ascii/binary/codec spelling"; imported
        # lazily -- mcp/_tools.py imports the package, so a module-scope
        # import here would be a cycle. (It never imports the mcp SDK.)
        from .mcp._tools import _variant_kwargs

        write_kwargs = _variant_kwargs(
            out_fmt, "auto" if encoding == "default" else encoding, codec
        )
    if out.get("FloatFormat"):
        write_kwargs["float_fmt"] = out["FloatFormat"]
    return write_kwargs


def _resolve_settings(settings):
    """A dict is used as-is; a string/PathLike is JSON text when it starts
    with '{' (after whitespace), else a settings-file path."""
    if isinstance(settings, dict):
        return settings
    if isinstance(settings, os.PathLike):
        settings = os.fspath(settings)
    if not isinstance(settings, str):
        raise TypeError(
            "meshio++: pipeline: settings must be a dict, JSON text or a path"
        )
    if settings.lstrip().startswith("{"):
        return json.loads(settings)
    with open(settings, "r", encoding="utf-8") as handle:
        return json.load(handle)


class _Model:
    """The CAD model a run carries (roadmap 6.1.1): ``kind`` is ``"brep"`` or
    ``"nurbs"`` and ``data`` the model dict, or both ``None`` and ``why`` the
    reason a model step raises. The C++ engine's ``PipeModel``."""

    def __init__(self, why):
        self.kind = None
        self.data = None
        self.why = why

    def has(self):
        return self.kind is not None

    def set(self, kind, data):
        self.kind, self.data = kind, data

    def drop(self, why):
        self.kind = self.data = None
        self.why = why


def _sniffed_format(path):
    """The format the C++ engine resolves an ambiguous extension to (``.mesh``
    is Medit or MFEM by its content), or ``""`` without the core. Called after
    the Input was read, so the file exists."""
    try:
        from . import _core
    except ImportError:  # pragma: no cover - the numpy-only install
        return ""
    return _core.sniff_format(str(path))


def _iga_triple(step, key, default):
    v = step.get(key)
    if v is None:
        return [default] * 3
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return [int(_number(step, key, default))] * 3
    if (
        not isinstance(v, (list, tuple))
        or not 1 <= len(v) <= 3
        or any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in v)
        or any(float(x) != int(x) for x in v)
    ):
        raise _err(step["Op"], f"'{key}' must be a number or 1 to 3 integers")
    return [int(x) for x in v] + [default] * (3 - len(v))


def _index_size_pairs(step, key):
    """A flat ``[index, size, index, size, ...]`` list as pairs (the C++
    engine's settings values are flat lists, so both read the same)."""
    flat = step.get(key)
    if flat is None:
        return None
    if not isinstance(flat, (list, tuple)) or len(flat) % 2:
        raise _err(step["Op"], f"'{key}' must list index, size pairs")
    pairs = []
    for index, size in zip(flat[0::2], flat[1::2]):
        if isinstance(index, bool) or float(index) != int(index):
            raise _err(step["Op"], f"'{key}': an index must be an integer")
        pairs.append((int(index), float(size)))
    return pairs


def _apply_model_step(step, model, steps, warnings):
    """A model step on the run's model; returns the mesh it makes."""
    from ._iga import iga_refine, iga_sample
    from ._mesh_cad import mesh_cad_surface

    op = step["Op"]
    entry = {"op": op}
    if not model.has():
        raise _err(op, "needs a CAD model, and the pipeline carries none: " + model.why)
    if op == "MeshCadSurface":
        if model.kind != "brep":
            raise _err(
                op, "meshes a B-rep (CAD JSON, STEP or IGES), not an MFEM NURBS model"
            )
        if "TargetSize" not in step:
            raise _err(op, "'TargetSize' is required")
        method = _text(step, "SmoothMethod", "laplacian")
        if method not in ("laplacian", "odt"):
            raise _err(op, "'SmoothMethod' must be 'laplacian' or 'odt'")
        metric = _text(step, "Metric", "mean")
        if metric not in ("mean", "varying"):
            raise _err(op, "'Metric' must be 'mean' or 'varying'")
        background = None
        if "Background" in step:
            background = read(_text(step, "Background", ""))
        metric_background = None
        if "MetricBackground" in step:
            metric_background = read(_text(step, "MetricBackground", ""))
        mesh, report = mesh_cad_surface(
            model.data,
            _number(step, "TargetSize", 0.0),
            max_chordal_deviation=_number(step, "MaxChordalDeviation", 0.0),
            min_points_per_edge=int(_number(step, "MinPointsPerEdge", 0)),
            tolerance=_number(step, "Tolerance", 0.0),
            record_face_ids=_flag(step, "RecordFaceIds", True),
            edges=_flag(step, "Edges", True),
            match_unreferenced_trims=_flag(step, "MatchUnreferencedTrims", True),
            strict=_flag(step, "Strict", False),
            max_points=int(_number(step, "MaxPoints", 50_000_000)),
            min_angle=_number(step, "MinAngle", 0.0),
            smooth_iterations=int(_number(step, "SmoothIterations", 0)),
            smooth_method=method,
            order=int(_number(step, "Order", 1)),
            curvature_deviation=_number(step, "CurvatureDeviation", 0.0),
            cells_across_gaps=_number(step, "CellsAcrossGaps", 0.0),
            grading=_number(step, "Grading", 0.0),
            pole_rings=int(_number(step, "PoleRings", 0)),
            face_sizes=_index_size_pairs(step, "FaceSizes"),
            edge_sizes=_index_size_pairs(step, "EdgeSizes"),
            background=background,
            background_array=_text(step, "BackgroundArray", "size"),
            metric=metric,
            virtual_edges=_flag(step, "VirtualEdges", False),
            reparam_condition_number=_number(step, "ReparamConditionNumber", 0.0),
            model_tolerances=_flag(step, "ModelTolerances", False),
            virtual_faces=_flag(step, "VirtualFaces", False),
            lop_corners=_flag(step, "LopCorners", False),
            metric_background=metric_background,
            metric_array=_text(step, "MetricArray", ""),
            metric_error=_number(step, "MetricError", 0.0),
            untangle=_flag(step, "Untangle", False),
            quads=_text(step, "Quads", "none"),
            quad_pairs=_text(step, "QuadPairs", "matching"),
            return_report=True,
        )
        warnings.extend("MeshCadSurface: " + w for w in report["warnings"])
        entry.update(
            NumPoints=float(report["num_points"]),
            NumTriangles=float(report["num_triangles"]),
            NumFaces=float(report["num_faces"]),
            NumFailed=float(report["num_failed"]),
            NumOpenEdges=float(report["num_open_edges"]),
            NumNonManifoldEdges=float(report["num_non_manifold_edges"]),
            NumRefinementPoints=float(report["num_refinement_points"]),
            NumCurvedStraightened=float(report["num_curved_straightened"]),
            NumReparametrised=float(sum(r["taken"] for r in report["reparametrised"])),
            NumMergedFaces=float(sum(m["merged"] for m in report["merged_faces"])),
            NumLoppedCorners=float(len(report["lopped_corners"])),
            NumUntangled=float(report["num_untangled"]),
            NumQuads=float(report["num_quads"]),
            MinAngleDeg=float(report["min_angle_deg"]),
            MaxAngleDeg=float(report["max_angle_deg"]),
        )
    else:
        # IgaSample, IgaRefine: the model is sampled (refined first), and the
        # run carries on the model the sampling describes.
        sample = dict(
            trim_samples=int(_number(step, "TrimSamples", 16)),
            edges=_flag(step, "Edges", True),
            trim=_flag(step, "Trim", True),
        )
        if op == "IgaRefine":
            mesh = iga_refine(
                model.data,
                subdivisions=_iga_triple(step, "Subdivisions", 1),
                elevate=_iga_triple(step, "Elevate", 0),
                sample_subdivisions=int(_number(step, "SampleSubdivisions", 0)),
                **sample,
            )
        else:
            mesh = iga_sample(
                model.data, subdivisions=int(_number(step, "Subdivisions", 0)), **sample
            )
        model.set(model.kind, mesh.info[model.kind])
        ncells = sum(len(c) for c in mesh.cells)
        entry.update(
            NumPatches=float(len(model.data["patches"])),
            NumPoints=float(len(mesh.points)),
            NumCells=float(ncells),
        )
    steps.append(entry)
    _prov_add_operation(_render_op(step))
    return mesh


def _apply_v2(mesh, step, steps, warnings):
    op = step["Op"]
    entry = {"op": op}
    if op == "Merge":
        extras = [read(path) for path in step["Inputs"]]
        mesh = merge(
            [mesh] + extras,
            weld=_flag(step, "Weld", False),
            atol=_number(step, "Atol", 1e-8),
            source_tag=_flag(step, "SourceTag", True),
            data_policy=_text(step, "DataPolicy", "intersection"),
            drop_duplicate_cells=_flag(step, "DropDuplicateCells", False),
        )
        entry["NumInputs"] = float(1 + len(extras))
    elif op == "Interpolate":
        mesh = interpolate(
            read(step["Inputs"][0]),
            mesh,
            method=_text(step, "Method", "nearest"),
            arrays=_svec(step, "Arrays") or None,
            extrapolate=_flag(step, "Extrapolate", False),
            default_value=_number(step, "DefaultValue", 0.0),
            on_conflict=_text(step, "OnConflict", "error"),
        )
    elif op == "UndoGreen":
        mesh, report = undo_green(read(step["Inputs"][0]), mesh, return_report=True)
        entry.update(
            NumGroupsUndone=float(report["num_groups_undone"]),
            NumCellsRemoved=float(report["num_cells_removed"]),
        )
    else:
        return _apply_step(mesh, step, steps, warnings)
    steps.append(entry)
    _prov_add_operation(_render_op(step))
    return mesh


def _piece_path(pattern, key, part=False):
    token = "{part}" if part else "{key}"
    key = "".join(
        (
            chr(c)
            if (65 <= c <= 90 or 97 <= c <= 122 or 48 <= c <= 57 or c in (45, 46, 95))
            else "_"
        )
        for c in str(key).encode("utf-8")
    )
    if key in ("", ".", ".."):
        key = "_"
    path = str(pattern).replace(token, key)
    if "{" in path or "}" in path:
        raise ValueError("meshio++: pipeline: unsupported Output.Pattern token")
    return path


def _fanout(mesh, step, pattern, out, inp, steps_spec):
    part = step["Op"] == "Partition"
    if part:
        pieces = list(
            enumerate(
                partition(
                    mesh,
                    int(_number(step, "Nparts", 2)),
                    method=_text(step, "Method", "auto"),
                    imbalance=_number(step, "Imbalance", 0.03),
                    mode=_text(step, "Mode", "eco"),
                    seed=int(_number(step, "Seed", 0)),
                    weights=_text(step, "WeightsKey", "") or None,
                    record_ids=_flag(step, "RecordIds", False),
                    ghost_layers=int(_number(step, "GhostLayers", 0)),
                )
            )
        )
    else:
        pieces = list(
            split(
                mesh, by=_text(step, "By", "type"), tag=_text(step, "Tag", "") or None
            ).items()
        )
    inputs = {pathlib.Path(inp).resolve()}
    for spec in steps_spec:
        if spec["Op"] in ("Merge", "Interpolate", "UndoGreen"):
            inputs.update(pathlib.Path(p).resolve() for p in spec["Inputs"])
    paths, unique, kwargs = [], set(), []
    for key, _ in pieces:
        path = _piece_path(pattern, key, part)
        guard = _OUTPUT_PATH_GUARD.get()
        if guard is not None:
            path = guard(path)
        canonical = pathlib.Path(path).resolve()
        if canonical in unique or canonical in inputs:
            raise ValueError(
                f"meshio++: pipeline: Output.Pattern collides with another output or an input: {path}"
            )
        unique.add(canonical)
        paths.append(path)
        kwargs.append(_write_kwargs_from(out, path))
    _prov_add_operation(_render_op(step))
    for path, (_, piece), options in zip(paths, pieces, kwargs):
        write(path, piece, file_format=out.get("Format") or None, **options)
    return {"op": step["Op"], "NumPieces": float(len(pieces))}


def run_pipeline(settings, input_path=None, output_path=None):
    """Run a whole settings pipeline: read ``Input.Path``, apply
    ``Operations`` in order, write ``Output.Path``.

    :param settings: the settings document -- a dict, the JSON text itself, or
        the path of a ``settings.json`` (see ``doc/pipeline.md`` for the
        schema; ops and keys are PascalCase, enum values lowercase).
    :param input_path: overrides ``Input.Path`` (the CLI's ``--input``).
    :param output_path: overrides ``Output.Path`` (the CLI's ``--output``).
    :returns: the run report ``{"steps": [{"op", ...counters}],
        "warnings": [...]}`` -- counter keys are PascalCase, matching the C++
        engine's report exactly.
    :raises ValueError: on any schema error -- an unknown op, an unknown key
        anywhere, a mis-typed parameter -- always naming the offender, never
        silently ignoring it.
    """
    doc = _resolve_settings(settings)
    if not isinstance(doc, dict):
        raise ValueError("meshio++: pipeline: the settings document must be an object")

    # A transient document (a glob/list input, several steps, a {step} output)
    # is routed to the sequence driver rather than rejected, so one entry point
    # serves both and nobody has to know which kind of document they hold. A
    # plain single-file document never reaches that code path at all -- its
    # behaviour below is physically unchanged.
    from ._sequence import _is_sequence_document, run_sequence_pipeline

    if _is_sequence_document(doc, input_path, output_path):
        return run_sequence_pipeline(
            doc, input_path=input_path, output_path=output_path
        )

    _check_keys(
        doc, "the settings document", ("Version", "Input", "Operations", "Output")
    )
    version = doc.get("Version", 1)
    if (
        not isinstance(version, int)
        or isinstance(version, bool)
        or version not in (1, 2)
    ):
        raise ValueError(
            f"meshio++: pipeline: unsupported Version {version!r} (this build knows 1 and 2)"
        )

    ops_spec = doc.get("Operations")
    generated = (
        "Input" not in doc
        and isinstance(ops_spec, list)
        and bool(ops_spec)
        and isinstance(ops_spec[0], dict)
        and ops_spec[0].get("Op") == "Generate"
    )
    if generated:
        # The first step builds the mesh from parameters: nothing to read.
        doc = {**doc, "Input": {}}
    if "Input" not in doc or not isinstance(doc["Input"], dict):
        raise ValueError("meshio++: pipeline: Input is required and must be an object")
    if "Output" not in doc or not isinstance(doc["Output"], dict):
        raise ValueError("meshio++: pipeline: Output is required and must be an object")
    inp, out = doc["Input"], doc["Output"]
    _check_keys(inp, "Input", ("Path", "Format", "Options"))
    _check_keys(
        out,
        "Output",
        ("Path", "Format", "Encoding", "Codec", "FloatFormat")
        + (("Pattern",) if version == 2 else ()),
    )
    if "Path" in out and "Pattern" in out:
        raise ValueError(
            "meshio++: pipeline: Output.Path and Output.Pattern are mutually exclusive"
        )
    if "Pattern" in out and not any(t in out["Pattern"] for t in ("{key}", "{part}")):
        raise ValueError("meshio++: pipeline: Output.Pattern requires {key} or {part}")

    in_path = input_path if input_path is not None else inp.get("Path")
    out_path = (
        output_path if output_path is not None else out.get("Pattern", out.get("Path"))
    )
    if not in_path and not generated:
        raise ValueError("meshio++: pipeline: Input.Path is required")
    if not out_path:
        raise ValueError("meshio++: pipeline: Output.Path is required")

    warnings = []
    read_kwargs = _read_kwargs_from(inp, warnings)

    steps_spec = doc.get("Operations", [])
    if not isinstance(steps_spec, list):
        raise ValueError("meshio++: pipeline: Operations must be an array of steps")
    fanout = any(t in str(out_path) for t in ("{key}", "{part}"))
    if version == 2 and not fanout and any(t in str(out_path) for t in ("{", "}")):
        raise ValueError(
            "meshio++: pipeline: unsupported Version 2 output token; transient sequence schemas remain Version 1"
        )
    for i, step in enumerate(steps_spec):
        _validate_step(step, version)
        if step["Op"] == "Split" and (not fanout or i + 1 != len(steps_spec)):
            raise ValueError(
                "meshio++: pipeline: Split must be terminal and requires Output.Pattern with {key}"
            )
        if (
            step["Op"] == "Partition"
            and not fanout
            and ("RecordIds" in step or "GhostLayers" in step)
        ):
            raise ValueError(
                "meshio++: pipeline: RecordIds/GhostLayers require partition fan-out"
            )
    if fanout:
        if (
            version != 2
            or not steps_spec
            or steps_spec[-1]["Op"] not in ("Split", "Partition")
        ):
            raise ValueError(
                "meshio++: pipeline: Output.Pattern requires Version 2 and a terminal Split/Partition"
            )
        part = steps_spec[-1]["Op"] == "Partition"
        if ("{part}" if part else "{key}") not in str(out_path):
            raise ValueError(
                "meshio++: pipeline: Output.Pattern needs {key} for Split or {part} for Partition"
            )
        _piece_path(out_path, "probe", part)

    write_kwargs = _write_kwargs_from(out, out_path)

    # A no-op unless the caller opened a provenance scope around this call
    # (see meshioplusplus._provenance) -- the pipeline never opens one
    # itself, matching the C++ engine's own posture (opt-in everywhere).
    model = _Model("the pipeline starts from Generate, not a CAD read")
    if generated:
        mesh = Mesh(np.zeros((0, 3)), [])
    else:
        in_fmt = (
            inp.get("Format")
            or (_filetypes_from_path(pathlib.Path(str(in_path))) or [""])[0]
        )
        _prov_set_source(str(in_path), in_fmt)

        mesh = read(in_path, file_format=inp.get("Format") or None, **read_kwargs)
        info = getattr(mesh, "info", None) or {}
        model_fmt = inp.get("Format") or _sniffed_format(in_path) or in_fmt
        if model_fmt not in _MODEL_FORMATS:
            model.drop(
                f"Input format '{model_fmt}' carries none (read CAD JSON, STEP, IGES "
                "or an MFEM NURBS mesh)"
            )
        elif model_fmt == "mfem":
            if info.get("nurbs") is not None and len(info["nurbs"].get("patches", ())):
                model.set("nurbs", info["nurbs"])
            else:
                model.drop(f"Input '{in_path}' is an MFEM mesh with no NURBS patches")
        else:
            model.set("brep", info.get("brep"))

    steps = []
    for i, step in enumerate(steps_spec):
        if fanout and i + 1 == len(steps_spec):
            steps.append(_fanout(mesh, step, out_path, out, in_path, steps_spec))
            return {"steps": steps, "warnings": warnings}
        if step["Op"] in _MODEL_OPS:
            mesh = _apply_model_step(step, model, steps, warnings)
            continue
        if step["Op"] == "ProjectToCad":
            # Projected onto it, the mesh matches the model again: it stays.
            entry = {"op": "ProjectToCad"}
            mesh = _project_to_cad_step(mesh, step, entry, model)
            steps.append(entry)
            _prov_add_operation(_render_op(step))
            continue
        mesh = (_apply_v2 if version == 2 else _apply_step)(mesh, step, steps, warnings)
        if step["Op"] == "FitCad":
            # The fitted solid is the model the rest of the run carries, when
            # the backend gave one (the native legacy mode gives none).
            brep = (getattr(mesh, "info", None) or {}).get("brep")
            if brep is not None and len(brep.get("faces", ())):
                model.set("brep", brep)
            else:
                model.drop(
                    "'FitCad' gave no B-rep (its native legacy freeform mode fits none)"
                )
        elif model.has():
            model.drop(
                f"step {i + 1} '{step['Op']}' changed the mesh, which drops the model"
            )

    out_fmt = out.get("Format") or _write_format_for_path(pathlib.Path(str(out_path)))
    out_encoding = out.get("Encoding", "default")
    _prov_set_target(
        out_fmt,
        "" if out_encoding in (None, "default") else out_encoding,
        out.get("Codec") or "",
        out.get("FloatFormat") or "",
    )

    write(out_path, mesh, file_format=out.get("Format") or None, **write_kwargs)
    return {"steps": steps, "warnings": warnings}
