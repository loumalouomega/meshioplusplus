# The settings pipeline

Since v9.11.0 meshio++ can run a whole *chain* of operations declaratively: one `settings.json` document describes read → operations → write, and every surface executes it — both CLIs (the `pipeline` verb), Python (`meshioplusplus.run_pipeline`), C (`mio_pipeline_run_file`/`_json`), Fortran, Julia, R, WASM (`runPipeline`) and the MCP server (the `pipeline` tool).

```bash
meshioplusplus pipeline settings.json            # either CLI
meshioplusplus pipeline settings.json --input other.msh --output out.vtu
meshioplusplus pipeline settings.json --json     # machine-readable report
```

```json
{
  "Version": 1,
  "Input":  { "Path": "bracket.msh", "Format": "gmsh" },
  "Operations": [
    { "Op": "Transform",    "RotateAxis": [0, 0, 1], "RotateDegrees": 45 },
    { "Op": "Gradient",     "Array": "temperature", "Operator": "gradient" },
    { "Op": "ConvertCells", "Mode": "simplexify" },
    { "Op": "Refine",       "Array": "temperature:gradient", "Compare": ">",
                            "Value": 0.5, "Closure": "redgreen" },
    { "Op": "Clean" },
    { "Op": "Quality" }
  ],
  "Output": { "Path": "out.vtu", "Encoding": "binary", "Codec": "zlib" }
}
```

![A settings.json is parsed strictly, every step is validated before the input is read, the chain runs and writes through registry_write_ex, and every surface drives the same run_pipeline_steps](/diagrams/pipeline_flow.svg)

## The rules

- **Vocabulary is PascalCase** for op names and parameter keys (`"Op": "ConvertCells"`, `"RemoveOrphans": true`). Enum *values* keep the exact lowercase spellings the rest of meshio++ uses (`"simplexify"`, `"redgreen"`, `"cell"`, `"rcm"`). Refine's comparison key is `Compare`, never `Op` — `Op` is the step discriminant.
- **Parsing is strict.** An unknown op, an unknown key on a step, an unknown top-level key, or a mis-typed value is an error naming the offender — never silently ignored (the same rule `registry_write_ex` applies to `Output` options a format cannot honour).
- **Version 1 runs over one mesh at a time.** Multi-input/output operations are rejected, and `Partition` attaches labels. [Version 2](#version-2-spatial-multi-mesh-steps) adds file-backed inputs and terminal spatial fan-out without changing v1 semantics. `Diff`, `Shrinkwrap` and `ConservativeInterpolate` remain CLI/API-only multi-mesh operations.
- Steps are validated **before** the input is read — a typo in step 7 never costs reading a 10 GB mesh first.
- **The CAD model travels with the mesh as far as the steps that use it** (v17.17.0, roadmap 6.1.1). A CAD read (CAD JSON, STEP, IGES) returns a B-rep beside its sampled mesh, and an MFEM NURBS read a NURBS model; the model steps `MeshCadSurface`, `IgaSample` and `IgaRefine` act on it, and `FitCad` hands on the solid it fits. `ProjectToCad` keeps it too, since it puts the mesh back on the model. Every other step drops it, since nothing else keeps a model in step with the mesh it changes, so a model step after one fails naming the step that dropped it (`step 2 'Smooth' changed the mesh, which drops the model`); a model step on an Input that carries none fails naming its format. Output writes the mesh, never the model. A step applied on its own (`apply_pipeline_step`, `run_pipeline_steps`, a sequence, the WASM viewer's `convertSurfaceOps`) has no Input read behind it, so the model steps fail there by name.
- The run returns a **report**: `{"steps": [{"op", ...counters}], "warnings": [...]}` with PascalCase counter keys (`PointsWelded`, `SectionFaces`, `NumSkipped`, ...).

## Top-level schema

| Key | Required | Meaning |
| --- | --- | --- |
| `Version` | no (default 1) | schema version; `1` or `2` |
| `Input` | yes | `{Path, Format?, Options?}` — `Format` defaults from the extension, with the `sniff_format` read fallback |
| `Operations` | no (default `[]`) | the step array; empty = a plain convert |
| `Output` | yes | `{Path, Format?, Encoding?, Codec?, FloatFormat?}` |

`Input.Options` narrows the read (see [selective reads](selective_read.md)): `PointsOnly` (bool), `DataArrays` (string array; absent = all, `[]` = none), `TimeStep` (int), `Lenient` (bool), `Mmap` (`"auto" | "on" | "off"`).

`Output` maps onto `registry_write_ex`: `Encoding` (`"ascii" | "binary" | "raw_appended"`; `raw_appended` is VTU only, see [VTU](formats/vtu.md)), `Codec` (`"none" | "zlib" | "lz4" | "zstd"`, VTU/VTP block codecs), and `FloatFormat` (a printf-style float format for ASCII writers that take one). An option the output format cannot honour is an error.

## Operations

| `Op` | Parameters (defaults) | Notes / counters |
| --- | --- | --- |
| `Quality` | — | attaches `quality:*` cell data |
| `Clean` | `Weld` (false), `Atol` (1e-8), `RemoveOrphans` (true), `DropDegenerate` (true), `DropDuplicateCells` (true) | counters `PointsWelded`, `PointsRemovedOrphan`, `CellsDroppedDegenerate`, `CellsDroppedDuplicate` |
| `Smooth` | `Method` ("taubin"), `Iterations` (10), `Lambda` (method default), `Mu` (−0.34), `FixBoundary` (true), `PreserveFeatures` (true), `FeatureAngle` (30), `GuardInversion` (true) | counters `NumNodesMoved`, `MaxDisplacement`, `NumSkippedInversion` |
| `Refine` | `Levels` (1), `Cells`, `Region`, `Array` + `Compare` ("<") + `Value` (0), `Closure` ("redgreen"), `RecordLevels` (false), `RecordHierarchy` (false) | at most one selector; `Compare`, not `Op`; `RecordHierarchy` attaches `refine:cell_id`/`refine:parent_id` and forces `refine:entity` |
| `Decimate` | `Ratio` \| `TargetFaces` \| `MaxError` (none → `Ratio` 0.5), `Placement` ("optimal"), `PreserveBoundary` (true), `PreserveFeatures` (true), `FeatureAngle` (30) | counters `FacesRemoved`, `CollapsesRejected` |
| `Adapt` | `Hsiz` \| `Size` \| `Metric` (none: MMG sizes from the mesh), `Hmin`, `Hmax`, `Hausd`, `Hgrad`, `Hgradreq` (MMG's defaults), `AngleDetection` (45), `Optim`, `NoInsert`, `NoSwap`, `NoMove`, `NoSurf`, `Nreg`, `Opnbdy` (false), `Frozen` (a point region), `Required` (region names) | [MMG](/mmg) remeshing (a build with `MESHIOPLUSPLUS_WITH_MMG`); counters `NumPointsIn`, `NumPointsOut`, `NumElementsIn`, `NumElementsOut`, `LowFailure`, `QualityMin`, `QualityMean`; MMG's warnings join the run's |
| `DiscretizeLevelset` | `LevelSet` (a `point_data` name, required), `Isovalue` (0), `RemoveSmallComponents` (off), `Materials` (a flat `[Ref, Split, RefInside, RefOutside, ...]` list), and `Adapt`'s MMG parameters | [MMG](/mmg#discretize_levelset) embeds the isovalue conformally: `medit:ref`, the cell regions `inside`/`outside` and the side region `interface`; `Adapt`'s counters |
| `MeshDomain` | `Hsiz` \| `Size`, `Holes` (a flat `[x0, y0, x1, y1, ...]` list), and `Adapt`'s MMG parameters | [MMG](/mmg#mesh_domain)'s `MMG2D_mmg2dmesh` triangulates the planar domain the `line` cells bound, one cell region `subdomain:<k>` per enclosed region; `Adapt`'s counters plus `NumSubdomains` and `NumHoles` |
| `Move` | `Displacement` (a `point_data` vector, required), `Moving` (region names, required), `Lagrangian` (1), and `Adapt`'s MMG parameters | [MMG](/mmg#move)'s Lagrangian mode (an MMG built with ELAS) moves the `Moving` regions' boundary facets by the displacement and remeshes; `Adapt`'s counters |
| `DecimateVolume` | `Ratio` \| `TargetCells` \| `MaxError` (none → `Ratio` 0.5), `Placement` ("optimal"), `PreserveBoundary` (**false**, unlike `Decimate`), `PreserveFeatures` (true), `FeatureAngle` (30) | counters `TetsRemoved`, `CollapsesRejected` |
| `Partition` | `Nparts` (2), `Method` ("auto"), `Imbalance` (0.03), `Mode` ("eco"), `Seed` (0), `WeightsKey` | attaches `partition:part`; counter `Nparts` |
| `Slice` (alias `Section`) | `Point`, `Normal` (required), `RecordParentIds` (false) | counter `SectionFaces`; warns when the plane misses |
| `Gradient` | `Array` (required), `Operator` ("gradient"), `Method` ("green-gauss"), `Location` ("cell"), `Output`, `Component` | counters `NumSkipped`, `NumFallback`; warns on skipped cells |
| `Normals` | `PointNormals` (true), `CellNormals` (false), `Weight` ("angle" \| "area"), `SplitAngle` (absent → one smooth normal per point; a number of degrees in `[0, 180]` splits vertices at creases), `RecordParentIds` (false), `Region` | attaches `normals`; counters `NumIsolated`, `NumUndefined`, `NumDegenerate`, `NumSplitPoints`, `NumAddedPoints`; warns when the input winds inconsistently and no split was asked for; see [normals](/normals) |
| `FeatureEdges` | `FeatureAngle` (30), `Feature`, `Boundary`, `NonManifold`, `Inconsistent` (all true), `Region` | replaces the mesh by its feature edges, `line` cells with `feature:kind`/`feature:angle`; counters `NumFeature`, `NumBoundary`, `NumNonManifold`, `NumInconsistent`; see [feature edges](/feature_edges) |
| `Segment` | `FeatureAngle` (30), `NormalTolerance` (0), `CurvatureTolerance` (0), `MinPatchCells` (1), `Region`, `Prefix` (`patch`), `Method` (`angle`; `primitive` ends a patch where its [primitive](/segment#primitive-method) does), `Tolerance` (`1e-3` of the diagonal, `primitive` only) | a data step: adds `segment:patch`, `segment:adjacency` and one cell region per patch ([segment](/segment)); reports `NumPatches`, `NumMerged`, `NumFeatureEdges`, `NumBoundaryEdges` |
| `FitPrimitives` | `PatchData` (`segment:patch`), `Tolerance` (`1e-3`), `Relative` (true), `MinInlierFraction` (0.9), `Types` (all), `MaxIterations` (50) | a data step: adds `fit:*` cell data and `fit:patch_*` field data ([fit primitives](/fit_primitives)); reports `NumPatches`, `NumPlane` … `NumTorus`, `NumUnfitted` |
| `FitCad` | `Path` (none), `Tolerance` (0: automatic), `ExactEdges` (true), `Types` (all; `[]` facets everything), `MeshDeflection` (0: automatic), `MaxDemotionRounds` (8), `Freeform` (true), `FreeformDegree` (3), `FreeformMaxControl` (64), `FreeformSplit` (true), `Fillets` (true), `Threads` (true), `MergeFaces` (true), `Backend` (`auto`, `occt` or `native`) | the fitted skin becomes a CAD solid through the optional OpenCASCADE backend, or the kernel-free [native backend](/fit_cad#the-native-backend) without it ([fit CAD](/fit_cad)), written as STEP to `Path` when given; the mesh continues as the solid re-meshed with `fit_cad:*` cell data (native: the skin itself); reports `Native` (1 when the native backend built it), `Valid`, `VolumeRatio`, `MaxDeviation`, `NumFaces`, `NumBspline`, `NumSplitPatches`, `NumFillet`, `NumHelical`, `NumVariableFillet`, `NumBlend`, `NumRoundedThread`, `NumVertexBlend`, `NumChamfer`, `NumConicBlend`, `NumTaperedThread`, `NumVariablePitchThread`, `NumFacetedFaces`, `NumDemoted`. `Backend: "occt"` raises naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON` in a build without it. `RefineSegmentation` (false; v17.11.0) first refines the patches as the `RefinePatches` step does, with `MaxRefineIterations` (10), `RefineEpsilon` (`1e-3`), `NormalWeight` (1), `AllowMerge` (true), `AllowSplit` (true), `RansacSeed` (an integer in `[0, 2^53]`; default `0x6d657368696f2b2b`), `RefineTolerance` (`1e-3`) and `RelativeTolerance` (true); only then does the report add `Refined` (1), `NumRefineSplits`, `NumRefineMerges`, `ExplainedAreaFraction` and `RefinementKept` (1 when the solid was built from the refined patches), so a report without it is unchanged. `FreeformGeneral` (false) tries a patch the plane and cylinder freeform fit leaves unexplained as a general B-spline face, with `FreeformTolerance` (0: a quarter of the tolerance), `FreeformSmoothing` (`1e-12`), `MaxFreeformPoles` (40), `MaxFreeformIterations` (30) and `MaxFreeformSplits` (2), and adds `NumFreeformGeneral`, `NumFreeformRefinements` and `AvgFreeformPoles`; OpenCASCADE only, `SewTolerance` (0: the tolerance) sets the sewing tolerance, `ExtendSurfaces` (false; adds `NumExtendedEdges`) extends freeform surfaces past their patches, `SnapJunctions` (false; adds `NumSnappedJunctions`) places the corners where analytic faces meet by least squares, and `FaceDeviation` (false; adds `HausdorffMax` and `HausdorffRms`) measures each face against the skin as `fit_cad:face_hausdorff` / `fit_cad:face_rms` field data. Every run reports `Closed` (1 when every shell is closed) after `NumDemoted`. `FreeformMode` (`legacy`, so a report without it is unchanged; `adaptive`, `network`, `auto`; v17.12.0) picks the [freeform mode](/fit_cad#freeform-modes), with `FreeformParameterization` (`auto`), `ExtraordinaryVertex` (`g1x`), `LocalRefine` (`knots`), `FreeformG1ToleranceDeg` (1), `FreeformMaxControlPoints` (400), `LayoutTargetPatchSize` (0.1), `LayoutMaxValenceDeviation` (2), `ResegmentRounds` (2), `LayoutSnapToFeatures` (true), `ValidateBySurfaceSampling` (true) and `ConstrainPrimitives` (false, the [constrained primitives](/fit_cad#constrained-primitives)); only outside `legacy` does the report add `BrepChecked`, `FacetAreaFraction`, `NumFreeformAdaptive`, `NumNetworkPatches`, `NumNetworkEdges`, `NumExtraordinaryVertices`, `NumResegmented`, `NumEscalations` and `MaxG1AngleDeg`, only with `ConstrainPrimitives` `NumConstrainedPrimitives` and `NumConstraintGroups`, and the mesh continues with the `fit_cad:face_table` field data |
| `RefinePatches` | `Tolerance` (`1e-3`), `Relative` (true: a fraction of the diagonal), `NormalWeight` (1), `Epsilon` (`1e-3`), `FeatureAngle` (30), `MaxIterations` (10), `MinPatchCells` (3), `RansacSeed` (an integer in `[0, 2^53]`, the most a JSON number holds exactly; default `0x6d657368696f2b2b`), `Types` (all), `AllowSplit` (true), `AllowMerge` (true), `Prefix` (`patch`) | a data step (v17.11.0): refines the `segment:patch` labels of a prior `Segment` step by variational shape approximation (`refine_patches`) — boundaries move, freeform patches are split by seeded RANSAC and neighbours one primitive explains merged — and rewrites `segment:patch`, `segment:adjacency`, the patch regions and the `fit:*` data, adding `refine:*` field data; points and connectivity are untouched. Raises without `segment:patch`. Reports `NumPatches`, `NumIterations`, `NumSplits`, `NumMerges`, `NumPrimitive`, `NumFreeform`, `ExplainedAreaFraction`, `HausdorffMax` and `Stop`, why it stopped as a number (report counters are numbers): 0 `unchanged`, 1 `small_gain`, 2 `iteration_cap`, 3 `energy_rose` |
| `EditRegions` | `Edit` (`union` \| `intersection` \| `difference` \| `rename` \| `retag` \| `delete`), `Inputs` (region names), `Output`, `Kind` (pins the inputs' kind), `Dim`, `Tag`, `KeepInputs` (true) | one region edit, points/cells/data untouched — a list of edits is a list of steps; counter `NumRegions`; see [editing regions](/regions#editing-regions) |
| `RegionAdjacency` | `Regions` (optional ordered list of Cell-region names) | replaces the mesh with conforming shared facets; when fewer than two Cell regions exist, multiple cell blocks are treated as groups; counters `NumFacets` and `Measure`; see [region adjacency](/region_adjacency) |
| `FindInterface` | `RegionA`, `RegionB` (Cell-region names); `Mode` ("conforming" or "proximity"), `Master` ("a" or "b"), `GapTolerance` (0), `AngleTolerance` (30), `OverlapTolerance` (0) | replaces the mesh with matched master facets; counters `NumPairs`, `Area`, `MaxGap`, `UnmatchedA`, `UnmatchedB`; this pipeline form selects both parts from the current mesh |
| `SplitInterface` | `Region` (Side-region name), `AddCohesive` (false) | duplicates node fans along the Side region and optionally inserts cohesive cells; counters `NumDuplicatedPoints`, `NumCohesiveCells`; see [interfaces and contact](/region_adjacency) |
| `QualityGate` | `Require` (threshold texts, e.g. `["scaled_jacobian >= 0.2"]`), `MaxInverted` (0), `MaxDegenerate` (0; negative disables either) | a gate, not a transform: the mesh passes through untouched and a failed check stops the pipeline with the summary as its error; counters `NumChecks`, `NumFailed`; see [quality gate](/mesh_quality#quality-gate) |
| `CheckConformity` | `Tolerance` (1e-9, a fraction of the facets' bounding-box diagonal), `Require` (true), `AttachFlags` (false) | a gate by default: a mesh that is not conforming (a hanging node, T-junction, partly overlapping or unwelded facet pair) stops the pipeline, naming the counts; with `Require` false it only reports, and `AttachFlags` writes `conformity:hanging` onto the mesh that flows on; counters `NumHangingNodes`, `NumTJunctions`, `NumPartialOverlaps`, `NumUnweldedFacets`; see [conformity](/check_conformity) |
| `Snapshot` | `Path` (required: `.png`, `.txt`, `.ansi`, `.html` or `.cast`), `Width` (320), `Height` (240), `Supersample` (1), `Cols` (100), `Rows` (40), `Encoding` ("halfblock" \| "quadrant" \| "sextant" \| "braille" \| "ascii"), `ColorDepth` ("truecolor" \| "256" \| "16" \| "mono"), `PngCompress` (0), the camera `View`, `Azimuth`, `Elevation`, `Roll`, `Perspective`, `Fov`, `Zoom`, `Shading`, `Edges`, `FeatureAngle`, the colours `Background` and `FillColor` (`#rrggbb[aa]` or `none`), colouring by `ColorBy`, `Component`, `Cmap`, `VMin`, `VMax`, `Colorbar`, `Axes`, `ScaleBar`, and the field options `Reduce`, `Expr`, `ClipLow`, `ClipHigh`, `Symmetric`, `Scale`, `ScaleThreshold`, `Categorical`, `ColorRegions`, `CategoryEdges`, `Isolines`, `IsoLevels`, `Vectors`, `VectorCount`, `Warp`, `WarpScale`, `WarpOutline`, `Diagnostic`, `QualityMetric`, and the cut-away `Cutaway` (a flat list of numbers, six per plane, at most two: a point, then the normal of the side kept) with `CutawayTint`, and the streamlines `Streamlines` (a vector point array), `StreamSeeds` (40), `StreamLength` (0.5) and `StreamColor`, and the [theme](./tui#the-synthwave-theme) `Theme` ("none" \| "synthwave"), `GridFloor`, `Bloom`, `Fringe` and `Scanlines` | a side output, not a transform: the mesh passes through untouched and the [software rasterizer](/tui) writes the frame as it stands at this point, so a `Smooth` or `Decimate` can be shown before and after in a report; counters `Width`, `Height`. The parameters are the options of `snapshot`, with the same names in PascalCase |
| `EstimateError` | `Array` (required), `Method` ("zz"), `Marking` ("none" \| "absolute" \| "fraction" \| "dorfler"), `MarkingValue` (0), `Output`, `Marked` | attaches `error:zz` (and `error:marked` when `Marking` isn't "none"); counters `GlobalError`, `NumSkipped`, `NumMarked`; warns on skipped cells |
| `ErrorSize` | `TargetError` (required, positive), `Array` ("error:zz"), `Order` (1), `MinSize` (0: no bound), `MaxSize` (0: no bound), `Output` ("size") | a data step (roadmap 6.10.2): the per-cell indicator an `EstimateError` step attached becomes a per-point size array by the Zienkiewicz-Zhu adaptive procedure ([error estimation](/error)), what a later mesh_cad_surface's `background` reads; geometry untouched, an existing `Output` replaced. Counters `NumCells`, `NumSkipped`, `NumClamped`, `MinSize`, `MaxSize`; warns on skipped cells |
| `Isosurface` | `Array` (required), `Isovalue` (0) or `Isovalues`, `Component`, `RecordParentIds` (false) | counter `ContourCells`; warns when empty |
| `Transform` | exactly one of `Translate[3]`, `Scale` (number or `[3]`), `RotateAxis[3]`+`RotateDegrees`, `Matrix[16]` (row-major), `ScaleUnits` (a factor, e.g. `0.001` for mm→m); plus `RotateData` (false) | |
| `ConvertCells` | `Mode` ("linearize" \| "simplexify" \| "elevate"), `RecordParentIds` (false) | |
| `Subdivide` | `RecordParentIds` (false) | one polyhedral child per 3D cell face, connected to a new interior point; no per-type template table |
| `Tessellate` | `Levels` (2), `Curved` (true), `Fields` (true), `RecordStencil` (false) | isoparametric subdivision of the curved cell types (`triangle6`, `quad8`, `quad9`, `tetra10`, `hexahedron27`) onto a reference lattice into `triangle`/`tetra` cells, every other cell passing through, with the `tessellate:*` provenance arrays ([curved tessellation](/tessellation)); counters `NumCurvedSourceCells`, `NumPassThroughSourceCells`, `NumPoints`, `NumCells`, `NumFacetsWithBadCount` (zero for a conforming tessellation) |
| `Agglomerate` | `TargetGroupSize` (8), `MergeCoplanarFaces` (false), `CoplanarAngle` (1), `MinSphericity` (0) | greedy seed-and-grow over the shared-face dual, merging face-adjacent cells into one polyhedron per group; the many-to-one counterpart to `Subdivide`; counters `NumFacesMerged`, `NumRejected`; see [agglomerate](/agglomerate) |
| `Crop` | one of `Bbox[6]` (`[xmin,ymin,zmin,xmax,ymax,zmax]`), `Point[3]`+`Normal[3]`, or `Where` (a scalar `cell_data` name) + `Compare` ("<") + `Value` (0); `Mode` ("all" \| "any", **not** with `Where`), `RecordIds` (false) | counter `CellsKept` |
| `ExtractSurface` | `RecordParentIds` (false) | |
| `ExtractSkin` | `Linearize` (false) | |
| `Reorder` | `Method` ("rcm" \| "morton" \| "hilbert") | |
| `DataDrop` | `Point`/`Cell`/`Field` (string arrays), `IgnoreMissing` (false) | |
| `DataKeep` | `Point`/`Cell`/`Field` — only locations mentioned are touched; `[]` drops all there | |
| `DataRename` | `Point`/`Cell`/`Field`, entries `"OLD:NEW"` (split on the **last** colon — names carry colons: `gmsh:physical`) | |
| `SetsToData` | `Location` (`cell`; `point` or `cell`), optional `Name`, `Join` (`-`), `Order` (all distinct set names once; omitted/empty uses mesh set order) | Scalar Int64 labels; later overlaps win, uncovered rows are -1; removes converted regions |
| `DataToSets` | `Location` (`cell`; `point` or `cell`), required `Key` | Scalar integer field → region-backed sets; removes source field |
| `DataCalc` | `Expr` (`"NAME = EXPRESSION"`, split on the **first** `=`), `Location` ("point"), `Overwrite` (false) | |
| `DataCondition` | `Mode` ("clamp" \| "normalize" \| "standardize"), `Location` ("point"), `Names`, `Scope` ("component" \| "magnitude"), `Lo` (0), `Hi` (1), `NanPolicy` ("ignore"), `NanReplacement` (0), `Suffix` | |
| `TensorInvariants` | `Location` ("point"), `Names` (every 6- or 9-component array by default), `Outputs` ("mises,principal,hydrostatic,deviatoric" by default), `Prefix`, `Suffix`, `Overwrite` (true) | |
| `ToCell` / `ToPoint` | `Names`; `ToPoint` also `Weight` ("uniform" \| "measure") | |
| `MeshCadSurface` | `TargetSize` (required), `MaxChordalDeviation` (0), `MinPointsPerEdge` (0), `Tolerance` (0: automatic), `RecordFaceIds` (true), `Edges` (true), `MatchUnreferencedTrims` (true), `Strict` (false), `MaxPoints` (50000000), `MinAngle` (0), `SmoothIterations` (0), `SmoothMethod` ("laplacian" \| "odt"), `Order` (1), `CurvatureDeviation` (0), `CellsAcrossGaps` (0), `Grading` (0: 1.5), `PoleRings` (0), `FaceSizes` and `EdgeSizes` (flat `[index, size, …]` lists, 0-based), `Background` (a mesh path) and `BackgroundArray` ("size"), `Metric` ("mean" \| "varying"), `VirtualEdges` (false), `ReparamConditionNumber` (0), `ModelTolerances` (false), `VirtualFaces` (false), `LopCorners` (false), `MetricBackground` (a mesh path), `MetricArray` and `MetricError` (0), `Untangle` (false), `Quads` (`none`, `dominant` or `all`), `QuadPairs` (`matching` or `greedy`) | a model step: meshes the B-rep the Input read carries ([mesh generation](/mesh_cad)), the options of `mesh_cad_surface` in PascalCase; the mesh becomes the surface mesh. Counters `NumPoints`, `NumTriangles`, `NumFaces`, `NumFailed`, `NumOpenEdges`, `NumNonManifoldEdges`, `NumRefinementPoints`, `NumCurvedStraightened`, `MinAngleDeg`, `MaxAngleDeg`, `NumReparametrised`, `NumMergedFaces`, `NumLoppedCorners`, `NumUntangled`, `NumQuads`; the operation's warnings join the run's, prefixed `MeshCadSurface:` |
| `IgaSample` | `Subdivisions` (0: automatic), `TrimSamples` (16), `Trim` (true), `Edges` (true) | a model step: samples the B-rep or NURBS model into Lagrange cells ([isogeometric layer](/iga)); counters `NumPatches`, `NumPoints`, `NumCells` |
| `IgaRefine` | `Subdivisions` (1), `Elevate` (0) — each a number or one to three integers, the missing directions at the default — `SampleSubdivisions` (0), `TrimSamples` (16), `Trim` (true), `Edges` (true) | a model step: refines the model exactly and samples it again; the refined model is the one later steps see, so `IgaRefine` then `MeshCadSurface` meshes the refined B-rep. Counters as `IgaSample` |
| `ProjectToCad` | `Model` (a CAD JSON, STEP or IGES path; without it, the B-rep the Input read carries while no step has dropped it), `FaceArray` ("cad:face"), `EdgeArray` ("cad:edge"), `MaxDistance` (0: no limit), `TrimSamples` (16) | moves the mesh's points back onto the model by their cells' `cad:face` and `cad:edge` ([project to CAD](/project_to_cad)); the model stays with the run after it. Counters `NumProjected`, `NumOnEdges`, `NumOnCurves`, `NumOutsideTrim`, `MaxMove` |
| `BoundaryLayer` | `FirstHeight` (required), `NumLayers` (3), `Growth` (1.2), `Reverse` (false), `CollisionFraction` (0.45), `SmoothIterations` (5), `MaxShrink` (6), `Collapse` ("points" \| "layers"), `NormalSmoothIterations` (0), `FrontCheck` (true) | grows prismatic layers from the triangle and quad surface ([boundary layers](/boundary_layer)); the mesh becomes the layers, wall and cap. Counters `NumCells`, `NumShrunk`, `NumCollapsed`, `NumFewerLayers`, `NumFrontCrossings`, `MinScaledJacobian` |
| `Extrude` | `Direction` ([0, 0, 1], the whole displacement), `NumLayers` (1), `Ratio` (1), `Offsets` (none), `RecordLayers` (false), `MaxCells` (20000000) | sweeps the surface into a volume (or a curve into a surface); counters `NumCells`, `NumFlipped`, `MinScaledJacobian`; see [sweeps](/extrude) |
| `Revolve` | `Axis` ([0, 1, 0]), `AxisPoint` ([0, 0, 0]), `Angle` (360, degrees), `NumLayers` (16), `Tolerance` (1e-9), `RecordLayers` (false), `MaxCells` (20000000) | sweeps it around an axis; counters `NumCells`, `NumOnAxis`, `MinScaledJacobian`; see [sweeps](/extrude) |
| `Mirror` | `Point` ([0, 0, 0]), `Normal` ([1, 0, 0]), `KeepOriginal` (false), `Weld` (true), `Tolerance` (1e-9), `ReflectData` (false) | reflects across the plane with every cell right side out; counters `NumFlippedCells`, `NumPointsOnPlane`, `NumDroppedCells`; see [mirror](/mirror) |
| `Tetrahedralize` | `Backend` ("auto" \| "netgen" \| "octree"), `MaxSize` (0: none), `MaxVolume` (0: none), `Grading` (0.3), `OptimizeSteps` (3), `MaxCells` (50000000), `FaceArray` ("cad:face") | fills the closed triangle surface with tetrahedra, its triangles kept through Netgen or its boundary fitted through the native octree backend ([tetrahedralize](/tetrahedralize)); the mesh becomes the tetrahedra and the surface. Counters `NumTetrahedra`, `NumPoints`, `NumNewPoints`, `NumRegions`, `NumBoundaryTriangles`, `BoundaryPreserved` (1 or 0), `HausdorffToSurface`, `MinDihedralDeg`, `MaxDihedralDeg`, `Volume`; the operation's warnings join the run's, prefixed `Tetrahedralize:`. `Backend` "auto" takes Netgen when the build has it (`-DMESHIOPLUSPLUS_WITH_NETGEN=ON`) and "octree" otherwise; "netgen" without it fails naming the flag |

## Sequences (transient / multi-file runs)

Since v9.12.0 the same document can describe a whole **transient** run: a glob/list input, the chain applied per step, and a fan-out or fan-in output. Eight additional keys, all optional:

| Key | Where | Meaning |
| --- | --- | --- |
| `Mode` | top level | `"sequence"` / `"fan-in"` / `"fan-out"`; **asserts** the inferred shape rather than selecting it, and errors naming both on a mismatch |
| `Parallel` | top level | run the steps in a thread pool (**Python driver only**; an error for a fan-in, and the C++ engine warns and runs serially) |
| `Workers` | top level | worker count for `Parallel`; 0 means one per core |
| `Pattern` | `Input` | a glob (`*` and `?` only); mutually exclusive with `Path`/`Paths` |
| `Paths` | `Input` | an explicit, ordered list; not re-sorted |
| `Times` | `Input` | explicit per-step times; the count must match |
| `TimeFrom` | `Input` | `"auto"` (the documented precedence) / `"file"` / `"filename"` / `"index"` |
| `Resample` | top level | resample onto new times: `Times` (an array or `{Start, Stop, Step}`) or `TimesFrom` (another sequence's glob), `Method` (`"linear"` \| `"nearest"` \| `"previous"`), `Extrapolate` (`"error"` \| `"clamp"`), `BlendPoints` (false); at most two source meshes live — see [Resampling onto new times](sequences.md#resampling-onto-new-times) |

`Output.Path` may carry `{step}` or `{index}` to write one file per step.

```json
{
  "Version": 1,
  "Input": { "Pattern": "raw/out_*.vtu" },
  "Operations": [{ "Op": "Quality" }, { "Op": "Clean" }],
  "Output": { "Path": "post/out_{step}.vtu" }
}
```

Two rules worth stating explicitly:

- **A document using none of these keys behaves exactly as before** — the C++ engine literally delegates to the single-file `run_pipeline`, and Python's `run_pipeline` never enters the sequence code path.
- Conversely, the typed single-file parser (`parse_pipeline_json`) **rejects** a sequence key by name rather than ignoring it. Ignoring one would run a transient document as its first step, which is exactly the silent truncation this feature exists to prevent. The CLI verb, Python's `run_pipeline` and the MCP tool all route a sequence document to the right engine automatically, so nobody has to know which kind of document they hold.

See [sequences](sequences.md) for the ordering rule, the time-value precedence, the mode-inference table and the streaming guarantee.

## Where the engine lives

The engine is `operations/pipeline.{hpp,cpp}` in the C++ core, split in two layers:

- The **typed layer** (`PipelineStep`, `apply_pipeline_step`, `run_pipeline_steps`, `run_pipeline`) always compiles and is the **single owner of the step dispatch** — the browser viewer's [`convertSurfaceOps`](wasm.md) pipeline goes through the same code, so the viewer's op chips and a settings.json cannot drift apart (the viewer's op specs are the same words in camelCase; the two casings differ by exactly the first character).
- The **JSON front-end** (`parse_pipeline_*`, `run_pipeline_json/_file`) parses the document with [nlohmann/json](https://github.com/nlohmann/json) v3.12.0, vendored as a git submodule at `src/cpp/third_party/json` exactly like Eigen. When the submodule is absent or `-DMESHIOPLUSPLUS_WITH_JSON=OFF`, the entry points still exist and throw naming the flag (`pipeline_has_json()` / `mio_pipeline_has_json()` / `_core.__has_json__` report which build you have). Wheels, the release CLI binaries and the from-source builds all carry it; the conan/vcpkg packages currently do not (the submodule is not in a source export — the Eigen rule), so there the C ABI entry points fail by name.

**Python's `run_pipeline` is deliberately a pure-Python twin** dispatching over the public Python API rather than the C++ engine: it inherits every operation's C++/numpy parity contract *and* the per-format Python fallbacks (a gmsh `$Periodic` input still runs), and it works on an sdist install with no submodule. `_core.run_pipeline_file`/`_core.run_pipeline_json` expose the C++ engine too, and `tests/python/test_pipeline.py` pins the two engines' outputs — and the transcribed op/key table (`_core.pipeline_op_table()`) — against each other.

## Per-surface entry points

| Surface | Call |
| --- | --- |
| both CLIs | `pipeline SETTINGS.json [--input P] [--output P] [--json] [--quiet]` |
| Python | `run_pipeline(settings, input_path=None, output_path=None)` — settings = dict \| JSON text \| path |
| C | `mio_pipeline_run_file(path)`, `mio_pipeline_run_json(text)`, `mio_pipeline_has_json()` |
| Fortran | `mio_pipeline_run_file(path [, stat, errmsg])`, `mio_pipeline_run_json(...)`, `mio_pipeline_has_json()` |
| Julia | `run_pipeline_file(path)`, `run_pipeline_json(text)`, `pipeline_has_json()` |
| R | `mio_pipeline_run_file(path)`, `mio_pipeline_run_json(text)`, `mio_pipeline_has_json()` |
| WASM | `runPipeline(settings)` — object \| JSON text \| MEMFS `.json` path (no nlohmann in the wasm build; `JSON.parse` does the text forms) |
| MCP | tool `pipeline(settings_path, input_path?, output_path?)` — the sandbox covers the paths *inside* the document too |

### Structured reports on the flat ABI

The existing status-only entry points remain unchanged. C adds `mio_pipeline_run_file_report` / `mio_pipeline_run_json_report` (and `mio_sequence_pipeline_run_*_report` for sequence documents), returning an owned `mio_pipeline_report*` or NULL on failure with `mio_last_error()`. Query `mio_pipeline_report_json(report, NULL, 0)` for the byte length excluding NUL, allocate that length plus one, then call it again to copy JSON; free the handle with `mio_pipeline_report_free`. Accessors never rerun the pipeline, including a size query or truncated copy. JSON-enabled builds serialize the same `{"steps": [{"op": "Clean", ...PascalCase counters}], "warnings": [...]}` shape as Python/WASM, with non-finite counters as `null`. JSON-disabled builds fail by name, not missing symbol.

Fortran adds `call mio_pipeline_run_file_report(path, report, stat, errmsg)` / `mio_pipeline_run_json_report(text, report, ...)`, with an allocatable JSON string, and the analogous sequence subroutines. Julia adds `run_pipeline_file_report(path)` / `run_pipeline_json_report(text)` and `run_sequence_file_report` / `run_sequence_json_report`, returning JSON strings. R adds `mio_pipeline_run_file_report(path)` / `mio_pipeline_run_json_report(text)` and analogous `mio_sequence_pipeline_run_*_report`, returning JSON character scalars (parse with `jsonlite::fromJSON` if desired). These wrappers copy the report and release its native owner automatically, without requiring a JSON parser dependency in the binding package.

## Version 2 spatial multi-mesh steps

Set `"Version": 2` to extend the v1 vocabulary. Per-step `Inputs` is a non-empty array of file paths; paths are resolved relative to the process working directory, like `Input.Path`. Each auxiliary file is read with inferred/sniffed format and default reader options. For transient auxiliary files, this selects the default first step; select a frame beforehand if another is needed. Downstream ordinary steps operate on the one mesh returned by a multi-input step.

| Step | Role of current mesh / `Inputs` | Additional parameters |
| --- | --- | --- |
| `Merge` | current mesh first, followed by all `Inputs` in declared order | `Weld` (false), `Atol` (1e-8), `SourceTag` (true), `DataPolicy` ("intersection" or "fill"), `DropDuplicateCells` (false); counter `NumInputs` includes the current mesh |
| `Interpolate` | current mesh is target; exactly one `Inputs` source | `Method` ("nearest"), `Arrays` (all source point data), `Extrapolate` (false), `DefaultValue` (0), `OnConflict` ("error") |
| `UndoGreen` | current mesh is fine; exactly one `Inputs` coarse mesh | counters `NumGroupsUndone`, `NumCellsRemoved` |
| terminal `Split` | splits the current mesh; no extra inputs | `By` ("type", "component", "region"/"tag", "regions"), `Tag`; requires `Output.Pattern` with `{key}` |
| terminal `Partition` with a pattern | decomposes the current mesh into pieces | v1 partition parameters plus `RecordIds` (false), `GhostLayers` (0); requires `{part}` |

`Split` and partition-to-pieces must be the last step: no implicit downstream branch processing. Without a pattern, v2 `Partition` still attaches labels as in v1; `RecordIds`/`GhostLayers` are valid only for fan-out. `Output.Pattern` and `Output.Path` are mutually exclusive. Patterns replace all occurrences of `{key}` for Split or `{part}` for Partition (zero-based part ids); unknown/mixed tokens are errors. Split keys retain ASCII letters, digits, `.`, `_` and `-`; every other UTF-8 byte becomes `_`, and empty/`.`/`..` keys become `_`. Before writing any piece, both engines reject canonical-path collisions between outputs or with any input, including collisions caused by sanitization or existing symlinks. Existing non-input destinations are overwritten; writes are not an all-or-nothing filesystem transaction. Zero pieces write no files and report `NumPieces: 0`.

```json
{
  "Version": 2,
  "Input": {"Path": "left.vtu"},
  "Operations": [
    {"Op": "Merge", "Inputs": ["right.vtu"], "Weld": true},
    {"Op": "Partition", "Nparts": 2, "Method": "sfc", "RecordIds": true}
  ],
  "Output": {"Pattern": "part_{part}.vtu", "Codec": "none"}
}
```

Python, native C++, both CLIs, C/Fortran/Julia/R status/report entry points and WASM `runPipeline` accept the same document. Native typed callers use `Pipeline.mVersion = 2` and place the pattern in `PipelineOutput.mPath`; installed layouts are unchanged. `pipeline_v2_op_table()` / `_core.pipeline_v2_op_table()` expose the extended vocabulary; the existing per-mesh `apply_pipeline_step` and `pipeline_op_table` remain v1. MCP checks all auxiliary inputs and each expanded output against its configured sandbox, including output symlinks.

Spatial Version 2 is deliberately separate from transient sequence documents (`Mode`, `Input.Paths`/`Pattern`, `{step}` outputs, resampling): combining that vocabulary is rejected, never an implicit cross-product of branches. A multi-step `Input.Path` materializes one state, by default the first; use `Input.Options.TimeStep` to select another, including `-1` for the last. Use a v1 sequence pipeline and standalone v2 pipelines for the spatial stages.
- conan/vcpkg packages shipping the parser via a registry `nlohmann_json/3.12.0` dependency instead of the submodule.

## `Generate`

```json
{ "Version": 1, "Operations": [ { "Op": "Generate", "Shape": "box", "Dims": [8, 8, 8], "CellType": "tetra", "Order": 2 } ], "Output": { "Path": "block.vtu" } }
```

Builds a mesh from parameters ([`doc/primitives.md`](primitives.md)) and **replaces** the incoming geometry; it never reads the incoming mesh. When `Generate` is the first step of a single-mesh document, `Input` may be omitted. A document with any other first step still requires `Input`, and sequence documents (`Mode`, `Input.Pattern`, `Input.Paths`) always do.

`Shape` is one of `rectangle`, `box`, `disk`, `uv_sphere`, `icosa_sphere`, `cylinder`, `reference_triangle`, `ngon`, `annulus`, `lshape`, `plate_with_hole`, `honeycomb`, `reference_cells`, `line`, `circle`, `helix`, `sphere_points`, `box_points`, `tetra_sphere`, `octa_sphere`, `cube_sphere`, `geo_sphere`, `torus`, `moebius`, `tube`, `cone`, `capsule`, `box_surface`, `platonic`, `ball`, `spherical_shell`, `pipe`, `fichera` or `solid_torus`. The other parameters are the shape's arguments in PascalCase, with the defaults of the Python function: `Dims` (rectangle, box, honeycomb), `Size` (two values for rectangle and reference_triangle, three for box, a number for lshape and plate_with_hole), `Origin`, `Center`, `Radius`, `Height`, `N` (disk, reference_triangle, lshape), `Nu`, `Nv` (uv_sphere), `Subdivisions` (icosa_sphere), `Nxy`, `Nz` (cylinder), `Sides`, `Rings` (ngon), `InnerRadius`, `OuterRadius`, `Ntheta`, `Nr`, `Angle` (annulus), `Nt`, `Nr` (plate_with_hole), `Edge` (honeycomb), `Start`, `End` (line), `Pitch`, `Turns` (helix), `Seed` (box_points), `MajorRadius`, `MinorRadius` (torus, solid_torus), `Width`, `NumTwists` (moebius), `Ncap` (capsule), `Ratio` (spherical_shell, pipe), `Name` (platonic), `CellType`, `Split`, `Order` (up to 10), `Complete` (full-Lagrange order 2), `MaxCells`, and `Types`, `Pairs` (reference_cells). The step reports `NumPoints` and `NumCells`.

## `Voxelize`

```json
{ "Op": "Voxelize", "Resolution": [64, 64, 64], "Fill": "surface" }
```

Keys: `Resolution` (two entries select 2-D pixels, as `[nx, ny, 0]` does), `CellSize`, `Bounds`, `Padding`, `PaddingRelative`, `Fill`, `AttachOccupancy`, `MaxCells`, `Sign`, `Plane` (`auto`, `xy`, `xz` or `yz`, a 2-D skin's plane). Reports `NumOccupied`.

It and `ComputeSdf` are the only steps that **replace** their input's geometry rather than transforming it — read a skin, voxelize it, write a grid. See [`doc/voxelize.md`](voxelize.md).

## `Octree`

```json
{ "Op": "Octree", "MaxLevel": 6, "Balance": "face", "Fill": "inside" }
```

Keys: `RootCells`, `MinLevel`, `MaxLevel`, `Bounds`, `Padding`, `PaddingRelative`, `Dimension`, `SampleDistance`, `Balance`, `Fill`, `AttachDistance`, `EvenOdd`, `MaxCells`, `Sign`, `Transition`, `Snap`, `SurfaceLevel`, `CurvatureThreshold`, `FeatureAngle`, `SizeArray`, `Plane`. Reports `NumLeaves`, `NumHanging`, `NumSnapped` and `NumRefused`.

Like `Voxelize` it **replaces** its input's geometry. See [`doc/octree.md`](octree.md).

## `OctreeMesh`

```json
{ "Op": "OctreeMesh", "Output": "castellated", "MinSize": 0.01, "MaxSize": 0.2, "Fill": "inside" }
```

Keys: `Output`, `MinSize`, `MaxSize`, `MaxLevel`, `Bounds`, `Padding`, `PaddingRelative`, `CurvatureAngle`, `FeatureAngle`, `CornerAngle`, `FeatureSize`, `GapCells`, `GapOutside`, `SizeArray`, `SurfaceRegion`, `GradingLayers`, `SkinGrading`, `SurfaceSize`, `WindingThreshold`, `WindingBeta`, `Conformity`, `AlphaLong`, `AlphaShort`, `PreserveFeatures`, `FitIterations`, `DihedralFloorGuard`, `Fill`, `KeepOutside`, `AttachDistance`, `AttachQuality`, `Pillow`, `HexGrading`, `UntangleIterations`, `MaxCells`, `MaxTets`. Reports `NumLeaves`, `NumCells`, `NumInside`, `NumSurrogateFaces`, `NumHanging`, `NumInverted`, `MinDihedral` and `MinScaledJacobian`.

It **replaces** its input's geometry with a volume mesh of the skin's inside. See [`doc/octree_mesh.md`](octree_mesh.md).

## `ComputeSdf`

```json
{ "Op": "ComputeSdf", "Structure": "octree", "RootResolution": 8, "MaxDepth": 4 }
```

Keys: `Structure` ("voxel" | "octree"), `Resolution`, `CellSize`, `Bounds`, `Padding`, `PaddingRelative`, `RootResolution`, `MaxDepth`, `BandCells`, `RecordLevels`, `MaxCells`, `Sign`, `Location`, `Band`. Reports `MaxDepth` and `NumBanded`.

`Resolution`/`CellSize` size a **voxel** grid and are an error with `"Structure": "octree"`, whose finest cell is `RootResolution / 2^MaxDepth` and is therefore already determined. See [`doc/sdf.md`](sdf.md).
