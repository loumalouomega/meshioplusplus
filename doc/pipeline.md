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
| `DecimateVolume` | `Ratio` \| `TargetCells` \| `MaxError` (none → `Ratio` 0.5), `Placement` ("optimal"), `PreserveBoundary` (**false**, unlike `Decimate`), `PreserveFeatures` (true), `FeatureAngle` (30) | counters `TetsRemoved`, `CollapsesRejected` |
| `Partition` | `Nparts` (2), `Method` ("auto"), `Imbalance` (0.03), `Mode` ("eco"), `Seed` (0), `WeightsKey` | attaches `partition:part`; counter `Nparts` |
| `Slice` (alias `Section`) | `Point`, `Normal` (required), `RecordParentIds` (false) | counter `SectionFaces`; warns when the plane misses |
| `Gradient` | `Array` (required), `Operator` ("gradient"), `Method` ("green-gauss"), `Location` ("cell"), `Output`, `Component` | counters `NumSkipped`, `NumFallback`; warns on skipped cells |
| `Normals` | `PointNormals` (true), `CellNormals` (false), `Weight` ("angle" \| "area"), `SplitAngle` (absent → one smooth normal per point; a number of degrees in `[0, 180]` splits vertices at creases), `RecordParentIds` (false), `Region` | attaches `normals`; counters `NumIsolated`, `NumUndefined`, `NumDegenerate`, `NumSplitPoints`, `NumAddedPoints`; warns when the input winds inconsistently and no split was asked for; see [normals](/normals) |
| `FeatureEdges` | `FeatureAngle` (30), `Feature`, `Boundary`, `NonManifold`, `Inconsistent` (all true), `Region` | replaces the mesh by its feature edges, `line` cells with `feature:kind`/`feature:angle`; counters `NumFeature`, `NumBoundary`, `NumNonManifold`, `NumInconsistent`; see [feature edges](/feature_edges) |
| `EditRegions` | `Edit` (`union` \| `intersection` \| `difference` \| `rename` \| `retag` \| `delete`), `Inputs` (region names), `Output`, `Kind` (pins the inputs' kind), `Dim`, `Tag`, `KeepInputs` (true) | one region edit, points/cells/data untouched — a list of edits is a list of steps; counter `NumRegions`; see [editing regions](/regions#editing-regions) |
| `RegionAdjacency` | `Regions` (optional ordered list of Cell-region names) | replaces the mesh with conforming shared facets; when fewer than two Cell regions exist, multiple cell blocks are treated as groups; counters `NumFacets` and `Measure`; see [region adjacency](/region_adjacency) |
| `FindInterface` | `RegionA`, `RegionB` (Cell-region names); `Mode` ("conforming" or "proximity"), `Master` ("a" or "b"), `GapTolerance` (0), `AngleTolerance` (30), `OverlapTolerance` (0) | replaces the mesh with matched master facets; counters `NumPairs`, `Area`, `MaxGap`, `UnmatchedA`, `UnmatchedB`; this pipeline form selects both parts from the current mesh |
| `SplitInterface` | `Region` (Side-region name), `AddCohesive` (false) | duplicates node fans along the Side region and optionally inserts cohesive cells; counters `NumDuplicatedPoints`, `NumCohesiveCells`; see [interfaces and contact](/region_adjacency) |
| `QualityGate` | `Require` (threshold texts, e.g. `["scaled_jacobian >= 0.2"]`), `MaxInverted` (0), `MaxDegenerate` (0; negative disables either) | a gate, not a transform: the mesh passes through untouched and a failed check stops the pipeline with the summary as its error; counters `NumChecks`, `NumFailed`; see [quality gate](/mesh_quality#quality-gate) |
| `Snapshot` | `Path` (required: `.png`, `.txt`, `.ansi`, `.html` or `.cast`), `Width` (320), `Height` (240), `Supersample` (1), `Cols` (100), `Rows` (40), `Encoding` ("halfblock" \| "quadrant" \| "sextant" \| "braille" \| "ascii"), `ColorDepth` ("truecolor" \| "256" \| "16" \| "mono"), `PngCompress` (0), the camera `View`, `Azimuth`, `Elevation`, `Roll`, `Perspective`, `Fov`, `Zoom`, `Shading`, `Edges`, `FeatureAngle`, the colours `Background` and `FillColor` (`#rrggbb[aa]` or `none`), colouring by `ColorBy`, `Component`, `Cmap`, `VMin`, `VMax`, `Colorbar`, `Axes`, `ScaleBar`, and the field options `Reduce`, `Expr`, `ClipLow`, `ClipHigh`, `Symmetric`, `Scale`, `ScaleThreshold`, `Categorical`, `ColorRegions`, `CategoryEdges`, `Isolines`, `IsoLevels`, `Vectors`, `VectorCount`, `Warp`, `WarpScale`, `WarpOutline`, `Diagnostic`, `QualityMetric`, and the cut-away `Cutaway` (a flat list of numbers, six per plane, at most two: a point, then the normal of the side kept) with `CutawayTint`, and the streamlines `Streamlines` (a vector point array), `StreamSeeds` (40), `StreamLength` (0.5) and `StreamColor` | a side output, not a transform: the mesh passes through untouched and the [software rasterizer](/tui) writes the frame as it stands at this point, so a `Smooth` or `Decimate` can be shown before and after in a report; counters `Width`, `Height`. The parameters are the options of `snapshot`, with the same names in PascalCase |
| `EstimateError` | `Array` (required), `Method` ("zz"), `Marking` ("none" \| "absolute" \| "fraction" \| "dorfler"), `MarkingValue` (0), `Output`, `Marked` | attaches `error:zz` (and `error:marked` when `Marking` isn't "none"); counters `GlobalError`, `NumSkipped`, `NumMarked`; warns on skipped cells |
| `Isosurface` | `Array` (required), `Isovalue` (0) or `Isovalues`, `Component`, `RecordParentIds` (false) | counter `ContourCells`; warns when empty |
| `Transform` | exactly one of `Translate[3]`, `Scale` (number or `[3]`), `RotateAxis[3]`+`RotateDegrees`, `Matrix[16]` (row-major), `ScaleUnits` (a factor, e.g. `0.001` for mm→m); plus `RotateData` (false) | |
| `ConvertCells` | `Mode` ("linearize" \| "simplexify" \| "elevate"), `RecordParentIds` (false) | |
| `Subdivide` | `RecordParentIds` (false) | one polyhedral child per 3D cell face, connected to a new interior point; no per-type template table |
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

## `Voxelize`

```json
{ "Op": "Voxelize", "Resolution": [64, 64, 64], "Fill": "surface" }
```

Keys: `Resolution`, `CellSize`, `Bounds`, `Padding`, `PaddingRelative`, `Fill`, `AttachOccupancy`, `MaxCells`, `Sign`. Reports `NumOccupied`.

It and `ComputeSdf` are the only steps that **replace** their input's geometry rather than transforming it — read a skin, voxelize it, write a grid. See [`doc/voxelize.md`](voxelize.md).

## `ComputeSdf`

```json
{ "Op": "ComputeSdf", "Structure": "octree", "RootResolution": 8, "MaxDepth": 4 }
```

Keys: `Structure` ("voxel" | "octree"), `Resolution`, `CellSize`, `Bounds`, `Padding`, `PaddingRelative`, `RootResolution`, `MaxDepth`, `BandCells`, `RecordLevels`, `MaxCells`, `Sign`, `Location`, `Band`. Reports `MaxDepth` and `NumBanded`.

`Resolution`/`CellSize` size a **voxel** grid and are an error with `"Structure": "octree"`, whose finest cell is `RootResolution / 2^MaxDepth` and is therefore already determined. See [`doc/sdf.md`](sdf.md).
