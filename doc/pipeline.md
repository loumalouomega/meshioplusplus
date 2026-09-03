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
- **A chain runs over one mesh at a time.** `Merge`, `Interpolate`, `Split`, `Diff` and `UndoGreen` need extra inputs or produce extra outputs, and a step naming one errors pointing at the matching CLI verb (`undo-green` for `UndoGreen`, which needs a second coarse mesh exactly as `Interpolate` needs a second source mesh). `Partition` as a step attaches the `partition:part` labels (colour-by-part) rather than splitting into pieces.
- Steps are validated **before** the input is read — a typo in step 7 never costs reading a 10 GB mesh first.
- The run returns a **report**: `{"steps": [{"op", ...counters}], "warnings": [...]}` with PascalCase counter keys (`PointsWelded`, `SectionFaces`, `NumSkipped`, ...).

## Top-level schema

| Key | Required | Meaning |
| --- | --- | --- |
| `Version` | no (default 1) | schema version; this build knows `1` |
| `Input` | yes | `{Path, Format?, Options?}` — `Format` defaults from the extension, with the `sniff_format` read fallback |
| `Operations` | no (default `[]`) | the step array; empty = a plain convert |
| `Output` | yes | `{Path, Format?, Encoding?, Codec?, FloatFormat?}` |

`Input.Options` narrows the read (see [selective reads](selective_read.md)): `PointsOnly` (bool), `DataArrays` (string array; absent = all, `[]` = none), `TimeStep` (int), `Lenient` (bool), `Mmap` (`"auto" | "on" | "off"`).

`Output` maps onto `registry_write_ex`: `Encoding` (`"ascii" | "binary"`), `Codec` (`"none" | "zlib" | "lz4" | "zstd"`, VTU/VTP block codecs), and `FloatFormat` (a printf-style float format for ASCII writers that take one). An option the output format cannot honour is an error.

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
| `EstimateError` | `Array` (required), `Method` ("zz"), `Marking` ("none" \| "absolute" \| "fraction" \| "dorfler"), `MarkingValue` (0), `Output`, `Marked` | attaches `error:zz` (and `error:marked` when `Marking` isn't "none"); counters `GlobalError`, `NumSkipped`, `NumMarked`; warns on skipped cells |
| `Isosurface` | `Array` (required), `Isovalue` (0) or `Isovalues`, `Component`, `RecordParentIds` (false) | counter `ContourCells`; warns when empty |
| `Transform` | exactly one of `Translate[3]`, `Scale` (number or `[3]`), `RotateAxis[3]`+`RotateDegrees`, `Matrix[16]` (row-major), `ScaleUnits` (a factor, e.g. `0.001` for mm→m); plus `RotateData` (false) | |
| `ConvertCells` | `Mode` ("linearize" \| "simplexify" \| "elevate"), `RecordParentIds` (false) | |
| `Subdivide` | `RecordParentIds` (false) | one polyhedral child per 3D cell face, connected to a new interior point; no per-type template table |
| `Agglomerate` | `TargetGroupSize` (8) | greedy seed-and-grow over the shared-face dual, merging face-adjacent cells into one polyhedron per group; the many-to-one counterpart to `Subdivide` |
| `Crop` | one of `Bbox[6]` (`[xmin,ymin,zmin,xmax,ymax,zmax]`), `Point[3]`+`Normal[3]`, or `Where` (a scalar `cell_data` name) + `Compare` ("<") + `Value` (0); `Mode` ("all" \| "any", **not** with `Where`), `RecordIds` (false) | counter `CellsKept` |
| `ExtractSurface` | `RecordParentIds` (false) | |
| `ExtractSkin` | `Linearize` (false) | |
| `Reorder` | `Method` ("rcm" \| "morton" \| "hilbert") | |
| `DataDrop` | `Point`/`Cell`/`Field` (string arrays), `IgnoreMissing` (false) | |
| `DataKeep` | `Point`/`Cell`/`Field` — only locations mentioned are touched; `[]` drops all there | |
| `DataRename` | `Point`/`Cell`/`Field`, entries `"OLD:NEW"` (split on the **last** colon — names carry colons: `gmsh:physical`) | |
| `DataCalc` | `Expr` (`"NAME = EXPRESSION"`, split on the **first** `=`), `Location` ("point"), `Overwrite` (false) | |
| `DataCondition` | `Mode` ("clamp" \| "normalize" \| "standardize"), `Location` ("point"), `Names`, `Scope` ("component" \| "magnitude"), `Lo` (0), `Hi` (1), `NanPolicy` ("ignore"), `NanReplacement` (0), `Suffix` | |
| `ToCell` / `ToPoint` | `Names`; `ToPoint` also `Weight` ("uniform" \| "measure") | |

## Sequences (transient / multi-file runs)

Since v9.12.0 the same document can describe a whole **transient** run: a glob/list input, the chain applied per step, and a fan-out or fan-in output. Seven additional keys, all optional:

| Key | Where | Meaning |
| --- | --- | --- |
| `Mode` | top level | `"sequence"` / `"fan-in"` / `"fan-out"`; **asserts** the inferred shape rather than selecting it, and errors naming both on a mismatch |
| `Parallel` | top level | run the steps in a process pool (**Python driver only**; an error for a fan-in, and the C++ engine warns and runs serially) |
| `Workers` | top level | worker count for `Parallel`; 0 means one per core |
| `Pattern` | `Input` | a glob (`*` and `?` only); mutually exclusive with `Path`/`Paths` |
| `Paths` | `Input` | an explicit, ordered list; not re-sorted |
| `Times` | `Input` | explicit per-step times; the count must match |
| `TimeFrom` | `Input` | `"auto"` (the documented precedence) / `"file"` / `"filename"` / `"index"` |

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

The flat ABI (C/Fortran/Julia/R) carries **JSON text only** and reports status + `mio_last_error()`; the structured report is a recorded follow-up.

## Follow-ups (recorded, not implemented)

- **Multi-mesh steps**: `Merge`/`Interpolate`/`UndoGreen` would need per-step `Inputs: [paths]`, and `Split`/partition-to-pieces an `Output.Pattern` with `{key}`/`{part}` — the v2 schema sketch; today the CLI verbs cover these. (v9.12.0's [sequences](sequences.md) added the *input*-list and `{step}`-output halves of this for the transient case, but a step that consumes or produces several meshes at once is still out of scope.)
- A **C ABI report accessor** (caller-buffer JSON string of the run report).
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
