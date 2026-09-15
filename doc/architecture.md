# Architecture

meshio++ is one C++ core with six language surfaces on top of it and a set of tools built on those surfaces. This page is the map; every box below links to the page that owns it.

![meshio++ architecture: one C++ core holding the format registry, the operations layer and one of three mesh backends; the Python, C, Fortran, Julia, R, WebAssembly and C++ surfaces over it; and the CLIs, viewers, MCP server, Blender add-on, ParaView plugin and PhysicsNeMo adapter built on those](/diagrams/architecture.svg)

## One C++ core

Everything compiles into a single object library, `meshioplusplus_core_obj` (C++20), and every optional dependency compiles out when it is absent, so the core is dependency-free by default. Four things live in it.

- **The format registry** — [43 readable and 46 writable formats](./formats.md), each a reader and a writer over the uniform mesh API, registered in one dispatch table (`registry.cpp`: format name to reader and writer, file extension to default format). The Python binding keeps one function per format and dispatches in Python; every flat binding (C, Fortran, Julia, R, WebAssembly, the native CLI) dispatches through this table, which is why a new C++ format must be registered there to be reachable from them.
- **The operations layer** — 34 mesh operations and 5 data operations, described below.
- **The `detail/` kernels** — the shared machinery the formats and operations are built from: the face and edge tables (`cell_faces`, `cell_edges`), the single owners of cell indexing and region remapping (`cell_index`, `region_remap`), the marching-tetrahedra cutter, the spatial hash, the polyhedron kernel, the refinement templates and the provenance scope.
- **The mesh backend** — exactly one per build, selected by `MESHIOPLUSPLUS_MESH_BACKEND`: the meshio-mirroring `MESHIO` layout the Python wheel uses, the canonical `NATIVE` layout the WebAssembly build uses, or the Kratos-style `KRATOS` ModelPart. Formats, operations and bindings never touch a backend's own members; they go through the [uniform mesh API](./cpp_backends.md), which is what lets one implementation compile under all three.

The in-memory model those pieces share is the [Mesh data model](./mesh_data_model.md): points, a list of homogeneous cell blocks, data aligned to points or to blocks, and [named regions](./regions.md). Ragged cells (polygons and polyhedra) have their own [representation](./polyhedra.md) and cross every flat binding as CSR arrays.

## The operations layer

An operation is a computation *on* a mesh rather than a file format: it is written against the uniform mesh API only, parallelises its hot loop with `parallel_for`, and is exposed on every binding surface plus a CLI verb and an MCP tool. They group naturally.

| Group | Operations |
| --- | --- |
| Inspection and topology | [quality](./mesh_quality.md), [stats](./stats.md), [diff](./diff.md), [extract_surface](./extract_surface.md), [extract_skin](./extract_skin.md), [reorder](./reorder.md), [transform](./transform.md), [clean](./clean.md), [merge](./merge.md), [crop](./crop.md), [split](./split.md), [partition](./partition.md), [convert_cells](./convert_cells.md) |
| Refinement and coarsening | [refine](./refine.md), [undo_green](./undo_green.md), [subdivide](./subdivide.md), [agglomerate](./agglomerate.md), [decimate](./decimate.md), [decimate_volume](./decimate_volume.md) |
| Remeshing and smoothing | [remesh](./remesh.md), [remesh_volume](./remesh_volume.md), [optimize_volume](./optimize_volume.md), [smooth](./smooth.md) |
| Fields | [interpolate](./interpolate.md), [conservative_interpolate](./conservative_interpolate.md), [gradient](./gradient.md), [hessian](./hessian.md), [estimate_error](./error.md), [data_integrate](./field_integration.md) |
| Cutting, grids and distance | [slice](./slice.md), [isosurface](./isosurface.md), [grid and voxelize](./voxelize.md), [signed distance](./sdf.md) |
| Data operations | [manage, average, calc, condition, info](./data_operations.md) — the five that never touch geometry |

Two things tie the operations together. Chains of them are described declaratively by the [settings pipeline](./pipeline.md), whose typed layer (`run_pipeline_steps`) is the single owner of the step dispatch on every surface, and the [sequence driver](./sequences.md) runs such a chain once per step of a transient dataset with at most one mesh alive at a time. Every write can carry a [provenance record](./provenance.md) of where the mesh came from and what was done to it.

## The language surfaces

- **Python** — the `meshioplusplus._core` pybind11 extension with zero-copy numpy at the I/O boundary, wrapped by one shim per format that falls back to the pure-Python reference implementation on any exception. This is the surface the Python CLI, the MCP server and the integrations are built on.
- **C** — [`libmeshioplusplus`](./c_api.md), a pure C99 header with `SOVERSION 0` and append-only option structs, compiled under every mesh backend. [Fortran](./fortran.md), [Julia](./julia.md) and [R](./r.md) ride on it, each in its own idiom for handles, ownership and indexing.
- **WebAssembly** — the [`@meshioplusplus/wasm`](./wasm.md) npm package (embind over the NATIVE backend), in a sequential and a threaded build, working on a MEMFS virtual filesystem.
- **C++** — the [installable C++ API](./cpp_api.md), one real library per backend exported into the same CMake package as the C API, under a deliberate [ABI contract](./abi.md); or the [single-header amalgamation](./single_header.md) for a build with no CMake at all.

## The tools built on them

The Python [CLI](./cli.md) and the Python-free native CLI mirror each other verb for verb. The [MCP server](./mcp.md) exposes the whole Python surface to AI agents as stateless, file-path-based tools. The [browser viewer and dataset manager](./viewer.md) consume the published WebAssembly package, while the Polyscope viewer is a Python extra. The [Blender add-on](./blender.md) and [ParaView plugin](./paraview_plugin.md) bring meshio++'s formats into those applications, and the [interoperability](./interop.md), [GPU](./gpu.md), [machine-learning](./ml.md), [dataset](./datasets.md) and [PhysicsNeMo](./physicsnemo.md) layers hand meshes to the wider Python ecosystem without a file round-trip.

## Two patterns that recur

**A pure layer under gated wrappers.** The interoperability, GPU, MCP, Blender and viewer modules all split the same way: the bulk of the logic is a pure payload layer that imports no optional library, mutates nothing and is tested in the default CI matrix with none of the targets installed, and the public functions are thin wrappers that import their target inside the function and raise a named install hint when it is missing. A feature that needs a heavy or version-pinned dependency never makes the core, or its tests, depend on it.

**A single owner for anything two places could compute.** The block-major cell index (`detail/cell_index.hpp`), the remapping of regions through an operation (`detail/region_remap.hpp`), the refinement templates, the marching cutter, the provenance credit line and the pipeline step dispatch each live in exactly one place and are reused by every consumer, including the numpy twins on the Python side, which are pinned byte-for-byte against the C++ path by tests. Where a second implementation would have to reproduce a discrete branch on a sign (a winding repair, a flip acceptance, a warp), the operation is C++-core only and says so by name.

## Where to go next

Start with the [quickstart](./quickstart.md), then the [format table](./formats.md) and the [CLI reference](./cli.md). The [roadmap](./roadmap.md) lists what is not built yet.
