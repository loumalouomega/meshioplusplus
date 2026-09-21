# glTF (`.glb`, `.gltf`)

[glTF 2.0](https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html) output — a **write-only**, web-native skin export: the browser viewer, three.js, Blender, model viewers and dashboards load it directly, with no meshio++ in the loop. The **surface** of the mesh is exported; `.glb` writes the single-file binary container and `.gltf` writes JSON with a `.bin` sidecar beside it (the debuggable alternative).

| | |
|---|---|
| **Format name** | `gltf` |
| **Extensions** | `.glb`, `.gltf` |
| **Read / Write** | — / ✓ |
| **Extra dependencies** | — |

## Reading & writing

There is no reader — `register_format` is called with `read=None`, and reading glTF is out of scope. Full write signature:

```python
import meshioplusplus

meshioplusplus.gltf.write(
    "out.glb", mesh,
    split_angle=30.0,
    normal_weight="angle",
    normals=True,
    fields=True,
    color_by=None,
    component=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    nan_color="#808080",
    unlit=True,
    up_axis="auto",
    recenter=True,
    scale=1.0,
    by_region=True,
    container="auto",
)
```

- **`container`** — `"auto"` (default) follows the suffix: `.glb` is binary, anything else is JSON plus `<name>.bin`. `"glb"` and `"gltf"` override it. A file-like object always receives the binary container.
- **`split_angle`** / **`normal_weight`** / **`normals`** — see [Normals](#normals).
- **`fields`** — export the `point_data` arrays as `_NAME` attributes, see [Fields](#fields).
- **`color_by`** … **`unlit`** — bake a field into `COLOR_0`, see [Colour](#colour).
- **`up_axis`** / **`recenter`** / **`scale`** — see [Conventions](#conventions).
- **`by_region`** — one node per cell region (default) or a single node named `mesh`.

`meshioplusplus convert in.vtu out.glb` writes the defaults. The colour flags of `convert` (`--color-by`, `--component`, `--cmap`, `--vmin`, `--vmax`, `--nan-color`) apply to glTF output too, plus `--split-angle` and `--up-axis`; `--colorbar` is svg/tikz only. The MCP server exposes every option as `export_gltf`. The flat bindings (C, Fortran, Julia, R, WASM) write the defaults through `mio_write("x.glb")`/`writeMesh`; there is no options struct on the flat ABI.

## What is exported

Everything is reduced to a surface, in this order:

* **The skin of volume blocks** — the boundary facets of `tetra`, `hexahedron`, `wedge`, `pyramid` and their higher-order forms, linearised to triangles and quads (see [Skin extraction](../extract_skin.md)). Interior faces are not written.
* **2-D cells** — `triangle`, `quad` and `polygon` blocks (ragged or rectangular), triangulated with the same fan [`convert_cells(simplexify)`](../convert_cells.md) uses; higher-order cells contribute their corner nodes. **A 2-D cell that coincides with a skin facet wins**: the facet is dropped and the cell is exported instead, so a boundary patch such as a `wall` region keeps its own cells and does not z-fight with the skin.
* **`line` cells** as `LINES` primitives (a two-node segment per cell).
* **`vertex` cells** as `POINTS`; and a mesh with **no cells at all** — which is how a `.pcd` or `.xyz` [point cloud](./pcd.md) arrives — as one `POINTS` primitive over every point, so those formats get a web output for free.

Triangles with a repeated point id are dropped. Other cell types (for example the `VTK_LAGRANGE_*` families) have no glTF equivalent: they are skipped with a warning and a provenance note.

## Normals

glTF normals are per vertex, so a crease needs its position twice. Triangles are grouped into smooth fans exactly as [`compute_normals`](../normals.md) does — the same kernel, run over the whole soup so that smoothing crosses region borders — and each fan becomes its own vertex. `split_angle` is the largest dihedral angle, in degrees and in `[0, 180]`, still treated as smooth; boundary edges, non-manifold edges and edges whose two triangles are wound the same way always cut, and the triangles fanned from one polygon or one volume facet are never separated from each other. `normal_weight` is `"angle"` (default) or `"area"`. `normals=False` writes no `NORMAL` and shares vertices.

A group whose triangles have no area gets the normal `(0, 0, 1)`. A `POINTS` primitive gets a `NORMAL` from the `normals` point data when every exported row of it is finite and non-zero (it is normalised first), which is how a cloud written by `compute_normals` keeps its normals; otherwise `NORMAL` is omitted with a warning.

## Fields

glTF has no scalar-field concept. Every `point_data` array of **one to four components** is exported **raw** as an underscore-prefixed custom attribute — `temperature` becomes `_TEMPERATURE`, a 3-component `velocity` becomes a `VEC3` `_VELOCITY` — cast to `float32`, with the original name in the accessor's `name`. Names are upper-cased (ASCII), every other character becomes `_`, and two names that collapse onto one get a numeric suffix (`_T`, `_T_2`). There is no 32-bit integer vertex attribute type, so integer arrays are cast to float. An array with more than four components, and an array with a non-finite value on an exported point (glTF cannot hold one), is skipped with a warning and a provenance note. `normals` is never exported as a field.

## Colour

`color_by` names a `point_data` or `cell_data` array (a point array wins if both exist) to bake through a colormap into `COLOR_0`, with the raw values alongside as `_NAME`:

* **Point data** colours each vertex by its point's value; the renderer interpolates between the vertex colours. **Cell data** colours each triangle by its cell's value — through the parent volume cell for a skin facet — and vertices are split per cell so the colour stays flat.
* Multi-component arrays reduce to `component` or to their magnitude. The range is `vmin`..`vmax`, defaulting to the finite range of the **exported** vertices. Non-finite values take `nan_color`, given as `#rrggbb`. `cmap` is `viridis`, `coolwarm` or `turbo`, the tables the [SVG](./svg.md) and [TikZ](./tikz.md) writers use.
* `COLOR_0` is **linear** in glTF, while the colormap tables hold 8-bit sRGB display values, so each byte goes through the inverse sRGB transfer function first — otherwise every viewer would brighten the map. The conversion is a 256-entry table of `float32` bit patterns (generated by `tools/gen_srgb_linear.py`), not a `pow`, so both engines write the same bytes whatever libm they link.
* The material becomes **`KHR_materials_unlit`** with a white base colour, because shading distorts a colour reading; `unlit=False` keeps the standard PBR material. Without `color_by` no vertex colours are written and the material is a plain grey.

## Structure

One glTF node, and one mesh, per **cell region**, named by the region. A cell that is in several regions goes to the **smallest** one (ties go to the first in the canonical `(kind, name, dim, tag)` order), regions with the same name share a node, and a cell in no region goes to a node named `unassigned`; a note in the provenance record counts the overlaps. Only regions of kind `cell` matter — point and side regions are not exported. A mesh with no cell regions, or `by_region=False`, is one node named `mesh`. Nodes that would hold nothing are not written. Under one root node named `meshio++`, the scene has a single material shared by every primitive, which is `doubleSided` because the operation does not reorient the surface.

## Conventions

glTF is `float32`, right-handed, Y-up and metric. CAE data is usually **Z-up** and in its own unit, and often far from the origin, so the changes are carried by the **root node's transform** rather than written into the coordinates:

* **`up_axis`** names the source axis that points up: `"z"` is a rotation of −90° about X (a source `(x, y, z)` lands at `(x, z, −y)`), `"x"` a rotation of +90° about Z (`(−y, x, z)`), `"y"` none. `"auto"` (default) picks `y` for a flat mesh — 2-D points, or every exported `z` within `1e-14` of zero, the SVG writer's rule — and `z` otherwise.
* **`recenter`** subtracts the exported points' bounding-box centre from every position, so a model with coordinates in the millions keeps `float32` precision; the root node's `translation` is that centre, rotated and scaled, so the world position of every vertex is exactly the source position.
* **`scale`** is the factor from the source unit to metres (`0.001` for millimetres), carried as the root node's `scale`.

The rotation is a permutation with signs, so the transform is exactly reversible.

## Container

The writer follows the specification's rules for the binary container. A 12-byte header (magic `0x46546C67`, version 2, total length), then chunks of `length, type, data`, each 4-byte aligned: the JSON chunk padded with **spaces**, and the BIN chunk (omitted when there is no geometry) padded with **zeros**. Every accessor has its own `bufferView` with its `target`; indices are `uint32` (a primitive with more than `2^32 − 2` vertices, or a container beyond 4 GiB, is a `WriteError`); `POSITION` carries `min` and `max`, written as the exact `float32` values, which is what the validator checks them against. Numbers are written with `%.17g`, so the JSON is byte-exact between the engines. `asset.generator` is the [provenance](../provenance.md) tag; with a scope open, the full record goes in `asset.extras["meshioplusplus:provenance"]`.

An empty mesh writes an empty scene (`"scenes":[{}]`) with no buffer. A non-finite coordinate on an exported point is a `WriteError`.

## C++ vs Python

The C++ core (`write_gltf`, in `formats/gltf.hpp`) and the Python reference (`gltf/_gltf.py`) follow the same steps in the same order with the same arithmetic — the smooth-fan grouping is the bit-exact pair `detail::vertex_normal_groups` / `_normals.vertex_normal_groups` — and write **identical bytes** for `.glb` and for `.gltf` and its sidecar; `tests/python/test_gltf.py` pins that across cubes, volumes, a sphere, an open cylinder, regions, ragged polygons, lines and point clouds under six option sets. Set `MESHIOPLUSPLUS_STRICT_CORE=1` to prove the native path handles a file.

## Validation

The output passes the Khronos [glTF-Validator](https://github.com/KhronosGroup/glTF-Validator) with **zero errors, zero warnings and zero infos** on those fixtures, in both containers: `tests/python/test_gltf.py` runs it through `tests/python/gltf_validate.mjs` (the `gltf-validator` npm package) when `node` and the package are present, and the `gltf-validator` CI job requires it. It also loads in [trimesh](https://trimesh.org/) with the correct Y-up bounds. Loading in three.js and Blender, with the orientation and colours checked by eye, was not automated.

## Not supported

Reading glTF; Draco or meshopt compression; textures and UVs; transient results as morph targets (a plausible follow-up); `EXT_structural_metadata`, the principled carrier of typed fields, if a consumer asks; point and side regions.
