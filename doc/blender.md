# Blender add-on

Blender ships Python and reads almost no FEA formats: STL, OBJ, PLY, and essentially nothing else. meshio++ reads 43. The add-on is the bridge — install it and every one of them appears under `File > Import`.

```
File > Import > Mesh via meshio++     # .msh, .vtu, .med, .cgns, .inp, .e, .unv, ...
File > Export > Mesh via meshio++
```

It is a **Blender 4.2+ extension**: the meshio++ wheel travels inside the zip, so there is no pip step, no network access at install time, and nothing to configure.

## Installing

Download the zip for your platform from the [latest release](https://github.com/loumalouomega/meshioplusplus/releases) — `meshioplusplus-<version>-linux-x64.zip`, `-windows-x64`, `-macos-x64` or `-macos-arm64` — then either drag it into a Blender window, or use `Edit > Preferences > Add-ons > Install from Disk`, or run `blender --command extension install-file <zip>`.

To build one yourself, see [`src/blender/README.md`](https://github.com/loumalouomega/meshioplusplus/blob/main/src/blender/README.md).

## Using the bridge from a script

The add-on is a shell over two public functions, and they are usable directly from Blender's scripting console or any `bpy` session:

```python
import meshioplusplus as mio

mesh = mio.read("bracket.msh")
data = mio.to_blender(mesh, name="bracket")     # a bpy.types.Mesh datablock

import bpy
bpy.context.collection.objects.link(bpy.data.objects.new("bracket", data))

back = mio.from_blender(data)                   # and out again
mio.write("bracket.vtu", back)
```

`to_blender` returns the **datablock**, not an `Object`, and links nothing into a collection — the same separation as `to_pyvista` returning a grid rather than a plotter.

## There is deliberately no `[blender]` extra

`has_blender()` is the predicate, and calling either function without `bpy` raises:

```
ImportError: meshio++: to_blender: bpy is not installed. There is deliberately no
pip extra for it; install it directly with `run this inside Blender, or `pip install bpy`
for a headless build` (the PyPI wheel is pinned to one CPython version). See doc/blender.md.
```

There is no `pip install meshioplusplus[blender]` because that phrasing would be wrong twice over. The `bpy` wheel on PyPI is several hundred megabytes and built against exactly one CPython minor, so an extra pinning it would break for most people; and inside Blender — the only place this code normally runs — `bpy` is a builtin that pip must never touch, since installing it there would install a second Blender. This is the same reasoning [`doc/ml.md`](ml.md) records for `torch`, `jax` and CuPy, and it uses the same error shape.

## What arrives in Blender

| meshio++ | Blender |
| --- | --- |
| `points` | vertices (padded to 3 columns, stored float32) |
| `triangle`, `quad`, `polygon` | polygons — **n-gons are kept, not triangulated** |
| `line` | edges |
| `vertex` | nothing of its own; every point is imported regardless |
| `point_data` | attributes on the `POINT` domain |
| `cell_data` | attributes on the `FACE` domain |
| `Point` / `Cell` regions | `BOOLEAN` attributes named `region:<name>` |
| `field_data` | custom properties on the mesh datablock, `meshioplusplus:field:<name>` |

Scalars become `FLOAT`, `INT` or `BOOLEAN`; a 3-component float array becomes `FLOAT_VECTOR` and a 2-component one `FLOAT2`. Anything wider, and any non-float vector, expands into one scalar attribute per component named `v_0`, `v_1`, … — the same suffix rule [`to_pandas`](interop.md) uses, which is the one component-flattening convention in the library. A name Blender owns (`position`, `material_index`, `id`, …) or that starts with a dot is renamed with a note rather than colliding.

Region `dim` and `tag` cannot live on a boolean mask, so they ride a JSON sidecar in the mesh's custom properties under `meshioplusplus:regions` — the same key and the same convention [`to_pyvista`](interop.md) uses, which is what lets a region survive a trip out through one bridge and back through the other. `Side` regions are dropped by name: Blender has no `(cell, local facet)` concept.

**n-gons are kept, and that is the point.** Blender's mesh is vertices + loops + polygons, so it holds triangles, quads and n-gons natively. `to_trimesh` triangulates because a `Trimesh` holds nothing else; doing that here would destroy exactly the quad topology someone importing a hex-dominant mesh came for.

## Volume meshes arrive as their boundary

Blender has no tetrahedron, hexahedron, wedge or pyramid — there is no representation to map them onto — so a volume mesh becomes its boundary surface. Each `cell_data` array is carried through to the faces of the cell that owned it, so a solid can still be coloured by its per-cell material or tag; and higher-order cells are linearized, because a triangle renderer has no concept of a mid-side node.

This is the *same* reduction the [browser viewer](viewer.md) uses, deliberately: a mesh should look the same whether it reached you through Blender, through `view()`, or through the hosted demo.

## Exporting from Blender

`from_blender` accepts a `bpy.types.Mesh` or a `bpy.types.Object`, and with `apply_modifiers=True` evaluates the object's modifier stack first. The result holds `triangle`, `quad`, `polygon<n>` and `line` blocks — always dimension 2 or below, because that is everything Blender has. `POINT` and `FACE` attributes come back as `point_data` and `cell_data`, and the `region:` masks are rebuilt into regions using the sidecar.

The export operator bakes each object's own placement into the coordinates (through `transform`, since no format carries an object matrix separately) and merges several selected objects with `merge`.

## Architecture: the pure payload layer

`src/python/meshioplusplus/_blender.py` is split the way [`_interop`](interop.md#architecture-the-pure-payload-layer) and the [MCP server](mcp.md) are. All of the behaviour lives in `_to_blender_payload` and its inverse `_mesh_from_blender_arrays`, which import no `bpy` at all and speak in flat numpy arrays; `to_blender` and `from_blender` are thin wrappers that import `bpy` inside the function. That is what makes the subtle parts — the loop/polygon layout, the block-major indexing, the attribute typing — testable in the default CI matrix with Blender nowhere in sight.

Two reuses carry most of the weight. The surface reduction is `_viewer_browser._renderable_surface`, *called* rather than copied — its docstring binds it to the WASM `convertSurface` binding, and a second implementation would put that contract out of reach of the file documenting it. And the inverse routes through `_vtk_common.vtk_cells_from_data`, the same function `from_pyvista` uses, because Blender's loop layout and VTK 9's cell array are the same object under different names: `loop_vertices` is the connectivity, `loop_total` the per-cell count, `loop_start` the offset.

## Limitations

Several of these are permanent, not deferred.

- **Volume cells become their boundary and nothing else.** Re-exporting an imported `.vtu` of tetrahedra yields triangles. This is Blender's data model, not a gap to close.
- **Regions are lost on a volume input.** Extracting a boundary drops them; the parent gather recovers `cell_data` but not regions. Surface meshes keep everything. Same cause and same behaviour as `to_trimesh`.
- **A volume mesh's own supplied boundary blocks are discarded** in favour of the computed boundary — already the browser viewer's behaviour, and consistency with it is worth more than a special case.
- **Higher-order nodes are gone** to linearization: `triangle6` becomes `triangle`, `quad9` becomes `quad`, and a curved element renders flat.
- **`CORNER` and `EDGE` domain attributes survive in neither direction** — UV maps and per-corner normals among them — because meshio++ has point and cell data locations and no others.
- **Blender's integer attributes are 32-bit.** A value outside that range is stored as a float, with a note.
- **The file-dialog filter is trimmed.** meshio++ knows more extensions than Blender's `filter_glob` buffer holds, so the least common ones are not listed by default; the import operator's *Show all files* toggle clears the filter, and no format is ever unreachable.
- **Platforms:** `linux-x64`, `windows-x64`, `macos-x64`, `macos-arm64`. Not Windows on ARM (no such wheel is built) and not Linux on ARM (Blender ships no official build).
- **glibc floor.** The Linux wheel targets glibc 2.34 while Blender's own Linux builds target roughly 2.28, so a system old enough to be between the two can run Blender but cannot load this extension.
- **Blender and Python move together.** The compiled core is `cp311`, matching the CPython that Blender 4.2–4.5 embeds. A Blender that moves to 3.12 needs a rebuilt extension.
- **numpy is inherited, not bundled**, deliberately: a second numpy ahead of Blender's own is a known way to break `bpy`'s numpy interop.
