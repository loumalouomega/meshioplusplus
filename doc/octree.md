# Octrees and quadtrees from a skin — `octree`

`octree` fills a skin with an adaptive lattice: leaves are small where the skin is and large everywhere else. A surface gives an octree of `hexahedron` leaves, a planar closed polyline gives a quadtree of `quad` leaves. It is [`voxelize`](voxelize.md)'s adaptive sibling and, like it, replaces its input's geometry.

```python
import meshioplusplus as mp

skin = mp.read("bunny.stl")
tree = mp.octree(skin, max_level=6)                      # fine only near the surface
solid = mp.octree(skin, max_level=6, fill="inside")      # the interior, adaptively
band = mp.octree(skin, max_level=6, sample_distance=0.05)  # refine within 0.05 of the skin
mp.write("tree.vtu", tree)
```

A uniform tree (`min_level == max_level`) is exactly the lattice [`voxelize`](voxelize.md) builds, so the 3-D voxel grid is the degenerate octree.

![a quadtree around a circle, coloured by level, with the hanging nodes marked](/images/octree_quadtree.png)

## The output is an ordinary mesh

Leaves are one `hexahedron` block (`quad` in 2-D), so every writer, viewer and operation already works on them. Three arrays describe the tree:

| array | where | meaning |
|---|---|---|
| `octree:level` | `cell_data`, Int64 | 0 for a root cell, `max_level` for the finest leaves |
| `octree:cut` | `cell_data`, Int64 | 1 when the skin passes through the leaf |
| `octree:hanging` | `point_data`, Int64 | 1 for a node on the interior of a face or edge of a leaf rather than at its corners |

`attach_distance=True` adds the signed distance at leaf centres as `sdf:distance`.

## How the tree is built

- **Root.** The padded bounding box (or `bounds`) is covered by `root_cells` cubic root cells, one tree each; the default is a single root. The root edge is the largest of `extent / root_cells` over the active axes.
- **Refinement.** A leaf is split while its level is below `max_level` and either the skin cuts it (exact triangle–box overlap in 3-D, segment–box in 2-D) or lies within `sample_distance` of it; every leaf is split down to `min_level`.
- **Balance.** `balance="face"` (the default) keeps neighbours across faces (edges in 2-D) within one level of each other, `"full"` also across edges and corners, `"none"` leaves the tree as refined. Balancing ripples to a fixed point.
- **Fill.** `fill` keeps `voxelize`'s meaning: `"all"` leaves, only the leaves the skin `"surface"` cuts, or only the leaves whose centre is `"inside"`.

## Planar skins: the quadtree

A mesh of `line` cells (a closed polyline) selects the 2-D variant automatically; `dimension=2` forces it for a mesh of `triangle`/`quad` cells, whose boundary edges become the skin. The skin lies in the plane of its points (`z`), the leaves are `quad` cells at that `z`, and padding applies to the two in-plane axes. `fill="inside"` uses the winding number of the polyline (`even_odd=True` selects the even-odd rule for self-overlapping input) and refuses a polyline that does not close, naming its open ends.

## Transitions: `transition=`

Where a fine leaf meets a coarse one there are three outputs.

- **`"hanging"`** (the default): plain hexahedra (quads in 2-D) with hanging nodes, flagged in `octree:hanging` and counted in the report (`num_hanging`). This is what adaptive-mesh-refinement, lattice-Boltzmann, immersed-boundary and cut-cell solvers and octree learning models consume. A node is flagged relative to the leaves actually kept, so a `surface` or `inside` fill does not report nodes that a dropped neighbour would have hung.
- **`"polyhedra"`**: conforming with **no new points**. A leaf next to a finer one becomes a `polyhedron` (a `polygon` in 2-D) whose faces are split into the fine faces, the split-hex representation snappyHexMesh produces; every other leaf stays a `hexahedron`/`quad`. It is readable by every format and finite-volume solver that takes [polyhedra](polyhedra.md) (OpenFOAM, Fluent, CGNS `NGON_n`, EnSight `nfaced`, VTU). A face with a hanging node on an edge gets that node as an extra vertex, so two leaves sharing a face always see the same polygon.
- **`"templates"`**: conforming and hex-dominant, for a finite-element code. A leaf with hanging nodes is coned to its centroid: plain quad faces give pyramids, and faces with hanging nodes are fanned from a new face-centre point into tetrahedra (triangles from the centroid in 2-D). New points (face centres, then leaf centroids) are appended after the lattice points; every other leaf stays a hexahedron/quad. Every cell has a positive volume and the volumes add up to the leaves' (the tests check this, and that each interior face is shared by exactly two cells, on random trees and on every neighbour pattern around a cell).

![the same quadtree with hanging nodes, as polygons, and as triangles](/images/octree_transitions.png)

Both conforming outputs need a balanced tree (`balance` `face` or `full`), carry `octree:leaf` (the leaf, in the `hanging` output's order, each cell was cut from) beside `octree:level`/`octree:cut` and `sdf:distance`, and have no `octree:hanging` array (`num_hanging` is 0). The cells are in blocks (the untouched leaves first, then the transition cells), in leaf order within a block. [`check_conformity`](check_conformity.md) holds them to it: on every neighbour pattern, the sphere and the quadtree it reports no hanging node and no partly overlapping face under `"polyhedra"` and `"templates"`, and under `"hanging"` it flags exactly the points `octree:hanging` marks.

## Snapping: `snap=True`

Moves the outer nodes of the kept leaves — the corners of faces with no kept leaf on the other side, except hanging nodes — onto the skin, for a body-fitted rather than stair-stepped boundary. In 3-D it is [`shrinkwrap`](shrinkwrap.md)'s projection; in 2-D the nearest point of the polyline. A leaf is never left inverted: any leaf with a non-positive corner Jacobian gets its moved nodes put back, repeatedly, and the report says how many were moved (`num_snapped`) and how many were put back (`num_refused`). It needs `fill` `surface` or `inside` (with `all` there is no outer boundary). ![a quadtree before and after snapping](/images/octree_snap.png)

Boundary-layer inflation stays a [non-goal](roadmap.md#non-goals-and-decisions-taken).

## Refinement criteria

By default every leaf the skin cuts is split down to `max_level`. `surface_level` stops cut leaves earlier, and the criteria take them further, up to `max_level`:

| option | refines where |
|---|---|
| `sample_distance` | the skin is within this distance of the leaf (also pulls in leaves it does not cut) |
| `curvature_threshold` | a triangle has a vertex whose \|`curvature:mean`\| exceeds this (3-D) |
| `feature_angle` | a sharp, open or non-manifold [feature edge](feature_edges.md) (dihedral above this many degrees) crosses the leaf (3-D) |
| `size_array` | the leaf edge exceeds the smallest value of this `point_data` array over the skin primitives reaching it |

## Planes: `plane=`

A 2-D skin lies in the plane of zero extent by default (`z`, then `y`, then `x`); `plane` `"xy"`, `"xz"` or `"yz"` chooses it. Leaves are emitted in that plane. The report's `origin` is in the skin's frame (the two in-plane axes first, the normal last). 2-D [`voxelize`](voxelize.md#2-d-pixels) takes the same `plane` and picks it the same way under `"auto"`.

## Determinism

Every point is an integer coordinate on the finest level's lattice (`2^max_level` per root cell), so shared corners are deduplicated arithmetically, with no hash and no tolerance, and coordinates are `origin + spacing * index`, never accumulated. Leaves come out in Morton order and points in lattice order, so the result is byte-identical across parallel backends, thread counts and the C++ core and the numpy reference implementation. `max_level` is at most 20, and `root_cells << max_level` must stay below 2²¹ on every axis. `max_cells` (default 20 000 000) refuses a tree with more leaves.

![an octree around a sphere, cut open](/images/octree_sphere.png)

## Report

`return_report=True` returns `dimension`, `origin`, `root_size`, `num_leaves`, `num_hanging`, `finest_level`, `num_snapped` and `num_refused`.

## Core only

`transition` other than `"hanging"`, `snap`, `surface_level`, `curvature_threshold`, `feature_angle` and `size_array` run in the compiled core; the numpy reference implements the hanging-node octree (and `plane` in 2-D) and says so by name rather than answering differently.

## Bindings

Python `mp.octree`, C `mio_octree`, Fortran `mesh%octree`, Julia `octree`, R `mio_octree`, WASM `octree`, the pipeline step `Octree`, the CLI verb `octree` and the MCP tool `octree` all take the same options; see [C API](c_api.md), [Fortran](fortran.md), [Julia](julia.md), [R](r.md), [WASM](wasm.md), [pipelines](pipeline.md), [CLI](cli.md) and [MCP](mcp.md).
