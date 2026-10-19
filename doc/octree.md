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

## Hanging nodes

Where a fine leaf meets a coarse one the output has hanging nodes, flagged in `octree:hanging` and counted in the report (`num_hanging`). This is what adaptive-mesh-refinement, lattice-Boltzmann, immersed-boundary and cut-cell solvers and octree learning models consume. A node is flagged relative to the leaves actually kept, so a `surface` or `inside` fill does not report nodes that a dropped neighbour would have hung.

## Determinism

Every point is an integer coordinate on the finest level's lattice (`2^max_level` per root cell), so shared corners are deduplicated arithmetically, with no hash and no tolerance, and coordinates are `origin + spacing * index`, never accumulated. Leaves come out in Morton order and points in lattice order, so the result is byte-identical across parallel backends, thread counts and the C++ core and the numpy reference implementation. `max_level` is at most 20, and `root_cells << max_level` must stay below 2²¹ on every axis. `max_cells` (default 20 000 000) refuses a tree with more leaves.

![an octree around a sphere, cut open](/images/octree_sphere.png)

## Report

`return_report=True` returns `dimension`, `origin`, `root_size`, `num_leaves`, `num_hanging` and `finest_level`.

## Not yet

Conforming transitions (`transition="polyhedra"` and `"templates"`), snapping the outer nodes to the skin, and the curvature, feature-edge and size-field refinement criteria are tracked in the [roadmap](roadmap.md#3-operations).

## Bindings

Python `mp.octree`, C `mio_octree`, Fortran `mesh%octree`, Julia `octree`, R `mio_octree`, WASM `octree`, the pipeline step `Octree`, the CLI verb `octree` and the MCP tool `octree` all take the same options; see [C API](c_api.md), [Fortran](fortran.md), [Julia](julia.md), [R](r.md), [WASM](wasm.md), [pipelines](pipeline.md), [CLI](cli.md) and [MCP](mcp.md).
