# Boundary layers: plan

Roadmap 6.9's design pass: what a native boundary-layer mesher in the core needs, which of the established schemes it follows, what it leaves, and why it is native rather than a backend. The operation that came of it is [`boundary_layer`](./boundary_layer.md).

## The problem

A viscous flow needs thin cells along its walls: layers whose thickness grows from a first height set by the wall's `y+` to the size of the volume mesh, aligned with the wall so the gradient across it is resolved. The standard shape is a prism stack grown from the wall's surface mesh along its normals — wedges from triangles, hexahedra from quads — with the rest of the volume meshed against the last layer (Pirzadeh 1996, Garimella and Shephard 2000). Three things make it hard: what "the normal" is at a ridge or a corner, layers that run into each other or into another part of the surface, and cells that fold where the surface is concave.

## Choices

**The direction.** The normal at a point shared by faces meeting at an angle is not their average: averaged normals at a convex corner leave layers thinner along the faces than asked, and at a concave one they can point into a face. Aubry and Löhner (2009) define the *most normal* normal as the unit vector maximising the least dot product with the faces' normals — the direction that sees every face best — and show it is the right choice for layer growth at ridges and corners. It is found here by ascent from the angle-weighted mean toward the face that sees it least, with steps shrinking as `1 / (k + 2)`; for the cube's corner it reaches `(1, 1, 1) / √3`. Dividing the offsets by the least cosine (Garimella and Shephard's stretching, capped at four) keeps every face's layers their thickness; the cube's cap then lies exactly on the offset planes. A point with no direction seeing all its faces (the least cosine under 0.05) cannot grow without folding a cell, and does not.

**Collisions.** A layer stack must stop before it meets the surface on the other side of a gap, or another stack growing toward it. A ray per point along its direction against the surface's triangles (the core's bounding-volume hierarchy and Möller–Trumbore test, already used by `octree_mesh`'s gap sizing) gives the free distance; taking at most 45% of it leaves both stacks across a gap short of the middle. The share each point keeps is smoothed over its neighbours by taking, each pass, the lesser of its own share and its mean with theirs, so a share is never raised past what its own ray allows and the layers taper over a few cells instead of stepping. Checking layer against layer (the advancing-layer method's front intersection test, Pirzadeh 1996) is not done: the ray and the validity check below catch the cases that matter for a mesh read back, and the front test needs a triangle–triangle kernel the core does not have.

**Folding cells.** Where the surface is concave, a stack grown at full height folds its cells over one another. Every cell is judged by its least corner Jacobian; the points of a folding cell have their stacks halved, and after six halvings collapsed. Collapsing rather than stopping keeps the mesh conforming without a front: a collapsed point stays on the surface in every layer, a wedge with one collapsed corner becomes a pyramid on its opposite side face and with two a tetrahedron, and a quad with a collapsed corner is split into two triangles' columns in every layer, so the faces the cells share agree. Garimella and Shephard grow fewer layers there and close the stack with pyramids and tetrahedra on top; collapsing to zero is the simpler end of the same idea and is what is done.

**Sweep.** Copying each point at L levels and joining the levels into cells is `extrude`'s operation too ([roadmap 3.1.2](./roadmap.md#_3-operations)). It is written once, in `detail/layer_sweep`, taking a direction and offsets per point, so `extrude` is a uniform direction and offsets on the same kernel.

## Native or a backend

The well-known open-source boundary-layer meshers are GPL (cfMesh, snappyHexMesh, Gmsh) and excluded from the core; Netgen's layers (LGPL-2.1) would be a bring-your-own backend. The native route needs only the core's own kernels — the normals, the ray cast, the validity measure and the sweep — at a few hundred lines, so it is native, which narrows the [non-goal](./roadmap.md#non-goals-and-decisions-taken) that left boundary-layer inflation to backends: layers grown from a surface are the core's, the volume behind them stays with roadmap 6.7's tetrahedral backends (or `octree_mesh`, fitted rather than conforming).

## What is left

- The volume behind the last layer, conforming to it: roadmap 6.7.
- Layer-against-layer intersection, and growing fewer layers rather than collapsing a point.
- Normals smoothed over the surface for very coarse meshes (Garimella and Shephard smooth them by Laplacian passes; the most normal normal is computed point by point).
- Curved (`triangle6`, `quad9`) surfaces.

## References

- A. Aubry and R. Löhner, "On the 'most normal' normal", *Communications in Numerical Methods in Engineering* 25 (2009).
- S. Pirzadeh, "Three-dimensional unstructured viscous grids by the advancing-layers method", *AIAA Journal* 34 (1996).
- R. V. Garimella and M. S. Shephard, "Boundary layer mesh generation for viscous flow simulations", *International Journal for Numerical Methods in Engineering* 49 (2000).
