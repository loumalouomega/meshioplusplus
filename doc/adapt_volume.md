# Size-field tetrahedral remeshing (`adapt_volume`)

`adapt_volume` turns a valid but over-refined tetrahedral mesh, typically the graded output of [`octree_mesh`](octree_mesh.md), into one with far fewer and better-shaped tetrahedra. It works towards an isotropic size field, with a gradation and a bound on how far the boundary may move, and it keeps corners, ridges, materials, regions and boundary tags. It is native, MIT-licensed and has no dependency: an alternative to post-processing with the LGPL [MMG backend](mmg.md).

This page is being written alongside the implementation; the algorithm, options, guarantees and measured tables follow.
