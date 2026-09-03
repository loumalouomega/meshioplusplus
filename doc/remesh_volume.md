# Volume retetrahedralization (isosurface stuffing)

`meshioplusplus.remesh_volume(mesh, cell_size=..., ...)` retetrahedralizes a **volume** mesh (or a closed surface) at a caller-chosen resolution by isosurface stuffing, discarding the input's own tetrahedra entirely and generating a fresh mesh from a body-centered cubic (BCC) lattice. It is the volumetric sibling of [`remesh`](/remesh) — the same relationship [`decimate_volume`](/decimate_volume) has to [`decimate`](/decimate). Nothing else in this repo can *raise* the quality of a tetrahedral mesh at a chosen target resolution: [`refine`](/refine) subdivides the input's own (possibly bad) cells, `decimate_volume` can only remove elements, [`smooth`](/smooth) (see "ODT smoothing" below) only moves points. `remesh_volume` instead generates a fresh lattice-based mesh, so output element quality is a property of the lattice construction, not of the input.

```python
import meshioplusplus

mesh = meshioplusplus.read("part.vtu")  # a volume mesh, or a closed surface

out = meshioplusplus.remesh_volume(mesh, cell_size=0.5)
out = meshioplusplus.remesh_volume(mesh, resolution=(64, 64, 64))

out, report = meshioplusplus.remesh_volume(mesh, cell_size=0.5, return_report=True)
print(report["num_tets"], report["num_vertices_warped"], report["num_non_manifold_edges"])
```

Both CLIs expose it as a verb:

```bash
meshioplusplus remesh-volume part.vtu out.vtu --cell-size 0.5
meshioplusplus remesh-volume part.vtu out.vtu --resolution 64,64,64 --warp-fraction 0.2
```

## Why a lattice, not Delaunay / CVD / ODT remeshing

[`doc/roadmap.md`](/roadmap)'s "Volumetric CVD/ODT" bullet names the literal 3D counterpart of `remesh`'s own surface clustering — optimal Delaunay triangulation / centroidal Voronoi tessellation in 3D. That method needs robust orientation and in-sphere predicates, exactly where this project's dependency-free, MIT posture stops paying: the permissive predicate library (Geogram, BSD-3) is a large new dependency, and the smaller, purpose-built alternatives are copyleft (TetGen is AGPL-3.0, CGAL's `Mesh_3` is GPL) — a poor fit even as an optional dependency for an MIT project.

Isosurface stuffing needs no predicate at all — only a signed-distance field, which this repo already has ([`sample_distance`](/sdf)) — and every uncut lattice tetrahedron has a dihedral angle from a small, mesh-size-independent fixed set, rather than an unbounded one. It closes the bullet's *capability* (retetrahedralize at a chosen resolution with bounded element quality), though it is honestly not the named method. [`SmoothMethod::Odt`](/smooth#odt-smoothing) (below) closes the "ODT" half of the bullet's name separately and honestly: ODT *smoothing* on an existing tet mesh's fixed connectivity, not ODT *remeshing*.

### Attribution and licence

Implemented from the **published description only** of Labelle & Shewchuk, "Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral Angles" (SIGGRAPH 2007). Neither the paper's own reference implementation (Stellar, non-commercial licence) nor TetGen (AGPL-3.0, also by Shewchuk) is read or vendored here — the same clean-room posture [`doc/remesh.md`](/remesh) documents for its ACVD attribution, applied to a second reference. See `CITATION.cff`.

## The algorithm

1. **Lattice.** A body-centered cubic lattice — two interleaved simple cubic lattices, the second offset by half a cell on every axis — is generated over the (padded) bounding box, sized with the exact same vocabulary [`compute_sdf`](/sdf)/`VoxelOptions` use (`resolution` / `cell_size` / `bounds` / `padding` / `padding_relative`), reused verbatim rather than reinvented. A per-axis `resolution` that would resolve to non-cubic cells is refused by name.

   Its natural tetrahedralization — 12 congruent tets per cell, each cell's body-center connected across one fixed diagonal of each of the cell's 6 faces — gives every **uncut** tet a dihedral angle from the fixed set **{45°, 60°, 90°, 120°}**, independent of cell size (every cell is a translate of the same unit construction). This is this project's own derivation, not a transcription of Labelle & Shewchuk's own tetrahedralization (never read): the paper reports a tighter {60°, 90°, 120°} set for its own construction; this one is a deliberately simpler, independently-derived and numerically-verified alternative — the 12-tet table tiles the unit cell with no gap or overlap, every tet is positively oriented, and the tets' volumes sum to the cell's own — that still delivers a small, mesh-size-independent bound, reported honestly as what it is rather than as a reproduction of the paper's exact angles.
2. **Classify.** Every lattice vertex gets a signed distance to the surface and an inside/outside label from that sign, fixed for the rest of the algorithm.
3. **Warp.** A lattice vertex within `warp_fraction * h` of the surface (`h` the lattice's own cell size) is moved — its label unchanged — onto the nearest point on the surface. This is what removes the arbitrarily thin slivers a naive cut would otherwise leave near the boundary; warping is **load-bearing, not an optimization** (see "The warp/quality tradeoff, measured" below).
4. **Cut.** A tet whose four vertices share one label is kept whole or discarded. A straddling tet is cut by a small sign-mask case table — the same idea marching tetrahedra uses for a level *surface* ([`detail/marching.hpp`](/slice)), here applied to fill the interior volume instead. Each crossing edge gets one cut point, reusing an incident vertex's warped position when one of the edge's two endpoints was warped, falling back to linear interpolation otherwise.

Determinism follows the repo's standard phase split: signed distances and warps are computed in parallel into disjoint per-vertex slots; cut points are deduplicated by a serial, ascending-edge-key pass (the `refine.cpp`/`surface.cpp` idiom), so two lattice tets sharing a face agree bit-for-bit on that face's cut points.

![A body-centred cubic lattice cell with its tetrahedra, then the classify, warp and cut steps of isosurface stuffing on a 2-D cut](/diagrams/remesh_volume_bcc.svg)

## The warp/quality tradeoff, measured

Every claim below is a measurement on this repo's own test fixtures (a 12×24 or 16×32 UV sphere), not an assertion carried over from the paper — the `RemeshOptions::max_anisotropy` precedent from the surface `remesh` operation, applied here.

Warping genuinely improves boundary-tet quality, even at a small fraction — on a 16×32 sphere at `cell_size=0.3`, the worst (minimum) dihedral angle in the mesh rises from about 0.1° unwarped to several degrees at even a small `warp_fraction`, and keeps improving as the fraction grows toward its default:

| `warp_fraction` | vertices warped | worst min dihedral |
| --- | --- | --- |
| 0.00 (off) | 0 | ~0.1° |
| 0.01 | 14 | ~0.9° |
| 0.05 | 47 | ~4.6° |
| 0.10 | 48 | ~10.0° |

Reusing a warped vertex's position for every cut edge it touches — the mechanism that carries this quality benefit into the cut region — can also leave a **small number of the output's own boundary edges non-manifold**, on some inputs, at `warp_fraction > 0` only. This is reported honestly as `RemeshVolumeResult.num_non_manifold_edges` rather than hidden. On a 12×24 sphere at `cell_size=0.25` (the fixture `tests/cpp/test_remesh_volume.cpp` uses):

| `warp_fraction` | vertices warped | non-manifold edges | boundary triangles | fraction |
| --- | --- | --- | --- | --- |
| 0.00 (off) | 0 | 0 (exact) | 1360 | 0.0% |
| 0.01 | 2 | 0 | — | — |
| 0.05 | 39 | 5 | — | — |
| 0.15 | 103 | 17 | — | — |
| 0.25 | 199 | 64 | — | — |
| 0.35 (default) | 273 | 78 | 1360 | 3.8% |

In every measured case, **orientation stays exact** (`num_inverted == 0`) and the boundary stays **hole-free** (`boundary_edges == 0`, `inconsistent_pairs == 0`) — the defect is confined to a small fraction of edges being shared by more than two boundary triangles, never a hole or a winding disagreement. A caller needing an exactly watertight boundary should either set `warp_fraction=0` (mathematically exact, at the cost of the quality improvement above) or post-process the result with [`clean(weld=True)`](/clean).

## What does NOT survive

Identical to `remesh`: the output has entirely new points and new connectivity, so there is no point or cell map, and `point_data`/`cell_data`/regions are dropped (a warning is logged, never silent); `field_data` carries through. Compose with [`interpolate`](/interpolate)/[`conservative_interpolate`](/conservative_interpolate) to transfer a field onto the result.

## Scope

Accepts a **volume** mesh (its boundary is extracted internally via [`extract_surface`](/extract_surface), which is closed by construction) **or** a **closed surface** directly (triangle, quad, rectangular polygon) — unlike plain `remesh`, which only ever works on a surface and throws on a volume input. The boundary of the result is only ever as good as the lattice resolution and the warp step can make it: sharp input edges and corners round off (output boundary vertices are warped lattice vertices), and features thinner than a cell merge or vanish. A caller wanting a genuinely high-quality boundary composes with `extract_surface` + `remesh` on the result rather than expecting this operation to build one in — this is a deliberate composition, not a built-in.

C++-core only, with **no numpy fallback at all** — the `subdivide`/`agglomerate`/`decimate_volume`/`remesh` precedent: warp and cut are discrete branches on a sign near a threshold, so a second, independently-written implementation could land on the other side of a near-tie and diverge into a different mesh, not a last-ulp difference. `meshioplusplus.remesh_volume` raises `NotImplementedError` naming the reason when `_core` is unavailable, for any input.
