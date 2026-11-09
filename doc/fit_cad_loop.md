# Closing the reverse-engineering loop

[`fit_cad`](./fit_cad.md) turns a skin into a solid in one pass and reports its Hausdorff deviation. Two operations close the loop around it (roadmap 7.7.1, 7.7.2): `face_support` says how much of each face of the solid the scan actually supports, and `fit_cad_iterate` refits towards a target, re-segmenting the faces that fall short.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

# A scan missing its bottom, the hole closed by repair (which marks the fill).
closed = meshioplusplus.repair(scan, max_hole_edges=0, record_provenance=True)
skin = meshioplusplus.fit_primitives(meshioplusplus.segment(closed, method="primitive"))
solid = meshioplusplus.fit_cad(skin)
marked, report = meshioplusplus.face_support(solid, closed, return_report=True)
report["num_measured"], report["num_inferred"]  # 5, 1: the bottom was guessed

best, rounds = meshioplusplus.fit_cad_iterate(closed, "part.step", part_class="cast", return_report=True)
[r["share_within"] for r in rounds["rounds"]], rounds["stop_reason"]
```

## Support per face

`face_support(solid, scan)` samples each face of a `fit_cad` re-mesh (its `fit_cad:face` cell data; each triangle's four sub-triangle centroids, weighted by a quarter of its area) and measures every sample's distance to the scan. A face is **measured** when at least `measured_share` (0.9) of its area lies within `tolerance` of the scan, **inferred** below `inferred_share` (0.1): no scan near it (a hidden interior, the far side of a through-hole, a face the solid closed by construction), and **partial** between. Scan cells whose `hole_data` cell data (`repair:hole`, `repair`'s provenance with `record_provenance=True`) is at least 0 are the fills `repair` added, and count as unmeasured, so a hole `repair` closed before `fit_cad` reads as inferred. `tolerance` defaults to `1e-3` of the scan's diagonal.

The solid gains the field data `fit_cad:face_support` (Int64 per face: 0 measured, 1 partial, 2 inferred) and `fit_cad:face_measured` (Float64, the share), and the cell data `fit_cad:support` (each triangle's face's class). The report gives `face_share`, `face_class`, the three counts, `share_within` (the share of the scan's points, on measured cells, within the tolerance of the solid) and `tolerance`.

A box scanned without its bottom, closed by `repair` and fitted, has five measured faces and an inferred bottom (`tests/cpp/test_fit_cad_loop.cpp`).

## Targets and the iterating driver

`fit_cad_iterate(scan)` segments the scan with [`segment`](./segment.md)'s primitive method, fits with [`fit_primitives`](./fit_primitives.md) and builds with `fit_cad` (`face_deviation` on), in rounds. Each round measures `share_within` (as `face_support` does, at `fit_cad`'s tolerance) and `explained` (the share of the solid's area on analytic or feature faces, `fit_cad:type` other than 0), then:

1. stops at the **target**: `share_within >= target_within` and `explained >= target_explained` (0: off);
2. stops at **no improvement**: a round no better than the best so far (by `share_within`, then `explained`);
3. stops after `max_rounds` (4);
4. otherwise re-segments the patches of the `worst_faces` (3) worst faces (those whose face deviation exceeds the tolerance, largest first, then the largest faceted faces) with the segment tolerance times `refine_factor` (0.25), splicing their new labels in after the others; with nothing to refine it stops.

The best round's solid is kept (and written to `step_path` as STEP); each round reports `share_within`, `explained`, `hausdorff_max`, `segment_tolerance`, `num_patches`, `num_faces` and the patches it re-segmented. `part_class` sets the defaults a part's process suggests:

| `part_class` | `target_within` | tolerance (of the diagonal) |
|---|---|---|
| `any` (default), `analytic` | 0.90 | `1e-3` |
| `cast` | 0.85 | `2e-3` |
| `freeform` | 0.75 | `4e-3` |

`segment_tolerance` and `fit_tolerance` override the class's tolerance for the first segmentation and for `fit_primitives`. Deterministic: the same input and options give the same rounds. A box whose top bends in its middle, segmented coarsely, merges the two halves into one patch no plane explains; the second round separates them and `explained` rises (`tests/cpp/test_fit_cad_loop.cpp`).

## Surfaces

C++ `face_support`, `fit_cad_iterate` and `fit_cad_part_class_from_name` (`operations/fit_cad_loop.hpp`); Python `meshioplusplus.face_support` and `fit_cad_iterate` (C++ core only); both CLIs (`face-support SOLID SCAN [OUT] [--tolerance T] [--measured-share S] [--inferred-share S] [--hole-data NAME] [--json]`, `fit-cad-iterate SCAN OUT [--step PATH] [--part-class C] [--target-within S] [--target-explained S] [--segment-tolerance T] [--fit-tolerance T] [--max-rounds N] [--worst-faces N] [--refine-factor F] [--backend B] [--json]`); the `FaceSupport` (`Scan`) and `FitCadIterate` [pipeline](./pipeline.md) steps; C `mio_face_support` (`mio_face_support_opts`, `mio_face_support_report`, 64 bytes each) and `mio_fit_cad_iterate` (`mio_fit_cad_iterate_opts`, 112 bytes; `mio_fit_cad_iterate_report`, 96; `mio_fit_cad_iterate_round`); Fortran, Julia, R, Octave/MATLAB, WASM `faceSupport` and `fitCadIterate`; the MCP tools `face_support` and `fit_cad_iterate`.
