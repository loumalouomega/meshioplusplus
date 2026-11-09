# Clipping and capping

`clip(mesh, origin, normal, keep="below", cap=True)` cuts a surface by a plane, keeps one side with the cells that cross it split at the plane, and caps the cut so a closed part stays closed (roadmap 7.8.3). [`crop`](./crop.md) keeps or drops whole cells and [`slice`](./slice.md) returns only the section; this is the operation between them. It is named `clip` because `cut` is already the [level-set cutting](./cut.md) operation and CLI verb.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

lower, report = meshioplusplus.clip(part, (0, 0, 10), (0, 0, 1), return_report=True)
report["num_loops"], report["num_cap_triangles"]
upper = meshioplusplus.clip(part, (0, 0, 10), (0, 0, 1), keep="above")
```

| Option | Default | Meaning |
|---|---|---|
| `origin`, `normal` | | the plane; `d = (x - origin) . normal`, a zero normal is refused |
| `keep` | `"below"` | `"below"` keeps `d < 0`, `"above"` keeps `d >= 0` |
| `cap` | `true` | triangulate the cut's closed loops |

The input is a surface: triangles, quads and polygons, simplexified as `slice` simplexifies them; a volume block is refused by name. A triangle wholly on the kept side stays, one across the plane is split at the crossings `slice` computes (the sorted-edge formula, so the cut boundary is `slice`'s section bit for bit), and the rest go. The cut's crossings are chained ([chained sections](./slice.md#chained-sections)); each closed loop is triangulated by the planar CDT on the crossings themselves, with no new points, nested loops by the even-odd rule, the cap facing out of the kept side. Open chains (a cut through an open surface's border) are counted, not capped.

The result is one `triangle` block. Point data is interpolated at the crossings (Float64) and kept elsewhere; cell data is the parent cell's, NaN (floating) or -1 (integer) on the cap, which the Int64 cell data `clip:cap` marks with 1. The report gives `num_kept`, `num_split`, `num_cap_triangles`, `num_loops`, `num_open_chains` and `warnings`.

A closed sphere clipped at its equator gives a closed hemisphere whose volume is half the sphere's, and the two halves' volumes add up to the sphere's to `1e-12` (`tests/cpp/test_clip.cpp`, `tests/python/test_scan_edit.py`).

## Surfaces

C++ `clip` (`operations/clip.hpp`), Python `meshioplusplus.clip` (C++ core only), both CLIs (`clip IN OUT [--origin x,y,z] [--normal x,y,z] [--keep below|above] [--no-cap] [--json]`), the `Clip` [pipeline](./pipeline.md) step (`Origin`, `Normal`, `Keep`, `Cap`), C `mio_clip(mesh, origin, normal, keep_below, cap, &report)` (`mio_clip_report`, 72 bytes; `mio_clip_warnings`), Fortran, Julia, R, Octave/MATLAB, WASM `clip` and the MCP tool `clip`.
