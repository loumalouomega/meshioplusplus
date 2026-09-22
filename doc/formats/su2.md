# SU2 (`.su2`)

The [SU2](https://su2code.github.io/docs_v7/Mesh-File/) mesh format: an ASCII format with `NDIME`, `NPOIN`, `NELEM` volume cells and `NMARK` boundary markers, using VTK-style numeric type codes. A single file can also hold several independent [zones](https://su2code.github.io/docs_v7/Multizone/) (`NZONE=`/`IZONE=`, v15.5.0).

| | |
|---|---|
| **Format name** | `su2` |
| **Extensions** | `.su2` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.su2")
meshioplusplus.su2.write("out.su2", mesh)
```

`write` takes no keyword arguments.

## File structure

Line-oriented `KEY= value` records (`%` = comment; blank lines skipped; malformed lines without `=` just warn and are skipped).

- `NDIME= 2|3` — dimension.
- `NPOIN= n [extra columns]` — the first point row is read separately to auto-detect whether an extra trailing global-index column is present (some SU2 files append one); the remaining `n-1` rows follow with the same column count, then the extra column(s) are stripped from all rows.
- `NELEM= n` or `MARKER_ELEMS= n` — reads exactly `n` element lines; the first line is peeked to determine the node count (from the VTK type code) and whether an extra trailing column is present, then all `n` lines are parsed as one integer block and binned by VTK type code (elements of different types within one block are separated into distinct cell blocks). `NELEM` cells get `su2:tag = 0`; `MARKER_ELEMS` cells get the current marker's tag id.
- `NMARK= n` — expected marker count (soft-checked; a mismatch only warns).
- `MARKER_TAG= <value>` — sets the tag used for the following `MARKER_ELEMS` block. A numeric value is used as `su2:tag` directly; a non-numeric one gets an auto-incrementing integer `su2:tag` **and** its own text is kept as the marker's name (see below) instead of being discarded.
- Boundary cell blocks of the **same type** are merged into one block per type after the full-file scan (across every zone, for a multizone file), concatenating their `su2:tag` (and, for a multizone file, `su2:zone`) values in the same order.

Write: `NDIME=` from `points.shape[1]`; `NPOIN=` + coordinates; volume cells (`triangle`/`quad` for 2D, `tetra`/`hexahedron`/`wedge`/`pyramid` for 3D) under `NELEM=`, each row prefixed by its type code; boundary markers grouped by the first integer-typed `cell_data` array found (via the same "first-int-array" convention used by several other formats, and excluding `su2:zone` from the candidates), one `MARKER_TAG=`/`MARKER_ELEMS=` pair per distinct tag value. A marker's `MARKER_TAG=` text is its matching named region (see below) when one exists, else the bare integer tag.

### Multizone (`NZONE=`)

A single file can hold several independent zones: `NZONE= n` (always the first key) followed by `n` `IZONE= i` sections, each a complete standalone mesh (its own `NDIME`/`NPOIN`/`NELEM`/`NMARK`). Zones are **never welded** — each keeps its own point numbering, offset when concatenated into one mesh — matching SU2's own semantics, where `IZONE`s are read by separate solver instances with no shared degrees of freedom. `cell_data["su2:zone"]` (only present for a multizone file) names each cell's zone; `su2:tag` is scoped per zone (a marker's auto-incremented id starts over at each `IZONE`).

Write: a mesh whose `su2:zone` carries more than one distinct value writes `NZONE=`/`IZONE=` sections, one per distinct zone value, each with its own point subset (only the points its cells actually reference, renumbered densely from 0) and its own `NMARK`/markers. A mesh with no `su2:zone`, or only one distinct value, writes the single-zone layout (no `NZONE=` line at all).

### Marker names as regions

A marker whose `MARKER_TAG` is a **name** (not a plain integer) gets a [`Region`](../regions.md) of kind `cell` over its boundary cells, so the name survives a round trip instead of collapsing into an anonymous `su2:tag` id: `"<name>"` for a single-zone file, `"zone_<i>/<name>"` for a multizone one. A multizone file also gets one `"zone_<i>"` region per zone, covering every cell (volume and boundary) of that zone. A purely numeric `MARKER_TAG` gets no region — `su2:tag` already carries it losslessly.

On write, **any** Cell region whose cells share one `su2:tag` value supplies that marker's name (for a multizone write, only a region already named `"zone_<i>/..."` is read back for zone `i`) — this is deliberately permissive, the same "best-effort, first match" spirit as the tag-source convention above, not a strict "only round-tripped regions" check. There is no per-volume-cell tag in the SU2 format itself (every `NELEM` row is implicitly `su2:tag = 0`), so a Cell region over *volume* cells has nothing to attach a marker to and is silently dropped on write.

## Cell types

| code | nodes | meshio++ type |
|---|---|---|
| 3 | 2 | `line` |
| 5 | 3 | `triangle` |
| 9 | 4 | `quad` |
| 10 | 4 | `tetra` |
| 12 | 8 | `hexahedron` |
| 13 | 6 | `wedge` |
| 14 | 5 | `pyramid` |

## Data mapping

- `cell_data["su2:tag"]` — volume cells always get tag `0`; boundary (`MARKER_ELEMS`) cells get their marker's tag id (an auto-incrementing int starting from 1 for the first non-numeric string tag encountered in that zone, if any).
- `cell_data["su2:zone"]` — which `IZONE` a cell came from; present only when the file is multizone.
- Named [regions](../regions.md) — see "Marker names as regions" above.

## Quirks & limitations

- Unsupported cell types on write trigger a warning, but a latent message bug means the warning text shows Python's `type` builtin rather than the actual cell type name — cosmetic only, doesn't affect behavior.
- Only **one** integer cell-data array can be used as the boundary-marker tag source on write (`su2:zone` is excluded from the candidates, since it drives zone splitting, not markers); if other candidates exist, they're dropped with a warning.
- Point coordinates use `np.savetxt`'s default format in the Python writer (`%.18e`) vs. `%.16e` in the C++ writer — a formatting difference, not a round-trip correctness issue.
- A multizone write always duplicates any point shared between zones (each zone's own point subset is written independently) — this matches SU2's own zone independence, not a limitation of the mapping.

## Notes

- `tests/python/meshes/su2/square.su2` (official SU2-docs example: structured 2D quad mesh, 9 points/8 quads/markers) — checked for 16 total cells, 9 points, 4 unique tags summing to 20.
- `tests/python/meshes/su2/mixgrid.su2` (mixed-cell 3D grid) — checked for 30 cells, 16 points, 6 unique tags summing to 62.
- Multizone reading, writing and marker-name regions are covered by hand-written fixtures in `tests/python/test_su2.py` and `tests/cpp/test_misc_formats.cpp` — there is no small, MIT-compatible real-world multizone `.su2` file to check against, unlike `square.su2`/`mixgrid.su2`.
- Fully handled by the C++ core.
