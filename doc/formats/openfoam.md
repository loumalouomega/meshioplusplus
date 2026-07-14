# OpenFOAM polyMesh (read-only)

A reader for [OpenFOAM](https://www.openfoam.com/)'s native `polyMesh`
representation: an unstructured, face-based mesh described by 4-5 sibling
files (`points`, `faces`, `owner`, `neighbour`, `boundary`) rather than a
single file. Supports both ASCII and binary (little-endian, `label=32/64`,
`scalar=32/64`) encodings, and reconstructs general polyhedra as well as
tetra/pyramid/wedge/hexahedron.

| | |
|---|---|
| **Format name** | `openfoam` |
| **Extensions** | `.foam` |
| **Read / Write** | ✓ / — |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("case.foam")                 # a <case>/case.foam marker file
mesh = meshio.openfoam.read("/path/to/case")     # or the case directory directly
mesh = meshio.openfoam.read("/path/to/constant/polyMesh")  # or polyMesh directly
```

There is no writer — `register_format("openfoam", [".foam"], read, {})` is
called with an empty writer map. `read(filename)` takes no keyword
arguments; the `polyMesh` directory is located from whichever of the three
input forms above is given (`_resolve_polymesh`): a `.foam` suffix looks for
`<parent>/constant/polyMesh`; a directory literally named `polyMesh` is used
as-is; any other directory is checked for `constant/polyMesh` then
`polyMesh` as subdirectories. A `FileNotFoundError` is raised if none match.

## File structure

Each `polyMesh` file (`points`, `faces`, `owner`, `neighbour`, `boundary`) is
an OpenFOAM "FoamFile": a header block declaring `format` (`ascii`/`binary`)
and (for binary) an `arch` string encoding `label=32|64`/`scalar=32|64`,
followed by a count and a parenthesized list.

```
FoamFile
{
    ...
    format      binary;
    arch        "LSB;label=32;scalar=64;";
}
// ...
<N>
(
<binary or ascii data>
)
```

- **Header detection** (`_detect_format`) scans line-by-line for a
  `format ...;` line and an `arch "...";` line, extracting `label=` /
  `scalar=` byte widths; defaults to ascii/8-byte if absent.
- **Binary `points`** (`vectorField`): `N (` followed directly by
  `N*3*scalar_bytes` raw floats, no per-row framing — read via a single
  `np.frombuffer`.
- **Binary `owner`/`neighbour`** (`labelList`): `N (` followed by `N*label_bytes`
  raw ints, same direct-buffer read.
- **Binary `faces`** (`faceList`, non-contiguous): each face is its own
  `labelList` — `<count> ( <count*label_bytes bytes> )` repeated `N` times.
  Read in two passes: a sequential ASCII scan locating each face's byte
  offset and node count (cheap, since `find(b"(")` only ever scans the short
  ASCII gap between faces, never binary data that might coincidentally equal
  `'('`), then a single vectorized gather of every face's binary blob via a
  cumulative byte mask — bounded peak memory even for tens of millions of
  faces (see `_RaggedArray`, a CSR-style `(conn, offsets)` pair standing in
  for a `list[list[int]]`).
- **ASCII** variants use simple regex/line-based parsing
  (`_parse_points_ascii`, `_parse_faces_ascii`, `_parse_int_list_ascii`) after
  comment-stripping (`/* */` and `//`) and header-skipping.
- **`boundary`**: a dict of `patch_name → {type, nFaces, startFace}`, parsed
  via a brace-matching regex over the whole (header/comment-stripped) text.
- **Cell reconstruction**: cell↔face topology is built once as a CSR
  `_RaggedArray` (`_cell_faces_csr`, vectorized via `argsort`+`bincount`) from
  `owner`/`neighbour`; per cell, each face is oriented outward (reversed if
  the cell is that face's neighbour, since the stored normal points
  owner→neighbour) and classified by `(n_faces, n_points)`:
  `(4,4)→tetra`, `(5,5)→pyramid`, `(5,6)→wedge`, `(6,8)→hexahedron`, anything
  else → a general `polyhedron` (kept as outward-oriented face lists). Each
  of the 4 named types has a dedicated orientation-fixing builder
  (`_build_tetra`/`_build_pyramid`/`_build_wedge`/`_build_hexahedron`) that
  computes a scalar triple product and flips the node order if it comes out
  negative, guaranteeing positive-volume connectivity regardless of the
  source mesh's face-normal convention.

## Cell types

Volume cells: `tetra`, `pyramid`, `wedge`, `hexahedron`, and general
`polyhedron<N>` (grouped by **unique node count** `N`, ragged per-cell face
lists stored as an `object`-dtype array — one `CellBlock` per distinct `N`).

Boundary (patch) faces: `triangle`, `quad`, and `polygon<N>` for `N > 4`
(grouped by vertex count `N`, one `CellBlock` per size — via
`_build_boundary_polygons`).

## Data mapping

- `cell_data["cell_tags"]` — per-cell-block tag array: `0` for every volume
  cell block, and a distinct negative "MED-style family id" `-(patch_index+1)`
  per boundary patch's face blocks (so a triangle patch and a quad patch on
  the *same* physical boundary would currently get *different* tag values —
  see Quirks).
- `mesh.cell_tags` — mesh-level attribute (not `cell_data`),
  `{family_id: [patch_name]}`, letting a MED write bridge these patch names
  through the same mechanism used for Gmsh physical groups (see
  [`med.md`](med.md)).
- `mesh.point_tags` — always set to `{}` (present for interface symmetry
  with the MED-derived tag convention; OpenFOAM has no point-tag concept).
- No point_data or field_data (OpenFOAM field files like `U`, `p`, `T` in the
  case's time directories are not read by this module — only the mesh
  topology under `constant/polyMesh`).

## Quirks & limitations

- **Read-only.** There is no `write` function at all.
- Degenerate volume cells that match a named type's `(n_faces, n_points)`
  signature but whose topology doesn't resolve cleanly (`_match_top` finds
  more or less than one vertical neighbour per base node) are **silently
  skipped** and logged as a warning count, rather than falling back to a
  general polyhedron.
- Boundary patches are tagged by **patch index**, not patch identity across
  face-size groups — if one named patch contributes both triangles and
  quads, its triangle `CellBlock` and quad `CellBlock` get the *same*
  `cell_tags` id (assigned once per patch, reused across whichever
  size-buckets that patch's faces fall into), but two *different* named
  patches always get distinct ids.
- All binary reads assume little-endian (`LSB`) — the format's own `arch`
  string is trusted for label/scalar width but not for byte order.
- Read goes through the C++ core (`meshio._core.openfoam_read`, using
  `std::filesystem` for the polyMesh directory), with the Python reference as
  an automatic fallback. General polyhedra cross the C++↔Python boundary via
  the ragged `polyhedron<N>` cell representation (a copied list of face
  arrays); boundary patch names travel through an `OpenFoamInfo` side-channel
  struct as `mesh.cell_tags`.

## Notes

- No `tests/meshes/` reference fixture (no case directory is checked in);
  `tests/test_openfoam.py` builds small ASCII/binary `polyMesh` file sets
  inline under `tmp_path`, covering ASCII and binary variants, all 4 named
  volume cell types, general polyhedra, boundary polygon grouping, and the
  `.foam`/case-dir/`polyMesh`-dir path-resolution forms. Most tests import
  the internal Python functions directly (keeping the Python reference
  exercised); the public-API test drives the C++ path.
- Ported from [Simvia's meshlane fork](https://github.com/simvia-tech/meshlane)
  (see `CHANGELOG.md`) — this format did not exist upstream before that.
