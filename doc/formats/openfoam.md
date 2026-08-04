# OpenFOAM polyMesh

A reader and writer for [OpenFOAM](https://www.openfoam.com/)'s native `polyMesh` representation: an unstructured, face-based mesh described by 4-5 sibling files (`points`, `faces`, `owner`, `neighbour`, `boundary`) rather than a single file. Supports both ASCII and binary (little-endian, `label=32/64`, `scalar=32/64`) encodings, and reconstructs general polyhedra as well as tetra/pyramid/wedge/hexahedron.

| | |
|---|---|
| **Format name** | `openfoam` |
| **Extensions** | `.foam` |
| **Read / Write** | ✓ / ✓ (ASCII) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("case.foam")                 # a <case>/case.foam marker file
mesh = meshioplusplus.openfoam.read("/path/to/case")     # or the case directory directly
mesh = meshioplusplus.openfoam.read("/path/to/constant/polyMesh")  # or polyMesh directly
```

```python
meshioplusplus.write("out/case.foam", mesh)              # infers the format from `.foam`
meshioplusplus.openfoam.write("/path/to/case", mesh)     # a case root (no extension)
```

`read(filename)` takes no keyword arguments; the `polyMesh` directory is located from whichever of the three input forms above is given (`_resolve_polymesh`): a `.foam` suffix looks for `<parent>/constant/polyMesh`; a directory literally named `polyMesh` is used as-is; any other directory is checked for `constant/polyMesh` then `polyMesh` as subdirectories. A `FileNotFoundError` is raised if none match.

## File structure

Each `polyMesh` file (`points`, `faces`, `owner`, `neighbour`, `boundary`) is an OpenFOAM "FoamFile": a header block declaring `format` (`ascii`/`binary`) and (for binary) an `arch` string encoding `label=32|64`/`scalar=32|64`, followed by a count and a parenthesized list.

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

- **Header detection** (`_detect_format`) scans line-by-line for a `format ...;` line and an `arch "...";` line, extracting `label=` / `scalar=` byte widths; defaults to ascii/8-byte if absent.
- **Binary `points`** (`vectorField`): `N (` followed directly by `N*3*scalar_bytes` raw floats, no per-row framing — read via a single `np.frombuffer`.
- **Binary `owner`/`neighbour`** (`labelList`): `N (` followed by `N*label_bytes` raw ints, same direct-buffer read.
- **Binary `faces`** (`faceList`, non-contiguous): each face is its own `labelList` — `<count> ( <count*label_bytes bytes> )` repeated `N` times. Read in two passes: a sequential ASCII scan locating each face's byte offset and node count (cheap, since `find(b"(")` only ever scans the short ASCII gap between faces, never binary data that might coincidentally equal `'('`), then a single vectorized gather of every face's binary blob via a cumulative byte mask — bounded peak memory even for tens of millions of faces (see `_RaggedArray`, a CSR-style `(conn, offsets)` pair standing in for a `list[list[int]]`).
- **ASCII** variants use simple regex/line-based parsing (`_parse_points_ascii`, `_parse_faces_ascii`, `_parse_int_list_ascii`) after comment-stripping (`/* */` and `//`) and header-skipping.
- **`boundary`**: a dict of `patch_name → {type, nFaces, startFace}`, parsed via a brace-matching regex over the whole (header/comment-stripped) text.
- **Cell reconstruction**: cell↔face topology is built once as a CSR `_RaggedArray` (`_cell_faces_csr`, vectorized via `argsort`+`bincount`) from `owner`/`neighbour`; per cell, each face is oriented outward (reversed if the cell is that face's neighbour, since the stored normal points owner→neighbour) and classified by `(n_faces, n_points)`: `(4,4)→tetra`, `(5,5)→pyramid`, `(5,6)→wedge`, `(6,8)→hexahedron`, anything else → a general `polyhedron` (kept as outward-oriented face lists). Each of the 4 named types has a dedicated orientation-fixing builder (`_build_tetra`/`_build_pyramid`/`_build_wedge`/`_build_hexahedron`) that computes a scalar triple product and flips the node order if it comes out negative, guaranteeing positive-volume connectivity regardless of the source mesh's face-normal convention.

## Cell types

Volume cells: `tetra`, `pyramid`, `wedge`, `hexahedron`, and general `polyhedron<N>` (grouped by **unique node count** `N`, ragged per-cell face lists stored as an `object`-dtype array — one `CellBlock` per distinct `N`).

Boundary (patch) faces: `triangle`, `quad`, and `polygon<N>` for `N > 4` (grouped by vertex count `N`, one `CellBlock` per size — via `_build_boundary_polygons`).

## Data mapping

- `cell_data["cell_tags"]` — per-cell-block tag array: `0` for every volume cell block, and a distinct negative "MED-style family id" `-(patch_index+1)` per boundary patch's face blocks (so a triangle patch and a quad patch on the *same* physical boundary would currently get *different* tag values — see Quirks).
- `mesh.cell_tags` — mesh-level attribute (not `cell_data`), `{family_id: [patch_name]}`, letting a MED write bridge these patch names through the same mechanism used for Gmsh physical groups (see [`med.md`](med.md)).
- `mesh.point_tags` — always set to `{}` (present for interface symmetry with the MED-derived tag convention; OpenFOAM has no point-tag concept).
- No point_data or field_data (OpenFOAM field files like `U`, `p`, `T` in the case's time directories are not read by this module — only the mesh topology under `constant/polyMesh`).

## Quirks & limitations

- **The only meshio++ writer that creates a directory.** `write` resolves its path exactly as `read` does and creates `<case>/constant/polyMesh/` as needed; a `.foam` target also gets its (empty) marker file written, which is what makes the case openable by ParaView. Because a case *directory* has no extension, that form needs an explicit `file_format="openfoam"` — `resolve_format` is a pure string function and deliberately does not stat the filesystem.
- **Write is ASCII only.** A binary polyMesh is a documented follow-up, so an explicit binary request fails by name rather than silently writing ASCII.
- **There is no Python fallback writer.** A twin would have to re-implement the per-cell winding repair — a discrete branch on the sign of an enclosed volume — and two implementations of such a branch can land on opposite sides for a near-degenerate cell, the same reasoning that keeps `smooth`'s inversion guard out of its numpy fallback. It would also be dead code: `openfoam_write` needs no optional dependency, so it ships in every wheel. The writer is registered only when the compiled core provides it.
- **Face ids are not preserved by a round trip.** They come from meshio++'s own cell-major deduplication, not from the source file's face order. The invariant the writer guarantees is topological (the ordering contract below), not byte-level.
- **A patch's `type` is downgraded when it needs companion entries meshio++ does not carry** (`cyclic`, `cyclicAMI`, `processor`, `mapped*`, …): those declare keys like `neighbourPatch` or `myProcNo`, and OpenFOAM refuses to *load* a case whose patch declares such a type without them. A downgraded case loads and solves with boundary conditions you can see and fix. An unknown type is written as `patch`, never `wall` — `wall` selects wall functions, so guessing it would silently change a solve's physics.
- **A boundary cell coinciding with an *internal* face is dropped with a warning.** OpenFOAM cannot put a face shared by two cells on a patch, so writing it would produce a case that does not load.
- Degenerate volume cells that match a named type's `(n_faces, n_points)` signature but whose topology doesn't resolve cleanly (`_match_top` finds more or less than one vertical neighbour per base node) are **silently skipped** and logged as a warning count, rather than falling back to a general polyhedron.
- Boundary patches are tagged by **patch index**, not patch identity across face-size groups — if one named patch contributes both triangles and quads, its triangle `CellBlock` and quad `CellBlock` get the *same* `cell_tags` id (assigned once per patch, reused across whichever size-buckets that patch's faces fall into), but two *different* named patches always get distinct ids.
- All binary reads assume little-endian (`LSB`) — the format's own `arch` string is trusted for label/scalar width but not for byte order.
- Read goes through the C++ core (`meshioplusplus._core.openfoam_read`, using `std::filesystem` for the polyMesh directory), with the Python reference as an automatic fallback. General polyhedra cross the C++↔Python boundary via the ragged `polyhedron<N>` cell representation (a copied list of face arrays); boundary patch names travel through an `OpenFoamInfo` side-channel struct as `mesh.cell_tags`.

## The ordering contract (write)

A polyMesh is not just a face list: OpenFOAM requires an ordering, and breaking it produces a mesh `checkMesh` rejects rather than an error at write time. All of it is enforced and then **re-validated** before any file is opened, with one check per clause naming the clause it broke. The validator runs in release builds too — release is where large cases get written, and its failure means a corrupt mesh was about to reach a solver.

1. `neighbour` holds **only** internal faces. (meshio++'s own reader also accepts a `-1`-padded full-length list; OpenFOAM does not — which is exactly why round-tripping through our reader is a weak oracle for this writer.)
2. `owner[i] < neighbour[i]` for every internal face.
3. Internal faces occupy `[0, nInternalFaces)`; boundary faces follow.
4. Internal faces are sorted by `(owner, neighbour)` — strictly upper-triangular.
5. Each internal face's normal points from its owner toward its neighbour; each boundary face's points out of the domain.
6. Patches partition the boundary range exactly, and every face in a patch's range really belongs to that patch.

Winding is **repaired**, not required: an inverted (negative-Jacobian) cell is written correctly oriented rather than rejected, because `detail/cell_faces.hpp`'s tables are outward on the *reference* element and would otherwise emit inward normals for every one of its faces.

## Patches on write

Patch names and types are recovered from `mesh.cell_tags`' negative family ids together with `mesh.openfoam_patch_types` — both of which `read` sets, so an OpenFOAM → OpenFOAM round trip preserves the `boundary` file's patch order.

A mesh carrying neither, which is anything converted from another format, gets a single patch named **`defaultFaces`** of type `patch` — OpenFOAM's own name for exactly this, and what `blockMesh` produces. Patches are never synthesized from geometry (connected boundary components, cell-block membership, …): that would invent structure the mesh never stated and make the patch list depend on the mesh's shape.

## Notes

- No `tests/python/meshes/` reference fixture (no case directory is checked in); `tests/python/test_openfoam.py` builds small ASCII/binary `polyMesh` file sets inline under `tmp_path`, covering ASCII and binary variants, all 4 named volume cell types, general polyhedra, boundary polygon grouping, and the `.foam`/case-dir/`polyMesh`-dir path-resolution forms. Most tests import the internal Python functions directly (keeping the Python reference exercised); the public-API test drives the C++ path.
- Ported from [Simvia's meshlane fork](https://github.com/simvia-tech/meshlane) (see `CHANGELOG.md`) — this format did not exist upstream before that.
