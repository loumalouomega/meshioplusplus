# Data array management (rename / drop / keep)

`data_rename`, `data_drop` and `data_keep` rewrite **which** data arrays a mesh carries and under what names. Values, dtypes and shapes are copied verbatim, and the geometry comes through bit-identical. They are [data operations](/data_operations), not file formats.

```python
import meshioplusplus as mp

mesh = mp.read("part.vtu")

mesh = mp.data_rename(mesh, "point", "T", "temperature")   # rename one array
mesh = mp.data_drop(mesh, "point", ["scratch", "tmp"])     # remove by name
mesh = mp.data_keep(mesh, "cell", ["material"])            # keep only these
```

For several changes at once, `data_manage` applies them in one pass and reports what it did:

```python
result = mp.data_manage(
    mesh,
    keep=[("point", "T"), ("point", "p")],
    drop=[("cell", "scratch")],
    rename=[("point", "T", "temperature")],
)
result["mesh"]      # the rewritten mesh
result["dropped"]   # ["cell_data:scratch", "point_data:u", ...] (sorted)
result["renamed"]   # [("point_data:T", "point_data:temperature")]
```

## Phase order

The three phases apply in a fixed, documented order: **keep**, then **drop**, then **rename**.

`keep` is a whitelist that only affects the locations it actually mentions — so keeping some `point_data` arrays leaves `cell_data` and `field_data` completely untouched. `data_keep(mesh, "point", [])` is meaningful: it drops every `point_data` array while leaving the other locations alone.

## Errors

An unknown key raises `ValueError` listing every key that *is* present:

```
meshio++: no point_data array named 'foo' (available: T, p, u)
```

Pass `ignore_missing=True` to skip absent keys instead. A rename whose target already exists also raises, unless that target is itself renamed away or dropped in the same pass — which makes swapping two names legal:

```python
# Legal: neither target survives under its old name.
mp.data_manage(mesh, rename=[("point", "a", "b"), ("point", "b", "a")])
```

Two renames targeting the same name, or two renames of the same source, always raise.

## CLI

```bash
meshioplusplus data rename in.vtu out.vtu --point T:temperature
meshioplusplus data drop   in.vtu out.vtu --point a,b --cell c
meshioplusplus data keep   in.vtu out.vtu --point T,p --cell mat
```

`--point`, `--cell` and `--field` select the location. For `rename` the `OLD:NEW` value is split on the **last** colon, because data names routinely contain colons — so `--point gmsh:physical:tag` renames `gmsh:physical` to `tag`. `--point` is repeatable for `rename`; for `drop`/`keep` it takes a comma-separated list. See the [CLI reference](/cli#meshioplusplus-data).

## Sets ↔ integer data {#sets--integer-data}

Python's mutating `mesh.point_sets_to_data(join_char="-")` and `mesh.cell_sets_to_data(data_name=None)` turn region-backed membership into scalar integer labels; `mesh.point_data_to_sets(key)` and `mesh.cell_data_to_sets(key)` reverse it. Native dispatch preserves Python's set insertion order and keeps geometry objects, `mesh.info`, non-numeric field data and other Python-only metadata in place. Non-membership cell-set metadata (for example `gmsh:bounding_entities`) stays on the reference path.

Each set receives its zero-based position as its label; later overlapping sets overwrite earlier ones, and uncovered entities receive -1 with a warning. The output field replaces any same-name field and converted point/cell regions are removed; side regions survive. Reversing requires a one-dimensional integer field, sorts its distinct tags numerically (including -1), removes duplicate names from splitting the key on `-`, and uses those names if their count matches the tags. Otherwise it uses `set-<key>-<tag>` for cells and the historical `set-key-<tag>` for points. It removes the source field and retains dimension/tag metadata when replacing an existing same-name region. Empty sets cannot be reconstructed from labels alone.

The native core exposes `sets_to_data(mesh, location, name, join, order)` and `data_to_sets(mesh, location, key)` in `operations/data_manage.hpp`; both return a new mesh and preserve geometry, unrelated data and property sets bit-identically. Only point/cell locations are accepted. Omit `order` for native region-name order, or supply every distinct set name exactly once to choose label order; `name` is optional (including an explicitly empty name) and defaults to joining set names with `join`. Integer labels are compared without converting through floating point; the usual mesh-backend dtype canonicalization still applies.

The matching flat APIs are `mio_sets_to_data` / `mio_data_to_sets` in C, `m%sets_to_data` / `m%data_to_sets` in Fortran, `sets_to_data` / `data_to_sets` in Julia, `mio_sets_to_data` / `mio_data_to_sets` in R, and `setsToData` / `dataToSets` in WASM. Labels are values, not indices: Fortran/Julia/R do not shift 0 or -1. Both CLIs accept `convert -s` / `convert -d` for single meshes; `-d` converts every point/cell field and fails on non-integer fields. Neither flag combines with selective reads or sequences; use the pipeline `SetsToData` / `DataToSets` steps for a sequence. MCP's `sets_data` tool uses the same semantics on sandboxed file paths.

## Other languages

- **C API** — `mio_data_drop`, `mio_data_keep` and `mio_data_rename` each return a new `mio_mesh*`. Name lists cross as `const char* const* names` plus an explicit `int64_t count`. The combined `data_manage` is deliberately not exposed; the three primitives compose. See the [C API reference](/c_api).
- **Fortran** — `m%data_drop(MIO_DATA_POINT, ["T"])`, `m%data_keep(...)`, `m%data_rename(MIO_DATA_POINT, "T", "temperature")`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataDrop(mesh, "point", ["T"], false)`, `dataKeep(...)`, `dataRename(mesh, "point", "T", "temperature")`. See the [WebAssembly reference](/wasm).
