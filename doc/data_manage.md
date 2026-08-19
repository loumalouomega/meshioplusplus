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

## Other languages

- **C API** — `mio_data_drop`, `mio_data_keep` and `mio_data_rename` each return a new `mio_mesh*`. Name lists cross as `const char* const* names` plus an explicit `int64_t count`. The combined `data_manage` is deliberately not exposed; the three primitives compose. See the [C API reference](/c_api).
- **Fortran** — `m%data_drop(MIO_DATA_POINT, ["T"])`, `m%data_keep(...)`, `m%data_rename(MIO_DATA_POINT, "T", "temperature")`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataDrop(mesh, "point", ["T"], false)`, `dataKeep(...)`, `dataRename(mesh, "point", "T", "temperature")`. See the [WebAssembly reference](/wasm).
