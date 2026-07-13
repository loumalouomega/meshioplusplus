# Abaqus (`.inp`)

The [Abaqus](http://abaqus.software.polimi.it/v6.14/index.html) input-deck
format: keyword-driven ASCII (`*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`, `*INCLUDE`,
…).

| | |
|---|---|
| **Format name** | `abaqus` |
| **Extensions** | `.inp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.inp")
meshio.abaqus.write("out.inp", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format.
- **`translate_cell_names`** — map Abaqus element names to human-readable ones.

## File structure

`*NODE` gives `id, x, y, z`; `*ELEMENT, TYPE=<abaqus type>` gives
`id, n1, n2, …`. Node/element ids are remapped to contiguous 0-based indices.

## Cell types

The Abaqus element-type maps (`T2D2`, `CPS3`, `C3D8`, `C3D20`, `C3D10`, …) ↔
meshio cell types.

## Data mapping

- `*NSET` / `*ELSET` → `point_sets` / `cell_sets`.

## Notes

- The C++ core handles `*NODE` + `*ELEMENT`. Files using `*NSET`/`*ELSET`/
  `*INCLUDE`, or meshes carrying `point_sets`/`cell_sets`, fall back to the
  Python implementation.
