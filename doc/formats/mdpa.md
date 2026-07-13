# Kratos / MDPA (`.mdpa`)

The [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data)
model-part data format (`.mdpa`): block-structured ASCII (`Begin … End …`).

| | |
|---|---|
| **Format name** | `mdpa` |
| **Extensions** | `.mdpa` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.mdpa")
meshio.mdpa.write("out.mdpa", mesh, float_fmt=".16e", binary=False)
```

- **`float_fmt`** — coordinate format.
- **`binary`** — currently ASCII only.

## File structure

`Begin ModelPartData … End`, `Begin Nodes … End`, `Begin Elements <Type> … End`,
`Begin Conditions … End` and `Begin SubModelPart … End` blocks.

## Cell types

The Kratos ↔ meshio element/condition type maps.

## Data mapping

- Sub-model-parts → `cell_sets`.

## Notes

- Implemented in pure Python (no C++ core path).
