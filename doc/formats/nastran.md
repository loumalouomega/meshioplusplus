# Nastran (`.bdf`, `.fem`, `.nas`)

The [MSC/NX Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00)
bulk-data format: fixed-width card entries (`GRID`, `CTRIA3`, `CTETRA`, `CHEXA`,
…) in small, large or free field format.

| | |
|---|---|
| **Format name** | `nastran` |
| **Extensions** | `.bdf`, `.fem`, `.nas` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.bdf")
meshio.nastran.write("out.bdf", mesh,
    point_format="fixed-large",   # "fixed-small", "fixed-large", or "free"
    cell_format="fixed-small",
)
```

- **`point_format`** / **`cell_format`** — field layout for `GRID` and element
  cards.

## Cell types

The Nastran ↔ meshio type maps, including second-order elements and the
`CHEXA`/`CPENTA` node reorders.

## Data mapping

- Property/material references → `cell_data["nastran:ref"]`.

## Notes

- The C++ core writes meshio's default layout (fixed-large `GRID*` with the
  16-char Nastran float encoding, fixed-small element cards with continuations).
  Its reader is **sentinel-gated** — it only parses meshio-written files. Reading
  arbitrary/real-world `.fem` files, buffer or long-format reads, non-default
  write formats and `nastran:ref` data use the robust Python implementation.
