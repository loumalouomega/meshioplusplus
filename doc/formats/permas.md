# PERMAS (`.post`, `.dato`)

The [PERMAS](https://www.intes.de) data-file format: `$`-delimited keyword
sections in plain text, optionally gzip-compressed.

| | |
|---|---|
| **Format name** | `permas` |
| **Extensions** | `.post`, `.post.gz`, `.dato`, `.dato.gz` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.post")
meshio.permas.write("out.post", mesh)
```

`write` takes no keyword arguments.

## File structure

`$COOR` gives node coordinates; `$ELEMENT TYPE=<permas type>` blocks give the
connectivity. `$NSET`/`$ESET` and other keywords are read but not mapped into the
mesh (meshio drops them).

## Cell types

The full PERMAS ↔ meshio type table, including second-order elements
(`triangle6`, `tetra10`, `quad9`, `wedge15`) with the appropriate node reorders.

## Notes

- The C++ core handles the plain-text `.post`/`.dato` files. The gzip `.gz`
  containers fall back to the Python implementation.
