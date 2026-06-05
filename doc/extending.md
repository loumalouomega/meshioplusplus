# Extending meshio

## Registering a custom format at runtime

Use `meshio.register_format` to add a format from outside the meshio package — for example, in application code or a third-party plugin.

```python
import meshio

def my_read(filename):
    # parse the file and return a meshio.Mesh
    ...

def my_write(filename, mesh, **kwargs):
    # serialize mesh to the file
    ...

meshio.register_format(
    "myformat",               # format name used in file_format=
    [".myfmt"],               # file extensions (lowercase, with leading dot)
    my_read,                  # reader function, or None if write-only
    {"myformat": my_write},   # dict mapping format name(s) to writer function(s)
)
```

After calling `register_format`, the format is immediately available through `meshio.read`, `meshio.write`, and the CLI.

A format can expose multiple writer variants under different names:

```python
meshio.register_format(
    "myformat",
    [".myfmt"],
    my_read,
    {
        "myformat":     my_write_v2,   # default
        "myformat-v1":  my_write_v1,   # legacy variant
    },
)
```

## Deregistering a format

```python
meshio.deregister_format("myformat")
```

This removes the format from all internal maps. Useful for overriding a built-in format or in tests.

---

## Adding a new built-in format

Follow the existing module layout under `src/meshio/`:

1. **Create the module directory** `src/meshio/<format>/` with an `__init__.py` that exports `read` and `write`.

2. **Implement `read(filename)`** — return a `meshio.Mesh`.

3. **Implement `write(filename, mesh, **kwargs)`** — serialize the mesh.

4. **Add cell type mappings** if the format uses its own element names. Convention: name them `_<format>_to_meshio_type` and `_meshio_to_<format>_type` in a `common.py` inside the module.

5. **Call `register_format` at module level** (bottom of the main implementation file):

   ```python
   from .._helpers import register_format
   register_format("myformat", [".myfmt"], read, {"myformat": write})
   ```

6. **Import the module in `src/meshio/__init__.py`** — add it to both the import list and `__all__`.

7. **Add `tests/test_<format>.py`** using `helpers.write_read`:

   ```python
   import pytest
   import meshio
   from . import helpers

   @pytest.mark.parametrize("mesh", [helpers.tri_mesh, helpers.tet_mesh])
   def test_myformat(mesh, tmp_path):
       helpers.write_read(tmp_path, meshio.myformat.write, meshio.myformat.read, mesh, atol=1e-15)
   ```

---

## Reader and writer function signatures

```python
def read(filename: str) -> meshio.Mesh:
    ...

def write(filename: str, mesh: meshio.Mesh, **kwargs) -> None:
    ...
```

Readers should return writeable numpy arrays (set `flags["WRITEABLE"] = True` if needed). Writers should not mutate the mesh object.

## Buffers

If your format can operate on open file buffers (not just paths), both `read` and `write` can accept buffer objects. Use `meshio._files.is_buffer(obj, mode)` to detect them:

```python
from meshio._files import is_buffer

def read(filename):
    if is_buffer(filename, "r"):
        return _read_buffer(filename)
    with open(filename, "rb") as f:
        return _read_buffer(f)
```

Note: formats that span multiple files (like TetGen) cannot support buffers.
