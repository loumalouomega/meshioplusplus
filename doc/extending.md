# Extending meshio++

## Registering a custom format at runtime

Use `meshioplusplus.register_format` to add a format from outside the meshio++ package — for example, in application code or a third-party plugin.

```python
import meshioplusplus

def my_read(filename):
    # parse the file and return a meshioplusplus.Mesh
    ...

def my_write(filename, mesh, **kwargs):
    # serialize mesh to the file
    ...

meshioplusplus.register_format(
    "myformat",               # format name used in file_format=
    [".myfmt"],               # file extensions (lowercase, with leading dot)
    my_read,                  # reader function, or None if write-only
    {"myformat": my_write},   # dict mapping format name(s) to writer function(s)
)
```

After calling `register_format`, the format is immediately available through `meshioplusplus.read`, `meshioplusplus.write`, and the CLI.

A format can expose multiple writer variants under different names:

```python
meshioplusplus.register_format(
    "myformat",
    [".myfmt"],
    my_read,
    {
        "myformat":     my_write_v2,   # default
        "myformat-v1":  my_write_v1,   # legacy variant
    },
)
```

## Shipping a format as a package

A separate package adds its format to every program that uses meshio++ by declaring a `meshioplusplus.formats` entry point; installing the package is all it takes. The entry point names a callable taking no arguments that calls `register_format` (or a module that calls it when imported):

```toml
# pyproject.toml of the package `meshioplusplus-foo`
[project.entry-points."meshioplusplus.formats"]
foo = "meshioplusplus_foo:register"
```

```python
# meshioplusplus_foo/__init__.py
import meshioplusplus


def register():
    meshioplusplus.register_format("foo", [".foo"], read, {"foo": write})
```

The entry points are read once, on the first use of the format registry (a `read`, a `write`, `formats()` or a CLI's format list), so `import meshioplusplus` stays as cheap as before and `read("x.foo")` just works. A plugin that fails to load is skipped with a `RuntimeWarning` naming it, and the others still load; `meshioplusplus._plugins.format_plugin_errors()` says why each one failed. `MESHIOPLUSPLUS_NO_PLUGINS=1` skips discovery, as it does for the add-on wheels. Unlike those add-ons (the `meshioplusplus.plugins` group, pinned to the release), a format package is not tied to a meshio++ version: `register_format` is the contract. Such a format is Python-only — the native CLI, the C API and the flat bindings read the formats of the C++ registry.

## Deregistering a format

```python
meshioplusplus.deregister_format("myformat")
```

This removes the format from all internal maps. Useful for overriding a built-in format or in tests.

---

## Adding a new built-in format

Follow the existing module layout under `src/python/meshioplusplus/`:

1. **Create the module directory** `src/python/meshioplusplus/<format>/` with an `__init__.py` that exports `read` and `write`.

2. **Implement `read(filename)`** — return a `meshioplusplus.Mesh`.

3. **Implement `write(filename, mesh, **kwargs)`** — serialize the mesh.

4. **Add cell type mappings** if the format uses its own element names. Convention: name them `_<format>_to_meshio_type` and `_meshio_to_<format>_type` in a `common.py` inside the module.

5. **Call `register_format` at module level** (bottom of the main implementation file):

   ```python
   from .._helpers import register_format
   register_format("myformat", [".myfmt"], read, {"myformat": write})
   ```

6. **Name the module in `src/python/meshioplusplus/__init__.py`** — add it to `__all__` and to the `TYPE_CHECKING` import list. The package does not import format modules itself: `read` and `write` load the one a path or format name needs from a generated table, so run `python tools/gen_format_table.py` to add the new format to `src/python/meshioplusplus/_format_table.py` (`tests/python/test_format_table.py` fails until you do).

7. **Add `tests/python/test_<format>.py`** using `helpers.write_read`:

   ```python
   import pytest
   import meshioplusplus
   from . import helpers

   @pytest.mark.parametrize("mesh", [helpers.tri_mesh, helpers.tet_mesh])
   def test_myformat(mesh, tmp_path):
       helpers.write_read(tmp_path, meshioplusplus.myformat.write, meshioplusplus.myformat.read, mesh, atol=1e-15)
   ```

---

## Reader and writer function signatures

```python
def read(filename: str) -> meshioplusplus.Mesh:
    ...

def write(filename: str, mesh: meshioplusplus.Mesh, **kwargs) -> None:
    ...
```

Readers should return writeable numpy arrays (set `flags["WRITEABLE"] = True` if needed). Writers should not mutate the mesh object.

## Buffers

If your format can operate on open file buffers (not just paths), both `read` and `write` can accept buffer objects. Use `meshioplusplus._files.is_buffer(obj, mode)` to detect them:

```python
from meshioplusplus._files import is_buffer

def read(filename):
    if is_buffer(filename, "r"):
        return _read_buffer(filename)
    with open(filename, "rb") as f:
        return _read_buffer(f)
```

Note: formats that span multiple files (like TetGen) cannot support buffers.
