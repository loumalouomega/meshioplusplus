# Vendor-runtime routes

Scripts that run inside a vendor's own runtime and write a file meshio++ reads natively, for formats readable only by the vendor's software. See [doc/vendor_routes.md](../doc/vendor_routes.md) for the convention and one page per route.

| Directory | Source | Runs in | Writes |
|---|---|---|---|
| `abaqus_odb/` | Abaqus `.odb` | `abaqus python` | `.pvd` + `.vtu` |
| `marc_t16/` | MSC Marc `.t16` | Marc/Mentat's Python (PyPost) | `.pvd` + `.vtu` |
| `ansys_dpf/` | Ansys results via DPF | Python with `ansys-dpf-core` | `.vtkhdf` |
| `femap/` | Femap `.modfem` | Femap's COM API (Windows) | neutral `.neu` |
| `tecplot/` | Tecplot `.szplt` | PyTecplot or a Tecplot macro | `.plt` |
| `common/` | `mio_vtu_series.py`, the numpy-only series writer the `.pvd` scripts use: copy it next to them | | |

These files are not installed with the Python package. They are tested in `tests/python/test_contrib_routes.py` against stand-ins of each vendor API.
