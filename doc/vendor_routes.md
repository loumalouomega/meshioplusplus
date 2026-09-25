# Vendor-runtime routes

Some result and model files have no public specification and are readable only by the vendor's own software: Abaqus `.odb`, MSC Marc `.t16`, Femap `.modfem`, Tecplot `.szplt`, and the Ansys results meshio++'s native reader does not cover. meshio++ does not reverse-engineer them and never redistributes vendor libraries. The route instead is a small script that runs inside the vendor's own runtime — the Python a licensed install ships, or its SDK — and writes a file meshio++ reads natively.

| Source | Runs in | Writes | Script | Page |
|---|---|---|---|---|
| Abaqus `.odb` | `abaqus python` (ODB API) | `.pvd` + `.vtu` | `contrib/abaqus_odb/odb_to_vtu.py` | [Abaqus `.odb`](./routes/abaqus_odb.md) |
| MSC Marc `.t16` | Marc/Mentat's Python (PyPost) | `.pvd` + `.vtu` | `contrib/marc_t16/t16_to_vtu.py` | [Marc `.t16`](./routes/marc_t16.md) |
| Ansys results via DPF | Python with `ansys-dpf-core` and a DPF Server | `.vtkhdf` | `contrib/ansys_dpf/rst_to_vtkhdf.py` | [Ansys via PyDPF](./routes/ansys_dpf.md) |
| Femap `.modfem` | Femap's COM API (Windows, `pywin32`) | neutral `.neu` | `contrib/femap/export_neutral.py` | [Femap `.modfem`](./routes/femap_modfem.md) |
| Tecplot `.szplt` | TecIO (a native build), PyTecplot or a Tecplot macro | read directly, or `.plt` | `contrib/tecplot/` | [Tecplot `.szplt`](./routes/tecplot_szplt.md) |

## The `contrib/` convention

- One directory per vendor under [`contrib/`](https://github.com/loumalouomega/meshioplusplus/tree/main/contrib), each script a single file with its usage in its docstring and an MIT header.
- The vendor Pythons (`abaqus python`, Mentat's) have numpy but neither h5py nor meshio++, so those scripts need numpy only, run on Python 2.7 and 3, and write a `.pvd` index with one `.vtu` per step and part through `contrib/common/mio_vtu_series.py` — copy it next to the script. meshio++ reads the result as a time sequence (`meshioplusplus.read_sequence("job.pvd")`) and ParaView opens it directly.
- A route that runs in an ordinary Python (PyDPF) uses meshio++ itself and writes VTKHDF.
- The scripts are not part of the installed package. They are tested in `tests/python/test_contrib_routes.py` against stand-ins of each vendor API — the calls each script makes, as the vendor documents them — which proves the scripts' own logic, not the vendor APIs. The run against each real tool, and what to install for it, is listed in the [roadmap](./roadmap.md#awaiting-a-licensed-run); each page records the vendor version once it has been run.
- An optional C++ plugin against a vendor SDK (the Abaqus ODB C++ API) comes only if a user asks. TecIO, being freely available, is the one SDK meshio++ links today ([`szplt`](./formats/szplt.md)).
