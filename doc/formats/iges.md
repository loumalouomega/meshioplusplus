# IGES (`.igs`, `.iges`)

The Initial Graphics Exchange Specification, read through OpenCASCADE (the `meshioplusplus-occt` add-on wheel, `pip install "meshioplusplus[occt]"`, and the release binaries): the twin of the [STEP](/formats/step#reading) reader for the older CAD exchange format many tools still write. Read-only.

| | |
|---|---|
| **Format name** | `iges` |
| **Extensions** | `.igs`, `.iges` |
| **Read / Write** | ✓ (with OpenCASCADE) / — |
| **Extra dependencies** | OpenCASCADE 7.6 or later (LGPL-2.1 with the OCCT exception; the `meshioplusplus-occt` wheel and the release binaries, on by default from source) |

```python
import meshioplusplus as mp

mesh = mp.read("part.igs")                         # triangles at 1e-3 of the diagonal
mesh = mp.iges.read("part.igs", deflection=0.01)       # a chordal deflection of your own
mp.write("part.cad.json", mesh)                    # the NURBS faces as Kratos CAD JSON
```

`IGESControl_Reader` transfers the file's shape, which is then handled exactly as a STEP file's is: converted to NURBS throughout, triangulated by `BRepMesh` into one `triangle` block with the `cad:*` cell data, and kept as faces, trimming loops and edges in `mesh.info["brep"]` (C++ `BrepInfo`); see [STEP: reading](/formats/step#reading) for the details and the options. Without the backend, reading raises `NotImplementedError` (C++ `Unsupported`) naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON`; there is no Python twin of a CAD kernel. In WASM the [full artifact](/wasm#lean-and-full-builds) (v17.2.0) links the backend and reads IGES; the lean one does not list it among its readers and names the full artifact. Writing IGES is not planned: [STEP](/formats/step) is the exchange format meshio++ writes (see [conformance](../conformance.md#iges)).

As for STEP, the NURBS B-rep in `mesh.info["brep"]` can be refined, evaluated and sampled exactly with the [isogeometric layer](../iga.md), and triangulated at a target size with [`mesh_cad_surface`](../mesh_cad.md). Each edge's `extra` holds OpenCASCADE's tolerances as for STEP (since v17.20.0), which `mesh_cad_surface(model_tolerances=True)` accepts.
