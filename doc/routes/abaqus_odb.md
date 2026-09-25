# Abaqus `.odb`

The Abaqus output database has no public specification and is version-locked (`abaqus upgrade -odb` moves an old one forward); only the ODB API of a licensed Abaqus reads it. The route is `contrib/abaqus_odb/odb_to_vtu.py`, run by Abaqus's own Python. Where no licence is available, request the `.fil` results file in the input deck (`*FILE FORMAT`, `*NODE FILE`, `*EL FILE`): meshio++ reads it natively ([Abaqus `.fil`](../formats/abaqus_fil.md)).

## Running it

Copy `odb_to_vtu.py` and `contrib/common/mio_vtu_series.py` into one directory, then:

```bash
abaqus python odb_to_vtu.py job.odb                    # job.pvd + job/
abaqus python odb_to_vtu.py job.odb out.pvd --fields U,S,PEEQ --steps Step-1
```

```python
import meshioplusplus

for time, mesh in meshioplusplus.read_sequence("job.pvd"):
    mesh.point_data["U"], mesh.cell_data["S"]
```

It runs on the Python 2.7 of Abaqus 2023 and before and the Python 3 of Abaqus 2024 and later, with the numpy Abaqus bundles.

## What it writes

- **Instances → parts.** Each part instance is a `.pvd` part named after it, so a read merges them with one region per instance (`part=` selects one, see [PVD](../formats/pvd.md#two-axes-time-and-part)).
- **Frames → steps.** A time-domain frame's time is its step's total time at the start plus the frame value; frame 0 of a later step, the previous step's last state, is skipped. A frequency or modal frame (whose frame value is a frequency or an eigenvalue) comes one after the previous time, so the times keep ascending.
- **Field outputs by position.** `NODAL` output (U, RF, NT…) is point data; integration-point, centroidal and whole-element output (S, E, PEEQ, ELEN…) is cell data from Abaqus's own `CENTROID` (or `WHOLE_ELEMENT`) subset, averaged over section points; output with only `ELEMENT_NODAL` values is averaged over the elements sharing a node. Components keep Abaqus's order (`S11, S22, S33, S12, S13, S23`).
- **Labels and sets.** `abaqus:node_label` and `abaqus:element_label` carry the labels; every node and element set, of the instance or of the assembly, is a `set:<name>` 0/1 mask on every part.
- **Elements.** Continuum (`C3D`, `DC3D`, `AC3D`, `COH3D`, `SC`), planar and axisymmetric (`CPE`, `CPS`, `CAX`, `CGAX`…), shell and membrane (`S4R`, `STRI3`, `M3D`, `R3D`, `SFM3D`…) and line (`T2D`, `T3D`, `B2`, `B3`, `PIPE`…) elements map by family and node count as the native `.inp` and `.fil` readers map them (Abaqus lists every such shape's nodes in meshio++'s order); connectors, springs, masses and other elements with no cell are skipped with a count.

## Checked against

The script's logic is tested against a stand-in of the ODB API (`tests/python/test_contrib_routes.py`). The run against Abaqus itself is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run). The references are [ODB2VTK](https://github.com/Arris-Composites/ODB2VTK) and the Abaqus Scripting Reference's `odbAccess` chapter.
